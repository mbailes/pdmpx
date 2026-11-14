# pdmpx

**Fast period/DM optimisation and candidate diagnostics for folded pulsar archives.**

pdmpx does the job of PSRCHIVE's `pdmp` (find the period and DM that maximise a folded candidate's S/N), but
thousands of times faster, and then goes further: for every archive it writes one CSV row of numbers that
separate pulsars from RFI and noise, and optionally a compact plot file and a diagnostic page.
It is built for pipelines that fold thousands of FFT-search candidates and need to triage them.

- **Standalone**: no PSRCHIVE needed. Reads **Timer** archives (dspsr's default `.ar`) and **fold-mode PSRFITS**.
- **Fast**: about 2 candidates per second per core for a 24 × 1024 × 128 fold (131,584 period/DM trials),
  ~20 per second on a 16-core machine. A 1024-channel archive that pdmp took 2.5 hours on takes about 1 s.
- **Batch mode**: `pdmpx -outdir results *.ar` - one CSV row per archive, parallel over archives, resumable,
  with a progress meter.
- **Diagnostics** for each candidate above an S/N threshold: S/N at DM 0 vs a real dispersed pulse,
  measured vs predicted noise, single-pulse (per-subint) significance, drift-scan beam transit, spectrum
  concentration, and flags for broad periodic RFI.

Full option and output reference: **[docs/HELP.md](docs/HELP.md)**.

## Building

You need a C++17 compiler (GCC ≥ 8 or clang), **FFTW3 in single precision** (`libfftw3f`), **cfitsio**,
**zlib**, and OpenMP (comes with GCC; optional, but without it pdmpx uses one core). PSRCHIVE is *not* needed.

```sh
git clone https://github.com/mbailes/pdmpx.git
cd pdmpx
make                     # builds ./pdmpx
./pdmpx -h               # options
make install             # optional: copies pdmpx and tools/*.py to ~/.local/bin (PREFIX=... to change)
```

### Ubuntu / Debian (your own machine)
```sh
sudo apt install g++ make pkg-config libfftw3-dev libcfitsio-dev zlib1g-dev
make
```

### OzSTAR / NT (Swinburne)
NT runs AlmaLinux 9 with Lmod/EasyBuild modules (and tcsh as the default shell - the commands below work in
tcsh and bash). The same FFTW and cfitsio as dspsr:
```sh
module load gcc/11.3.0 fftw/3.3.10 cfitsio/4.2.0
git clone https://github.com/mbailes/pdmpx.git
cd pdmpx
make
./pdmpx -h
```
The Makefile picks the libraries up from the modules (`$EBROOTFFTW`, `$EBROOTCFITSIO`) and records their
paths in the binary, so `pdmpx` still runs later without the modules loaded. For the Python tools:
`module load python-scientific/3.10.4-foss-2022a` (numpy + matplotlib).

### Another cluster without root
Check whether the development headers are installed system-wide (then a plain `make` works):
```sh
ls /usr/include/fftw3.h /usr/include/fitsio.h
```
If not, load FFTW and cfitsio from the module system (`module spider fftw`, `module spider cfitsio`; if you run
dspsr there, `ldd $(which dspsr)` shows which ones it uses), or give their install directories:
```sh
make FFTW_DIR=/path/to/fftw CFITSIO_DIR=/path/to/cfitsio
```
or install them yourself with conda (no root needed):
```sh
conda create -n pdmpx -c conda-forge fftw cfitsio zlib compilers make
conda activate pdmpx
make FFTW_DIR=$CONDA_PREFIX CFITSIO_DIR=$CONDA_PREFIX
```

### Red Hat family (AlmaLinux / Rocky / RHEL 9, with root)
```sh
sudo dnf install epel-release && sudo dnf --enablerepo=crb install gcc-c++ make fftw-devel cfitsio-devel zlib-devel
make
```

### macOS (Homebrew)
```sh
brew install fftw cfitsio libomp
make
```

### Notes
- FFTW must include the **single-precision** library (`libfftw3f`). Some installations are double only.
- The default flags are portable (`-O3`). To tune for one CPU type: `make CXXFLAGS="-O3 -march=native -std=c++17"`,
  but then build on the same kind of node you run on (a login node's CPU may differ from the compute nodes').
- `make OPENMP=0` builds without OpenMP; `make CXX=g++` picks the compiler; `make clean` starts again.
- The Python tools in `tools/` need Python 3 with numpy and matplotlib (`drift_psrs.py` also uses `psrcat`).

Tested: macOS 26 (Apple clang + libomp), Ubuntu 24.04 (GCC 13) and AlmaLinux 9.8 (GCC 11), identical results.

## Quick start

One archive (a pdmp-style report):
```sh
./pdmpx candidate.ar
```
```
pdmpx: J0405+3347 nsub=22 nchan=1024 (searched 32) nbin=128  DM: 512  P0: 257  = 131584 trials on 16 threads
grid: sensible  DM 24.8201..82.2834 step 0.112453   P 63.86697143..64.05282857 ms step 0.726004 us
Best S/N = 22.17  width = 3 bins   (coho's definition: S/N 18.269, width 3 bins)
Best TC Period (ms) = 63.95772199  Correction (ms) = -0.002178013291  Error (ms) = 0.0001444
Best DM = 53.6080  Correction = 0.0000  Error = 0.04746
Noise: rms predicted 0.00034834, measured 0.00039585 (ratio 1.136), chi2_red 4.964
Dispersion: S/N at DM 0 = 3.432 (a real pulse: 2.820), dispersed at 13.30 sigma [256 channels]
Subints: brightest S/N 13.727 (#17, 12.93 sigma), modulation index 1.119, t_frac90 0.455, f_frac90 0.422
Drift transit: peak t = 34.426 +- 0.285 s (MJD 59685.30706445), FWHM 14.583 s
time: load 0.002 s, prepare 0.007 s, search 0.060 s, total 0.149 s
```

Thousands of archives (one CSV row each; plot files for S/N ≥ 7):
```sh
./pdmpx -snr 7 -outdir results /data/folds/*.ar
python3 tools/pdmpx_summary.py results/pdmpx.csv      # ranked table, flags, RFI period families, histograms
python3 tools/pdmpx_plot.py results/*.pdx             # diagnostic PNG pages (or add -png to the pdmpx run)
```
The run can be interrupted and restarted: archives already in `results/pdmpx.csv` are skipped. On a cluster,
run it inside a job with as many cores as you like (`-j N` or `OMP_NUM_THREADS`); with many archives each
core takes one archive at a time.

## What it does

1. **Load** the archive (the whole file is mapped at once and subints are decoded in parallel), form total
   intensity, dedisperse at the archive DM at full frequency resolution.
2. **Scrunch** to at most 32 channels (enough to see the DM curve), keeping a ≥128-channel copy for the
   dispersion test.
3. **Search** period and DM over **±1 turn of drift** (DM across the band, period over the observation), never
   below DM 0 and never more than 100 above the input DM (`-dmup`). Each profile is Fourier-transformed once;
   every trial is a phase rotation and a sum. S/N is pdmp's (PSRCHIVE's `smooth_mmw` boxcar matched filter with
   pdmp's noise estimate), so results compare directly with pdmp (`-natural -binrot` reproduces pdmp exactly).
4. **Diagnose** (above `-snr`): S/N vs DM down to 0 on the ≥128-channel copy, S/N vs period and width,
   per-subint and per-channel S/N, drift-scan transit fit and the RFI flags; write the CSV row and plot file.

See [docs/HELP.md](docs/HELP.md) for every option, every CSV column and how to read the diagnostic page, and
[BENCHMARK.md](BENCHMARK.md) for speed comparisons with pdmp.

## Tools

| | |
|---|---|
| `tools/pdmpx_summary.py` | text summary of a `pdmpx.csv`: ranked table, triage flags, harmonic period families (RFI), histograms |
| `tools/pdmpx_plot.py` | render plot files (`.pdx`) as diagnostic PNG pages (in parallel) |
| `tools/fil_window.py` | cut a time window (by MJD or seconds) out of a SIGPROC filterbank, with a corrected header |
| `tools/drift_psrs.py` | which catalogued pulsars drift through which beams of a drift scan, and when (from the headers) |

## Validation

- Reads Timer and PSRFITS identically, and matches PSRCHIVE's `pdv`, including 4-pol (AA+BB) total intensity.
- With `-natural -binrot`, reproduces pdmp's S/N on every trial (pdmp rounds all shifts to whole bins).
- The default fractional shifts agree with PSRCHIVE's `pam -d DM -D`.

## Licence

MIT - see [LICENSE](LICENSE).
