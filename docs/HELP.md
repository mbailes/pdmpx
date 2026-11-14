# pdmpx help

```
pdmpx [options] archive [archive ...]
```
One archive and no `-outdir`: a pdmp-style report on the terminal. Several archives, or `-outdir`: batch mode,
one CSV row per archive.

## Options

### Search
| option | default | meaning |
|---|---|---|
| `-mc N` | ≤32 | scrunch to N channels before searching (must divide nchan). Default: the largest divisor of nchan ≤ 32 |
| `-ms N` | all | scrunch to N subints |
| `-mb N` | all | scrunch to N phase bins (useful for 1024-bin folds: the trial grid scales as nbin²) |
| `-turns T` | 1 | half-range of the search in turns of drift: DM across the band, period over the observation |
| `-dmup X` | 100 | never search more than X above the archive DM (`-1`: no limit). Never below DM 0 |
| `-ds/-dr/-do X` | | DM step / half-range / offset (override the defaults) |
| `-ps/-pr/-po X` | | period step / half-range / offset, in microseconds |
| `-mw S` | 0 | minimum boxcar width in seconds (default 1 bin) |
| `-natural` | | pdmp's own "natural" grid and no automatic scrunch (for comparison with pdmp) |
| `-binrot` | | round every phase shift to whole bins, exactly as pdmp does (default: exact fractional shifts) |
| `-P S` | | folding period in seconds, for a PSRFITS file that carries none |

### Batch and diagnostics
| option | default | meaning |
|---|---|---|
| `-outdir D` | | write `D/pdmpx.csv` (appended; archives already in it are skipped) and plot files `D/<archive>.pdx` |
| `-snr X` | 0 | full diagnostics and a plot file only when S/N ≥ X. Below it the CSV row has the search results only |
| `-png` | | after the run, render a PNG page from each plot file (`tools/pdmpx_plot.py`, needs matplotlib) |
| `-dispch N` | 128 | minimum channels kept for the DM-0 / dispersion test (raised automatically, see below) |
| `-j N` | all cores | threads. With several archives each thread takes one archive at a time |

### Single-archive extras
| option | meaning |
|---|---|
| `-grid F` | write the full S/N grid as text: `period_us DM S/N`, one DM per block |
| `-profile F` | write the best profile: `bin amplitude` |
| `-dump s c` | print one raw profile (subint s, channel c) as loaded, then exit (debugging the readers) |
| `-v` | verbose (archive header details) |

### Environment variables
| variable | meaning |
|---|---|
| `OMP_NUM_THREADS` | threads, if `-j` is not given |
| `PDMPX_PLOT` | path to `pdmpx_plot.py` for `-png` (default: `tools/` next to the binary) |
| `PDMPX_TRACE=file` | log each archive as it starts and ends (to find an archive that crashes a batch) |
| `PDMPX_IO=read` | read files with one `read()` instead of `mmap` |
| `PDMPX_CFITSIO=1` | read PSRFITS through cfitsio row by row instead of directly (slower; for checking) |

## The search

1. The archive is dedispersed at its own DM at full frequency resolution, then scrunched to ≤32 channels.
2. Trials: DM from (input DM − 1 turn) to (input DM + 1 turn) in steps of ½ bin of delay across the band,
   clipped to [0, input DM + 100]; period ±1 turn of drift over the observation in 1-bin steps. For
   128 bins that is up to 512 × 257 = 131,584 trials.
   For slow pulsars one turn of DM across the band can be hundreds of DM units - hence the `-dmup` cap.
3. S/N of every trial is pdmp's: the PSRCHIVE `smooth_mmw` boxcar matched filter (baseline = the lowest
   running mean over 40% of the turn; boxcars from 1 bin growing ×1.5 to nbin/2) with pdmp's noise estimate
   (the rms of each profile's quietest half, propagated with the weights). Channel and subint weights are used
   throughout; zero-weight data are skipped.

## The CSV (`pdmpx.csv`)

One row per archive; columns are fixed and named with their units. Empty = not computed (below threshold,
or not applicable).

**Identity and data**
| column | meaning |
|---|---|
| `file`, `status` | path; `ok`, `below_threshold` or `error: ...` (every input gets a row) |
| `source`, `mjd`, `tobs_s` | source name, start MJD, length (s) |
| `fc_MHz`, `bw_MHz` | centre frequency, bandwidth |
| `nsub`, `nchan`, `nbin` | archive dimensions |
| `nchan_search`, `nchan_disp` | channels searched; channels used for the dispersion test |
| `zap_frac` | fraction of (subint, channel) weights that are zero |

**The fit**
| column | meaning |
|---|---|
| `P_fold_ms`, `DM_fold` | the folding period and DM (the input) |
| `P_ms`, `P_err_ms`, `DM`, `DM_err` | best period and DM; errors from the curvature of S/N at the peak |
| `dP_ms`, `dDM` | best − folded. Large moves (especially in DM) at low S/N are usually noise |
| `edge` | `0`, or `DM`/`P`/`DM+P` if the best trial is within 2% of the grid edge (the true peak may be outside) |
| `n_trials` | DM × period trials searched |

**The pulse**
| column | meaning |
|---|---|
| `snr`, `width_bins`, `width_ms`, `duty` | pdmp-definition S/N and its boxcar width |
| `snr_coho`, `width_coho` | the same best profile scored with coho's definition (median baseline, 1.4826×MAD over all bins, boxcar centred on the peak) - typically ~15-20% lower than `snr` |
| `w_maxed` | 1 if S/N is still rising at the widest boxcar: a sinusoid-like or baseline signal, not a pulse (but some wide-profile MSPs look like this) |

**Noise quality** - the best RFI detectors
| column | meaning |
|---|---|
| `baseline`, `rms_pred`, `rms_meas` | baseline level; the predicted rms (pdmp's) and the rms measured in the profile's baseline window |
| `rms_ratio` | measured/predicted. ≈1 for clean data; ≫1 when the "baseline" is full of other spikes (e.g. periodic RFI folded at a multiple of its period) |
| `chi2_red` | reduced χ² of the profile against a flat line |

**Dispersion** (above threshold; computed on ≥`-dispch` channels)
| column | meaning |
|---|---|
| `snr_dm0` | S/N at DM 0, same period and width |
| `snr_dm0_model` | what a real pulse of this width and spectrum would give at DM 0 (plus the noise floor) |
| `disp_sigma` | (S/N at best DM − S/N at DM 0)/√2: ≫3 for a dispersed source, ≤0 for terrestrial RFI |
| `dm_peak2` | second separate peak of S/N vs DM relative to the highest. Near 1 = a comb of DM peaks: a broad periodic signal re-aligning every turn (RFI), not a pulsar |

The dispersion test uses at least 128 channels, raised until the archive-DM sweep within one channel is
below half the pulse width. With fewer channels, zero-DM RFI stays smeared inside each channel and can look
"dispersed".

**Time and frequency behaviour** (above threshold)
| column | meaning |
|---|---|
| `sp_max`, `sp_max_sub` | S/N of the brightest single subint in the on-pulse window, and which subint |
| `sp_sig` | its significance against same-width off-pulse windows, allowing for picking the brightest of nsub (empty for pulses too wide to leave off-pulse windows) |
| `mod_index` | std/mean of the per-subint S/N |
| `t_frac90` | fraction of subints holding 90% of the signal: ≈1 steady; ≪1 bursty / one transit / impulsive RFI |
| `f_frac90` | fraction of channels holding 90% of the signal: ≈1 broadband; ≪1 narrowband RFI (or one scintle) |
| `drift_t_s`, `drift_t_err_s`, `drift_fwhm_s`, `drift_mjd_peak` | Gaussian beam-transit fit to S/N vs time, when significant (drift scans) |

**Reserved, empty for now**: `n_eff`, `fap`, `sigma` (noise-search false-alarm calibration), `fft_snr`,
`fold_gain` (FFT-search S/N, to be joined from the candidate list).

**Bookkeeping**: `plotfile`, `runtime_s`, `version`, `options`.

## Triage flags (`tools/pdmpx_summary.py`)

Heuristics, for sorting - not verdicts:

| flag | rule |
|---|---|
| `RFI?` | not dispersed: S/N at DM 0 ≥ 0.8 × best, or `disp_sigma` < 1.5 with best DM < 10 |
| `RMS` | `rms_ratio` > 1.3 |
| `NB` | `f_frac90` < 0.25 (narrowband) |
| `BURSTY` | `t_frac90` < 0.25 |
| `WIDE` | duty cycle > 30% |
| `SINE` | `w_maxed` = 1 |
| `DMCOMB` | `dm_peak2` > 0.7 |
| `EDGE` | best trial at the grid edge |
| `noise` | S/N < 7 and not dispersed above 3σ |
| `LOOK` | S/N ≥ 7, dispersed ≥ 3σ, `rms_ratio` < 1.3, none of NB/RFI?/SINE/DMCOMB |

It also groups candidates whose periods are harmonically related (n/1, 1/n, or ratios of numbers ≤ 4):
large families at 50/100 Hz, 10/20/30 Hz etc. are RFI.

## Reading the diagnostic page (`tools/pdmpx_plot.py`)

All phase-resolved panels are aligned at the best period and DM, so a pulsar is **vertical**; 1.5 turns are
shown (phase −0.25 to 1.25, pulse centred at 0.5; the repeated quarter-turns are greyed).

- **Profile**: black = best, grey = the same data at DM 0. Blue shading = the baseline window; blue band = ±1
  predicted rms; red dashed = ±1 measured rms; orange = the bins summed for the S/N.
- **Phase vs time**, **phase vs frequency** (at the best DM and at DM 0, with the predicted dispersion sweep):
  rows scrunched to S/N ≈ 3 each. A pulsar is a vertical stripe at the best DM and follows the red sweep at DM 0.
- **Per-subint S/N** with the cumulative S/N and its √t expectation, and the beam-transit fit if any.
- **Single pulses vs noise**: per-subint on-pulse S/N against off-pulse windows of the same width.
- **Spectrum**: on-pulse S/N per channel.
- **Bullseye**: S/N over the whole period-DM grid. **S/N vs DM** (from DM 0, with the curve a real pulse would
  follow), **S/N vs period**, **S/N vs boxcar width** (a pulse has a single peak; still rising at the widest
  boxcar = not a pulse).

## Plot file format (`.pdx`)

gzip-compressed. A text header, then named float32 arrays:
```
PDMPX-BUNDLE 1
key=value            (every CSV column, plus grid parameters, baseline window, boxcar position, ...)
...
ARRAYS
name ny nx           followed by ny*nx little-endian float32
...
```
Arrays: `profile`, `profile_dm0`, `sub` (nsub × nbin), `chan`, `chan_dm0` (nchan × nbin), `sub_wt`, `chan_wt`,
`chan_freq`, `epoch`, `dur`, `grid` (nDM × nP), `dmcurve` (DM, S/N best width, S/N fixed width, model),
`pcurve`, `wcurve`, `sp_on`, `sp_off`, `chan_snr`. `tools/pdmpx_plot.py` has a reader (`read_pdx`).

## Other tools

```
tools/fil_window.py IN.fil OUT.fil --mjd MJD --half SECONDS     # cut a window out of a filterbank
tools/fil_window.py IN.fil --info                               # header, duration, MJD range
tools/drift_psrs.py *.fil                                        # known pulsars drifting through each beam
tools/drift_psrs.py *.hdr --sizes ls-l.txt                       # same, from headers + `ls -l` only
```
