# pdmpx vs pdmp benchmarks

Apple M-series, 16 cores. pdmp = PSRCHIVE pdmp (patched, middle-subint period range), `-g /null`.
Archives: `cand.ar` = a noise-level FAST drift-scan candidate fold, 25 subints x 1024 ch x 128 bins, P ~ 1 s, DM ~ 360;
`psr_M05.ar` = PSR J0405+3347 (CRAFTS beam M05), 22 x 1024 x 128, P = 64 ms, DM 53.6.
Times are the search plus preparation (load is ~0.005 s). pdmpx `-binrot` reproduces pdmp's best trial exactly.

## 1. pdmp's own "natural" grid (`pdmpx -natural`)

| archive | -mc | trials (DM x P) | pdmp | pdmpx 1 thread | pdmpx 16 threads | speed-up |
|---|---|---|---|---|---|---|
| cand.ar | 32 | 130 x 51 | ~2 s | 0.07 s | 0.008 s | ~250x |
| cand.ar | 128 | 514 x 51 | 42.2 s | 0.47 s | 0.043 s | ~980x |
| cand.ar | 1024 | 4098 x 51 | **9246 s (2 h 34 m)** | 18.4 s | 1.6 s | ~5800x |
| psr_M05.ar | 128 | 514 x 46 | 32.2 s | 0.43 s | 0.038 s | ~850x |
| psr_M05.ar | 1024 | 4098 x 46 | (running, >1.5 h) | 16.5 s | 1.1 s | |

The natural DM half-range is the DM giving one bin of delay across a *single channel*, i.e. ~nchan bins
(many turns) across the band, and it grows with nchan. pdmp picked DM = -24 (at -mc 128) and DM = -2047
(at full resolution) for that candidate (input DM ~ 360).

## 2. A +-0.5 turn grid (the first sensible default; now the default is +-1 turn, see 3)

Scrunch to <= 32 channels; DM half-range = 0.5 turn of delay across the band, step 1/2 bin across the band
(2 x nbin trials, independent of nchan); period half-range = min(pdmp's, 0.5 turn over the observation),
step 1 bin over it. Negative DMs are dropped. pdmp was run with the same grid via the `-dr/-ds/-pr/-ps`
values pdmpx prints.

| archive | trials | pdmp | pdmpx 1 thread | pdmpx 16 threads | best (both, -binrot) |
|---|---|---|---|---|---|
| cand.ar | 256 x 51 | 8.6 s | 0.15 s | 0.024 s | S/N 5.14 (noise level; same best trial in both) |
| psr_M05.ar | 256 x 46 | 6.4 s | 0.14 s | 0.022 s | S/N 21.68, DM 53.38 (fractional shifts: 22.17, DM 53.61) |

Recovery of J0405+3347 with the search centred away from the true DM: `-do 10` -> 53.60;
`-do 20` is outside +-0.5 turn (+-14.4 DM here) and fails; `-do 20 -turns 1` -> 53.59; `-do 40 -turns 2` -> 53.57 (0.05 s).

## 3. Default grid: +-1 turn in DM and period, <= 32 channels

512 DM x 257 P = 131,584 trials on psr_M05.ar: pdmp 561 s; pdmpx 0.063 s search (16 threads), 0.75 s (1 thread);
load 0.005 s, prepare (dedisperse 1024 ch, scrunch to 32) 0.006 s. pdmp's FortranSNR baseline re-summed a
51-bin window at every bin (O(nbin^2), 58% of the search); it is now a running sum (identical S/N on all trials).
Recovers DM 53.6 with the search centred 20 DM units off.

## 4. Reading a big archive (MPTA-sized): slurp vs per-profile I/O

Synthetic archive made with PSRCHIVE: 450 subints x 4 pol (AABBCRCI) x 1024 ch x 1024 bins, MeerKAT L-band,
a 5.7575 ms pulsar at DM 25, every 97th channel zapped. Timer 3.84 GB; PSRFITS 3.80 GB (written by `pam -a PSRFITS`).
Page cache warm (48 GB machine; cold reads need `purge`/sudo, not run). PSRCHIVE time = `Archive::load` plus touching
every profile (PSRFITS loads each subint lazily).

| reader | Timer | PSRFITS |
|---|---|---|
| PSRCHIVE `Archive::load` | 5.38 s | 3.28 s |
| old pdmpx (stdio freads per profile / cfitsio per row & column) | 0.98 s | 1.36 s |
| **new pdmpx: mmap, subints decoded in parallel (16 threads)** | **0.21 s** | **0.20 s** |
| new pdmpx: mmap, 1 thread | 0.31 s | 0.29 s |
| new pdmpx: one read() into a buffer, then decode | 0.54 s | 0.51 s |

Whole job with `-mc 32 -mb 128` (16 threads): load 0.21 s, prepare 0.43 s (dedisperse 1024 ch x 1024 bins, scrunch),
search 0.79 s (512 DM x 257 P): 1.4 s total; finds S/N 129.7 at DM 25.0000, P = 5.7575 ms.
Without `-mb`, the default +-1 turn grid at 1024 bins is 4096 DM x 2049 P = 8.4 M trials (374 s): too fine; see README.

Also found while building this: the Timer reader summed AA+BB+CR for 4-pol Coherence data (now uses the header's
`banda.correlator_mode`, as PSRCHIVE; matches pdv AA+BB to 4e-6; old reader gave S/N 106 instead of 130 here), and
the byte-order test misread big-endian nbin = 1024 as 262144.
