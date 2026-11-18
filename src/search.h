#pragma once
#include "archive.h"
#include <vector>

struct Grid {
    double dm0 = 0, ddm = 0; int ndm = 1;     // trial DM_i = dm0 + i*ddm (absolute)
    double p0 = 0, dp = 0; int np = 1;        // trial topocentric period P_j = p0 + j*dp (seconds)
};

struct Result {
    int idm = 0, ip = 0, width = 0;
    float snr = 0;
    double dm = 0, period = 0;
    std::vector<float> snr_grid;              // [idm*np + ip]
    std::vector<float> snr_grid_w;            // same, at the fixed width SearchConfig::fixw (if > 0)
    std::vector<float> best_profile;
};

struct SearchConfig {
    int minwidth = 1, maxwidth = 0;           // bins (maxwidth 0 -> nbin/2)
    double rms = 1;
    bool keep_grid = false;
    int fixw = 0;                             // > 0 (with keep_grid): also the S/N at this one width, same pass
};

// Brute-force grid in the Fourier domain: each profile is transformed once; a trial is a phase
// rotation of the per-channel spectra (DM), then of the per-subint spectra (period), one inverse
// FFT and pdmp's S/N. Same trials and same S/N definition as pdmp; parallel over DM.
Result search_grid(const Archive& a, const Grid& g, const SearchConfig& cfg);
