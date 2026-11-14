// pdmpx - fast period/DM optimiser for folded pulsar archives (standalone: no PSRCHIVE).
// The in-memory archive: total intensity only, data[isub][ichan][ibin].
#pragma once
#include <string>
#include <vector>

struct Archive {
    std::string filename, source, format;
    int nsub = 0, nchan = 0, nbin = 0;
    double dm = 0;              // the DM the data are (or will be) dedispersed at
    bool dedispersed = false;   // channels already aligned at dm
    double centre_freq = 0;     // MHz: the dedispersion reference frequency (as PSRCHIVE)
    double bandwidth = 0;       // MHz (signed)
    std::vector<double> freq;   // per channel, MHz (first subint)
    std::vector<double> sfreq;  // per subint and channel [isub*nchan+ichan], MHz: after fscrunch the weighted mean, as PSRCHIVE
    double f(int s, int c) const { return sfreq.empty() ? freq[c] : sfreq[(size_t)s * nchan + c]; }
    std::vector<double> epoch;  // per subint: mid-time, seconds since the first subint's start
    std::vector<double> dur;    // per subint duration, s
    std::vector<double> pfold;  // per subint folding period, s
    double mjd0 = 0;            // MJD of the archive start (for reporting)
    std::vector<float> wts;     // [isub*nchan + ichan]
    std::vector<float> data;    // [(isub*nchan + ichan)*nbin + ibin]

    float* prof(int s, int c) { return &data[((size_t)s * nchan + c) * nbin]; }
    const float* prof(int s, int c) const { return &data[((size_t)s * nchan + c) * nbin]; }
    float& wt(int s, int c) { return wts[(size_t)s * nchan + c]; }
    float wt(int s, int c) const { return wts[(size_t)s * nchan + c]; }
    double start_time() const { return epoch.front() - 0.5 * dur.front(); }
    double end_time() const { return epoch.back() + 0.5 * dur.back(); }
};

Archive load_archive(const std::string& filename, bool verbose);
extern double g_force_period;   // -P: folding period (s) to assume when the file has none
Archive load_timer(const std::string& filename, bool verbose);
Archive load_psrfits(const std::string& filename, bool verbose);

// preparation (prep.cpp) - each mirrors the PSRCHIVE operation pdmp performs
void rotate_phase(float* amps, int nbin, double phase);           // PSRCHIVE convention: +phase = earlier
extern bool g_binrot;   // round every rotation to whole bins, as pdmp (Profile::rotate_in_phase_domain = true)
void dedisperse(Archive& a);                                       // at a.dm about a.centre_freq
void tscrunch_to_nsub(Archive& a, int nsub);
void fscrunch_to_nchan(Archive& a, int nchan);
void bscrunch_to_nbin(Archive& a, int nbin);
double dispersion_delay(double dm, double fref, double f);         // s, PSRCHIVE's 1/2.41e-4
