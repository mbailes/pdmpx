// Archive preparation, mirroring the PSRCHIVE operations pdmp applies before its search.
#include "archive.h"
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <fftw3.h>
#include "fftplan.h"

double dispersion_delay(double dm, double fref, double f)
{
    return (dm / 2.41e-4) * (1.0 / (f * f) - 1.0 / (fref * fref));   // PSRCHIVE Physical.C
}

// PSRCHIVE Profile::rotate_phase -> FTransform::shift: X_h *= exp(+2 pi i h phase) for 0 < h < nbin/2
// (DC and the Nyquist term untouched): positive phase moves features EARLIER.
bool g_binrot = false;

void rotate_phase(float* amps, int nbin, double phase)
{
    phase -= std::floor(phase);
    if (phase == 0) return;
    if (g_binrot) {   // PSRCHIVE ::shift(nbin, binshift, amps): bin j takes bin j+binshift (earlier)
        unsigned k = unsigned(phase * double(nbin) + 0.5) % nbin;
        if (!k) return;
        std::vector<float> t(amps, amps + nbin);
        for (int j = 0; j < nbin; j++) amps[j] = t[(j + k) % nbin];
        return;
    }
    int nh = nbin / 2 + 1;
    // the shared plans for this length (made once; executing them on new arrays is thread-safe)
    const FFTPlans& pl = fft_plans(nbin);
    thread_local std::vector<float> in;
    thread_local std::vector<std::complex<float>> X;
    if ((int)in.size() != nbin) { in.assign(nbin, 0.f); X.assign(nh, 0.f); }
    std::copy(amps, amps + nbin, in.begin());
    fftwf_execute_dft_r2c(pl.r2c, in.data(), reinterpret_cast<fftwf_complex*>(X.data()));
    for (int h = 1; h < nbin / 2; h++) {         // h < nbin/2 (integer), exactly as PSRCHIVE
        double ph = 2 * M_PI * h * phase;
        X[h] *= std::complex<float>((float)std::cos(ph), (float)std::sin(ph));
    }
    fftwf_execute_dft_c2r(pl.c2r, reinterpret_cast<fftwf_complex*>(X.data()), in.data());
    for (int b = 0; b < nbin; b++) amps[b] = in[b] / nbin;
}

void dedisperse(Archive& a)
{
    if (a.dedispersed || a.dm == 0) { a.dedispersed = true; return; }
    #pragma omp parallel for collapse(2) schedule(static)
    for (int s = 0; s < a.nsub; s++)
        for (int c = 0; c < a.nchan; c++)
            rotate_phase(a.prof(s, c), a.nbin, dispersion_delay(a.dm, a.centre_freq, a.f(s, c)) / a.pfold[s]);
    a.dedispersed = true;
}

// Weighted means, summed weights - as PSRCHIVE's scrunching.
void tscrunch_to_nsub(Archive& a, int new_nsub)
{
    if (new_nsub <= 0 || new_nsub >= a.nsub) return;
    int factor = a.nsub / new_nsub + (a.nsub % new_nsub ? 1 : 0);   // Archive::tscrunch_to_nsub
    int out = (a.nsub + factor - 1) / factor;
    Archive b = a;
    b.nsub = out; b.epoch.assign(out, 0); b.dur.assign(out, 0); b.pfold.assign(out, 0);
    b.wts.assign((size_t)out * a.nchan, 0); b.data.assign((size_t)out * a.nchan * a.nbin, 0);
    if (!a.sfreq.empty()) b.sfreq.assign((size_t)out * a.nchan, 0);
    for (int o = 0; o < out; o++) {
        int s0 = o * factor, s1 = std::min(a.nsub, s0 + factor);
        double t0 = a.epoch[s0] - 0.5 * a.dur[s0], t1 = a.epoch[s1 - 1] + 0.5 * a.dur[s1 - 1];
        b.epoch[o] = 0.5 * (t0 + t1);
        double ps = 0;
        for (int s = s0; s < s1; s++) { b.dur[o] += a.dur[s]; ps += a.pfold[s]; }
        b.pfold[o] = ps / (s1 - s0);
        for (int c = 0; c < a.nchan; c++) {
            double W = 0; float* d = b.prof(o, c);
            for (int s = s0; s < s1; s++) {
                float w = a.wt(s, c); if (w <= 0) continue; W += w;
                const float* x = a.prof(s, c);
                for (int i = 0; i < a.nbin; i++) d[i] += w * x[i];
            }
            if (W > 0) for (int i = 0; i < a.nbin; i++) d[i] /= (float)W;
            b.wt(o, c) = (float)W;
            if (!a.sfreq.empty()) {
                double Wf = 0, Ws = 0;
                for (int s = s0; s < s1; s++) { double w = std::max(a.wt(s, c), 0.f); Wf += w * a.f(s, c); Ws += w; }
                b.sfreq[(size_t)o * a.nchan + c] = Ws > 0 ? Wf / Ws : a.f(s0, c);
            }
        }
    }
    a = std::move(b);
}

void fscrunch_to_nchan(Archive& a, int new_nchan)
{
    if (new_nchan <= 0 || new_nchan >= a.nchan) return;
    if (a.nchan % new_nchan)
        throw std::runtime_error("fscrunch_to_nchan: requested nchan does not divide nchan (as PSRCHIVE)");
    int f = a.nchan / new_nchan;
    Archive b = a;
    b.nchan = new_nchan; b.freq.assign(new_nchan, 0);
    b.wts.assign((size_t)a.nsub * new_nchan, 0); b.data.assign((size_t)a.nsub * new_nchan * a.nbin, 0);
    for (int o = 0; o < new_nchan; o++) {
        double fs = 0; for (int c = o * f; c < (o + 1) * f; c++) fs += a.freq[c];
        b.freq[o] = fs / f;
    }
    b.sfreq.assign((size_t)a.nsub * new_nchan, 0);
    for (int s = 0; s < a.nsub; s++)
        for (int o = 0; o < new_nchan; o++) {
            double W = 0, Wf = 0, fs = 0; float* d = b.prof(s, o);
            for (int c = o * f; c < (o + 1) * f; c++) {
                fs += a.f(s, c);
                float w = a.wt(s, c); if (w <= 0) continue; W += w; Wf += w * a.f(s, c);
                const float* x = a.prof(s, c);
                for (int i = 0; i < a.nbin; i++) d[i] += w * x[i];
            }
            if (W > 0) for (int i = 0; i < a.nbin; i++) d[i] /= (float)W;
            b.wt(s, o) = (float)W;
            // PSRCHIVE: the combined channel's frequency is the weighted mean of its members
            b.sfreq[(size_t)s * new_nchan + o] = W > 0 ? Wf / W : fs / f;
        }
    a = std::move(b);
}

void bscrunch_to_nbin(Archive& a, int new_nbin)
{
    if (new_nbin <= 0 || new_nbin >= a.nbin) return;
    int f = a.nbin / new_nbin;
    if (f < 2) return;
    int nb = a.nbin / f;
    std::vector<float> d((size_t)a.nsub * a.nchan * nb, 0);
    for (size_t p = 0; p < (size_t)a.nsub * a.nchan; p++)
        for (int i = 0; i < nb; i++) {
            float s = 0; for (int j = 0; j < f; j++) s += a.data[p * a.nbin + i * f + j];
            d[p * nb + i] = s / f;
        }
    a.data.swap(d); a.nbin = nb;
}
