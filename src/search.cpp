#include "search.h"
#include "snr.h"
#include <cmath>
#include <complex>
#include <vector>
#include <fftw3.h>
#include "fftplan.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using cf = std::complex<float>;
using cd = std::complex<double>;

Result search_grid(const Archive& a, const Grid& g, const SearchConfig& cfg)
{
    const int nsub = a.nsub, nchan = a.nchan, nbin = a.nbin, nh = nbin / 2 + 1;
    const int hrot = nbin / 2;                 // harmonics 1..hrot-1 are rotated; DC and Nyquist are not (PSRCHIVE)
    const int maxw = cfg.maxwidth > 0 ? cfg.maxwidth : nbin / 2;

    // 1. every (subint, channel) profile -> weighted spectrum, once
    std::vector<cf> X((size_t)nsub * nchan * nh);
    double W = 0;
    {
        std::vector<float> in(nbin);
        const fftwf_plan p = fft_plans(nbin).r2c;
        for (int s = 0; s < nsub; s++)
            for (int c = 0; c < nchan; c++) {
                float w = a.wt(s, c);
                cf* out = &X[((size_t)s * nchan + c) * nh];
                if (w <= 0) { std::fill(out, out + nh, cf(0)); continue; }
                const float* x = a.prof(s, c);
                for (int b = 0; b < nbin; b++) in[b] = w * x[b];
                fftwf_execute_dft_r2c(p, in.data(), reinterpret_cast<fftwf_complex*>(out));   // thread-safe
                W += w;
            }
    }
    const float invWN = (float)(1.0 / (W * nbin));   // weighted mean, and FFTW's unnormalised inverse

    // the subint time offsets from the mid-observation reference (pdmp: start + (end-start)/2)
    const double tref = 0.5 * (a.start_time() + a.end_time());
    const double Pmid = a.pfold[nsub / 2];
    std::vector<double> dt(nsub);
    for (int s = 0; s < nsub; s++) dt[s] = a.epoch[s] - tref;

    Result R;
    if (cfg.keep_grid) R.snr_grid.assign((size_t)g.ndm * g.np, 0.f);
    if (cfg.keep_grid && cfg.fixw > 0) R.snr_grid_w.assign((size_t)g.ndm * g.np, 0.f);
    R.snr = -1;

    #pragma omp parallel
    {
        std::vector<cf> Z((size_t)nsub * nh), Y(nh);
        std::vector<float> prof(nbin), best_local;
        std::vector<cd> ph(nh);
        const fftwf_plan inv = fft_plans(nbin).c2r;
        float bsnr = -1; int bi = 0, bj = 0, bw = 0;

        // rotation of a spectrum by 'phase' turns, PSRCHIVE sense (+ = earlier): X_h *= exp(+2 pi i h phase)
        auto make_phasors = [&](double phase) {
            double f = phase - std::floor(phase);
            if (g_binrot) {
                unsigned k = unsigned(f * double(nbin) + 0.5) % nbin;
                double r = 2 * M_PI * double(k) / nbin;
                cd z1(std::cos(r), std::sin(r)), z(1, 0);
                for (int h = 0; h < nh; h++) { ph[h] = z; z *= z1; }
                return;
            }
            double r = 2 * M_PI * f;
            cd z1(std::cos(r), std::sin(r)), z(1, 0);
            for (int h = 0; h < nh; h++) { ph[h] = (h >= 1 && h < hrot) ? z : cd(1, 0); z *= z1; }
        };

        #pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < g.ndm; i++) {
            const double ddm = g.dm0 + i * g.ddm - a.dm;     // relative to the DM the data are aligned at
            // 2. dedisperse at this trial and sum channels: Z[s][h] = sum_c X[s][c][h] * rot(delay/P_s)
            std::fill(Z.begin(), Z.end(), cf(0));
            for (int s = 0; s < nsub; s++) {
                cf* z = &Z[(size_t)s * nh];
                for (int c = 0; c < nchan; c++) {
                    if (a.wt(s, c) <= 0) continue;
                    const cf* x = &X[((size_t)s * nchan + c) * nh];
                    double phase = dispersion_delay(ddm, a.centre_freq, a.f(s, c)) / a.pfold[s];
                    if (phase == 0) { for (int h = 0; h < nh; h++) z[h] += x[h]; continue; }
                    if (g_binrot) {   // whole-bin cyclic shift: every harmonic, Nyquist included
                        make_phasors(phase);
                        for (int h = 0; h < nh; h++) z[h] += x[h] * cf((float)ph[h].real(), (float)ph[h].imag());
                        continue;
                    }
                    double r = 2 * M_PI * (phase - std::floor(phase));
                    cd z1(std::cos(r), std::sin(r)), zz(1, 0);
                    z[0] += x[0];
                    zz = z1;
                    for (int h = 1; h < hrot; h++) { z[h] += x[h] * cf((float)zz.real(), (float)zz.imag()); zz *= z1; }
                    for (int h = std::max(hrot, 1); h < nh; h++) z[h] += x[h];
                }
            }
            // 3. each trial period: rotate each subint by -dnu*dt (pdmp counter_drift), sum, invert, S/N
            for (int j = 0; j < g.np; j++) {
                const double P = g.p0 + j * g.dp;
                const double dnu = 1.0 / P - 1.0 / Pmid;
                std::fill(Y.begin(), Y.end(), cf(0));
                for (int s = 0; s < nsub; s++) {
                    make_phasors(-dnu * dt[s]);
                    const cf* z = &Z[(size_t)s * nh];
                    for (int h = 0; h < nh; h++) Y[h] += z[h] * cf((float)ph[h].real(), (float)ph[h].imag());
                }
                fftwf_execute_dft_c2r(inv, reinterpret_cast<fftwf_complex*>(Y.data()), prof.data());
                for (int b = 0; b < nbin; b++) prof[b] *= invWN;
                std::vector<float> keep;
                int kw = 0;
                float sn;
                {
                    std::vector<float> tmp(prof);
                    sn = smooth_mmw(tmp.data(), nbin, cfg.minwidth, maxw, (float)cfg.rms, &kw);
                }
                if (cfg.keep_grid) R.snr_grid[(size_t)i * g.np + j] = sn;
                if (cfg.keep_grid && cfg.fixw > 0) {     // the fixed-width curve from the same dedispersed profile
                    std::vector<float> tmp(prof);
                    R.snr_grid_w[(size_t)i * g.np + j] = smooth_mmw(tmp.data(), nbin, cfg.fixw, cfg.fixw, (float)cfg.rms, nullptr);
                }
                if (sn > bsnr) { bsnr = sn; bi = i; bj = j; bw = kw; best_local = prof; }
            }
        }
        #pragma omp critical(best)
        {
            // ties resolved towards pdmp's scan order (first DM, then first period)
            if (bsnr > R.snr || (bsnr == R.snr && (bi < R.idm || (bi == R.idm && bj < R.ip)))) {
                R.snr = bsnr; R.idm = bi; R.ip = bj; R.width = bw; R.best_profile = best_local;
            }
        }
    }
    R.dm = g.dm0 + R.idm * g.ddm;
    R.period = g.p0 + R.ip * g.dp;
    return R;
}
