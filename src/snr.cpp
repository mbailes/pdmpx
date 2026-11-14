// pdmp's S/N: PSRCHIVE FortranSNR -> Util/genutil/smooth.f smooth_mmw, ported line for line,
// and pdmp's getRMS (the noise of the weighted-mean profile from each profile's quietest half).
#include "snr.h"
#include <cmath>
#include <vector>

float smooth_mmw(float* pr, int nbin, int minwidth, int maxwidth, float rmsp, int* kwmax_out)
{
    return smooth_mmw_ex(pr, nbin, minwidth, maxwidth, rmsp, kwmax_out, nullptr);
}

float smooth_mmw_ex(float* pr, int nbin, int minwidth, int maxwidth, float rmsp, int* kwmax_out, MMWInfo* info)
{
    int kwmax = 0;
    if (minwidth < 1) minwidth = 1;
    if (minwidth > maxwidth) minwidth = maxwidth;          // (Fortran stops; never hit by pdmp's clamps)
    float snrmax = 0.f;
    int ksm = (int)(nbin / 2.5f + 0.5f);
    float smax = 1.e30f;
    int bstart = 0, kstart = 0;
    {                                                      // baseline: the lowest ksm-bin running mean
        double s = 0;                                      // O(nbin) running sum (Fortran re-sums each window)
        for (int k = 0; k < ksm; k++) s += pr[k % nbin];
        for (int j = 0; j < nbin; j++) {
            if ((float)s < smax) { smax = (float)s; bstart = j; }
            s += pr[(j + ksm) % nbin] - pr[j];
        }
    }
    smax /= ksm;
    for (int j = 0; j < nbin; j++) pr[j] -= smax;

    int nvalid = 0, kw = minwidth;
    for (int nn = 0; nn < 1000; nn++) {
        if (kw >= minwidth && kw <= maxwidth) nvalid++;
        kw = (int)(kw * 1.5f); if (kw == 1) kw = 2;
        if (kw > 1000000) break;
    }
    bool bail;
    if (nvalid == 0) { minwidth = maxwidth = nbin / 2; bail = true; } else bail = false;

    std::vector<float> wrk(nbin + 1);
    kw = minwidth;
    for (int nn = 0; nn < 1000; nn++) {
        if (nvalid != 0) {
            if (kw < minwidth) goto next;                  // (unreachable: kw starts at minwidth)
            if (kw > maxwidth) break;
        } else { kw = nbin / 2; bail = true; }
        if (kw > nbin * 0.9) break;                        // Fortran 'return'
        {
            float s = 0.f;
            for (int k = 0; k < kw; k++) { s += pr[k]; wrk[k] = pr[k]; }
            int ja = 0;                                    // 1-based in Fortran: ja runs 1..kw
            float smx = s; int st = 0;
            for (int j = 2; j <= nbin; j++) {
                ja++; if (ja > kw) ja -= kw;
                float al = wrk[ja - 1];
                int jj = j + kw - 1; if (jj > nbin) jj -= nbin;
                float an = pr[jj - 1];
                s = s + an - al;
                wrk[ja - 1] = an;
                if (s > smx) { smx = s; st = j - 1; }   // the window now starts at bin j-1 (0-based)
            }
            float sn = smx / (rmsp * std::sqrt(kw * (1.f + (float)kw / nbin)));
            if (sn > snrmax) { snrmax = sn; kwmax = kw; kstart = st; }
        }
        if (bail) break;
    next:
        kw = (int)(kw * 1.5f); if (kw == 1) kw = 2;
    }
    if (kwmax_out) *kwmax_out = kwmax;
    if (info) { info->base = smax; info->bstart = bstart; info->ksm = ksm; info->start = kstart; info->width = kwmax; }
    return snrmax;
}

double pdmp_rms(const Archive& a)
{
    int nbin = a.nbin, half = nbin / 2;
    double wtv = 0, tw = 0;
    for (int s = 0; s < a.nsub; s++)
        for (int c = 0; c < a.nchan; c++) {
            const float* x = a.prof(s, c);
            double smin = -1; int itmin = 0;
            for (int ib = 0; ib < nbin; ib++) {
                double sum = 0;
                for (int j = ib; j < ib + half; j++) sum += x[j % nbin];
                if (sum < smin || smin < 0) { smin = sum; itmin = ib; }
            }
            double mean = smin / half, ss = 0;
            for (int i = itmin; i < itmin + half; i++) { double d = x[i % nbin] - mean; ss += d * d; }
            double var = ss / (half - 1);
            double w = a.wt(s, c);
            wtv += w * w * var; tw += w;
        }
    return tw > 0 ? std::sqrt(wtv / (tw * tw)) : 1.0;
}
