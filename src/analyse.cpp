#include "analyse.h"
#include "snr.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <random>
#include <zlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef PDMPX_VERSION
#define PDMPX_VERSION "0.4"
#endif

// ------------------------------------------------------------------------------------------------ CSV
const std::vector<std::string>& csv_columns()
{
    static const std::vector<std::string> c = {
        "file", "status", "source", "mjd", "tobs_s", "fc_MHz", "bw_MHz", "nsub", "nchan", "nbin",
        "nchan_search", "nchan_disp", "zap_frac",
        "P_fold_ms", "DM_fold", "P_ms", "P_err_ms", "DM", "DM_err", "dP_ms", "dDM", "edge",
        "snr", "width_bins", "width_ms", "duty", "snr_coho", "width_coho", "baseline", "rms_pred", "rms_meas", "rms_ratio", "chi2_red",
        "n_trials", "n_eff", "fap", "sigma",
        "snr_dm0", "snr_dm0_model", "disp_sigma", "dm_peak2", "w_maxed",
        "sp_max", "sp_max_sub", "sp_sig", "mod_index", "t_frac90", "f_frac90",
        "drift_t_s", "drift_t_err_s", "drift_fwhm_s", "drift_mjd_peak",
        "fft_snr", "fold_gain",
        "plotfile", "runtime_s", "version", "options"};
    return c;
}
std::string csv_header()
{
    std::string s;
    for (auto& c : csv_columns()) s += (s.empty() ? "" : ",") + c;
    return s;
}
std::string csv_line(const Row& r)
{
    std::string s;
    bool first = true;
    for (auto& c : csv_columns()) {
        std::string v = r.get(c);
        if (v.find_first_of(",\"\n") != std::string::npos) {   // quote, as CSV
            std::string q = "\"";
            for (char ch : v) { if (ch == '"') q += '"'; q += ch == '\n' ? ' ' : ch; }
            v = q + "\"";
        }
        s += (first ? "" : ",") + v;
        first = false;
    }
    return s;
}
void Row::set(const std::string& k, double x, const char* fmt)
{
    if (!std::isfinite(x)) { v[k] = ""; return; }
    char b[64]; snprintf(b, sizeof b, fmt, x); v[k] = b;
}

// ------------------------------------------------------------------------------------------------ grid
static const double MICROSEC = 1e6;

static Grid make_grid(const Archive& a, const Options& o, double& refP, double& refDM, std::string& kind)
{
    const double DMCONST = 4.15e6, MILLISEC = 1e3;
    double ds = o.ds, dr = o.dr, ps = o.ps, pr = o.pr;
    const bool user_dr = dr >= 0, user_ds = ds >= 0, user_pr = pr >= 0, user_ps = ps >= 0;
    double P0 = a.pfold[0];
    if (!o.natural) {   // sensible grid: +-turns of drift in DM across the band; period steps of one bin over the observation
        double fmin = 1e30, fmax = -1e30;
        for (int s = 0; s < a.nsub; s++)
            for (int c = 0; c < a.nchan; c++)
                if (a.wt(s, c) > 0) { fmin = std::min(fmin, a.f(s, c)); fmax = std::max(fmax, a.f(s, c)); }
        double D1 = fmax > fmin ? std::fabs(dispersion_delay(1.0, a.centre_freq, fmin) - dispersion_delay(1.0, a.centre_freq, fmax)) : 0;
        double Pm = a.pfold[a.nsub / 2], tspan = a.end_time() - a.start_time();
        if (D1 > 0) {
            if (!user_dr) dr = o.turns * Pm / D1;
            if (!user_ds) ds = 0.5 * (Pm / a.nbin) / D1;
        }
        if (!user_ps) ps = MICROSEC * Pm * (Pm / a.nbin) / tspan;
        // period half-range: pdmp's own (below: one bin of drift within the middle subint); -turns sets the DM range only
    }
    {   // pdmp's natural values for anything still unset (getNatural..., setSensibleStepSizes)
        double cp_ms = P0 * MILLISEC, tbin = cp_ms / a.nbin, bwband = std::fabs(a.bandwidth / a.nchan), fc = a.centre_freq;
        double nat_dmhr = tbin / DMCONST / (std::pow(1.0 / (fc - 0.5 * bwband), 2) - std::pow(1.0 / (fc + 0.5 * bwband), 2));
        if (ds < 0) ds = (dr < 0 ? nat_dmhr : dr) / (2.0 * a.nchan + 1);
        if (dr < 0) dr = nat_dmhr;
        double tspan = a.end_time() - a.start_time();
        if (ps < 0) ps = MICROSEC * P0 * (P0 / a.nbin) / tspan;
        if (pr < 0) {   // patched pdmp (ece2d87): the middle subint's length
            double cpu = P0 * MICROSEC, tsub = a.dur[a.nsub / 2] * MICROSEC;
            pr = cpu * (cpu / a.nbin) / tsub;
        }
        if (a.nchan == 1 || ds <= 0) dr = ds / 2;
        if (a.nsub == 1 || ps <= 0) pr = ps / 2;
    }
    int dmBins = (int)std::ceil(std::fabs(dr) * 2 / ds);
    int periodBins = (int)std::ceil(std::fabs(pr) * 2 / ps);
    double refP_us = a.pfold[a.nsub / 2] * MICROSEC;
    refP = refP_us * 1e-6; refDM = a.dm;
    Grid g;
    g.ndm = dmBins; g.ddm = ds; g.dm0 = refDM + o.doff - std::fabs(ds * std::floor(dmBins / 2.0));
    g.np = periodBins + 1; g.dp = ps * 1e-6;
    g.p0 = (refP_us + o.poff - std::fabs(ps * std::floor(periodBins / 2.0))) * 1e-6;
    if (!o.natural && g.dm0 < 0) {   // physical DMs only
        int k = (int)std::ceil(-g.dm0 / g.ddm - 1e-9);
        g.dm0 += k * g.ddm; g.ndm -= k;
    }
    if (o.dmup >= 0) {               // never more than dmup above the input DM
        double top = refDM + o.dmup;
        if (g.dm0 + (g.ndm - 1) * g.ddm > top) g.ndm = std::max(1, (int)std::floor((top - g.dm0) / g.ddm + 1e-9) + 1);
    }
    if (g.ndm < 1) g.ndm = 1;
    kind = o.natural ? "pdmp natural" : "sensible";
    return g;
}

// ------------------------------------------------------------------------------------------------ helpers
// every profile rotated to the trial (DM, P) exactly as the search does; weighted-mean sums per subint,
// per channel and overall (each normalised by its weight)
struct Aligned {
    std::vector<float> sub, chan, prof;       // [nsub*nbin], [nchan*nbin], [nbin]
    std::vector<double> subw, chanw;
};
static Aligned align(const Archive& a, double dm, double P)
{
    const int nsub = a.nsub, nchan = a.nchan, nbin = a.nbin;
    const double tref = 0.5 * (a.start_time() + a.end_time()), Pmid = a.pfold[nsub / 2];
    const double dnu = 1.0 / P - 1.0 / Pmid;
    Aligned r;
    r.sub.assign((size_t)nsub * nbin, 0.f); r.chan.assign((size_t)nchan * nbin, 0.f); r.prof.assign(nbin, 0.f);
    r.subw.assign(nsub, 0); r.chanw.assign(nchan, 0);
    std::vector<float> rot((size_t)nsub * nchan * nbin);
    #pragma omp parallel for collapse(2) schedule(static)
    for (int s = 0; s < nsub; s++)
        for (int c = 0; c < nchan; c++) {
            float* d = &rot[((size_t)s * nchan + c) * nbin];
            if (a.wt(s, c) <= 0) { std::fill(d, d + nbin, 0.f); continue; }
            std::copy(a.prof(s, c), a.prof(s, c) + nbin, d);
            rotate_phase(d, nbin, dispersion_delay(dm - a.dm, a.centre_freq, a.f(s, c)) / a.pfold[s] - dnu * (a.epoch[s] - tref));
        }
    double W = 0;
    for (int s = 0; s < nsub; s++)
        for (int c = 0; c < nchan; c++) {
            float w = a.wt(s, c);
            if (w <= 0) continue;
            const float* d = &rot[((size_t)s * nchan + c) * nbin];
            for (int b = 0; b < nbin; b++) { r.sub[(size_t)s * nbin + b] += w * d[b]; r.chan[(size_t)c * nbin + b] += w * d[b]; r.prof[b] += w * d[b]; }
            r.subw[s] += w; r.chanw[c] += w; W += w;
        }
    for (int s = 0; s < nsub; s++) if (r.subw[s] > 0) for (int b = 0; b < nbin; b++) r.sub[(size_t)s * nbin + b] /= r.subw[s];
    for (int c = 0; c < nchan; c++) if (r.chanw[c] > 0) for (int b = 0; b < nbin; b++) r.chan[(size_t)c * nbin + b] /= r.chanw[c];
    if (W > 0) for (int b = 0; b < nbin; b++) r.prof[b] /= W;
    return r;
}

// S/N of a W-bin window starting at 'start', against the mean and rms of the off-pulse bins
static double window_snr(const float* p, int nbin, int start, int W, const std::vector<char>& offmask, int skip_start = -1)
{
    double s = 0, ss = 0; int n = 0;
    for (int b = 0; b < nbin; b++) {
        if (!offmask[b]) continue;
        if (skip_start >= 0 && ((b - skip_start + nbin) % nbin) < W) continue;
        s += p[b]; ss += (double)p[b] * p[b]; n++;
    }
    if (n < 3) return NAN;
    double mean = s / n, var = (ss - n * mean * mean) / (n - 1);
    if (!(var > 0)) return NAN;
    double on = 0;
    for (int k = 0; k < W; k++) on += p[(start + k) % nbin] - mean;
    return on / std::sqrt(var * W);
}

// one-sided Gaussian-equivalent significance of a tail probability p (Acklam's inverse-normal approximation)
static double p_to_sigma(double p)
{
    if (!(p > 0)) return INFINITY;
    if (p >= 0.5) return 0.0;
    static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00};
    double q, x;
    if (p < 0.02425) {
        q = std::sqrt(-2 * std::log(p));
        x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    } else {
        q = p - 0.5; double r = q * q;
        x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
    }
    return -x;   // x = Phi^-1(p) < 0 for p < 0.5
}

// fraction of the entries holding 90% of the (positive) signal
static double frac90(std::vector<double> v)
{
    for (auto& x : v) x = std::isfinite(x) ? std::max(x, 0.0) : 0.0;
    double tot = std::accumulate(v.begin(), v.end(), 0.0);
    if (!(tot > 0) || v.empty()) return NAN;
    std::sort(v.rbegin(), v.rend());
    double acc = 0;
    for (size_t i = 0; i < v.size(); i++) { acc += v[i]; if (acc >= 0.9 * tot) return double(i + 1) / v.size(); }
    return 1.0;
}

// E[max over phase] of the W-bin boxcar S/N for pure noise (the noise floor of a max-over-phase statistic)
static double noise_floor(int nbin, int W)
{
    std::mt19937 rng(12345);
    std::normal_distribution<double> nd;
    const int N = 300;
    double acc = 0;
    std::vector<double> x(nbin);
    for (int it = 0; it < N; it++) {
        for (auto& v : x) v = nd(rng);
        double s = 0, best = -1e30;
        for (int k = 0; k < W; k++) s += x[k % nbin];
        for (int j = 0; j < nbin; j++) { best = std::max(best, s); s += x[(j + W) % nbin] - x[j]; }
        acc += best / std::sqrt(W * (1.0 + double(W) / nbin));
    }
    return acc / N;
}

// fraction of a top-hat pulse's boxcar S/N that survives a DM error ddm, with the measured per-channel signal
static double dm_loss(double ddm, const std::vector<double>& freq, const std::vector<double>& wsig, double fc, double P, int nbin, int W)
{
    std::vector<double> acc(nbin, 0.0);
    double wsum = 0;
    for (size_t c = 0; c < freq.size(); c++) {
        if (!(wsig[c] > 0)) continue;
        double sh = dispersion_delay(ddm, fc, freq[c]) / P * nbin;   // bins
        sh -= std::floor(sh / nbin) * nbin;
        int i0 = (int)std::floor(sh); double fr = sh - i0;
        for (int k = 0; k < W; k++) {                               // top-hat at bins 0..W-1, shifted by sh
            acc[(k + i0) % nbin] += wsig[c] * (1 - fr);
            acc[(k + i0 + 1) % nbin] += wsig[c] * fr;
        }
        wsum += wsig[c];
    }
    if (!(wsum > 0)) return NAN;
    // subtract the baseline as smooth_mmw does (lowest ksm-bin running mean): a pulse smeared over many turns
    // becomes a flat offset and is removed, exactly as in the data
    const int ksm = (int)(nbin / 2.5 + 0.5);
    double bs = 0, bmin = 1e300;
    for (int k = 0; k < ksm; k++) bs += acc[k % nbin];
    for (int j = 0; j < nbin; j++) { bmin = std::min(bmin, bs); bs += acc[(j + ksm) % nbin] - acc[j]; }
    for (auto& x : acc) x -= bmin / ksm;
    double s = 0, best = 0;
    for (int k = 0; k < W; k++) s += acc[k % nbin];
    for (int j = 0; j < nbin; j++) { best = std::max(best, s); s += acc[(j + W) % nbin] - acc[j]; }
    return best / (wsum * W);
}

// Gaussian beam transit + constant fitted to y(t): grid in (t0, sigma), linear least squares in (A, c)
struct Drift { bool ok = false; double t0 = NAN, t0err = NAN, fwhm = NAN, amp = NAN; };
static Drift fit_drift(const std::vector<double>& t, const std::vector<double>& y, double span)
{
    Drift d;
    const int n = (int)t.size();
    if (n < 6) return d;
    double ym = 0; for (double v : y) ym += v; ym /= n;
    double chi_c = 0; for (double v : y) chi_c += (v - ym) * (v - ym);
    auto fit = [&](double t0, double sg, double& A, double& c) {
        double Sg = 0, Sgg = 0, Sy = 0, Sgy = 0;
        for (int i = 0; i < n; i++) { double gi = std::exp(-0.5 * std::pow((t[i] - t0) / sg, 2)); Sg += gi; Sgg += gi * gi; Sy += y[i]; Sgy += gi * y[i]; }
        double det = n * Sgg - Sg * Sg;
        if (std::fabs(det) < 1e-12) { A = 0; c = Sy / n; } else { A = (n * Sgy - Sg * Sy) / det; c = (Sy - A * Sg) / n; }
        double chi = 0;
        for (int i = 0; i < n; i++) { double r = y[i] - (A * std::exp(-0.5 * std::pow((t[i] - t0) / sg, 2)) + c); chi += r * r; }
        return chi;
    };
    double best = 1e300, bt = 0, bs = 0, bA = 0;
    double dt = (t.back() - t.front()) / (n - 1);
    for (int is = 0; is < 40; is++) {
        double sg = dt * std::pow(span / dt, is / 39.0) / 2.0;
        for (int it = 0; it <= 4 * n; it++) {
            double t0 = t.front() + (t.back() - t.front()) * it / (4.0 * n);
            double A, c, chi = fit(t0, sg, A, c);
            if (A > 0 && chi < best) { best = chi; bt = t0; bs = sg; bA = A; }
        }
    }
    if (best >= 1e300) return d;
    double s2 = best / std::max(1, n - 4);                       // residual variance
    // t0 error from the curvature of chi2 in t0 (Delta chi2 = s2)
    double h = std::max(dt * 0.05, 1e-6), A, c;
    double cp = fit(bt + h, bs, A, c), cm = fit(bt - h, bs, A, c);
    double curv = (cp + cm - 2 * best) / (h * h);
    d.t0 = bt; d.fwhm = 2.3548 * bs; d.amp = bA;
    d.t0err = curv > 0 ? std::sqrt(2 * s2 / curv) : NAN;
    double dchi = (chi_c - best) / std::max(s2, 1e-12);
    d.ok = bt > t.front() && bt < t.back() && d.fwhm < 0.8 * span && dchi > 25 && bA > 5 * std::sqrt(s2);
    return d;
}

// error on a parameter from the S/N curvature at the peak: sigma^2 = 1 / (-S * S'') (matched-filter likelihood)
static double peak_error(const float* v, int n, int i, double step)
{
    if (n < 3 || i <= 0 || i >= n - 1) return NAN;
    double S = v[i], d2 = (v[i + 1] + v[i - 1] - 2.0 * v[i]) / (step * step);
    return d2 < 0 ? 1.0 / std::sqrt(-S * d2) : NAN;
}

// ------------------------------------------------------------------------------------------------ bundle
struct Bundle {
    std::string meta;
    struct Arr { std::string name; int ny, nx; std::vector<float> d; };
    std::vector<Arr> arrays;
    void kv(const std::string& k, const std::string& v) { meta += k + "=" + v + "\n"; }
    void kv(const std::string& k, double x) { char b[64]; snprintf(b, sizeof b, "%.12g", x); kv(k, std::string(b)); }
    template <class V> void add(const std::string& name, int ny, int nx, const V& v)
    {
        Arr a{name, ny, nx, std::vector<float>(v.begin(), v.begin() + (size_t)ny * nx)};
        arrays.push_back(std::move(a));
    }
    bool write(const std::string& fn) const
    {
        gzFile f = gzopen(fn.c_str(), "wb6");
        if (!f) return false;
        std::string head = "PDMPX-BUNDLE 1\n" + meta + "ARRAYS\n";
        gzwrite(f, head.data(), (unsigned)head.size());
        for (auto& a : arrays) {
            char h[256]; int n = snprintf(h, sizeof h, "%s %d %d\n", a.name.c_str(), a.ny, a.nx);
            gzwrite(f, h, n);
            gzwrite(f, a.d.data(), (unsigned)(a.d.size() * sizeof(float)));   // little-endian float32
        }
        gzclose(f);
        return true;
    }
};

static std::string basename_noext(const std::string& p)
{
    size_t s = p.find_last_of('/');
    return s == std::string::npos ? p : p.substr(s + 1);
}

// coho's folded-boxcar S/N (coho.cpp foldBox), applied to a profile: baseline = median of the bins,
// rms = 1.4826 x MAD over ALL bins (no off-pulse mask), peak = max of a 3-bin running mean, boxcars of
// coho's widths {1,2,3,4,6,8,12,16,24,36,54} bins centred on that peak, S/N = sum / (rms sqrt(W)).
static double coho_snr(const std::vector<float>& prof, int* width)
{
    static const int WIDTHS[] = {1, 2, 3, 4, 6, 8, 12, 16, 24, 36, 54};
    const int NB = (int)prof.size();
    std::vector<double> srt(prof.begin(), prof.end());
    std::sort(srt.begin(), srt.end());
    const double bs = srt[NB >> 1];
    std::vector<double> dv(NB);
    for (int b = 0; b < NB; b++) dv[b] = std::fabs(prof[b] - bs);
    std::sort(dv.begin(), dv.end());
    double orms = 1.4826 * dv[NB >> 1]; if (orms <= 0) orms = 1e-12;
    int pk = 0; double pv = -INFINITY;
    for (int b = 0; b < NB; b++) {
        double v = (prof[(b - 1 + NB) % NB] + prof[b] + prof[(b + 1) % NB]) / 3;
        if (v > pv) { pv = v; pk = b; }
    }
    double sn = -INFINITY; int bw = 1; const int wmax = std::max(1, NB / 2);
    for (int W : WIDTHS) {
        if (W > wmax) break;
        int st = ((pk - (W >> 1)) % NB + NB) % NB;
        double acc = 0; for (int k = 0; k < W; k++) acc += prof[(st + k) % NB] - bs;
        double v = acc / (orms * std::sqrt((double)W));
        if (v > sn) { sn = v; bw = W; }
    }
    if (width) *width = bw;
    return sn;
}

// ------------------------------------------------------------------------------------------------ analysis
Row analyse_file(const std::string& file, const Options& o, Extra* ex)
{
    using clk = std::chrono::steady_clock;
    auto sec = [](auto d) { return std::chrono::duration<double>(d).count(); };
    auto t0 = clk::now();
    Row row;
    row.set("file", file);
    row.set("version", std::string(PDMPX_VERSION));
    row.set("options", o.optstring);
    try {
        Archive a = load_archive(file, o.verbose);
        auto t1 = clk::now();
        row.set("source", a.source);
        row.set("mjd", a.mjd0, "%.8f");
        row.set("tobs_s", a.end_time() - a.start_time(), "%.4f");
        row.set("fc_MHz", a.centre_freq, "%.4f");
        row.set("bw_MHz", a.bandwidth, "%.4f");
        row.set("nsub", a.nsub, "%.0f"); row.set("nchan", a.nchan, "%.0f"); row.set("nbin", a.nbin, "%.0f");
        row.set("P_fold_ms", a.pfold[a.nsub / 2] * 1e3, "%.10g");
        row.set("DM_fold", a.dm, "%.5f");
        {
            size_t z = 0; for (float w : a.wts) z += (w <= 0);
            row.set("zap_frac", double(z) / a.wts.size(), "%.4f");
        }

        // as pdmp: dedisperse at the archive DM (full resolution), then the partial scrunches
        dedisperse(a);
        if (o.ms > 1) tscrunch_to_nsub(a, o.ms);
        if (o.mb > 1) bscrunch_to_nbin(a, o.mb);
        Archive full = a;                       // kept for the dispersion analysis (dedispersed at the archive DM)
        int mc = o.mc;
        if (mc < 0 && !o.natural && a.nchan > 32) { mc = 32; while (a.nchan % mc) mc--; }
        if (mc > 1 && mc < a.nchan) fscrunch_to_nchan(a, mc);
        row.set("nchan_search", a.nchan, "%.0f");

        double refP, refDM; std::string kind;
        Grid g = make_grid(a, o, refP, refDM, kind);
        SearchConfig cfg;
        int minw = (int)(o.mw / refP * a.nbin);
        minw = std::max(1, std::min(minw, a.nbin / 2));
        cfg.minwidth = minw; cfg.maxwidth = a.nbin / 2;
        cfg.rms = pdmp_rms(a);
        cfg.keep_grid = true;
        auto t2 = clk::now();
        Result R = search_grid(a, g, cfg);
        auto t3 = clk::now();

        const int nbin = a.nbin, W = std::max(1, R.width);
        const double P = R.period, DM = R.dm;
        row.set("P_ms", P * 1e3, "%.10g"); row.set("DM", DM, "%.5f");
        row.set("dP_ms", (P - refP) * 1e3, "%.6g"); row.set("dDM", DM - refDM, "%.5g");
        row.set("snr", R.snr, "%.3f");
        row.set("width_bins", W, "%.0f"); row.set("width_ms", W * P / nbin * 1e3, "%.5g"); row.set("duty", double(W) / nbin, "%.4f");
        row.set("n_trials", (double)g.ndm * g.np, "%.0f");
        { int wc = 0; double sc = coho_snr(R.best_profile, &wc); row.set("snr_coho", sc, "%.3f"); row.set("width_coho", wc, "%.0f"); }
        {
            int e = 0;
            const int mdm = std::max(2, g.ndm / 50), mp = std::max(2, g.np / 50);   // within 2% of an edge
            if (g.ndm > 2 * mdm && (R.idm < mdm || R.idm >= g.ndm - mdm)) e |= 1;
            if (g.np > 2 * mp && (R.ip < mp || R.ip >= g.np - mp)) e |= 2;
            row.set("edge", std::string(e == 0 ? "0" : e == 1 ? "DM" : e == 2 ? "P" : "DM+P"));
            std::vector<float> dmcut(g.ndm), pcut(g.np);
            for (int i = 0; i < g.ndm; i++) dmcut[i] = R.snr_grid[(size_t)i * g.np + R.ip];
            for (int j = 0; j < g.np; j++) pcut[j] = R.snr_grid[(size_t)R.idm * g.np + j];
            row.set("DM_err", peak_error(dmcut.data(), g.ndm, R.idm, g.ddm), "%.4g");
            row.set("P_err_ms", peak_error(pcut.data(), g.np, R.ip, g.dp) * 1e3, "%.4g");
        }
        // the profile, its baseline window and boxcar, the measured rms
        MMWInfo mi;
        {
            std::vector<float> tmp(R.best_profile);
            smooth_mmw_ex(tmp.data(), nbin, cfg.minwidth, cfg.maxwidth, (float)cfg.rms, nullptr, &mi);
        }
        double rms_meas = NAN, chi2 = 0;
        {
            double s = 0, ss = 0;
            for (int k = 0; k < mi.ksm; k++) { double v = R.best_profile[(mi.bstart + k) % nbin]; s += v; ss += v * v; }
            double m = s / mi.ksm; rms_meas = std::sqrt(std::max(0.0, (ss - mi.ksm * m * m) / (mi.ksm - 1)));
            double pm = 0; for (float v : R.best_profile) pm += v; pm /= nbin;
            for (float v : R.best_profile) chi2 += std::pow((v - pm) / cfg.rms, 2);
            chi2 /= (nbin - 1);
        }
        {
            std::vector<float> t1(R.best_profile), t2(R.best_profile);
            int wmax = nbin / 2, wprev = cfg.minwidth;
            for (int kw = cfg.minwidth; kw <= wmax; ) { int nk = (int)(kw * 1.5f); if (nk == 1) nk = 2; if (nk > wmax) break; wprev = kw; kw = nk; }
            int wlast = wprev; for (int kw = cfg.minwidth; kw <= wmax; ) { wlast = kw; int nk = (int)(kw * 1.5f); if (nk == 1) nk = 2; kw = nk; }
            row.set("w_maxed", std::string(W >= wlast ? "1" : "0"));
        }
        row.set("baseline", mi.base, "%.6g"); row.set("rms_pred", cfg.rms, "%.5g"); row.set("rms_meas", rms_meas, "%.5g");
        row.set("rms_ratio", rms_meas / cfg.rms, "%.3f"); row.set("chi2_red", chi2, "%.3f");

        if (ex) {
            ex->g = g; ex->R = R; ex->refP = refP; ex->refDM = refDM; ex->grid_kind = kind; ex->rms = cfg.rms;
            ex->t_load = sec(t1 - t0); ex->t_prep = sec(t2 - t1); ex->t_search = sec(t3 - t2);
        }

        if (!(R.snr >= o.threshold)) {
            row.set("status", std::string("below_threshold"));
            row.set("runtime_s", sec(clk::now() - t0), "%.4f");
            if (ex) ex->t_total = sec(clk::now() - t0);
            return row;
        }

        // ---------------------------------------------------------------- above threshold: the diagnostics
        // dispersion analysis on enough channels that the archive-DM sweep within one channel is < W/2 bins
        int nfull = full.nchan, nd = std::min(nfull, std::max(o.disp_nchan, a.nchan));
        {
            double fmin = 1e30, fmax = -1e30;
            for (int c = 0; c < nfull; c++) { fmin = std::min(fmin, full.freq[c]); fmax = std::max(fmax, full.freq[c]); }
            double sweep = std::fabs(dispersion_delay(full.dm, full.centre_freq, fmin) - dispersion_delay(full.dm, full.centre_freq, fmax)) / P * nbin;
            int need = (int)std::ceil(2.0 * sweep / W);
            nd = std::min(nfull, std::max(nd, need));
            while (nd < nfull && nfull % nd) nd++;
        }
        Archive ad = full;
        if (nd < ad.nchan) fscrunch_to_nchan(ad, nd);
        row.set("nchan_disp", ad.nchan, "%.0f");

        // aligned at the best trial and at DM 0
        Aligned A = align(ad, DM, P), A0 = align(ad, 0.0, P);
        const int nsub = ad.nsub, nch = ad.nchan, on0 = mi.start;
        // off-pulse = everything outside the boxcar plus a guard of up to one width each side, always leaving
        // at least a quarter of the turn
        std::vector<char> offmask(nbin, 1);
        const int guard = std::max(0, std::min(W, (3 * nbin / 4 - W) / 2));
        for (int k = -guard; k < W + guard; k++) offmask[((on0 + k) % nbin + nbin) % nbin] = 0;
        // per-subint (single-pulse) S/N, and the same statistic for off-pulse windows
        std::vector<double> sp(nsub, NAN), spoff;
        for (int s = 0; s < nsub; s++) {
            if (A.subw[s] <= 0) continue;
            const float* p = &A.sub[(size_t)s * nbin];
            sp[s] = window_snr(p, nbin, on0, W, offmask);
            for (int j = 0; j < nbin; j += std::max(1, W / 4)) {   // overlapping windows: fine for the mean and rms
                bool inoff = true;
                for (int k = 0; k < W; k++) inoff &= offmask[(j + k) % nbin] != 0;
                if (inoff) { double v = window_snr(p, nbin, j, W, offmask, j); if (std::isfinite(v)) spoff.push_back(v); }
            }
        }
        std::vector<double> chs(nch, NAN);
        for (int c = 0; c < nch; c++) if (A.chanw[c] > 0) chs[c] = window_snr(&A.chan[(size_t)c * nbin], nbin, on0, W, offmask);
        {
            int imax = -1; double vmax = -1e30, s1 = 0, s2 = 0; int n = 0;
            for (int s = 0; s < nsub; s++) if (std::isfinite(sp[s])) { if (sp[s] > vmax) { vmax = sp[s]; imax = s; } s1 += sp[s]; s2 += sp[s] * sp[s]; n++; }
            if (imax >= 0) {
                row.set("sp_max", vmax, "%.3f"); row.set("sp_max_sub", imax, "%.0f");
                // the noise distribution of this statistic, calibrated on the off-pulse windows (Gaussian tail beyond them)
                double om = 0, os = 0;
                for (double v : spoff) { om += v; os += v * v; }
                double sig = NAN;
                if (spoff.size() > 10) {
                    om /= spoff.size(); os = std::sqrt(std::max(1e-12, os / spoff.size() - om * om));
                    double z = (vmax - om) / os;
                    double p1 = 0.5 * std::erfc(z / std::sqrt(2.0));               // one window
                    if (p1 > 1e-300) sig = p_to_sigma(-std::expm1(n * std::log1p(-p1)));   // the brightest of n subints
                    else sig = std::sqrt(std::max(0.0, z * z - 2 * std::log((double)n)));
                }
                row.set("sp_sig", sig, "%.2f");
                double mean = s1 / n, sd = std::sqrt(std::max(0.0, s2 / n - mean * mean));
                row.set("mod_index", mean > 0 ? sd / mean : NAN, "%.3f");
            }
            row.set("t_frac90", frac90(sp), "%.3f");
            row.set("f_frac90", frac90(chs), "%.3f");
        }
        // drift-scan transit
        std::vector<double> tt, yy;
        for (int s = 0; s < nsub; s++) if (std::isfinite(sp[s])) { tt.push_back(ad.epoch[s]); yy.push_back(sp[s]); }
        Drift dr = fit_drift(tt, yy, ad.end_time() - ad.start_time());
        if (dr.ok) {
            row.set("drift_t_s", dr.t0, "%.3f"); row.set("drift_t_err_s", dr.t0err, "%.3f"); row.set("drift_fwhm_s", dr.fwhm, "%.3f");
            row.set("drift_mjd_peak", ad.mjd0 + dr.t0 / 86400.0, "%.8f");
        }

        // S/N vs DM from 0 (on the dispersion channels), best width and fixed width
        SearchConfig cd = cfg; cd.rms = pdmp_rms(ad); cd.keep_grid = true;
        SearchConfig cw = cd; cw.minwidth = cw.maxwidth = W;
        double dmmax = std::max(2.0 * DM, g.dm0 + (g.ndm - 1) * g.ddm);
        Grid gd; gd.ndm = 401; gd.dm0 = 0; gd.ddm = dmmax > 0 ? dmmax / (gd.ndm - 1) : 1; gd.np = 1; gd.p0 = P; gd.dp = 0;
        // The S/N-vs-DM curve (and the period curve below) run on the SEARCH channels (<= 32) at the best period, best
        // width and fixed width W in one pass. On the nchan_disp channels they took ~90% of a single-threaded run.
        // Far from the best DM the curve then carries intra-subband smearing; the model curve does not, so compare
        // the two near the peak. The DM-0 interference test stays on the dispersion channels: one trial, and a
        // zero-DM signal smeared within 32 subbands formed at the archive DM would look falsely dispersed.
        SearchConfig ca = cfg; ca.keep_grid = true; ca.fixw = W;
        Result Rd = search_grid(a, gd, ca);
        const std::vector<float>& Rdw_grid = Rd.snr_grid_w;
        Grid g1; g1.ndm = 1; g1.dm0 = DM; g1.np = 1; g1.p0 = P;
        Result Rbw = search_grid(ad, g1, cw);                                   // fixed width at the best trial
        Grid g0 = g1; g0.dm0 = 0;
        Result R0 = search_grid(ad, g0, cw);                                    // fixed width at DM 0
        double snr_dm0 = R0.snr, snr_best_w = Rbw.snr;
        std::vector<double> wsig(nch), fr(nch);
        for (int c = 0; c < nch; c++) { fr[c] = ad.f(nsub / 2, c); wsig[c] = std::isfinite(chs[c]) ? std::max(0.0, chs[c]) : 0.0; }
        if (std::accumulate(wsig.begin(), wsig.end(), 0.0) <= 0) for (int c = 0; c < nch; c++) wsig[c] = A.chanw[c];
        const double floor_ = noise_floor(nbin, W);
        auto model = [&](double ddm) { double l = dm_loss(ddm, fr, wsig, ad.centre_freq, P, nbin, W); return std::sqrt(std::pow(snr_best_w * l, 2) + floor_ * floor_); };
        {
            const std::vector<float>& c = Rdw_grid;
            const int n = (int)c.size();
            int imax = (int)(std::max_element(c.begin(), c.end()) - c.begin());
            // the main lobe: walk down from the peak until the curve stops falling or halves
            int lo = imax, hi = imax;
            while (lo > 0 && c[lo - 1] <= c[lo] && c[lo - 1] > 0.5f * c[imax]) lo--;
            while (hi < n - 1 && c[hi + 1] <= c[hi] && c[hi + 1] > 0.5f * c[imax]) hi++;
            float second = 0;
            for (int i = 1; i < n - 1; i++)
                if ((i < lo - 1 || i > hi + 1) && c[i] >= c[i - 1] && c[i] >= c[i + 1]) second = std::max(second, c[i]);
            row.set("dm_peak2", c[imax] > 0 ? second / c[imax] : NAN, "%.3f");
        }
        row.set("snr_dm0", snr_dm0, "%.3f");
        row.set("snr_dm0_model", model(-DM), "%.3f");
        row.set("disp_sigma", (snr_best_w - snr_dm0) / std::sqrt(2.0), "%.2f");
        // S/N vs period at the best DM, and S/N vs width
        Grid gp; gp.ndm = 1; gp.dm0 = DM; gp.np = std::max(2, 4 * g.np); gp.p0 = g.p0; gp.dp = g.np > 1 ? g.dp * (g.np - 1) / (gp.np - 1) : 0;
        Result Rp = search_grid(a, gp, ca);
        std::vector<float> wc;
        for (int kw = cfg.minwidth; kw <= nbin / 2 && kw <= 0.9 * nbin; ) {
            std::vector<float> tmp(R.best_profile);
            float sn = smooth_mmw(tmp.data(), nbin, kw, kw, (float)cfg.rms, nullptr);
            wc.push_back((float)kw); wc.push_back(sn);
            kw = (int)(kw * 1.5f); if (kw == 1) kw = 2;          // smooth_mmw's width sequence
        }

        // ---------------------------------------------------------------- the plot file
        row.set("status", std::string("ok"));
        if (!o.outdir.empty()) {
            Bundle B;
            for (auto& c : csv_columns()) if (c != "plotfile" && c != "runtime_s") B.kv(c, row.get(c));
            B.kv("grid_kind", kind);
            B.kv("base", mi.base); B.kv("bstart", mi.bstart); B.kv("ksm", mi.ksm); B.kv("on_start", mi.start); B.kv("width", W);
            B.kv("grid_dm0", g.dm0); B.kv("grid_ddm", g.ddm); B.kv("grid_ndm", g.ndm);
            B.kv("grid_p0", g.p0); B.kv("grid_dp", g.dp); B.kv("grid_np", g.np);
            B.kv("best_period", P); B.kv("best_dm", DM); B.kv("rms_disp", cd.rms); B.kv("noise_floor", floor_);
            B.kv("snr_best_fixed_width", snr_best_w); B.kv("mjd0", ad.mjd0); B.kv("centre_freq", ad.centre_freq);
            B.add("profile", 1, nbin, R.best_profile);
            B.add("profile_dm0", 1, nbin, A0.prof);
            B.add("sub", nsub, nbin, A.sub);
            B.add("chan", nch, nbin, A.chan);
            B.add("chan_dm0", nch, nbin, A0.chan);
            std::vector<float> v;
            v.assign(A.subw.begin(), A.subw.end()); B.add("sub_wt", 1, nsub, v);
            v.assign(A.chanw.begin(), A.chanw.end()); B.add("chan_wt", 1, nch, v);
            v.assign(fr.begin(), fr.end()); B.add("chan_freq", 1, nch, v);
            v.assign(ad.epoch.begin(), ad.epoch.end()); B.add("epoch", 1, nsub, v);
            v.assign(ad.dur.begin(), ad.dur.end()); B.add("dur", 1, nsub, v);
            B.add("grid", g.ndm, g.np, R.snr_grid);
            v.clear();
            for (int i = 0; i < gd.ndm; i++) { double d = gd.dm0 + i * gd.ddm; v.push_back((float)d); v.push_back(Rd.snr_grid[i]); v.push_back(Rdw_grid[i]); v.push_back((float)model(d - DM)); }
            B.add("dmcurve", gd.ndm, 4, v);
            v.clear();
            for (int j = 0; j < gp.np; j++) { v.push_back((float)((gp.p0 + j * gp.dp - P) * 1e6)); v.push_back(Rp.snr_grid[j]); v.push_back(Rp.snr_grid_w[j]); }
            B.add("pcurve", gp.np, 3, v);
            B.add("wcurve", (int)wc.size() / 2, 2, wc);
            v.assign(sp.begin(), sp.end()); B.add("sp_on", 1, nsub, v);
            v.assign(spoff.begin(), spoff.end()); B.add("sp_off", 1, (int)spoff.size(), v);
            v.assign(chs.begin(), chs.end()); B.add("chan_snr", 1, nch, v);
            if (dr.ok) { B.kv("drift_amp", dr.amp); }
            std::string fn = o.outdir + "/" + basename_noext(file) + ".pdx";
            if (B.write(fn)) row.set("plotfile", fn);
        }
        row.set("status", std::string("ok"));
        row.set("runtime_s", sec(clk::now() - t0), "%.4f");
        if (ex) { ex->t_total = sec(clk::now() - t0); }
        return row;
    } catch (std::exception& e) {
        row.set("status", std::string("error: ") + e.what());
        row.set("runtime_s", sec(clk::now() - t0), "%.4f");
        return row;
    }
}
