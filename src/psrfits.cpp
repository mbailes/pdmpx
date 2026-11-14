// Reader for fold-mode PSRFITS (cfitsio). Total intensity from POL_TYPE AABBCRCI / AA+BB / IQUV / INTEN.
// Folding period per subint: the SUBINT PERIOD column if present, else the POLYCO table (nearest set,
// f = REF_F0 + (1/60) sum i c_i dt^(i-1), dt in minutes). A tempo2 T2PREDICT-only file needs -P.
#include "archive.h"
#include "slurp.h"
#include <cstdint>
#include <cstdlib>
#include <fitsio.h>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <mutex>

namespace {
struct Fits {
    fitsfile* f = nullptr; int st = 0;
    void chk(const char* what) {
        if (st) { char m[FLEN_STATUS]; fits_get_errstatus(st, m); throw std::runtime_error(std::string("PSRFITS ") + what + ": " + m); }
    }
    bool has_col(const char* name, int& col) { int s = 0; fits_get_colnum(f, CASEINSEN, (char*)name, &col, &s); return s == 0; }
    bool key_d(const char* k, double& v) { int s = 0; fits_read_key(f, TDOUBLE, k, &v, nullptr, &s); return s == 0; }
    bool key_i(const char* k, long& v) { int s = 0; fits_read_key(f, TLONG, k, &v, nullptr, &s); return s == 0; }
    bool key_s(const char* k, std::string& v) { char b[FLEN_VALUE]; int s = 0; fits_read_key(f, TSTRING, k, b, nullptr, &s); if (s) return false; v = b; return true; }
    bool goto_ext(const char* name) { int s = 0; fits_movnam_hdu(f, BINARY_TBL, (char*)name, 0, &s); return s == 0; }
};
}

// --- direct table access: the SUBINT data are read straight from the slurped file -------------------------
// cfitsio parses the header and supplies each column's byte offset within the row, its type and repeat;
// the (big-endian) bytes are then decoded in place, all rows in parallel. Returns false to fall back.
static inline uint64_t be_u64(const unsigned char* p) { uint64_t v = 0; for (int i = 0; i < 8; i++) v = (v << 8) | p[i]; return v; }
static inline uint32_t be_u32(const unsigned char* p) { return (uint32_t)p[0] << 24 | (uint32_t)p[1] << 16 | (uint32_t)p[2] << 8 | p[3]; }
static inline double cell(const unsigned char* p, int type, long i)
{
    switch (type) {
    case TDOUBLE: { uint64_t u = be_u64(p + 8 * i); double d; memcpy(&d, &u, 8); return d; }
    case TFLOAT:  { uint32_t u = be_u32(p + 4 * i); float f; memcpy(&f, &u, 4); return f; }
    case TSHORT:  return (int16_t)(p[2 * i] << 8 | p[2 * i + 1]);
    case TLONG: case TINT: return (int32_t)be_u32(p + 4 * i);
    case TBYTE:   return p[i];
    default: return NAN;
    }
}

static bool read_subints_direct(fitsfile* ff, const std::string& fn, Archive& a, std::vector<double>& offs,
                                int c_tsub, int c_offs, int c_freq, int c_wts, int c_scl, int c_off, int c_data,
                                int c_per, long npol, int npol_use)
{
    int st = 0;
    LONGLONG head0, data0, end0;
    if (fits_get_hduaddrll(ff, &head0, &data0, &end0, &st)) return false;
    long rowlen = 0;
    if (fits_read_key(ff, TLONG, "NAXIS1", &rowlen, nullptr, &st)) return false;
    tcolumn* tc = ff->Fptr->tableptr;
    if (!tc) return false;
    const int nchan = a.nchan, nbin = a.nbin, nsub = a.nsub;
    auto ok = [&](int c, long need) { return c > 0 && tc[c - 1].tdatatype > 0 && tc[c - 1].trepeat >= need; };
    if (!ok(c_tsub, 1) || !ok(c_offs, 1) || !ok(c_freq, nchan) || !ok(c_wts, nchan) || !ok(c_scl, nchan * npol) ||
        !ok(c_off, nchan * npol) || !ok(c_data, (long)nchan * npol * nbin) || (c_per && !ok(c_per, 1)))
        return false;
    for (int c : {c_tsub, c_offs, c_freq, c_wts, c_scl, c_off})
        if (tc[c - 1].tscale != 1 || tc[c - 1].tzero != 0) return false;
    const int dtype = tc[c_data - 1].tdatatype;
    if (dtype != TSHORT && dtype != TFLOAT && dtype != TBYTE) return false;
    const double dscale = tc[c_data - 1].tscale, dzero = tc[c_data - 1].tzero;

    Slurp file(fn);
    if ((size_t)data0 + (size_t)rowlen * nsub > file.n) return false;
    auto col = [&](const unsigned char* row, int c) { return row + tc[c - 1].tbcol; };

    #pragma omp parallel for schedule(dynamic, 1)
    for (int s = 0; s < nsub; s++) {
        const unsigned char* row = file.p + data0 + (size_t)rowlen * s;
        a.dur[s] = cell(col(row, c_tsub), tc[c_tsub - 1].tdatatype, 0);
        offs[s] = cell(col(row, c_offs), tc[c_offs - 1].tdatatype, 0);
        if (c_per) a.pfold[s] = cell(col(row, c_per), tc[c_per - 1].tdatatype, 0);
        const unsigned char *pf = col(row, c_freq), *pw = col(row, c_wts), *ps = col(row, c_scl), *po = col(row, c_off),
                            *pd = col(row, c_data);
        const int tf = tc[c_freq - 1].tdatatype, tw = tc[c_wts - 1].tdatatype, tsc = tc[c_scl - 1].tdatatype, to = tc[c_off - 1].tdatatype;
        for (int c = 0; c < nchan; c++) {
            if (s == 0) a.freq[c] = cell(pf, tf, c);
            a.wt(s, c) = (float)cell(pw, tw, c);
            float* d = a.prof(s, c);
            for (int p = 0; p < npol_use; p++) {
                const double scl = cell(ps, tsc, (long)p * nchan + c);
                const float sc = (float)(scl * dscale), of = (float)(cell(po, to, (long)p * nchan + c) + scl * dzero);
                const size_t i0 = ((size_t)p * nchan + c) * nbin;
                if (dtype == TSHORT) {
                    const unsigned char* q = pd + 2 * i0;
                    for (int b = 0; b < nbin; b++) d[b] += (float)(int16_t)(q[2 * b] << 8 | q[2 * b + 1]) * sc + of;
                } else
                    for (int b = 0; b < nbin; b++) d[b] += (float)cell(pd, dtype, (long)(i0 + b)) * sc + of;
            }
        }
    }
    return true;
}

double g_force_period = 0;

Archive load_psrfits(const std::string& fn, bool verbose)
{
    // cfitsio is only safe to call from several threads when built reentrant; otherwise one file at a time
    static std::mutex fits_mutex;
    std::unique_lock<std::mutex> lock(fits_mutex, std::defer_lock);
    if (!fits_is_reentrant()) lock.lock();
    Fits F; Archive a; a.filename = fn; a.format = "PSRFITS";
    fits_open_file(&F.f, fn.c_str(), READONLY, &F.st); F.chk("open");
    double imjd = 0, smjd = 0, soffs = 0, obsfreq = 0, obsbw = 0, chan_dm = 0;
    F.key_d("STT_IMJD", imjd); F.key_d("STT_SMJD", smjd); F.key_d("STT_OFFS", soffs);
    F.key_d("OBSFREQ", obsfreq); F.key_d("OBSBW", obsbw); F.key_d("CHAN_DM", chan_dm);
    F.key_s("SRC_NAME", a.source);
    a.mjd0 = imjd + (smjd + soffs) / 86400.0;
    a.centre_freq = obsfreq; a.bandwidth = obsbw; a.dm = chan_dm;

    // dedispersion state + DM from the last HISTORY row (PSRCHIVE writes DEDISP, DM_MODEL etc.)
    if (F.goto_ext("HISTORY")) {
        long nrows = 0; int s = 0; fits_get_num_rows(F.f, &nrows, &s);
        int col;
        if (nrows > 0 && F.has_col("DEDISP", col)) { int v = 0; fits_read_col(F.f, TINT, col, nrows, 1, 1, nullptr, &v, nullptr, &s); a.dedispersed = v != 0; }
        // PSRCHIVE takes the centre frequency and bandwidth from the last HISTORY row (dspsr leaves OBSFREQ = 0)
        double cf = 0, cbw = 0; int nch = 0;
        if (nrows > 0 && F.has_col("CTR_FREQ", col)) fits_read_col(F.f, TDOUBLE, col, nrows, 1, 1, nullptr, &cf, nullptr, &s);
        if (nrows > 0 && F.has_col("CHAN_BW", col)) fits_read_col(F.f, TDOUBLE, col, nrows, 1, 1, nullptr, &cbw, nullptr, &s);
        if (nrows > 0 && F.has_col("NCHAN", col)) fits_read_col(F.f, TINT, col, nrows, 1, 1, nullptr, &nch, nullptr, &s);
        if (s == 0 && cf > 0) a.centre_freq = cf;
        if (s == 0 && cbw != 0 && nch > 0) a.bandwidth = cbw * nch;
    }
    if (!F.goto_ext("SUBINT")) throw std::runtime_error("PSRFITS: no SUBINT table (search-mode or not PSRFITS?)");
    long nbin = 0, nchan = 0, npol = 0, nrows = 0;
    F.key_i("NBIN", nbin); F.key_i("NCHAN", nchan); F.key_i("NPOL", npol);
    fits_get_num_rows(F.f, &nrows, &F.st); F.chk("rows");
    double dmk; if (F.key_d("DM", dmk) && dmk != 0) a.dm = dmk;
    std::string pol_type; F.key_s("POL_TYPE", pol_type);
    a.nsub = (int)nrows; a.nchan = (int)nchan; a.nbin = (int)nbin;
    if (a.nsub <= 0 || a.nchan <= 0 || a.nbin <= 1) throw std::runtime_error("PSRFITS: not a fold-mode file");
    int npol_use = 1;  // AABBCRCI, AA+BB: sum first two; IQUV/INTEN: first only
    if ((pol_type == "AABBCRCI" || pol_type == "AABB" || pol_type == "XXYY" || pol_type == "RRLL") && npol >= 2) npol_use = 2;

    int c_tsub, c_offs, c_freq, c_wts, c_scl, c_off, c_data, c_per = 0;
    if (!F.has_col("TSUBINT", c_tsub) || !F.has_col("OFFS_SUB", c_offs) || !F.has_col("DAT_FREQ", c_freq) ||
        !F.has_col("DAT_WTS", c_wts) || !F.has_col("DAT_SCL", c_scl) || !F.has_col("DAT_OFFS", c_off) || !F.has_col("DATA", c_data))
        throw std::runtime_error("PSRFITS: SUBINT table lacks a required column");
    bool have_period = F.has_col("PERIOD", c_per);

    a.freq.assign(a.nchan, 0); a.wts.assign((size_t)a.nsub * a.nchan, 1.f);
    a.data.assign((size_t)a.nsub * a.nchan * a.nbin, 0.f);
    a.epoch.resize(a.nsub); a.dur.resize(a.nsub); a.pfold.assign(a.nsub, 0);
    std::vector<double> offs(a.nsub);
    if (!getenv("PDMPX_CFITSIO") && read_subints_direct(F.f, fn, a, offs, c_tsub, c_offs, c_freq, c_wts, c_scl, c_off, c_data,
                                                         have_period ? c_per : 0, npol, npol_use)) {
        // done: the whole table decoded from the mapped file
    } else {
    // fallback (unusual column types, or PDMPX_CFITSIO set): cfitsio, one call per column per row
    std::vector<double> fr(a.nchan); std::vector<float> w(a.nchan), scl(nchan * npol), off(nchan * npol);
    std::vector<float> raw((size_t)nbin * nchan * npol);
    for (long r = 1; r <= nrows; r++) {
        int s = (int)r - 1;
        fits_read_col(F.f, TDOUBLE, c_tsub, r, 1, 1, nullptr, &a.dur[s], nullptr, &F.st);
        fits_read_col(F.f, TDOUBLE, c_offs, r, 1, 1, nullptr, &offs[s], nullptr, &F.st);
        fits_read_col(F.f, TDOUBLE, c_freq, r, 1, nchan, nullptr, fr.data(), nullptr, &F.st);
        fits_read_col(F.f, TFLOAT, c_wts, r, 1, nchan, nullptr, w.data(), nullptr, &F.st);
        fits_read_col(F.f, TFLOAT, c_scl, r, 1, nchan * npol, nullptr, scl.data(), nullptr, &F.st);
        fits_read_col(F.f, TFLOAT, c_off, r, 1, nchan * npol, nullptr, off.data(), nullptr, &F.st);
        fits_read_col(F.f, TFLOAT, c_data, r, 1, (LONGLONG)raw.size(), nullptr, raw.data(), nullptr, &F.st);
        if (have_period) fits_read_col(F.f, TDOUBLE, c_per, r, 1, 1, nullptr, &a.pfold[s], nullptr, &F.st);
        F.chk("read row");
        if (s == 0) for (int c = 0; c < a.nchan; c++) a.freq[c] = fr[c];
        for (int c = 0; c < a.nchan; c++) {
            a.wt(s, c) = w[c];
            float* d = a.prof(s, c);
            for (int p = 0; p < npol_use; p++) {
                const float* src = &raw[((size_t)p * nchan + c) * nbin];
                float sc = scl[p * nchan + c], of = off[p * nchan + c];
                for (int b = 0; b < a.nbin; b++) d[b] += src[b] * sc + of;
            }
        }
    }
    }
    // last resort for the band: the channel frequencies themselves
    if (!(a.centre_freq > 0) || a.bandwidth == 0) {
        double lo = a.freq[0], hi = a.freq[0];
        for (double x : a.freq) { lo = std::min(lo, x); hi = std::max(hi, x); }
        double chbw = a.nchan > 1 ? (hi - lo) / (a.nchan - 1) : 0;
        if (!(a.centre_freq > 0)) a.centre_freq = 0.5 * (lo + hi);
        if (a.bandwidth == 0) a.bandwidth = (a.freq.back() < a.freq.front() ? -1 : 1) * chbw * a.nchan;
    }
    // a PERIOD column of zeros (PSRCHIVE writes one when it folded from a predictor) is no period at all
    if (have_period) for (int s = 0; s < a.nsub; s++) if (!(a.pfold[s] > 0)) { have_period = false; break; }
    // epochs relative to the first subint's start (OFFS_SUB is the subint centre after STT)
    double start = offs[0] - 0.5 * a.dur[0];
    for (int s = 0; s < a.nsub; s++) a.epoch[s] = offs[s] - start;
    a.mjd0 += start / 86400.0;

    if (!have_period && g_force_period > 0) { a.pfold.assign(a.nsub, g_force_period); have_period = true; }
    if (!have_period) {
        if (!F.goto_ext("POLYCO")) {
            fits_close_file(F.f, &F.st);
            throw std::runtime_error("PSRFITS: no PERIOD column and no POLYCO table (tempo2 predictor?) - give -P <period s>");
        }
        long np = 0; int st = 0; fits_get_num_rows(F.f, &np, &st);
        int c_mjd, c_f0, c_coef, c_nc;
        if (!F.has_col("REF_MJD", c_mjd) || !F.has_col("REF_F0", c_f0) || !F.has_col("COEFF", c_coef) || !F.has_col("NCOEF", c_nc))
            throw std::runtime_error("PSRFITS: POLYCO table lacks columns");
        struct Poly { double mjd, f0; std::vector<double> c; };
        std::vector<Poly> P(np);
        for (long r = 1; r <= np; r++) {
            Poly& q = P[r - 1]; int nc = 0;
            fits_read_col(F.f, TDOUBLE, c_mjd, r, 1, 1, nullptr, &q.mjd, nullptr, &st);
            fits_read_col(F.f, TDOUBLE, c_f0, r, 1, 1, nullptr, &q.f0, nullptr, &st);
            fits_read_col(F.f, TINT, c_nc, r, 1, 1, nullptr, &nc, nullptr, &st);
            q.c.resize(nc); fits_read_col(F.f, TDOUBLE, c_coef, r, 1, nc, nullptr, q.c.data(), nullptr, &st);
        }
        if (st || P.empty()) throw std::runtime_error("PSRFITS: could not read POLYCO");
        for (int s = 0; s < a.nsub; s++) {
            double t = a.mjd0 + a.epoch[s] / 86400.0;
            const Poly* best = &P[0];
            for (auto& q : P) if (std::fabs(q.mjd - t) < std::fabs(best->mjd - t)) best = &q;
            double dtmin = (t - best->mjd) * 1440.0, df = 0, x = 1;
            for (size_t i = 1; i < best->c.size(); i++) { df += i * best->c[i] * x; x *= dtmin; }
            a.pfold[s] = 1.0 / (best->f0 + df / 60.0);
        }
    }
    fits_close_file(F.f, &F.st);
    if (verbose)
        fprintf(stderr, "PSRFITS: %s nsub=%d nchan=%d nbin=%d npol=%ld (%s, using %d) dm=%.4f cf=%.3f bw=%.3f dedispersed=%d P=%.9f s (%s)\n",
                a.source.c_str(), a.nsub, a.nchan, a.nbin, npol, pol_type.c_str(), npol_use, a.dm, a.centre_freq, a.bandwidth,
                (int)a.dedispersed, a.pfold[a.nsub / 2], have_period ? "PERIOD column" : "POLYCO");
    return a;
}

Archive load_archive(const std::string& fn, bool verbose)
{
    FILE* f = fopen(fn.c_str(), "rb");
    if (!f) throw std::runtime_error("cannot open " + fn);
    char m[6] = {0}; size_t n = fread(m, 1, 6, f); fclose(f);
    if (n == 6 && memcmp(m, "SIMPLE", 6) == 0) return load_psrfits(fn, verbose);
    return load_timer(fn, verbose);
}
