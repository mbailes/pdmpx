// Reader for PSRCHIVE "Timer" archives (the format dspsr writes by default, .ar).
// Layout (from psrchive/Base/Formats/Timer): a 1024-byte 'timer' header, an optional backend block
// (be_data_size bytes when the backend code is recognised), the polyco and ephemeris text
// (nbytespoly + nbytesephem bytes), then nsub_int sub-integrations, each a 128-byte 'mini' header and
//   wts_and_bpass != 0: wts[nchan], med[npol][nchan], bpass[npol][nchan], then for each pol, chan a
//                        profile: cfreq(f4) nbin(i4) poln(i4) wt(f4) scale(f4) offset(f4) u16[nbin]
//                        with amps = u16/scale + offset;
//   wts_and_bpass == 0:  scale(f4) offset(f4) u16[npol*nchan*nbin], amps = offset + u16*scale.
// Endianness: archives are written big-endian; detected from the header's nbin.
#include "archive.h"
#include "slurp.h"
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace {
struct Reader {          // a cursor into the slurped file
    const unsigned char* b; size_t n, pos; bool swap;
    void raw(void* p, size_t k) { if (pos + k > n) throw std::runtime_error("Timer: unexpected end of file"); memcpy(p, b + pos, k); pos += k; }
    const unsigned char* here(size_t k) { if (pos + k > n) throw std::runtime_error("Timer: unexpected end of file"); auto q = b + pos; pos += k; return q; }
    template <class T> T get() {
        T v; raw(&v, sizeof v);
        if (swap) { auto* b = reinterpret_cast<unsigned char*>(&v); for (size_t i = 0; i < sizeof v / 2; i++) std::swap(b[i], b[sizeof v - 1 - i]); }
        return v;
    }
    void skip(long k) { if (k < 0 || pos + (size_t)k > n) throw std::runtime_error("Timer: unexpected end of file"); pos += (size_t)k; }
};
template <class T> T rd(const unsigned char* h, size_t off, bool swap) {
    T v; memcpy(&v, h + off, sizeof v);
    if (swap) { auto* b = reinterpret_cast<unsigned char*>(&v); for (size_t i = 0; i < sizeof v / 2; i++) std::swap(b[i], b[sizeof v - 1 - i]); }
    return v;
}
// byte offsets in struct timer (timer.h; 8-byte aligned Fortran-compatible layout)
// ram_boards 0..32, corr_boards 32..64, machine_id 64..72, version 72, minorversion 76, tape 80, file 84,
// utdate 88..104, fracmjd 104, mjd 112, ticks 116, offset 120, lst 128, coord 136..144, psrname 144..160,
// ra 160, dec 168, l 176, b 180, nominal_period 184, dm 192, ftr 196, npl 200, nps 204, nbin 208, tsmp 212,
// sub_int_time 216, ndump 220, narchive 224, junk 228, nsub_int 232, junk2 236, dump_time 240, nfreq 244,
// nsub_band 248, feedmode 252, tree 256..264, telid 264..280, tpover 280..288, nspan 288, ncoeff 292,
// nbytespoly 296, nbytesephem 300, banda 304.. (struct band), bandb, ...
const size_t O_FRACMJD = 104, O_MJD = 112, O_PSRNAME = 144, O_NOMP = 184, O_DM = 192, O_NBIN = 208,
             O_NSUB = 232, O_NCHAN = 248, O_NBPOLY = 296, O_NBEPH = 300, O_BANDA = 304;
// struct band: lo1 0, lo2 8, loUP 16, loDOWN 24, centrefreq 32, bw 40, flux_A 48, inv_mode 52, auto_atten 56,
// correlator_mode 60, f_atten_A 64, f_atten_B 68, polar 72, feed_offset 76, nlag 80, flux_B 84, flux_err 88, npol 92
const size_t B_CF = 32, B_BW = 40, B_NPOL = 92;
}

static size_t band_size() { return 96 + 32; }  // checked below against the 'corrected' field position

Archive load_timer(const std::string& fn, bool verbose)
{
    Slurp file(fn);                       // the whole archive in one gulp
    if (file.n < 1024) throw std::runtime_error("Timer: short header");
    const unsigned char* h = file.p;

    // endianness from nbin (a small positive integer)
    bool swap = false;
    {
        // nbin is a small positive number: read the wrong way round it becomes large (1024 -> 262144) or
        // negative, so take the byte order that gives the smaller positive value
        uint32_t nb; memcpy(&nb, h + O_NBIN, 4);
        uint32_t sw = (nb >> 24) | ((nb >> 8) & 0xff00) | ((nb << 8) & 0xff0000) | (nb << 24);
        auto ok = [](uint32_t v) { return v > 0 && v < (1u << 31); };
        swap = !ok(nb) || (ok(sw) && sw < nb);
    }
    Reader r{file.p, file.n, 1024, swap};
    Archive a;
    a.filename = fn; a.format = "Timer";
    char name[17] = {0}; memcpy(name, h + O_PSRNAME, 16); a.source = name;
    while (!a.source.empty() && (a.source.back() == ' ' || a.source.back() == 0)) a.source.pop_back();
    a.nbin = rd<int32_t>(h, O_NBIN, swap);
    a.nsub = rd<int32_t>(h, O_NSUB, swap);
    a.nchan = rd<int32_t>(h, O_NCHAN, swap);
    a.dm = rd<float>(h, O_DM, swap);
    a.centre_freq = rd<double>(h, O_BANDA + B_CF, swap);
    a.bandwidth = rd<double>(h, O_BANDA + B_BW, swap);
    int npol = rd<int32_t>(h, O_BANDA + B_NPOL, swap);
    // the polarisation state is banda.correlator_mode (PSRCHIVE TimerArchive::get_state):
    // 0 Intensity, 1 PPQQ, 2 Coherence (AA BB CR CI), 5 Invariant, 6 Stokes (IQUV).  Total intensity =
    // the first profile for Intensity/Invariant/Stokes, the sum of the first two for PPQQ/Coherence.
    const int corr_mode = rd<int32_t>(h, O_BANDA + 60, swap);
    const int npol_use = (npol == 2 || (npol == 4 && (corr_mode == 1 || corr_mode == 2))) ? 2 : 1;
    a.mjd0 = rd<int32_t>(h, O_MJD, swap) + rd<double>(h, O_FRACMJD, swap);
    int nbpoly = rd<int32_t>(h, O_NBPOLY, swap), nbeph = rd<int32_t>(h, O_NBEPH, swap);

    // the tail of the header: find 'corrected', wts_and_bpass, backend, be_data_size by scanning from the end:
    // ... wts_and_bpass(4) wtscheme(4) software(128) backend(8) be_data_size(4) rcvr_id(8) space(184) = 1024
    const size_t O_SPACE = 1024 - 184, O_RCVR = O_SPACE - 8, O_BESIZE = O_RCVR - 4, O_BACKEND = O_BESIZE - 8,
                 O_SOFT = O_BACKEND - 128, O_WTSCHEME = O_SOFT - 4, O_WTSBP = O_WTSCHEME - 4,
                 O_SCAL = O_WTSBP - 24, O_CALF = O_SCAL - 24, O_OBST = O_CALF - 4, O_CALIB = O_OBST - 4,
                 O_CORR = O_CALIB - 4;
    int corrected = rd<int32_t>(h, O_CORR, swap);
    int wts_and_bpass = rd<int32_t>(h, O_WTSBP, swap);
    char backend[9] = {0}; memcpy(backend, h + O_BACKEND, 8);
    uint32_t be_size = rd<uint32_t>(h, O_BESIZE, swap);
    a.dedispersed = (corrected & 0x0010) != 0;
    (void)band_size;

    // backend block: PSRCHIVE skips be_data_size only for recognised backend codes ("baseband"/"FAKE"/... ).
    // Treat any non-empty backend code with a plausible size as present.
    std::string be(backend); while (!be.empty() && (be.back() == ' ' || be.back() == 0)) be.pop_back();
    // PSRCHIVE (timer++.C) skips be_data_size only for the backends it recognises: "baseband", "scint"
    if (be == "baseband" || be == "scint") r.skip(be_size);
    // the predictor and ephemeris are only present for pulsar observations (obstype PULSAR == 0)
    int obstype = rd<int32_t>(h, O_OBST, swap);
    if (obstype == 0) r.skip((long)std::max(nbpoly, 0) + std::max(nbeph, 0));

    if (verbose)
        fprintf(stderr, "Timer: %s nsub=%d nchan=%d nbin=%d npol=%d dm=%.4f cf=%.3f bw=%.3f wts_and_bpass=%d corr_mode=%d "
                "backend='%s' be_size=%u poly=%d eph=%d dedispersed=%d swap=%d\n", a.source.c_str(), a.nsub, a.nchan,
                a.nbin, npol, a.dm, a.centre_freq, a.bandwidth, wts_and_bpass, corr_mode, be.c_str(), be_size, nbpoly, nbeph,
                (int)a.dedispersed, (int)swap);
    if (a.nsub <= 0 || a.nchan <= 0 || a.nbin <= 0 || npol <= 0 || npol > 4)
        throw std::runtime_error("Timer: implausible dimensions");

    a.freq.assign(a.nchan, 0); a.wts.assign((size_t)a.nsub * a.nchan, 1.f);
    a.data.assign((size_t)a.nsub * a.nchan * a.nbin, 0.f);
    a.epoch.resize(a.nsub); a.dur.resize(a.nsub); a.pfold.resize(a.nsub);
    std::vector<double> mjd(a.nsub);
    // every subint has the same size, so subint s starts at base + s*subsize: decode them in parallel
    const size_t prof_hdr = 24;
    const size_t subsize = 128 + (wts_and_bpass ? (size_t)a.nchan * 4 + (size_t)npol * a.nchan * 8 +
                                                  (size_t)npol * a.nchan * (prof_hdr + 2 * (size_t)a.nbin)
                                                : 8 + 2 * (size_t)npol * a.nchan * a.nbin);
    const size_t base = r.pos;
    if (base + subsize * a.nsub > file.n) throw std::runtime_error("Timer: file shorter than nsub_int sub-integrations");
    std::string err;
    #pragma omp parallel for schedule(dynamic, 1)
    for (int s = 0; s < a.nsub; s++) {
      try {
        Reader q{file.p, file.n, base + subsize * (size_t)s, swap};
        const unsigned char* m = q.here(128);
        // mini: junk 0, junk2 4, junk3 8, mjd 12, fracmjd 16, lst 24, pfold 32, tel_zen 40, tel_az 44, feed 48,
        // para 52, version 56, pulse_phase 60, flux_A 64, flux_B 68, integration 72, space 80..128
        mjd[s] = rd<int32_t>(m, 12, swap) + rd<double>(m, 16, swap);
        a.pfold[s] = rd<double>(m, 32, swap);
        a.dur[s] = rd<double>(m, 72, swap);
        if (wts_and_bpass) {
            q.skip((long)a.nchan * 4 + (long)npol * a.nchan * 4 * 2);   // wts, medians, bandpass
            for (int p = 0; p < npol; p++)
                for (int c = 0; c < a.nchan; c++) {
                    const unsigned char* ph = q.here(prof_hdr);
                    float cf = rd<float>(ph, 0, swap); int nb = rd<int32_t>(ph, 4, swap);
                    float wt = rd<float>(ph, 12, swap), scale = rd<float>(ph, 16, swap), offset = rd<float>(ph, 20, swap);
                    if (nb != a.nbin) throw std::runtime_error("Timer: profile nbin mismatch");
                    const unsigned char* pk = q.here(2 * (size_t)nb);
                    if (scale == 0 || !std::isfinite(scale)) throw std::runtime_error("Timer: bad profile scale");
                    if (s == 0 && p == 0) a.freq[c] = cf;
                    if (p == 0) a.wt(s, c) = wt;
                    if (p >= npol_use) continue;
                    float* d = a.prof(s, c);
                    const float inv = 1.f / scale;
                    if (swap) for (int b = 0; b < nb; b++) d[b] += (float)(uint16_t)((pk[2 * b] << 8) | pk[2 * b + 1]) * inv + offset;
                    else      for (int b = 0; b < nb; b++) d[b] += (float)(uint16_t)(pk[2 * b] | (pk[2 * b + 1] << 8)) * inv + offset;
                }
        } else {
            const unsigned char* sh = q.here(8);
            float scale = rd<float>(sh, 0, swap), offset = rd<float>(sh, 4, swap);
            const unsigned char* pk = q.here(2 * (size_t)npol * a.nchan * a.nbin);
            size_t idx = 0;
            for (int p = 0; p < npol; p++)
                for (int c = 0; c < a.nchan; c++) {
                    float* d = a.prof(s, c);
                    for (int b = 0; b < a.nbin; b++, idx++) {
                        uint16_t v = swap ? (uint16_t)((pk[2 * idx] << 8) | pk[2 * idx + 1]) : (uint16_t)(pk[2 * idx] | (pk[2 * idx + 1] << 8));
                        if (p < npol_use) d[b] += offset + v * scale;
                    }
                }
            if (s == 0)
                for (int c = 0; c < a.nchan; c++)
                    a.freq[c] = a.centre_freq - a.bandwidth / 2 + (c + 0.5) * a.bandwidth / a.nchan;
        }
      } catch (std::exception& e) {
        #pragma omp critical(timer_err)
        if (err.empty()) err = e.what();
      }
    }
    if (!err.empty()) throw std::runtime_error(err);
    // epochs: mini MJD is the subint epoch (mid-time, PSRCHIVE Integration::get_epoch)
    double start = (mjd[0] - a.mjd0) * 86400.0 - 0.5 * a.dur[0];
    for (int s = 0; s < a.nsub; s++) a.epoch[s] = (mjd[s] - a.mjd0) * 86400.0 - start;
    a.mjd0 += start / 86400.0;
    return a;
}
