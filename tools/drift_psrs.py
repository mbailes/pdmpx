#!/usr/bin/env python3
"""Which catalogued pulsars drift through which FAST beams, and when?

usage: drift_psrs.py FILE_OR_HDR [...] [--sizes ls-l.txt] [--radius ARCMIN] [--lon 106.8566]
  FILE_OR_HDR  SIGPROC .fil files, or their first ~1 KB saved as *.hdr (then give --sizes: the output of
               `ls -l` on the real files, so the duration is known)
For a drift (transit/arcdrift) scan the telescope is fixed, so each beam's RA advances at the sidereal rate
from its header position (src_raj, src_dej at tstart). Lists pulsars (psrcat) passing within --radius of a
beam centre, with the closest-approach time (s from the file start, and MJD) and offset.
"""
import argparse, os, re, subprocess, sys
import math

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fil_window import read_header

SIDEREAL = 1.00273790935     # sidereal seconds per solar second
FAST_FWHM_ARCMIN = 2.9       # at 1.25 GHz


def sigproc_angle(v, hours):
    """SIGPROC packs RA as hhmmss.s and Dec as ddmmss.s"""
    s = -1 if v < 0 else 1
    v = abs(v)
    d = int(v // 10000); m = int((v - d * 10000) // 100); x = v - d * 10000 - m * 100
    return s * (d + m / 60 + x / 3600) * (15 if hours else 1)   # degrees


def psrcat_positions():
    out = subprocess.run(["psrcat", "-o", "short", "-nohead", "-nonumber", "-c", "jname raj decj p0 dm"],
                         capture_output=True, text=True).stdout
    res = []
    for line in out.splitlines():
        f = line.split()
        if len(f) < 3 or not re.match(r"^[JB]\d", f[0]):
            continue
        try:
            rh = [float(x) for x in f[1].split(":")]; dd = f[2].split(":")
            ra = 15 * (rh[0] + (rh[1] if len(rh) > 1 else 0) / 60 + (rh[2] if len(rh) > 2 else 0) / 3600)
            sg = -1 if dd[0].startswith("-") else 1
            dv = [abs(float(dd[0]))] + [float(x) for x in dd[1:]]
            dec = sg * (dv[0] + (dv[1] if len(dv) > 1 else 0) / 60 + (dv[2] if len(dv) > 2 else 0) / 3600)
        except ValueError:
            continue
        res.append((f[0], ra, dec, f[3] if len(f) > 3 else "*", f[4] if len(f) > 4 else "*"))
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--sizes", help="`ls -l` listing of the real .fil files (for .hdr inputs)")
    ap.add_argument("--radius", type=float, default=3 * FAST_FWHM_ARCMIN, help="arcmin (default 3 FWHM)")
    a = ap.parse_args()
    sizes = {}
    if a.sizes:
        for line in open(a.sizes):
            f = line.split()
            if len(f) >= 9 and f[-1].endswith(".fil"):
                sizes[os.path.basename(f[-1])] = int(f[4])
    beams = []
    for fn in a.files:
        with open(fn, "rb") as f:
            h, _, hlen = read_header(f)
        base = os.path.basename(fn)
        real = base[:-4] if base.endswith(".hdr") else base
        size = sizes.get(real, os.path.getsize(fn) if fn.endswith(".fil") else None)
        if size is None:
            sys.exit(f"{fn}: need --sizes for the real file size")
        bps = h["nchans"] * h["nbits"] * h.get("nifs", 1) // 8
        T = (size - hlen) // bps * h["tsamp"]
        m = re.search(r"-(M\d\d)", real)
        beams.append(dict(name=m.group(1) if m else real, file=real, ra0=sigproc_angle(h["src_raj"], True),
                          dec=sigproc_angle(h["src_dej"], False), mjd0=h["tstart"], T=T))
    beams.sort(key=lambda b: b["name"])
    T = max(b["T"] for b in beams)
    print(f"{len(beams)} beams, {T:.0f} s ({T/3600:.2f} h) from MJD {beams[0]['mjd0']:.5f}; "
          f"RA drifts {T*SIDEREAL/240:.2f} deg; beam Decs {min(b['dec'] for b in beams):+.3f}..{max(b['dec'] for b in beams):+.3f}")
    hits = []
    for name, ra, dec, p0, dm in psrcat_positions():
        for b in beams:
            if abs(dec - b["dec"]) * 60 > a.radius:
                continue
            cosd = math.cos(math.radians(b["dec"]))
            dra = ((ra - b["ra0"] + 180) % 360) - 180             # degrees of RA the beam must drift
            t = dra * 240 / SIDEREAL                               # 1 deg of RA = 240 sidereal s
            if not (-60 <= t <= b["T"] + 60):
                continue
            off = abs(dec - b["dec"]) * 60                         # arcmin at closest approach
            hits.append((t, name, b["name"], off, p0, dm, b["mjd0"] + t / 86400, b["file"]))
    hits.sort()
    if not hits:
        print(f"no catalogued pulsar passes within {a.radius:.1f}' of any beam")
    print(f"{'t (s)':>8} {'MJD':>13} {'beam':>4} {'offset':>7} {'pulsar':12} {'P0 (s)':>12} {'DM':>8}  gain(FWHM {FAST_FWHM_ARCMIN}')")
    for t, name, bn, off, p0, dm, mjd, fn in hits:
        g = math.exp(-4 * math.log(2) * (off / FAST_FWHM_ARCMIN) ** 2)
        print(f"{t:8.1f} {mjd:13.6f} {bn:>4} {off:6.2f}' {name:12} {p0:>12} {dm:>8}  {g:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
