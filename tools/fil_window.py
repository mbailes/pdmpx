#!/usr/bin/env python3
"""Cut a time window out of a SIGPROC filterbank (.fil), keeping a valid header with the corrected tstart.

usage:  fil_window.py IN.fil OUT.fil --mjd MJD_CENTRE --half SECONDS
        fil_window.py IN.fil OUT.fil --start SECONDS --dur SECONDS      (seconds from the file start)
        fil_window.py IN.fil --info                                      (print the header and exit)

Plain Python 3, no dependencies: reads the header, then copies only the needed bytes (streams in 64 MB blocks),
so it is quick even on 45 GB files. Works for any nbits (whole samples are always whole bytes for
nchans*nbits a multiple of 8).
"""
import argparse, os, struct, sys

INT = {"telescope_id", "machine_id", "data_type", "barycentric", "pulsarcentric", "nbits", "nsamples",
       "nchans", "nifs", "nbeams", "ibeam", "sumifs"}
CHR = {"signed"}             # SIGPROC writes 'signed' as a single byte
DBL = {"az_start", "za_start", "src_raj", "src_dej", "tstart", "tsamp", "fch1", "foff", "refdm", "period"}
STR = {"source_name", "rawdatafile"}


def read_header(f):
    def rstr():
        n = struct.unpack("<i", f.read(4))[0]
        return f.read(n).decode()
    if rstr() != "HEADER_START":
        raise ValueError("not a SIGPROC filterbank")
    h, fields = {}, []
    while True:
        k = rstr()
        if k == "HEADER_END":
            break
        if k in INT:
            v = struct.unpack("<i", f.read(4))[0]
        elif k in DBL:
            v = struct.unpack("<d", f.read(8))[0]
        elif k in CHR:
            v = struct.unpack("<b", f.read(1))[0]
        elif k in STR:
            v = rstr()
        else:
            raise ValueError(f"unknown header keyword {k!r}")
        h[k] = v
        fields.append(k)
    return h, fields, f.tell()


def write_header(f, h, fields):
    def wstr(s):
        b = s.encode(); f.write(struct.pack("<i", len(b))); f.write(b)
    wstr("HEADER_START")
    for k in fields:
        wstr(k)
        v = h[k]
        if k in INT: f.write(struct.pack("<i", int(v)))
        elif k in DBL: f.write(struct.pack("<d", float(v)))
        elif k in CHR: f.write(struct.pack("<b", int(v)))
        else: wstr(v)
    wstr("HEADER_END")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inp"); ap.add_argument("out", nargs="?")
    ap.add_argument("--mjd", type=float, help="window centre (MJD)"); ap.add_argument("--half", type=float, help="half-width (s)")
    ap.add_argument("--start", type=float, help="window start, seconds from the file start"); ap.add_argument("--dur", type=float)
    ap.add_argument("--info", action="store_true")
    a = ap.parse_args()
    with open(a.inp, "rb") as f:
        h, fields, hlen = read_header(f)
    nchans, nbits, nifs = h["nchans"], h["nbits"], h.get("nifs", 1)
    bps = nchans * nbits * nifs // 8                       # bytes per time sample
    nsamp = (os.path.getsize(a.inp) - hlen) // bps
    T = nsamp * h["tsamp"]
    if a.info or not a.out:
        for k in fields: print(f"{k:14s} {h[k]}")
        print(f"header {hlen} bytes, {bps} bytes/sample, {nsamp} samples = {T:.1f} s = {T/3600:.2f} h, "
              f"MJD {h['tstart']:.6f} .. {h['tstart'] + T/86400:.6f}")
        return 0
    if a.mjd is not None:
        t0 = (a.mjd - h["tstart"]) * 86400 - a.half; t1 = t0 + 2 * a.half
    else:
        t0 = a.start; t1 = a.start + a.dur
    s0 = max(0, int(t0 / h["tsamp"])); s1 = min(nsamp, int(t1 / h["tsamp"]))
    if s1 <= s0:
        sys.exit(f"window {t0:.1f}..{t1:.1f} s is outside the file (0..{T:.1f} s)")
    h["tstart"] += s0 * h["tsamp"] / 86400.0
    if "nsamples" in h: h["nsamples"] = s1 - s0
    with open(a.inp, "rb") as fi, open(a.out, "wb") as fo:
        write_header(fo, h, fields)
        fi.seek(hlen + s0 * bps)
        left = (s1 - s0) * bps
        while left > 0:
            b = fi.read(min(left, 64 << 20))
            if not b: break
            fo.write(b); left -= len(b)
    print(f"{a.out}: samples {s0}..{s1} ({(s1-s0)*h['tsamp']:.1f} s, {(s1-s0)*bps/1e9:.2f} GB), tstart MJD {h['tstart']:.8f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
