#!/usr/bin/env python3
"""Quick text summary of a pdmpx.csv: ranked table, simple flags, period families, distributions.

usage: pdmpx_summary.py results/pdmpx.csv [-n 30]
Flags (heuristics, for triage only):
  RFI?    not dispersed: S/N at DM 0 >= 0.8 x best, or disp_sigma < 1.5 with best DM < 10
  NB      narrowband: 90% of the on-pulse signal in < 25% of the channels
  BURSTY  90% of the signal in < 25% of the subints (or brightest subint > 5 sigma with low t_frac90)
  WIDE    duty cycle > 30%
  SINE    S/N still rising at the widest boxcar (sinusoid-like, not a pulse)
  DMCOMB  second separate peak in S/N vs DM > 0.7 of the highest (broad periodic signal)
  EDGE    best trial on the edge of the search grid
  noise   S/N < 7 and not dispersed above 3 sigma
  LOOK    S/N >= 7, dispersed >= 3 sigma, rms ratio < 1.3, not narrowband
"""
import argparse, csv, math, os, re
from collections import Counter


def fv(r, k):
    try:
        return float(r[k])
    except (KeyError, ValueError, TypeError):
        return math.nan


def flags(r):
    f = []
    snr, d0, ds, dm = fv(r, "snr"), fv(r, "snr_dm0"), fv(r, "disp_sigma"), fv(r, "DM")
    if (d0 >= 0.8 * snr) or (ds < 1.5 and dm < 10):
        f.append("RFI?")
    if fv(r, "f_frac90") < 0.25:
        f.append("NB")
    if fv(r, "t_frac90") < 0.25:
        f.append("BURSTY")
    if fv(r, "duty") > 0.30:
        f.append("WIDE")
    if r.get("w_maxed") == "1":
        f.append("SINE")          # S/N still rising at the widest boxcar: sinusoid / baseline structure, not a pulse
    if fv(r, "dm_peak2") > 0.7:
        f.append("DMCOMB")        # several comparable DM peaks: a broad periodic signal re-aligning every turn
    if r.get("edge", "0") not in ("0", ""):
        f.append("EDGE")
    if fv(r, "rms_ratio") > 1.3:
        f.append("RMS")
    if snr >= 7 and ds >= 3 and fv(r, "rms_ratio") < 1.3 and not {"NB", "RFI?", "SINE", "DMCOMB"} & set(f):
        f.append("LOOK")
    elif snr < 7 and not (ds >= 3):
        f.append("noise")
    return f


def short(path):
    b = os.path.basename(path)
    m = re.search(r"_(\d+)_candid_", b)
    return m.group(1) if m else b[:24]


def families(rows, tol=2e-4, maxh=16):
    """group candidates whose periods are in small-integer ratios (harmonics / RFI families)"""
    P = [(fv(r, "P_ms"), short(r["file"])) for r in rows if fv(r, "P_ms") > 0]
    P.sort()
    used, fams = set(), []
    for i, (p, n) in enumerate(P):
        if n in used:
            continue
        fam = [(p, n, "1")]
        for q, m in P:
            if m == n or m in used:
                continue
            for a in range(1, maxh + 1):
                for b in range(1, maxh + 1):
                    simple = min(a, b) == 1 or max(a, b) <= 4      # harmonics n/1, 1/n, or low-order ratios
                    if simple and abs(q / p - a / b) < tol * a / b and math.gcd(a, b) == 1:
                        fam.append((q, m, f"{a}/{b}"))
                        break
                else:
                    continue
                break
        if len(fam) > 1:
            fams.append(fam)
            used.update(x[1] for x in fam)
    return fams


def hist(vals, edges, label):
    c = [0] * (len(edges) - 1)
    for v in vals:
        for i in range(len(edges) - 1):
            if edges[i] <= v < edges[i + 1]:
                c[i] += 1
                break
    out = [f"{label}:"]
    scale = max(1.0, max(c) / 60.0)          # bars at most 60 characters wide
    for i, n in enumerate(c):
        out.append(f"  {edges[i]:>6g} – {edges[i+1]:<6g} {'#' * int(round(n / scale))} {n}")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("-n", type=int, default=1000, help="rows in the ranked table")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.csv)))
    ok = [r for r in rows if r["status"] in ("ok", "below_threshold")]
    bad = [r for r in rows if r not in ok]
    ok.sort(key=lambda r: -fv(r, "snr"))
    print(f"{len(rows)} archives: {sum(r['status']=='ok' for r in rows)} with diagnostics, "
          f"{sum(r['status']=='below_threshold' for r in rows)} below threshold, {len(bad)} errors;"
          f"  total runtime {sum(fv(r,'runtime_s') for r in rows if fv(r,'runtime_s')==fv(r,'runtime_s')):.1f} s (single-thread sum)\n")
    hdr = f"{'cand':>5} {'S/N':>6} {'P_ms':>12} {'DM':>8} {'dDM':>8} {'W':>3} {'duty':>5} {'W_ms':>7} {'rmsR':>5} {'S/N@0':>6} {'disp':>6} {'spmax':>6} {'spσ':>5} {'t90':>5} {'f90':>5}  flags"
    print(hdr); print("-" * len(hdr))
    for r in ok[: a.n]:
        print(f"{short(r['file']):>5} {fv(r,'snr'):6.2f} {fv(r,'P_ms'):12.6f} {fv(r,'DM'):8.2f} {fv(r,'dDM'):8.1f} {r['width_bins']:>3} "
              f"{fv(r,'duty'):5.3f} {fv(r,'width_ms'):7.3g} {fv(r,'rms_ratio'):5.2f} {fv(r,'snr_dm0'):6.2f} {fv(r,'disp_sigma'):6.2f} "
              f"{fv(r,'sp_max'):6.2f} {fv(r,'sp_sig'):5.1f} {fv(r,'t_frac90'):5.2f} {fv(r,'f_frac90'):5.2f}  {' '.join(flags(r))}")
    for r in bad:
        print(f"{short(r['file']):>5}  {r['status']}")

    print()
    fl = Counter(f for r in ok for f in flags(r))
    print("flags: " + ", ".join(f"{k} {v}" for k, v in fl.most_common()))
    look = [short(r["file"]) for r in ok if "LOOK" in flags(r)]
    print(f"LOOK (worth opening first): {', '.join(look) if look else 'none'}")
    print()
    fams = families(ok)
    if fams:
        fams.sort(key=lambda f: -len(f))
        print(f"period families (harmonics, or RFI such as 50 Hz mains): {len(fams)} families, "
              f"{sum(len(f) for f in fams)} candidates; the largest 12:")
        for fam in fams[:12]:
            print(f"  {len(fam):4d} members, fundamental {fam[0][0]:.4f} ms ({1000/fam[0][0]:.3f} Hz): "
                  + ", ".join(f"#{n}" for p, n, ratio in fam[:12]) + (" ..." if len(fam) > 12 else ""))
        print()
    print(hist([fv(r, "snr") for r in ok], [4, 5, 6, 7, 8, 10, 15, 30, 1e9], "S/N"))
    print(hist([int(r["width_bins"]) for r in ok if r["width_bins"]], [1, 2, 3, 5, 8, 13, 20, 30, 45, 65], "best width (bins of 128)"))
    print(hist([fv(r, "disp_sigma") for r in ok if not math.isnan(fv(r, 'disp_sigma'))], [-5, 0, 1, 2, 3, 5, 10, 1e9], "dispersed (sigma)"))


if __name__ == "__main__":
    main()
