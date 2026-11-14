#!/usr/bin/env python3
"""Render pdmpx plot files (.pdx) as diagnostic PNG pages.

usage: pdmpx_plot.py [-j N] [-cycles 1.5] [-list FILE] [plotfile.pdx ...]
Each FILE.pdx becomes FILE.png. The page is drawn only from the plot file - no archive needed.
"""
import argparse, gzip, os, sys
import numpy as np


# ------------------------------------------------------------------------------------------- the plot file
def read_pdx(path):
    """-> (meta dict of strings, arrays dict of float32 numpy arrays)"""
    with gzip.open(path, "rb") as f:
        if f.readline().strip() != b"PDMPX-BUNDLE 1":
            raise ValueError(f"{path}: not a pdmpx plot file")
        meta = {}
        while True:
            line = f.readline().decode()
            if line.strip() == "ARRAYS" or not line:
                break
            k, _, v = line.rstrip("\n").partition("=")
            meta[k] = v
        arrays = {}
        while True:
            h = f.readline()
            if not h:
                break
            name, ny, nx = h.decode().split()
            ny, nx = int(ny), int(nx)
            a = np.frombuffer(f.read(4 * ny * nx), dtype="<f4").astype(float)
            arrays[name] = a.reshape(ny, nx) if ny > 1 else a
    return meta, arrays


def fnum(meta, k, default=np.nan):
    try:
        return float(meta.get(k, ""))
    except ValueError:
        return default


# ------------------------------------------------------------------------------------------- the page
def render(path, cycles=1.5):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    m, A = read_pdx(path)
    nbin = A["profile"].size
    W, on0, ksm, bst = int(fnum(m, "width")), int(fnum(m, "on_start")), int(fnum(m, "ksm")), int(fnum(m, "bstart"))
    rms = fnum(m, "rms_pred"); P = fnum(m, "best_period"); DM = fnum(m, "best_dm"); snr = fnum(m, "snr")
    # centre the boxcar at phase 0.5: every phase-resolved array is rolled by the same number of bins
    shift = int(round(nbin / 2 - (on0 + W / 2.0)))
    roll = lambda x: np.roll(x, shift, axis=-1)
    lo = 0.5 - cycles / 2                                     # displayed phase range: [lo, lo + cycles]

    def tile(img):                                            # one stored cycle -> the displayed cycles
        img = np.atleast_2d(img)
        reps = int(np.ceil(cycles)) + 2
        big = np.concatenate([img] * reps, axis=-1)           # phases [-1, reps-1)
        ph = (np.arange(big.shape[-1]) + 0.5) / nbin - 1
        keep = (ph >= lo) & (ph < lo + cycles)
        return big[..., keep], ph[keep]

    def rowclean(img, okrows=None):                           # each row: subtract its mean, scale by its rms
        img = img - img.mean(axis=1, keepdims=True)
        sd = img.std(axis=1, keepdims=True); sd[sd == 0] = 1
        img = img / sd
        if okrows is not None:
            img[~okrows] = np.nan
        return img

    def scrunch_rows(img, wts, nrow):
        n = img.shape[0]
        if n <= nrow:
            return img, wts > 0, np.arange(n + 1)
        e = np.linspace(0, n, nrow + 1).astype(int)
        out = np.array([(img[a:b] * wts[a:b, None]).sum(0) / max(wts[a:b].sum(), 1e-30) for a, b in zip(e[:-1], e[1:])])
        ok = np.array([wts[a:b].sum() > 0 for a, b in zip(e[:-1], e[1:])])
        return out, ok, e

    nrow = lambda n: int(np.clip((snr / 3) ** 2, 8, max(8, n)))   # ~S/N 3 per row
    def turn_lines(ax):
        for x in (0, 1):
            if lo < x < lo + cycles:
                ax.axvline(x, color="tab:cyan", lw=0.6, alpha=0.7)
        for a, b in ((lo, 0), (1, lo + cycles)):              # grey the repeated parts
            if b > a:
                ax.axvspan(a, b, color="0.5", alpha=0.10, lw=0)

    fig = plt.figure(figsize=(21, 11))
    gs = GridSpec(3, 5, figure=fig, height_ratios=[1.05, 1, 0.9], hspace=0.42, wspace=0.34, left=0.04, right=0.985, top=0.92, bottom=0.05)
    fname = os.path.basename(m.get("file", path))
    fig.suptitle(f"{fname}\n{m.get('source','')}   S/N {snr:.2f}   P {P*1e3:.6f} ms   DM {DM:.3f}   status {m.get('status','')}",
                 fontsize=11, x=0.04, ha="left")

    # --- profile with what went into the S/N
    ax = fig.add_subplot(gs[0, 0:2])
    prof = roll(A["profile"]); prof0 = roll(A["profile_dm0"]); base = fnum(m, "base")
    y, ph = tile((prof - base) / rms); y = y[0]
    y0, _ = tile((prof0 - np.sort(prof0)[: ksm].mean()) / rms); y0 = y0[0]
    ax.step(ph, y0, where="mid", color="0.65", lw=0.8, label="same data at DM 0")
    ax.step(ph, y, where="mid", color="k", lw=1.0, label=f"best: DM {DM:.3f}, P {P*1e3:.6f} ms")
    on_c = (np.arange(on0, on0 + W) + shift) % nbin
    for b in (np.arange(bst, bst + ksm) + shift) % nbin:
        ax.axvspan(b / nbin, (b + 1) / nbin, color="tab:blue", alpha=0.10, lw=0)
    for b in on_c:
        ax.axvspan(b / nbin, (b + 1) / nbin, color="tab:orange", alpha=0.4, lw=0)
    rr = fnum(m, "rms_ratio")
    ax.axhline(0, color="tab:blue", lw=1); ax.axhspan(-1, 1, color="tab:blue", alpha=0.12)
    if np.isfinite(rr):
        ax.axhspan(-rr, rr, fill=False, ls="--", ec="tab:red", lw=0.8)
    turn_lines(ax); ax.set_xlim(lo, lo + cycles)
    ax.set_xlabel("phase (boxcar centred at 0.5)"); ax.set_ylabel("amplitude / predicted rms")
    ax.set_title(f"Profile   S/N {snr:.2f}, width {W} bins = {fnum(m,'width_ms'):.4g} ms ({100*W/nbin:.1f}%)", fontsize=10)
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    ax.text(0.01, 0.97, f"blue: baseline window ({ksm} bins), ±1 predicted rms  ·  red dashed: ±1 measured rms (ratio {rr:.2f})  ·  orange: the {W} bins summed",
            transform=ax.transAxes, va="top", fontsize=8)

    # --- phase vs time
    sub, sw, ep, du = roll(A["sub"]), A["sub_wt"], A["epoch"], A["dur"]
    sub = np.atleast_2d(sub)
    img, okt, et = scrunch_rows(sub, sw, nrow(sub.shape[0]))
    img, _ = tile(rowclean(img, okt))
    ax = fig.add_subplot(gs[0, 2])
    T = ep[-1] + du[-1] / 2
    ax.imshow(img, aspect="auto", origin="lower", cmap="Greys", interpolation="nearest", extent=[lo, lo + cycles, 0, T],
              vmin=np.nanpercentile(img, 1), vmax=np.nanpercentile(img, 99.7))
    turn_lines(ax); ax.set_xlabel("phase"); ax.set_ylabel("time (s)"); ax.set_title(f"Phase vs time ({img.shape[0]} rows)", fontsize=10)

    # --- per-subint S/N, cumulative, drift transit
    ax = fig.add_subplot(gs[0, 3])
    sp = A["sp_on"]
    ax.plot(ep, sp, "o-", ms=3, lw=0.8, label="on-pulse S/N per subint")
    ok = np.isfinite(sp)
    if ok.sum() > 1:
        cum = np.nancumsum(np.where(ok, sp, 0)) / np.sqrt(np.cumsum(ok))
        ax.plot(ep, cum, "-", color="tab:green", lw=1.2, label="cumulative (Σ S/N / √n)")
        ax.plot(ep, cum[-1] * np.sqrt(np.cumsum(ok) / ok.sum()), ":", color="tab:green", lw=1, label="√t expectation")
    if m.get("drift_t_s"):
        t0, fw, amp = fnum(m, "drift_t_s"), fnum(m, "drift_fwhm_s"), fnum(m, "drift_amp", np.nanmax(sp))
        tt = np.linspace(ep[0], ep[-1], 300); sg = fw / 2.3548
        base_d = np.nanmedian(sp[np.abs(ep - t0) > fw]) if np.any(np.abs(ep - t0) > fw) else 0
        ax.plot(tt, amp * np.exp(-0.5 * ((tt - t0) / sg) ** 2) + base_d, "r-", lw=1,
                label=f"beam transit: t={t0:.1f}±{fnum(m,'drift_t_err_s'):.1f} s, FWHM {fw:.1f} s")
    ax.axhline(0, color="0.7", lw=0.6); ax.set_xlabel("time (s)"); ax.set_ylabel("S/N")
    ax.set_title("Per-subint (single-pulse) S/N", fontsize=10); ax.legend(fontsize=7, frameon=False, loc="upper left")

    # --- single pulses vs noise
    ax = fig.add_subplot(gs[0, 4])
    off = A.get("sp_off", np.array([]))
    allv = np.concatenate([off[np.isfinite(off)], sp[ok]]) if off.size else sp[ok]
    if allv.size:
        hb = np.linspace(allv.min() - 0.5, allv.max() + 0.5, 40)
        if off.size:
            ax.hist(off, bins=hb, density=True, color="0.7", label="off-pulse windows (noise)")
        ax.hist(sp[ok], bins=hb, density=True, histtype="step", color="tab:orange", lw=1.5, label="on-pulse, each subint")
        xx = np.linspace(hb[0], hb[-1], 200); ax.plot(xx, np.exp(-xx ** 2 / 2) / np.sqrt(2 * np.pi), "k:", lw=0.8, label="N(0,1)")
    ax.set_xlabel("window S/N"); ax.legend(fontsize=7, frameon=False)
    ax.set_title(f"Single pulses: brightest {m.get('sp_max','')} = {m.get('sp_sig','?')} σ", fontsize=10)

    # --- stats
    ax = fig.add_subplot(gs[1, 0:2]); ax.axis("off")
    g = lambda k: m.get(k, "") or "–"
    lines = [
        f"source {g('source')}   MJD {g('mjd')}   T = {g('tobs_s')} s   fc {g('fc_MHz')} MHz   BW {g('bw_MHz')} MHz",
        f"data   {g('nsub')} sub × {g('nchan')} ch × {g('nbin')} bins   zapped {g('zap_frac')}   searched {g('nchan_search')} ch   dispersion test {g('nchan_disp')} ch",
        f"folded P {g('P_fold_ms')} ms   DM {g('DM_fold')}",
        f"best   P {g('P_ms')} ± {g('P_err_ms')} ms   DM {g('DM')} ± {g('DM_err')}   (ΔP {g('dP_ms')} ms, ΔDM {g('dDM')})   grid edge: {g('edge')}",
        f"S/N {g('snr')}   width {g('width_bins')} bins = {g('width_ms')} ms (duty {g('duty')})   trials {g('n_trials')}",
        f"noise  baseline {g('baseline')}   rms predicted {g('rms_pred')}  measured {g('rms_meas')}  ratio {g('rms_ratio')}   χ²_red {g('chi2_red')}",
        f"dispersion  S/N at DM 0 {g('snr_dm0')} (a real pulse: {g('snr_dm0_model')})   dispersed at {g('disp_sigma')} σ",
        f"subints  brightest {g('sp_max')} (#{g('sp_max_sub')}, {g('sp_sig')} σ)   modulation index {g('mod_index')}   t_frac90 {g('t_frac90')}   f_frac90 {g('f_frac90')}",
    ]
    if m.get("drift_t_s"):
        lines.append(f"drift  transit peak t = {g('drift_t_s')} ± {g('drift_t_err_s')} s (MJD {g('drift_mjd_peak')}), FWHM {g('drift_fwhm_s')} s")
    lines.append(f"significance (noise-search FAP): not yet calibrated      pdmpx {g('version')}   {g('options')}")
    ax.text(0, 1, "\n".join(lines), va="top", family="monospace", fontsize=8.3, transform=ax.transAxes)

    # --- phase vs frequency at the best DM and at DM 0
    cf, cw = A["chan_freq"], A["chan_wt"]
    for k, (key, title) in enumerate([("chan", f"Phase vs freq at best DM {DM:.2f}"), ("chan_dm0", "Phase vs freq at DM 0 (+ predicted sweep)")]):
        ch = np.atleast_2d(roll(A[key]))
        img, okf, ef = scrunch_rows(ch, cw, nrow(ch.shape[0]))
        fr = np.array([cf[a:b][cw[a:b] > 0].mean() if cw[a:b].sum() > 0 else np.nan for a, b in zip(ef[:-1], ef[1:])])
        img, _ = tile(rowclean(img, okf))
        ax = fig.add_subplot(gs[1, 2 + k])
        ax.imshow(img, aspect="auto", origin="lower", cmap="Greys", interpolation="nearest", extent=[lo, lo + cycles, 0, img.shape[0]],
                  vmin=np.nanpercentile(img, 1), vmax=np.nanpercentile(img, 99.7))
        yl = np.arange(img.shape[0]) + 0.5
        if k == 1:
            fc = fnum(m, "centre_freq")
            sweep = 0.5 + DM / 2.41e-4 * (fr ** -2 - fc ** -2) / P
            for wrap in range(-3, 4):
                ax.plot(sweep + wrap + 0.06, yl, "r-", lw=0.7, alpha=0.6)
            ax.set_xlim(lo, lo + cycles)
        step = max(1, len(yl) // 6)
        ax.set_yticks(yl[::step]); ax.set_yticklabels([f"{f:.0f}" for f in fr[::step]])
        turn_lines(ax); ax.set_xlabel("phase"); ax.set_ylabel("freq (MHz)"); ax.set_title(title, fontsize=10)

    # --- spectrum
    ax = fig.add_subplot(gs[1, 4])
    ax.plot(A["chan_snr"], cf, ".", ms=4, color="tab:purple")
    ax.axvline(0, color="0.7", lw=0.6); ax.set_ylabel("freq (MHz)"); ax.set_xlabel("on-pulse S/N per channel")
    ax.set_title("Spectrum", fontsize=10)

    # --- bullseye
    ax = fig.add_subplot(gs[2, 0])
    grid = np.atleast_2d(A["grid"])
    dms = fnum(m, "grid_dm0") + np.arange(grid.shape[0]) * fnum(m, "grid_ddm")
    pus = (fnum(m, "grid_p0") + np.arange(grid.shape[1]) * fnum(m, "grid_dp") - P) * 1e6
    ax.imshow(grid, aspect="auto", origin="lower", cmap="viridis", extent=[pus[0], pus[-1], dms[0], dms[-1]])
    if grid.shape[0] > 1 and grid.shape[1] > 1:
        ax.contour(pus, dms, grid, levels=[0.5 * grid.max(), 0.8 * grid.max(), 0.95 * grid.max()], colors="w", linewidths=0.6)
    ax.plot(0, DM, "r+", ms=12); ax.set_xlabel("P − P_best (µs)"); ax.set_ylabel("DM"); ax.set_title("P–DM S/N (bullseye)", fontsize=10)

    # --- S/N vs DM (from 0), S/N vs P, S/N vs width
    ax = fig.add_subplot(gs[2, 1]); d = A["dmcurve"]
    ax.plot(d[:, 0], d[:, 1], "k-", lw=1, label="best width per trial")
    ax.plot(d[:, 0], d[:, 2], "-", color="tab:orange", lw=1, label=f"fixed width {W}")
    ax.plot(d[:, 0], d[:, 3], "r--", lw=1, label="real pulse (measured spectrum) + floor")
    ax.axvline(DM, color="0.6", lw=0.6); ax.axvline(0, color="0.6", lw=0.6)
    ax.set_xlabel("DM"); ax.set_ylabel("S/N"); ax.set_title("S/N vs DM (from 0)", fontsize=10); ax.legend(fontsize=7, frameon=False)
    ax = fig.add_subplot(gs[2, 2]); p = A["pcurve"]
    ax.plot(p[:, 0], p[:, 1], "k-", lw=1, label="best width"); ax.plot(p[:, 0], p[:, 2], "-", color="tab:orange", lw=1, label=f"fixed width {W}")
    ax.set_xlabel("P − P_best (µs)"); ax.set_ylabel("S/N"); ax.set_title("S/N vs period", fontsize=10); ax.legend(fontsize=7, frameon=False)
    ax = fig.add_subplot(gs[2, 3]); w = np.atleast_2d(A["wcurve"])
    ax.plot(w[:, 0], w[:, 1], "o-", ms=4, color="k"); ax.axvline(W, color="tab:orange", lw=0.8)
    ax.set_xscale("log"); ax.set_xlabel("boxcar width (bins)"); ax.set_ylabel("S/N"); ax.set_title("S/N vs boxcar width", fontsize=10)
    ax = fig.add_subplot(gs[2, 4]); ax.axis("off")
    ax.text(0.5, 0.5, "False-alarm panel\n(noise-search calibration)\n\nnot yet implemented", ha="center", va="center", fontsize=9,
            color="0.45", transform=ax.transAxes, bbox=dict(ec="0.7", fc="none", ls="--"))

    out = os.path.splitext(path)[0] + ".png"
    fig.savefig(out, dpi=100)
    plt.close(fig)
    return out


def _render(args):
    path, cycles = args
    try:
        return render(path, cycles)
    except Exception as e:  # one bad file must not stop the batch
        return f"ERROR {path}: {e}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*")
    ap.add_argument("-list", help="file with one plot-file path per line")
    ap.add_argument("-j", type=int, default=os.cpu_count())
    ap.add_argument("-cycles", type=float, default=1.5)
    a = ap.parse_args()
    files = list(a.files)
    if a.list:
        files += [l.strip() for l in open(a.list) if l.strip()]
    if not files:
        ap.error("no plot files")
    jobs = [(f, a.cycles) for f in files]
    if a.j > 1 and len(files) > 1:
        from multiprocessing import Pool
        with Pool(min(a.j, len(files))) as pool:
            res = pool.map(_render, jobs, chunksize=1)
    else:
        res = [_render(j) for j in jobs]
    bad = [r for r in res if r.startswith("ERROR")]
    for r in bad:
        print(r, file=sys.stderr)
    print(f"pdmpx_plot: {len(res) - len(bad)} PNGs written, {len(bad)} errors", file=sys.stderr)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
