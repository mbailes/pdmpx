// pdmpx: period/DM optimisation and diagnostics for folded pulsar archives - pdmp's search, done fast.
//   Standalone (no PSRCHIVE): reads Timer (.ar) and fold-mode PSRFITS archives.
//   One archive: pdmp-style report. Many archives: one CSV row each (+ plot files above the S/N threshold).
#include "analyse.h"
#include "archive.h"
#include "search.h"
#include "snr.h"
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <set>
#include <string>
#include <sys/stat.h>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

static void usage()
{
    puts("pdmpx - fast period/DM optimiser and candidate diagnostics for folded pulsar archives (Timer .ar or PSRFITS)\n"
         "usage: pdmpx [options] archive [archive ...]\n"
         "search:\n"
         "  -mc N / -ms N / -mb N   scrunch to N channels / subints / bins first (default: <=32 channels)\n"
         "  -ds/-dr/-do  DM step / half-range / offset   -ps/-pr/-po  period step / half-range / offset (us)\n"
         "  -turns T  DM half-range in turns of delay across the band (default 1; period range = pdmp's)   -natural  pdmp's own grid, no automatic scrunch\n"
         "  -dmup X   never search more than X above the input DM (default 100; -1: no limit); never below DM 0\n"
         "  -mw S     minimum pulse width (s)   -binrot  whole-bin shifts, exactly as pdmp   -P S  period for PSRFITS without one\n"
         "batch / diagnostics:\n"
         "  -snr X    full diagnostics + plot file only when S/N >= X (default 0: always)\n"
         "  -outdir D write D/pdmpx.csv (one row per archive, resumable) and the plot files D/<archive>.pdx\n"
         "  -png      also render D/<archive>.png from each plot file (tools/pdmpx_plot.py, matplotlib)\n"
         "  -dispch N minimum channels kept for the DM-0 / dispersion analysis (default 128)\n"
         "  -j N      threads (default: all cores; with several archives, one archive per thread)\n"
         "single archive extras:  -grid F (S/N grid as text)  -profile F (best profile)  -v (verbose)");
}

static bool file_exists(const std::string& p) { struct stat st; return stat(p.c_str(), &st) == 0; }

static std::string exe_dir(const char* argv0)
{
    char buf[PATH_MAX];
    if (realpath(argv0, buf)) { std::string s(buf); size_t k = s.find_last_of('/'); return k == std::string::npos ? "." : s.substr(0, k); }
    return ".";
}

int main(int argc, char** argv)
{
    Options o;
    int threads = 0, dump_s = -1, dump_c = -1;
    bool png = false;
    std::vector<std::string> files;
    std::string gridfile, proffile;
    for (int i = 1; i < argc; i++) {
        std::string op = argv[i];
        auto nx = [&]() -> double { if (i + 1 >= argc) { usage(); exit(1); } return atof(argv[++i]); };
        auto sx = [&]() -> std::string { if (i + 1 >= argc) { usage(); exit(1); } return argv[++i]; };
        int before = i;
        if (op == "-mc") o.mc = (int)nx(); else if (op == "-ms") o.ms = (int)nx(); else if (op == "-mb") o.mb = (int)nx();
        else if (op == "-ds") o.ds = nx(); else if (op == "-dr") o.dr = nx(); else if (op == "-do") o.doff = nx();
        else if (op == "-ps") o.ps = nx(); else if (op == "-pr") o.pr = nx(); else if (op == "-po") o.poff = nx();
        else if (op == "-mw") o.mw = nx(); else if (op == "-turns") o.turns = nx();
        else if (op == "-dmup") o.dmup = nx();
        else if (op == "-natural") o.natural = true;
        else if (op == "-binrot") g_binrot = true;
        else if (op == "-P") g_force_period = nx();
        else if (op == "-snr") o.threshold = nx();
        else if (op == "-outdir") o.outdir = sx();
        else if (op == "-png") png = true;
        else if (op == "-dispch") o.disp_nchan = (int)nx();
        else if (op == "-j") threads = (int)nx();
        else if (op == "-grid") gridfile = sx();
        else if (op == "-profile") proffile = sx();
        else if (op == "-dump") { dump_s = (int)nx(); dump_c = (int)nx(); }
        else if (op == "-v") o.verbose = true;
        else if (op == "-h" || op == "--help") { usage(); return 0; }
        else if (op[0] == '-') { fprintf(stderr, "pdmpx: unknown option %s\n", op.c_str()); return 1; }
        else { files.push_back(op); continue; }
        if (op != "-outdir" && op != "-grid" && op != "-profile" && op != "-j" && op != "-png" && op != "-v") {
            o.optstring += (o.optstring.empty() ? "" : " ") + op;   // the options that change results, for the record
            for (int k = before + 1; k <= i; k++) o.optstring += std::string(" ") + argv[k];
        }
    }
    if (files.empty()) { usage(); return 1; }
#ifdef _OPENMP
    if (threads > 0) omp_set_num_threads(threads);
    int nthr = omp_get_max_threads();
#else
    int nthr = 1;
#endif

    if (dump_s >= 0) {   // debugging: one raw profile as loaded
        Archive a;
        try { a = load_archive(files[0], o.verbose); } catch (std::exception& e) { fprintf(stderr, "pdmpx: %s\n", e.what()); return 1; }
        printf("# sub %d chan %d freq %.6f wt %g epoch %.9f dur %.9f pfold %.12f mjd0 %.10f\n", dump_s, dump_c,
               a.freq[dump_c], a.wt(dump_s, dump_c), a.epoch[dump_s], a.dur[dump_s], a.pfold[dump_s], a.mjd0);
        for (int b = 0; b < a.nbin; b++) printf("%d %.7g\n", b, a.prof(dump_s, dump_c)[b]);
        return 0;
    }

    // ------------------------------------------------------------ one archive, no output directory: pdmp-style report
    if (files.size() == 1 && o.outdir.empty()) {
        Extra ex;
        Row r = analyse_file(files[0], o, &ex);
        if (r.get("status").rfind("error", 0) == 0) { fprintf(stderr, "pdmpx: %s\n", r.get("status").c_str()); return 1; }
        const Grid& g = ex.g; const Result& R = ex.R;
        printf("pdmpx: %s nsub=%s nchan=%s (searched %s) nbin=%s  DM: %d  P0: %d  = %lld trials on %d threads\n",
               r.get("source").c_str(), r.get("nsub").c_str(), r.get("nchan").c_str(), r.get("nchan_search").c_str(),
               r.get("nbin").c_str(), g.ndm, g.np, (long long)g.ndm * g.np, nthr);
        printf("grid: %s  DM %.4f..%.4f step %.6g   P %.10g..%.10g ms step %.6g us\n", ex.grid_kind.c_str(),
               g.dm0, g.dm0 + (g.ndm - 1) * g.ddm, g.ddm, g.p0 * 1e3, (g.p0 + (g.np - 1) * g.dp) * 1e3, g.dp * 1e6);
        printf("Best S/N = %.2f  width = %d bins   (coho's definition: S/N %s, width %s bins)\n", R.snr, R.width,
               r.get("snr_coho").c_str(), r.get("width_coho").c_str());
        printf("Best TC Period (ms) = %.10g  Correction (ms) = %.10g  Error (ms) = %s\n", R.period * 1e3, (R.period - ex.refP) * 1e3, r.get("P_err_ms").c_str());
        printf("Best DM = %.4f  Correction = %.4f  Error = %s\n", R.dm, R.dm - ex.refDM, r.get("DM_err").c_str());
        printf("Noise: rms predicted %s, measured %s (ratio %s), chi2_red %s\n", r.get("rms_pred").c_str(), r.get("rms_meas").c_str(),
               r.get("rms_ratio").c_str(), r.get("chi2_red").c_str());
        if (r.get("status") == "ok") {
            printf("Dispersion: S/N at DM 0 = %s (a real pulse: %s), dispersed at %s sigma [%s channels]\n", r.get("snr_dm0").c_str(),
                   r.get("snr_dm0_model").c_str(), r.get("disp_sigma").c_str(), r.get("nchan_disp").c_str());
            printf("Subints: brightest S/N %s (#%s, %s sigma), modulation index %s, t_frac90 %s, f_frac90 %s\n", r.get("sp_max").c_str(),
                   r.get("sp_max_sub").c_str(), r.get("sp_sig").c_str(), r.get("mod_index").c_str(), r.get("t_frac90").c_str(), r.get("f_frac90").c_str());
            if (!r.get("drift_t_s").empty())
                printf("Drift transit: peak t = %s +- %s s (MJD %s), FWHM %s s\n", r.get("drift_t_s").c_str(), r.get("drift_t_err_s").c_str(),
                       r.get("drift_mjd_peak").c_str(), r.get("drift_fwhm_s").c_str());
        }
        printf("time: load %.3f s, prepare %.3f s, search %.3f s, total %.3f s\n", ex.t_load, ex.t_prep, ex.t_search, ex.t_total);
        if (!proffile.empty()) {
            FILE* f = fopen(proffile.c_str(), "w");
            for (int b = 0; b < (int)R.best_profile.size(); b++) fprintf(f, "%d %.7g\n", b, R.best_profile[b]);
            fclose(f);
        }
        if (!gridfile.empty()) {
            FILE* f = fopen(gridfile.c_str(), "w");
            for (int i = 0; i < g.ndm; i++) {
                for (int j = 0; j < g.np; j++)
                    fprintf(f, "%.10g %.10g %.6g\n", (g.p0 + j * g.dp) * 1e6, g.dm0 + i * g.ddm, R.snr_grid[(size_t)i * g.np + j]);
                fprintf(f, "\n");
            }
            fclose(f);
        }
        return 0;
    }

    // ------------------------------------------------------------ batch: one CSV row per archive
    FILE* csv = stdout;
    std::set<std::string> done;
    std::string csvpath;
    if (!o.outdir.empty()) {
        mkdir(o.outdir.c_str(), 0755);
        csvpath = o.outdir + "/pdmpx.csv";
        bool exists = file_exists(csvpath);
        if (exists) {   // resume: skip archives already in the table
            std::ifstream in(csvpath);
            std::string line;
            std::getline(in, line);
            while (std::getline(in, line)) {
                std::string f;
                if (!line.empty() && line[0] == '"') { size_t e = line.find("\",", 1); f = line.substr(1, e == std::string::npos ? std::string::npos : e - 1); }
                else f = line.substr(0, line.find(','));
                done.insert(f);
            }
        }
        csv = fopen(csvpath.c_str(), "a");
        if (!csv) { fprintf(stderr, "pdmpx: cannot write %s\n", csvpath.c_str()); return 1; }
        if (!exists) fprintf(csv, "%s\n", csv_header().c_str());
    } else {
        fprintf(csv, "%s\n", csv_header().c_str());
    }
    std::vector<std::string> todo;
    for (auto& f : files) if (!done.count(f)) todo.push_back(f);
    if (!o.outdir.empty())
        fprintf(stderr, "pdmpx: %zu archives to do (%zu already in %s), %d threads, S/N threshold %g\n", todo.size(),
                files.size() - todo.size(), csvpath.c_str(), nthr, o.threshold);

    std::vector<std::string> plots;
    int nok = 0, nbelow = 0, nerr = 0;
    // progress meter on stderr: one line, refreshed about once a second (and at the end)
    const bool meter = !o.outdir.empty() && todo.size() > 1;
    const auto tstart = std::chrono::steady_clock::now();
    double last_print = -1;
    auto progress = [&](bool final) {
        double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - tstart).count();
        if (!final && el - last_print < 1.0) return;
        last_print = el;
        size_t done_n = nok + nbelow + nerr;
        double rate = el > 0 ? done_n / el : 0, eta = rate > 0 ? (todo.size() - done_n) / rate : 0;
        fprintf(stderr, "\r[%zu/%zu] %5.1f%%  %.1f/s  elapsed %dm%02ds  left %dm%02ds   above %d  below %d  errors %d   ",
                done_n, todo.size(), 100.0 * done_n / todo.size(), rate, (int)el / 60, (int)el % 60, (int)eta / 60, (int)eta % 60,
                nok, nbelow, nerr);
        if (final) fprintf(stderr, "\n");
        fflush(stderr);
    };
    FILE* trace = getenv("PDMPX_TRACE") ? fopen(getenv("PDMPX_TRACE"), "w") : nullptr;
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < todo.size(); i++) {
        if (trace) {   // PDMPX_TRACE=file: record each archive as it starts (to find one that crashes the run)
            #pragma omp critical(trace)
            { fprintf(trace, "start %s\n", todo[i].c_str()); fflush(trace); }
        }
        Row r = analyse_file(todo[i], o);   // searches inside run single-threaded here (no nested parallelism)
        if (trace) {
            #pragma omp critical(trace)
            { fprintf(trace, "end %s\n", todo[i].c_str()); fflush(trace); }
        }
        #pragma omp critical(csv_out)
        {
            fprintf(csv, "%s\n", csv_line(r).c_str());
            fflush(csv);
            std::string st = r.get("status");
            if (st == "ok") nok++; else if (st == "below_threshold") nbelow++; else { nerr++; fprintf(stderr, "\npdmpx: %s: %s\n", todo[i].c_str(), st.c_str()); }
            if (!r.get("plotfile").empty()) plots.push_back(r.get("plotfile"));
            if (meter) progress(false);
        }
    }
    if (meter) progress(true);
    if (csv != stdout) fclose(csv);
    if (!o.outdir.empty()) fprintf(stderr, "pdmpx: done: %d above threshold, %d below, %d errors\n", nok, nbelow, nerr);

    if (png && !plots.empty()) {
        std::string script = getenv("PDMPX_PLOT") ? getenv("PDMPX_PLOT") : exe_dir(argv[0]) + "/tools/pdmpx_plot.py";
        std::string list = o.outdir + "/.pdmpx_png_list";
        FILE* f = fopen(list.c_str(), "w");
        for (auto& p : plots) fprintf(f, "%s\n", p.c_str());
        fclose(f);
        std::string cmd = "python3 \"" + script + "\" -j " + std::to_string(nthr) + " -list \"" + list + "\"";
        fprintf(stderr, "pdmpx: rendering %zu PNGs\n", plots.size());
        int rc = system(cmd.c_str());
        if (rc != 0) fprintf(stderr, "pdmpx: PNG rendering failed (%d)\n", rc);
        remove(list.c_str());
    }
    return 0;
}
