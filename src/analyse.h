// pdmpx: the per-archive analysis - search, the one-line summary (CSV row) and, above the S/N threshold,
// the diagnostics and the plot file (bundle).
#pragma once
#include "archive.h"
#include "search.h"
#include <map>
#include <string>
#include <vector>

struct Options {
    int mc = -1, ms = -1, mb = -1;                  // scrunch the searched data to (default: <=32 channels)
    double ds = -1, dr = -1, doff = 0;              // DM step / half-range / offset
    double ps = -1, pr = -1, poff = 0;              // period step / half-range / offset (us)
    double mw = 0, turns = 1;                       // min width (s); half-range in turns of drift
    double dmup = 100;                              // never search above the input DM + dmup (<0: no limit); never below 0
    bool natural = false, verbose = false;
    double threshold = 0;                           // S/N above which the full diagnostics + plot file are made
    int disp_nchan = 128;                           // minimum channels kept for the dispersion (DM 0) analysis
    std::string outdir;                             // where plot files go ("" = none)
    std::string optstring;                          // the command-line options, for the record
};

// one CSV row: column -> value (empty when not computed)
struct Row {
    std::map<std::string, std::string> v;
    void set(const std::string& k, const std::string& s) { v[k] = s; }
    void set(const std::string& k, double x, const char* fmt = "%.6g");
    std::string get(const std::string& k) const { auto i = v.find(k); return i == v.end() ? "" : i->second; }
};
const std::vector<std::string>& csv_columns();
std::string csv_header();
std::string csv_line(const Row& r);

// what a single-file run also wants to print or write
struct Extra {
    Grid g;
    Result R;
    double refP = 0, refDM = 0;
    double t_load = 0, t_prep = 0, t_search = 0, t_total = 0;
    int nthreads = 1;
    std::string grid_kind;
    double rms = 0;
    std::vector<float> profile;
};

// The whole analysis of one archive. Never throws: errors go into the row's status.
Row analyse_file(const std::string& file, const Options& o, Extra* extra = nullptr);
