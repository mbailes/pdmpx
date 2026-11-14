#pragma once
#include "archive.h"

// pdmp's matched-filter S/N (PSRCHIVE smooth_mmw). pr is modified (baseline removed).
float smooth_mmw(float* pr, int nbin, int minwidth, int maxwidth, float rms, int* kwmax);
// the same, also reporting where the baseline window and the winning boxcar are
struct MMWInfo { float base = 0; int bstart = 0, ksm = 0, start = 0, width = 0; };
float smooth_mmw_ex(float* pr, int nbin, int minwidth, int maxwidth, float rms, int* kwmax, MMWInfo* info);
// pdmp's getRMS: noise rms of the weighted-mean profile
double pdmp_rms(const Archive& a);
