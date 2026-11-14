// One FFTW plan pair per transform length, made once and shared by every thread.
// fftwf_execute_dft_r2c / _c2r on new arrays is thread-safe; only the planner is not, so it is
// called here alone, under a mutex. FFTW_UNALIGNED: the plans are executed on arrays (and offsets
// into larger buffers) other than the ones they were planned with.
#pragma once
#include <fftw3.h>

struct FFTPlans { fftwf_plan r2c = nullptr, c2r = nullptr; };
const FFTPlans& fft_plans(int n);
