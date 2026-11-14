#include "fftplan.h"
#include <map>
#include <memory>
#include <mutex>

const FFTPlans& fft_plans(int n)
{
    static std::mutex m;
    static std::map<int, std::unique_ptr<FFTPlans>> cache;   // never freed: plans live for the whole run
    std::lock_guard<std::mutex> lock(m);
    auto& p = cache[n];
    if (!p) {
        p.reset(new FFTPlans);
        float* in = fftwf_alloc_real(n);
        fftwf_complex* out = fftwf_alloc_complex(n / 2 + 1);
        p->r2c = fftwf_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE | FFTW_UNALIGNED);
        p->c2r = fftwf_plan_dft_c2r_1d(n, out, in, FFTW_ESTIMATE | FFTW_UNALIGNED);
        fftwf_free(in); fftwf_free(out);
    }
    return *p;
}
