// Whole-file access in one gulp: the file is mapped (or, with PDMPX_IO=read, read with one read() call)
// and profiles are then addressed by byte offsets into it - no per-profile I/O calls.
#pragma once
#include <cstddef>
#include <string>
#include <vector>

struct Slurp {
    const unsigned char* p = nullptr;
    size_t n = 0;
    explicit Slurp(const std::string& filename);
    ~Slurp();
    Slurp(const Slurp&) = delete;
    Slurp& operator=(const Slurp&) = delete;
private:
    void* map = nullptr;
    std::vector<unsigned char> buf;
};
