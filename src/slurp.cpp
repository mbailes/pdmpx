#include "slurp.h"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

Slurp::Slurp(const std::string& fn)
{
    int fd = open(fn.c_str(), O_RDONLY);
    if (fd < 0) throw std::runtime_error("cannot open " + fn);
    struct stat st;
    if (fstat(fd, &st) != 0) { close(fd); throw std::runtime_error("cannot stat " + fn); }
    n = (size_t)st.st_size;
    const char* io = getenv("PDMPX_IO");
#ifdef F_NOCACHE
    if (getenv("PDMPX_NOCACHE")) fcntl(fd, F_NOCACHE, 1);   // benchmarking: bypass the page cache (macOS)
#endif
    if (io && !strcmp(io, "read")) {
        buf.resize(n);
        size_t got = 0;
        while (got < n) {
            ssize_t k = ::read(fd, buf.data() + got, std::min(n - got, (size_t)1 << 30));   // macOS: <= INT_MAX per call
            if (k <= 0) { close(fd); throw std::runtime_error("read failed: " + fn); }
            got += (size_t)k;
        }
        p = buf.data();
    } else if (n > 0) {
        map = mmap(nullptr, n, PROT_READ, MAP_PRIVATE, fd, 0);
        if (map == MAP_FAILED) { map = nullptr; close(fd); throw std::runtime_error("mmap failed: " + fn); }
        madvise(map, n, MADV_SEQUENTIAL);
        madvise(map, n, MADV_WILLNEED);
        p = static_cast<const unsigned char*>(map);
    }
    close(fd);
}

Slurp::~Slurp()
{
    if (map) munmap(map, n);
}
