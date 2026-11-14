# pdmpx - build with:   make            (or: make CXX=g++, make OPENMP=0, make FFTW_DIR=... CFITSIO_DIR=...,
#                                         make CXXFLAGS="-O3 -march=native -std=c++17" for this CPU only)
#
# Needs: a C++17 compiler, FFTW3 (single precision, libfftw3f), cfitsio, zlib, and OpenMP (optional but
# strongly recommended - without it pdmpx runs on one core).
#   macOS (Homebrew):  brew install fftw cfitsio libomp
#   Debian/Ubuntu:     sudo apt install g++ make pkg-config libfftw3-dev libcfitsio-dev zlib1g-dev
#   OzSTAR / NT:       module load gcc/11.3.0 fftw/3.3.10 cfitsio/4.2.0   (then just: make)
#   other HPC:         module load gcc fftw cfitsio (names vary - see README), or set FFTW_DIR / CFITSIO_DIR

UNAME   := $(shell uname -s)
OPENMP  ?= 1
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wno-unused-variable
PKGCONF := $(shell command -v pkg-config 2>/dev/null)

# --- libraries: explicit directories win, then EasyBuild modules ($EBROOT...), then pkg-config, then the usual places
# (on OzSTAR/NT: module load gcc/11.3.0 fftw/3.3.10 cfitsio/4.2.0 sets EBROOTFFTW and EBROOTCFITSIO)
ifndef FFTW_DIR
  ifdef EBROOTFFTW
    FFTW_DIR := $(EBROOTFFTW)
  endif
endif
ifndef CFITSIO_DIR
  ifdef EBROOTCFITSIO
    CFITSIO_DIR := $(EBROOTCFITSIO)
  endif
endif
ifdef FFTW_DIR
  CPPFLAGS += -I$(FFTW_DIR)/include
  LDLIBS   += -L$(FFTW_DIR)/lib -L$(FFTW_DIR)/lib64 -Wl,-rpath,$(FFTW_DIR)/lib -Wl,-rpath,$(FFTW_DIR)/lib64
endif
ifdef CFITSIO_DIR
  CPPFLAGS += -I$(CFITSIO_DIR)/include
  LDLIBS   += -L$(CFITSIO_DIR)/lib -L$(CFITSIO_DIR)/lib64 -Wl,-rpath,$(CFITSIO_DIR)/lib -Wl,-rpath,$(CFITSIO_DIR)/lib64
endif
ifneq ($(PKGCONF),)
  PCFLAGS := $(shell pkg-config --cflags fftw3f cfitsio zlib 2>/dev/null)
  PCLIBS  := $(shell pkg-config --libs fftw3f cfitsio zlib 2>/dev/null)
endif
ifeq ($(UNAME),Darwin)
  BREW := $(shell brew --prefix 2>/dev/null)
  ifneq ($(BREW),)
    CPPFLAGS += -I$(BREW)/include
    LDLIBS   += -L$(BREW)/lib
  endif
endif
CPPFLAGS += $(PCFLAGS)
LDLIBS   += $(if $(PCLIBS),$(PCLIBS),-lfftw3f -lcfitsio -lz) -lm

# --- OpenMP: Apple clang needs libomp and -Xpreprocessor; GCC / LLVM clang take -fopenmp -------------------
ifeq ($(OPENMP),1)
  ifeq ($(UNAME),Darwin)
    ifneq ($(shell $(CXX) --version 2>/dev/null | grep -c 'Apple clang'),0)
      OMP := $(shell brew --prefix libomp 2>/dev/null)
      CPPFLAGS += -Xpreprocessor -fopenmp -I$(OMP)/include
      LDLIBS   += -L$(OMP)/lib -lomp
    else
      CXXFLAGS += -fopenmp
      LDFLAGS  += -fopenmp
    endif
  else
    CXXFLAGS += -fopenmp
    LDFLAGS  += -fopenmp
  endif
endif

SRC := $(wildcard src/*.cpp)
OBJ := $(SRC:.cpp=.o)

pdmpx: $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

src/%.o: src/%.cpp src/*.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f src/*.o pdmpx

# install the binary and the Python tools (PREFIX=$HOME/.local by default)
PREFIX ?= $(HOME)/.local
install: pdmpx
	mkdir -p $(PREFIX)/bin
	cp pdmpx $(PREFIX)/bin/
	for t in tools/*.py; do cp $$t $(PREFIX)/bin/; chmod +x $(PREFIX)/bin/$$(basename $$t); done
	@echo "installed to $(PREFIX)/bin (make sure it is on your PATH)"

.PHONY: clean install
