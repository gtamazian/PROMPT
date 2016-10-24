FC = gfortran

FCFLAGS = -O3

GFORT_DFLAGS = -g -O0 -Wall -Wextra -Wconversion \
	-fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all \
	-fbackslash -ffpe-trap=zero,overflow,underflow -finit-real=nan

IFORT_DFLAGS = -g -O0 -check all -fpe0 -warn -traceback -debug extended

ifeq ($(FC), ifort)
	DFLAGS = $(IFORT_DFLAGS)
else
	DFLAGS = $(GFORT_DFLAGS)
endif

ifndef MATLAB_ROOT
	MATLAB_ROOT := $(shell which matlab | sed 'sA/bin/matlabAA')
	MATLAB_ROOT := $(strip $(MATLAB_ROOT))
endif

TARGETOS := $(shell uname -s)

ifndef BIN_SUFFIX
ifeq ($(TARGETOS), Darwin)
	BIN_SUFFIX = maci64
else
	BIN_SUFFIX = glnxa64
endif
endif

ifndef MEX_SUFFIX
ifeq ($(TARGETOS), Darwin)
	MEX_SUFFIX = maci64
else
	MEX_SUFFIX = a64
endif
endif

all: trmobjfunc.mex$(MEX_SUFFIX) trmoptim.mex$(MEX_SUFFIX)

debug: FCFLAGS += $(DFLAGS)
debug: all

prompt.o: prompt.f90
	$(FC) $(FCFLAGS) -c $< -fPIC

trmobjfunc.o: trmobjfunc.F prompt.o
	$(FC) $(FCFLAGS) -c -DMX_COMPAT_32 \
		-I"$(MATLAB_ROOT)/extern/include/" -fexceptions \
		-fPIC -fno-omit-frame-pointer -o $@ -O $<

trmoptim.o: trmoptim.F prompt.o
	$(FC) $(FCFLAGS) -c -DMX_COMPAT_32 \
		-I"$(MATLAB_ROOT)/extern/include/" -fexceptions \
		-fPIC -fno-omit-frame-pointer -o $@ -O $<

trmobjfunc.mex$(MEX_SUFFIX): trmobjfunc.o prompt.o
	$(FC) $(FCFLAGS) -pthread -shared -O $^ \
		-Wl,-rpath,$(MATLAB_ROOT)/bin/$(BIN_SUFFIX) \
		-L"$(MATLAB_ROOT)/bin/$(BIN_SUFFIX)" -lmx -lmex -lmat -lm \
		-o $@

trmoptim.mex$(MEX_SUFFIX): trmoptim.o prompt.o
	$(FC) $(FCFLAGS) -pthread -shared -O $^ \
		-Wl,-rpath,$(MATLAB_ROOT)/bin/$(BIN_SUFFIX) \
		-L"$(MATLAB_ROOT)/bin/$(BIN_SUFFIX)" -lmx -lmex -lmat -lm \
		-o $@

.PHONY: clean

clean:
	rm -f *.o *.mod *.mex$(MEX_SUFFIX)

