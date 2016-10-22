FC = gfortran

FCFLAGS = -g -O3 -fimplicit-none

# The following values are given for the OS X systm. For Linux, they
# should be the following:
#   BIN_SUFFIX = glnxa64
#   MEX_SUFFIX = a64
# Also one should specify the appropriate MATLAB location in the
# MATLAB_ROOT variable.

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

debug: FCFLAGS += -O0 -fbounds-check -Wall -Wextra -pedantic \
    -fbacktrace
debug: all

prompt.o: prompt.f90
	$(FC) $(FCFLAGS) -c $< -fPIC

trmobjfunc.o: trmobjfunc.F prompt.o
	$(FC) $(FCFLAGS) -c -DMX_COMPAT_32 \
		-I"$(MATLAB_ROOT)/extern/include/" -fexceptions \
		-fbackslash -fPIC -fno-omit-frame-pointer -o $@ -O $<

trmoptim.o: trmoptim.F prompt.o
	$(FC) $(FCFLAGS) -c -DMX_COMPAT_32 \
		-I"$(MATLAB_ROOT)/extern/include/" -fexceptions \
		-fbackslash -fPIC -fno-omit-frame-pointer -o $@ -O $<

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
	rm *.o *.mod *.mex$(MEX_SUFFIX)

