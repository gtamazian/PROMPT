FC = gfortran

FCFLAGS = -g -O3 -fimplicit-none

MATLAB_ROOT = /Applications/MATLAB_R2016b.app

all: trmobjfunc.mexmaci64 trmoptim.mexmaci64

debug: FCFLAGS += -Og -fbounds-check -Wall -Wextra -pedantic \
    -fbacktrace
debug: all

prompt.o: prompt.f90
	$(FC) $(FCFLAGS) -c $<

trmobjfunc.o: trmobjfunc.F prompt.o
	$(FC) $(FCFLAGS) -c -DMX_COMPAT_32 \
		-I"$(MATLAB_ROOT)/extern/include/" -fexceptions \
		-fbackslash -fPIC -fno-omit-frame-pointer -o $@ -O $<

trmoptim.o: trmoptim.F prompt.o
	$(FC) $(FCFLAGS) -c -DMX_COMPAT_32 \
		-I"$(MATLAB_ROOT)/extern/include/" -fexceptions \
		-fbackslash -fPIC -fno-omit-frame-pointer -o $@ -O $<

trmobjfunc.mexmaci64: trmobjfunc.o prompt.o
	$(FC) $(FCFLAGS) -pthread -shared -O $^ \
		-Wl,-rpath,$(MATLAB_ROOT)/bin/maci64 \
		-L"$(MATLAB_ROOT)/bin/maci64" -lmx -lmex -lmat -lm \
        -L"/usr/local/lib/gcc/6/" \
		-lgfortran -o $@

trmoptim.mexmaci64: trmoptim.o prompt.o
	$(FC) $(FCFLAGS) -pthread -shared -O $^ \
		-Wl,-rpath,$(MATLAB_ROOT)/bin/maci64 \
		-L"$(MATLAB_ROOT)/bin/maci64" -lmx -lmex -lmat -lm \
        -L"/usr/local/lib/gcc/6/" \
		-lgfortran -o $@

.PHONY: clean

clean:
	rm *.o *.mod *.mexmaci64

