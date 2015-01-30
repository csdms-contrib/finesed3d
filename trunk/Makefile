#### Makefile for channel code #### 
EXEC = FineSed3D.x
TGZFILE = all.tgz
#
SRC = \
advance_compact.F  io.F  pstep.F  timers.F  viscxyz_mkl_compact.F  \
courant.F  main.F rhs.F compact.F initial_compact.F   \
divg.F partial_mkl_compact.F stats.F tt_rhs.F post_proc.F nltrms_up.F \
enlred.F fft2d_new_mkl.F mat_deriv.F idmax.F Helmholtz_compact.F
#
INC = common.inc fft.inc flags.inc global.inc timers.inc
#
OBJS = $(SRC:.F=.o)
#
#LIBS       = -lmkl_lapack -lmkl -lguide -lpthread -lscs
MKL_HOME   = /usr/local/intel/11.1/Compiler/11.1/075/mkl
MKLIBS     = -L $(MKL_HOME)/lib/intel64 -L $(MKL_HOME)/lib/em64t -L /usa/zcheng/FFTW/ -lmkl_intel_lp64 -lfftw3xf_intel -lmkl_intel_thread -lmkl_core -lpthread
#LIBS       = -lscs -lguide -lpthread
FLAGS      = -fpp -traceback -D COBALT
DEBUGFLAGS = -g -C -traceback
OPTFLAGS   = -O2
MPFLAGS    = -openmp 
FCMP       = ifort -I $(MKL_HOME)/include -I $(MKL_HOME)/include/fftw
FCSP       = ifort -I $(MKL_HOME)/include -I $(MKL_HOME)/include/fftw
#
  FC = $(FCMP)
  FLAGS += $(MPFLAGS)

%.o: %.F $(INC)
	$(FC) -c $(FLAGS) $<

$(EXEC): $(OBJS)
	$(FC) -o $@ $(FLAGS) $^ $(MKLIBS)

clean:
	-rm -f $(OBJS)
