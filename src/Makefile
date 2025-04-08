# -*- Makefile -*-
#
# Makefile for geoflac
#
# Author: Eh Tan <tan2@earth.sinica.edu.tw>
#

## Execute "make" if making production run. Or "make debug=1" for debugging run.
##
## debug = 0: optimized build; 1: debugging build
## omp = 1: enable openmp support
## acc = 1: enable openacc support (requiring PGI compiler)
## coverage = 1: enable code coverage analysis (used for debugging)
## gprof = 1: enable profiling (used for tuning optimization)


########################################################################
## Select default fortran 90 compiler
########################################################################

## GNU version >= 5, otherwise openmp and ieee_is_nan won't work
F90 = gfortran

## Intel, tested with v11.1
#F90 = ifort

## PGI group, tested with 2017, 2019, 2021
#F90 = pgf90

## Sun
#F90 = f90


########################################################################
## any customized flags
########################################################################
CFLAGS = -g
LDFLAGS =


########################################################################
## default command line options
########################################################################
debug = 0
omp = 1
acc = 0
coverage = 0
gprof = 0
nprof = 0
cuda_path = #/path/to/hpc_sdk/Linux_x86_64/25.1/cuda


########################################################################
## Select compiler and linker flags
## (Usually you won't need to modify anything below)
########################################################################

ifeq ($(nprof), 1)
	F90 = nvfortran
	CFLAGS += -cpp -DUSE_NVTX
else
	CFLAGS += -cpp
endif

## for GNU gfortran >= 5
ifeq ($(F90), gfortran)
	ifeq ($(omp), 1)
		CFLAGS += -fopenmp
		LDFLAGS += -fopenmp
	endif
	ifeq ($(debug), 0)
		CFLAGS += -O2
	else
		CFLAGS += -O0 -Wall -fbounds-check -ftrapv -fbacktrace -fdump-core \
			-Waliasing -Walign-commons -Wampersand -Wintrinsic-shadow \
			-Wline-truncation -Wsurprising -Wunderflow \
			-frange-check -Werror -Wno-unused-variable
	endif

	ifeq ($(coverage), 1)
		CFLAGS += -fprofile-arcs -ftest-coverage
		LDFLAGS += -fprofile-arcs -ftest-coverage
	endif

	ifeq ($(gprof), 1)
		CFLAGS += -pg
		LDFLAGS += -pg
	endif
endif

## for Intel v11.1
ifeq ($(F90), ifort)
	CFLAGS += -assume byterecl
	ifeq ($(omp), 1)
		CFLAGS += -fopenmp
		LDFLAGS += -fopenmp
	endif
	ifeq ($(debug), 0)
		CFLAGS += -fast
	else
	       CFLAGS += -O0 -no-opt-subscript-in-range -ftrapuv -check all -traceback
	       LDFLAGS +=
	endif
endif

## for PGI
ifeq ($(F90), pgf90)
	ifeq ($(omp), 1)
		CFLAGS += -mp
		LDFLAGS += -mp
	endif
	ifeq ($(acc), 1)
		CFLAGS += -acc -ta=tesla:cc70,managed -Mcuda -Minfo=accel
		LDFLAGS += -acc -ta=tesla:cc70,managed -Mcuda
	endif
        ifeq ($(debug), 0)
                CFLAGS += -O2 -Minfo=all
        else
               CFLAGS += -O0 -Mbounds -Mchkstk -traceback
               LDFLAGS +=
        endif
endif

## for Solaris
ifeq ($(F90), f90)
	ifeq ($(omp), 1)
		CFLAGS += -xopenmp=parallel
		LDFLAGS += -xopenmp
	endif

# Explanation of flags --
# -r8const: Promote single precision const to double precision.
# -C: Enable runtime subscript range checking
# -fpover=yes: Provide run-time overflow check during READ
# -xcheck=%all: Enable all run-time check
	ifeq ($(debug), 0)
		CFLAGS += -O2 -r8const -xloopinfo -s -xipo
	else
		CFLAGS += -w3 -r8const -C -fpover=yes -xcheck=%all
	endif
endif

ifeq ($(F90), nvfortran)
	ifeq ($(debug), 0)
		CFLAGS += -O
		ifeq ($(omp), 1)
			CFLAGS += -fopenmp
			LDFLAGS += -fopenmp
		endif
	else
		CFLAGS += -O0 -traceback
	endif

	CFLAGS += -L$(cuda_path)/lib64 -lnvToolsExt
	LDFLAGS += -L$(cuda_path)/lib64 -Wl,-rpath,$(cuda_path)/lib64 -lnvToolsExt
endif

########################################################################
## list of source files
########################################################################

# module files have to be compiled before other f90 files
MODULE_SRCS = \
	myrandom_mod.f90 \
	params.f90 \
	arrays.f90 \
	phases.f90 \
	marker_data.f90

ifeq ($(nprof), 1)
	MODULE_SRCS += nvtx_mod.f90
endif

SRCS =	\
	bar2euler.f90 \
	bc_update.f90 \
	change_phase.f90 \
	dt_mass.f90 \
	euler2bar.f90 \
	fl_move.f90 \
	fl_node.f90 \
	fl_rheol.f90 \
	fl_srate.f90 \
	fl_therm.f90 \
	flac.f90 \
	functions.f90 \
	init_areas.f90 \
	init_bc.f90 \
	init_cord.f90 \
	init_marker.f90 \
	init_phase.f90 \
	init_stress.f90 \
	init_temp.f90 \
	init_visc.f90 \
	lpeuler2bar.f90 \
	marker2elem.f90 \
	matprops.f90 \
	outflac.f90 \
	outmarker.f90 \
	par.f90 \
	read_params.f90 \
	rem_cord.f90 \
	remesh.f90 \
	rem_test.f90 \
	rh_elastic.f90 \
	rh_maxwell.f90 \
	rh_plastic.f90 \
	rmasses.f90 \
	rsflac.f90 \
	saveflac.f90 \
	setflac.f90 \
	user_luc.f90

SRCS2 = check_nan.f90

MODULES = $(MODULE_SRCS:.f90=.mod)
MODULE_OBJS = $(MODULE_SRCS:.f90=.o)
OBJS = $(SRCS:.f90=.o)
OBJS2 = $(SRCS2:.f90=.o)

EXE = flac

## Is git in the path?
HAS_GIT := $(shell git --version 2>/dev/null)
ifneq ($(HAS_GIT),)
        ## Is this a git repository?
        IS_REPO := $(shell git rev-parse --git-dir 2>/dev/null)
endif

########################################################################
## actions
########################################################################

.PHONY: all clean clean-data take-snapshot

all: $(EXE) take-snapshot

$(EXE): $(MODULE_OBJS) $(OBJS) $(OBJS2)
	$(F90) $(LDFLAGS) $(MODULE_OBJS) $(OBJS) $(OBJS2) -o $@

$(MODULE_OBJS): %.o : %.f90
	$(F90) $(CFLAGS) -c $<

$(OBJS): %.o : %.f90 $(MODULES)
	$(F90) $(CFLAGS) -c $<

# these subroutine cannot be compiled with high optimization
$(OBJS2): %.o : %.f90 $(MODULES)
	$(F90) -O2 -c $<

take-snapshot:
	@# snapshot of the code for building the executable
	@echo Flags used to compile the code: > snapshot.diff
	@echo '  ' F90=$(F90) debug=$(debug) openmp=$(omp) acc=$(acc) >> snapshot.diff
	@echo '  ' CFLAGS=$(CFLAGS) >> snapshot.diff
	@echo '  ' LDFLAGS=$(LDFLAGS) >> snapshot.diff
	@echo '  ' PATH= "$(PATH)" >> snapshot.diff
	@echo '  ' LD_LIBRARY_PATH=$(LD_LIBRARY_PATH) >> snapshot.diff
	@echo '  ' PWD=$(PWD) >> snapshot.diff
	@echo >> snapshot.diff
	@echo >> snapshot.diff
ifneq ($(HAS_GIT),)
ifneq ($(IS_REPO),)
	@echo '==== Summary of the code ====' >> snapshot.diff
	@git show -s >> snapshot.diff
	@echo >> snapshot.diff
	@git status >> snapshot.diff
	@echo >> snapshot.diff
	@echo >> snapshot.diff
	@echo '== Code modification (not checked-in) ==' >> snapshot.diff
	@git diff >> snapshot.diff
	@echo >> snapshot.diff
	@echo >> snapshot.diff
	@echo '== Code modification (checked-in but not in "origin") ==' >> snapshot.diff
	@git log --patch --no-merges origin..HEAD >> snapshot.diff
else
	@echo "Warning: Not a git repository. Cannot take code snapshot." | tee -a snapshot.diff
	@echo "Warning: Use 'git clone' to copy the code!" | tee -a snapshot.diff
endif
else
	@echo "'git' is not in path, cannot take code snapshot." >> snapshot.diff
endif


clean:
	@rm -f $(MODULES) $(MODULE_OBJS) $(OBJS) $(OBJS2) $(EXE)

clean-data:
	@rm -f *.0 *.rs sys.msg marker_init.asc vbc.s area.dat output.asc temp0.dat
