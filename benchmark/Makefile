
CASE := test-tiny-rift.inp
# CASE := test-tiny-long.inp
#CASE := test-big.inp
FRAME := 4
EXE := ../src/flac
OMP := 4
GPROF := 0
NPROF := 1
FLAG :=
STATS = 1


ifeq (${GPROF}, 1)
  RUN_PROF = gprof ../src/flac gmon.out -bp | head -20
endif

ifeq (${NPROF}, 1)
  	FLAGS = nsys profile -t nvtx,cuda --output=${CASE}.qdrep --force-overwrite true
  	ifeq (${STATS}, 1)
		FLAGS += --stats=true
	endif
endif

# Use single cpu for omp to avoid round-off errors
RUN_TEST = rm -f *.0 sys.msg output.asc; OMP_NUM_THREADS=${OMP} ${FLAGS} ${EXE} ${CASE}

ORIG = orig-${CASE}

all: cmp

set:
	# run test, store test output to old/
	${RUN_TEST}
	${MAKE} store

cmp:
	${RUN_TEST}
	python3 compare.py ${ORIG} ${FRAME}
	${RUN_PROF}

store:
	rm -rf ${ORIG}
	mkdir -p ${ORIG}
	mv *.0 *.rs sys.msg output.asc _contents.save vbc.s ${ORIG}
	cp ../src/snapshot.diff ${ORIG}

.PHONY: all set cmp store
