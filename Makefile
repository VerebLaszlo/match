# Makefile
#############################################################################
#############################################################################

CC=gcc
CONFUSE=0
include config.mk
ifdef BUILD
	ifeq (${BUILD}, prod)
		CFLAGS+=-O3
	else
		CFLAGS+=-Wall -Wextra -g3
		ifeq (${BUILD}, dev)
			RUN=valgrind --leak-check=full
		endif	
	endif
else
	CFLAGS+=-O3 -Wall -Wextra -g3 -ggdb3
endif

ifeq (${CONFUSE},1)
    CFLAGS+=-DCONFUSE $(shell pkg-config --cflags --libs libconfuse)
	NEWDEPS+=confuse-parser.o
endif

LAL_INC=$(shell pkg-config --cflags lalinspiral)
LAL_LIB=$(shell pkg-config --libs lalinspiral)

test: main_Test.c generator.o util_math.o detector.o match.o match_Multi.o variables.h match_Multi.h
	${CC} -o test main_Test.c generator.o util_math.o detector.o match.o match_Multi.o ${CFLAGS} ${LAL_INC} ${LAL_LIB}

testRun: test
	${RUN} ./test ${PAR}

OBJ=LALSQTPNWaveformInterface.o LALSQTPNWaveform.o LALSQTPNIntegrator.o
NEWDEPS=main_Match.c match.o
all: new

gen: main_Generator.c generator.o util_math.o detector.o match.o match_Multi.o variables.h match_Multi.h
	${CC} -o gen main_Generator.c generator.o util_math.o detector.o match.o match_Multi.o ${CFLAGS} ${LAL_INC} ${LAL_LIB}

main: main_Old_Match.c match.c detector.c
	${CC} -o main main_Old_Match.c match.c detector.c ${CFLAGS} ${LAL_INC} ${LAL_LIB} -lm

mainRun: main
	./main `head -n 1 input.data` twoPointFivePN ALL

new: $(NEWDEPS)
	${CC} ${CFLAGS} ${LAL_INC} ${LAL_LIB} -o new $(NEWDEPS)

%.o: %.c %.h
	${CC} -c ${CFLAGS} ${LAL_INC} $<

newRun: new
	${RUN} ./new $(shell head -n 1 input.data)

matchRun: match
	./match $(shell head -n 1 input.data)

match: main_Spin_Match.c match.h match.c
	${CC} ${CFLAGS} ${LAL_INC} ${LAL_LIB} -o match main_Spin_Match.c match.c


lal: LALSTPNWaveformTestMod.c
	${CC} ${CFLAGS} ${LAL_INC} ${LAL_LIB} -o lal LALSTPNWaveformTestMod.c -lm
	@echo ''
 
sqt: LALSQTPNWaveformTest.c
	${CC} ${CFLAGS} ${LAL_INC} ${LAL_LIB} -o sqt LALSQTPNWaveformTest.c -lm
	@echo ''

clean:
	rm -rf *.o *.out *.b
	@echo ''

cleanrun:
	rm -rf lal sqt overlap lal.out sqt.out main match new gen
	@echo ''

cleanall:
	make clean
	make cleanrun
	@echo ''

run: sqt lal
	@echo "`head -n 1 input.data`"
	./sqt `head -n 1 input.data`
	@echo "`tail -n 1 input.data`"
	./lal `tail -n 1 input.data`

val: sqt lal
	@echo "valgrind `head -n 1 input.data` sqt.out `tail -n 1 input.data`"
	valgrind ./sqt `head -n 1 input.data` sqt.out `tail -n 1 input.data`
	@echo "valgrind `head -n 1 input.data` lal.out"
	valgrind ./lal `head -n 1 input.data` lal.out
help :
	@echo 'all       : makes everything'
	@echo 'lal       : makes just the LALSTPNWaveform.c part'
	@echo 'sqt       : makes the whole sqt code'
	@echo 'clean     : deletes the object files'
	@echo 'cleanrun : deletes the exe files'
	@echo 'cleanall : invokes the "clean" and "cleanrun" commands'
	@echo 'run       : runs the two programs'
	@echo 'help      : prints this message'
	@echo ''

.PHONY: all clean cleanall cleanrun run help
