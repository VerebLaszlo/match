# Makefile
#############################################################################
#############################################################################

CC=gcc
CONFUSE=1
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

NEWDEPS=new_match.c match.o

LAL_INC=$(shell pkg-config --cflags lalinspiral)
LAL_LIB=$(shell pkg-config --libs lalinspiral)
OBJ=LALSQTPNWaveformInterface.o LALSQTPNWaveform.o LALSQTPNIntegrator.o

ifeq (${CONFUSE},1)
    CFLAGS+=-DCONFUSE $(shell pkg-config --cflags --libs libconfuse)
	NEWDEPS+=confuse-parser.o
endif

all: new

main: main.c match.c detector.c
	${CC} -o main main.c match.c detector.c ${CFLAGS} ${LAL_INC} ${LAL_LIB} -lm

mainRun: main
	./main `head -n 1 input.data`

new: $(NEWDEPS)
	${CC} ${CFLAGS} ${LAL_INC} ${LAL_LIB} -o new $(NEWDEPS)

%.o: %.c %.h
	${CC} -c ${CFLAGS} ${LAL_INC} $<

newRun: new
	${RUN} ./new $(shell head -n 1 input.data)

matchRun: match
	./match $(shell head -n 1 input.data)

match: sqt_match.c match.h match.c
	${CC} ${CFLAGS} ${LAL_INC} ${LAL_LIB} -o match sqt_match.c match.c


lal: LALSTPNWaveformTestMod.c
	${CC} ${CFLAGS} ${LAL_INC} ${LAL_LIB} -o lal LALSTPNWaveformTestMod.c -lm
	@echo ''
 
own: LALSQTPNWaveformTest.c
	${CC} ${CFLAGS} ${LAL_INC} ${LAL_LIB} -o own LALSQTPNWaveformTest.c -lm
	@echo ''

overlap: RandomInspiralSignalTest.c
	${CC} ${CFLAGS} ${LAL_INC} ${LAL_LIB} -o overlap RandomInspiralSignalTest.c

clean:
	rm -rf *.o *.out *.b
	@echo ''

cleanrun:
	rm -rf lal own overlap lal.out own.out main match new
	@echo ''

cleanall:
	make clean
	make cleanrun
	@echo ''

run: own lal
	@echo "`head -n 1 input.data` own.out `tail -n 1 input.data`"
	./own `head -n 1 input.data` own.out `tail -n 1 input.data`
	@echo "`head -n 1 input.data` lal.out"
	./lal `head -n 1 input.data` lal.out

val: own lal
	@echo "valgrind `head -n 1 input.data` own.out `tail -n 1 input.data`"
	valgrind ./own `head -n 1 input.data` own.out `tail -n 1 input.data`
	@echo "valgrind `head -n 1 input.data` lal.out"
	valgrind ./lal `head -n 1 input.data` lal.out
help :
	@echo 'all       : makes everything'
	@echo 'lal       : makes just the LALSTPNWaveform.c part'
	@echo 'own       : makes the whole own code'
	@echo 'clean     : deletes the object files'
	@echo 'cleanrun : deletes the exe files'
	@echo 'cleanall : invokes the "clean" and "cleanrun" commands'
	@echo 'run       : runs the two programs'
	@echo 'help      : prints this message'
	@echo ''

.PHONY: all clean cleanall cleanrun run help
