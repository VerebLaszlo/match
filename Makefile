# Makefile
#############################################################################
#		LALInspiral instalálása												#
#	1. instalálld fel a lal csomagot követve a képernyő és a				#
#	https://www.lsc-group.phys.uwm.edu/daswg/docs/howto/lal-install.html	#
#	oldalon lévő utasításokat.												#
#	2. tölts le és instalálld fel a metaio csomagot
#	3. instalálld fel a lalmetaio csomagot hasonlóan						#
#	4. végül a lalinspiral csomagot											#
#	5. ha nem találja a könyvtárakat, add hozzá az útvonalakat, az			#
#	/etc/ld.so.conf fájlhoz, és hajtsd végra a ldconfig parancsot.			#
#############################################################################

##
#	Módosított fájlok
#	LALSQTPN*
#	GenerateInspiral.*
#	LALInspiralWave.c
#	LALInspiral.h
#	LALRandomInspiralSignal.c
#	LALInspiralWaveOverlap.c
#	LIGOMetadataTables.h
#	GenerateInspiral.c
#	GenerateInspiral.h
###

RENORM=0 # 0,1
BUILD_TYPE=prod# debug, normal, prod
CFLAGS=-std=gnu99
ifeq (${BUILD_TYPE},debug)
	CFLAGS+=-Wall -W -g3
	DEBUG=1
else
ifeq (${BUILD_TYPE},prod)
	CFLAGS+=-O3
	DEBUG=0
else
	CFLAGS+=-Wall -W -g3
	DEBUG=0
endif
endif

CC=colorgcc -c
LAL_INC=$(shell pkg-config --cflags lalinspiral)
LAL_LIB=$(shell pkg-config --libs lalinspiral)
OBJ=LALSQTPNWaveformInterface.o LALSQTPNWaveform.o LALSQTPNIntegrator.o

all: new

main: main.c match.c detector.c
	gcc -o main main.c match.c detector.c ${CFLAGS} ${LAL_INC} ${LAL_LIB} -lm

mainRun: main
	./main `head -n 1 input.data`

new: new_match.c match.h match.c
	colorgcc ${CFLAGS} ${LAL_INC} ${LAL_LIB} -o new new_match.c match.c

match: sqt_match.c match.h match.c
	colorgcc ${CFLAGS} ${LAL_INC} ${LAL_LIB} -o match sqt_match.c match.c


lal: LALSTPNWaveformTestMod.c
	colorgcc ${CFLAGS} ${LAL_INC} ${LAL_LIB} -o lal LALSTPNWaveformTestMod.c -lm
	@echo ''
 
own: LALSQTPNWaveformTest.c
	colorgcc ${CFLAGS} ${LAL_INC} ${LAL_LIB} -o own LALSQTPNWaveformTest.c -lm
	@echo ''

overlap: RandomInspiralSignalTest.c
	colorgcc ${CFLAGS} ${LAL_INC} ${LAL_LIB} -o overlap RandomInspiralSignalTest.c

clean:
	rm -rf *.o *.out *.b
	@echo ''

cleanrun:
	rm -rf lal own overlap lal.out own.out
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
