#Makefile

INCL_LAL= $(env pkg-config --cflags lal)
LIBS_LAL= $(env pkg-config --libs lal)

all: match
.PHONY: all

# lal: LALSTPNWaveformTestMod.c
#	gcc -o lal LALSTPNWaveformTestMod.c $(INCL_LAL) $(LIBS_LAL)

match: main.o match.o ioHandling.o generator.o detector.o util.o
	gcc -o match main.o match.o ioHandling.o generator.o detector.o util.o -lm -lfftw3 -Wall

main.o: main.c match.o generator.o ioHandling.o util.c util.h
	gcc -c main.c match.c generator.c ioHandling.c util.c -Wall

util.o: util.c util.h
	gcc -c util.c util.h -Wall

detector.o: detector.c detector.h util.c util.h
	gcc -c detector.c util.c -Wall

generator.o: generator.c generator.h util.c util.h
	gcc -c generator.c util.c -Wall

ioHandling.o: ioHandling.c ioHandling.h generator.o util.c util.h
	gcc -c ioHandling.c generator.c util.c -Wall

match.o: match.c match.h generator.o detector.o ioHandling.o
	gcc -c match.c generator.c detector.c ioHandling.c -Wall

clean_all: clean_gcc clean_run

clean_gcc: clean
	rm -f match lal

clean: 
	rm -f main.o match.o ioHandling.o generator.o detector.o util.o

clean_run:
	rm -f own.match parameters.data to_Plot.data gen.out [0-9]* xparameters.data xto_Plot.data x[0-9]*

run: all
	time ./match data.init

.PHONY: clean clean_gcc clean_all clean_run run

#help:
#	echo all: make match and lal
#	echo match: make match
#	echo lal: make lal
#	echo clean_run: to clean after run
#	echo clean: to clean object files
#	echo clean_gcc: like clean and also delete match and lal
#	echo clean_all: all previus clean options
#	echo run: starts the program with data.init and watches the time
#.PHONY: help
