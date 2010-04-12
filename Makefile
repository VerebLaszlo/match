#Makefile

UTIL = util.c util.h

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

clean: main.o match.o ioHandling.o generator.o detector.o util.o
	rm match main.o match.o ioHandling.o generator.o detector.o util.o
.PHONY: clean

clean_other: own.match parameters.data to_Plot.data
	rm own.match parameters.data to_Plot.data
.PHONY: clear
