# Makefile
####################################################################################################
# Autómatikus változók:												\
$^ az összes feltétel listája belértve a mappákat is				\
$< az első feltétel													\
$@ a cél															\
$? változtatott fájlok												\
$* ??????															\
#####################################################################
# értékadás															\
:= értékadás														\
?= csak ha nincs definiálva											\
#####################################################################

include system.mk

DER=	# ha 'make DER=ok' ként hivom meg, akkor DER-nek "ok" lesz az értéke!!!!!

CC := gcc
CFLAGS := -std=gnu99 -O3 -ggdb3
include config.mk

errorExtraFlags := -Wshadow -Winit-self -Wunsafe-loop-optimizations -Wbad-function-cast
errorExtraFlags += -Wcast-qual -Wcast-align -Wwrite-strings -Wlogical-op -Waggregate-return
errorExtraFlags += -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wvla -Wdisabled-optimization
errorFlags := -Wall -Wextra -Wformat-nonliteral -Wformat-security -Wmissing-include-dirs
errorFlags += -Wswitch-default -Wswitch-enum -Wconversion -Wstrict-prototypes -Wold-style-definition $(errorExtraFlags)
CFLAGS += -march=$(shell arch) $(errorFlags)

srcdir := src
incdir := include
objdir := object_dir
objects := $(addprefix $(objdir)/,$(patsubst %.c,%.o,$(wildcard $(srcdir)/*.c)))
testdir := test

vpath %.c $(srcdir):$(testdir)
vpath %.h $(incdir)
# EZT MÉG LE KELL ELLENŐRIZNI
vpath lib%.so $(subst -L,,$(subst lib\ -L,lib:,$(shell pkg-config --libs-only-L lalinspiral)))
vpath lib%.a $(subst -L,,$(subst lib\ -L,lib:,$(shell pkg-config --libs-only-L lalinspiral)))

LAL_INC := $(shell pkg-config --cflags lalinspiral)
LAL_LIB := $(shell pkg-config --libs lalinspiral)

proba :
	$(info a fordító: $(CC))
	@echo "$(subst -L,,$(subst lib\ -L,lib:,$(shell pkg-config --libs-only-L lalinspiral)))"
	@echo "$(CC)"
	@echo "$(CFLAGS)"
	@echo "$(LAL_INC)"
	@echo "$(LAL_LIB)"
	@echo "$(objdir)"
	@echo "$(objects)"

all : main

release debug :
	@echo "$(CFLAGS)"

test : main

main : $(objects)
	$(CC) -o main $(objects)

$(objdir)/%.o : %.c
	@echo -e 'Building file: $<'
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -MMD -MF$(@:%.o=%.d) -MT$(@:%.o=%.d) $< -o $@	# implicit szabál objektumok létrehozására
	@echo -e 'Finished building: $<'

# kisegítők
$(objects): | $(objdir)		# az objektumok előfeltétele a mappájuk, de nem kell mindig frissíteni

$(objdir) :					# objektumok mappája
	mkdir $(objdir)

print : %.c					# kilistázza a változtatott fájlokat az utolsó print óta
	lpr -p $?
	touch $@

debug : 

release : CFLAGS += -O3

# parancsok
.PHONY : clean cleanobj cleanall all	# csak utasítás név, nem cél

clean : cleanobj

cleanall : cleanobj
	-rm main

cleanobj :
	-rm -f $(objects)

####################################################################################################
# régi

test: main_Test.c generator.o util_math.o detector.o match.o match_Multi.o variables.h match_Multi.h
	${CC} -o test main_Test.c generator.o util_math.o detector.o match.o match_Multi.o ${CFLAGS} ${LAL_INC} ${LAL_LIB}

testRun: test
	${RUN} ./test ${PAR}

OBJ=LALSQTPNWaveformInterface.o LALSQTPNWaveform.o LALSQTPNIntegrator.o
NEWDEPS=main_Match.c match.o
all: new

gen: main_Generator.c generator.o util_math.o detector.o match.o match_Multi.o variables.h match_Multi.h
	${CC} -o gen main_Generator.c generator.o util_math.o detector.o match.o match_Multi.o ${CFLAGS} ${LAL_INC} ${LAL_LIB}

#main: main_Old_Match.c match.c detector.c
#	${CC} -o main main_Old_Match.c match.c detector.c ${CFLAGS} ${LAL_INC} ${LAL_LIB} -lm

mainRun: main
	./main `head -n 1 input.data` twoPointFivePN ALL

new: $(NEWDEPS)
	${CC} ${CFLAGS} ${LAL_INC} ${LAL_LIB} -o new $(NEWDEPS)

%.o: %.c %.h
	${CC} -std=gnu99 -c ${CFLAGS} ${LAL_INC} $<

newRun: new
	${RUN} ./new $(shell head -n 1 input.data)

matchRun: match
	./match $(shell head -n 1 input.data)

match: main_Spin_Match.c match.h match.c
	${CC} ${CFLAGS} ${LAL_INC} ${LAL_LIB} -o match main_Spin_Match.c match.c


lal: LALSTPNWaveformTestMod.c
	${CC} ${CFLAGS} ${LAL_INC} ${LAL_LIB} -o lal LALSTPNWaveformTestMod.c -lm
	@echo ''
 
sqt: LALSTPNWaveformTest.c
	${CC} ${CFLAGS} ${LAL_INC} ${LAL_LIB} -std=gnu99 -o sqt LALSTPNWaveformTest.c -lm
	@echo ''

clean:
	rm -rf *.o *.out *.b
	@echo ''

cleanrun:
	rm -rf lal sqt overlap lal.out sqt.out main match new gen
	@echo ''

#cleanall:
#	make clean
#	make cleanrun
#	@echo ''

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
