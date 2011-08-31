# Makefile
####################################################################################################
# Autómatikus változók:												\
$^ az összes feltétel listája belértve a mappákat is				\
$< az első feltétel													\
$@ a cél															\
$? változtatott fájlok												\
$* ??????															\
#####################################################################
# ?= csak ha nincs definiálva										\
#####################################################################

include system.mk

DER=	# ha 'make DER=ok' ként hivom meg, akkor DER-nek "ok" lesz az értéke!!!!!

CC := gcc
CFLAGS := -std=gnu99 -O3 -ggdb3
include config.mk

errorFlags := -Wall -Wextra -Wformat-security -Wmissing-include-dirs -Wswitch-default
errorFlags += -Wstrict-prototypes -Wold-style-definition -Wno-aggregate-return

errorExtraFlags := -Wshadow -Winit-self -Wunsafe-loop-optimizations -Wcast-qual
errorExtraFlags += -Wcast-align -Wwrite-strings -Wlogical-op -Waggregate-return
errorExtraFlags += -Wmissing-prototypes -Wmissing-declarations -Wdisabled-optimization
errorExtraFlags += -Wconversion -Wswitch-enum

errorOptionalFlags := -Wformat-nonliteral -Wconversion -Wswitch-enum -Wbad-function-cast
errorOptionalFlags += -Wredundant-decls

CFLAGS += -march=$(shell arch) $(errorFlags)
srcdir := src
incdir := include
objdir := object_dir
objects := $(patsubst $(srcdir)/%.c,$(objdir)/%.o,$(wildcard $(srcdir)/*.c))
makes := $(patsubst $(srcdir)/%.c,$(objdir)/%.d,$(wildcard $(srcdir)/*.c))
testdir := test
objects += $(patsubst $(testdir)/%.c,$(objdir)/%.o,$(wildcard $(testdir)/*.c))
makes += $(patsubst $(testdir)/%.c,$(objdir)/%.d,$(wildcard $(testdir)/*.c))
includes := -I$(incdir) #-I/usr/include

objs_test := main_test.o signals.o detector.o binary_system.o binary_system_mass.o
objs_test += binary_system_spin.o util_math.o util_IO.o util.o test.o parameters.o lal_wrapper.o

vpath
vpath %.c $(srcdir)
vpath %.c $(testdir)
vpath %.h $(incdir)
vpath %.o $(objdir)
vpath %.d $(objdir)
# EZT MÉG LE KELL ELLENŐRIZNI
#vpath lib%.so $(subst -L,,$(subst lib\ -L,lib:,$(shell pkg-config --libs-only-L lalinspiral)))
#vpath lib%.a $(subst -L,,$(subst lib\ -L,lib:,$(shell pkg-config --libs-only-L lalinspiral)))

lal_includes := $(shell pkg-config --cflags lalinspiral)
includes += $(lal_includes)
lal_libraries := $(shell pkg-config --libs-only-l lalinspiral)
lal_libraries_path := $(shell pkg-config --libs-only-L lalinspiral)

all : test

test main : CFLAGS += $(errorExtraFlags) $(lal_libraries_path)
test : macros += -DTEST

test main : $(objects) -lfftw3 -lm
	@echo -e '\e[36mLinking: $@\e[0m'
	$(hide_echo)$(CC) $(CFLAGS) $(macros) $(lal_libraries) -o $@ $^
	@echo -e '\e[35mFinished linking: $@\e[0m'
	@echo ' '

$(objdir)/%.o : %.h

$(objdir)/%.o : %.c | $(objdir)
	@echo -e '\e[36mBuilding file: $<\e[0m'
	$(hide_echo)$(CC) $(CFLAGS) $(CPPFLAGS) $(includes) $(macros) -c -MMD -MF$(@:%.o=%.d) -MT$(@:%.o=%.d) $< -o $@
	@echo -e '\e[35mFinished building: $<\e[0m'
	@echo ' '

# kisegítők
#$(objdir)/%.o: | $(objdir)		# NEM MŰKÖDIK!!!!!!!!!!!!!!!!!!!!!!!!!
#	@echo 'obj'

#$(objdir)/%.d: | $(objdir)		# NEM MŰKÖDIK!!!!!!!!!!!!!!!!!!!!!!!!!
#	@echo 'makes'

$(objdir) :
	mkdir $(objdir)

print : %.c					# kilistázza a változtatott fájlokat az utolsó print óta
	lpr -p $?
	touch $@

debug : 

release : CFLAGS += -O3

# parancsok
.PHONY : doxy cleandoxy clean cleanobj cleanall all # csak utasítás név, nem cél

doxy : doxyutil
	doxygen Doxyfile$(TEST)

doxyutil :
	-mkdir doc
	-mkdir doc/html
	-mkdir doc/latex
	-rm src/*.h
	cp doc/style.sty doc/html/style.sty
	cp doc/style.sty doc/latex/style.sty
	clear scr

cleandoxy :
	-rm -R doc/*

clean : cleanobj

cleanall : cleanobj
	-rm $(objdir)/*.d
	-rm main test
	clear

cleanobj :
	-rm $(objdir)/*.o
	clear
