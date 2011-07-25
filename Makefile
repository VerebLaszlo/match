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
CFLAGS += -march=$(shell arch) #$(errorFlags)

srcdir := src
incdir := include
objdir := object_dir
objects := $(patsubst $(srcdir)/%.c,$(objdir)/%.o,$(wildcard $(srcdir)/*.c))
makes := $(patsubst $(srcdir)/%.c,$(objdir)/%.d,$(wildcard $(srcdir)/*.c))
testdir := test

vpath
vpath %.c $(srcdir)
vpath %.c $(testdir)
vpath %.h $(incdir)
vpath %.o $(objdir)
vpath %.d $(objdir)
# EZT MÉG LE KELL ELLENŐRIZNI
vpath lib%.so $(subst -L,,$(subst lib\ -L,lib:,$(shell pkg-config --libs-only-L lalinspiral)))
vpath lib%.a $(subst -L,,$(subst lib\ -L,lib:,$(shell pkg-config --libs-only-L lalinspiral)))

LAL_INC := $(shell pkg-config --cflags lalinspiral)
LAL_LIB := $(shell pkg-config --libs lalinspiral)

proba :
	@echo "$(objects)"
	@echo ''
	@echo "$(makes)"

all : release

release :
	@echo "$(CFLAGS)"

debug : CFLAGS += $(errorFlags)

debug :
	@echo "$(CFLAGS)"

test : main

main : $(objects)
	$(CC) -o main $(objects)

$(objdir)/%.o : %.c | $(objdir)
	@echo 'Building file: $<'
	ln -s ../$(incdir)/$(subst o,h,$(notdir $@)) $(srcdir)/$(subst o,h,$(notdir $@))
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -MMD -MF$(@:%.o=%.d) -MT$(@:%.o=%.d) $< -o $@
	@echo 'Finished building: $<'

$(objdir)/%.d : %.c | $(objdir)
	@set -e; rm -f $@; \
		$(CC) -M $(CPPFLAGS) $< > $@.$$$$; \
		set 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
		rm -f $@.$$$$
#	@set -e; rm -f $@; \
		$(CC) -M $(CPPFLAGS) $< > $@.$$$$; \
		set 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
		rm -f $@.$$$$

# kisegítők
$(objdir)/%.o: | $(objdir)		# NEM MŰKÖDIK!!!!!!!!!!!!!!!!!!!!!!!!!
	@echo 'obj'

$(objdir)/%.d: | $(objdir)		# NEM MŰKÖDIK!!!!!!!!!!!!!!!!!!!!!!!!!
	@echo 'makes'

$(objdir) :
	mkdir $(objdir)

print : %.c					# kilistázza a változtatott fájlokat az utolsó print óta
	lpr -p $?
	touch $@

debug : 

release : CFLAGS += -O3

# parancsok
.PHONY : doxy cleandoxy clean cleanobj cleanall all	# csak utasítás név, nem cél

doxy :
	-mkdir doc
	-mkdir doc/html
	-mkdir doc/latex
	-rm src/*.h
	cp doc/style.sty doc/html/style.sty
	cp doc/style.sty doc/latex/style.sty
	doxygen

cleandoxy :
	-rm -R doc/*

clean : cleanobj

cleanall : cleanobj
	-rm main

cleanobj :
	-rm $(objects)

