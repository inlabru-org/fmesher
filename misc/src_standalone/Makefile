SHELL = /bin/bash
MAKE = make

# USER OPTIONS

# the full path to where to install the program/library, headerfiles and docs.
PREFIX = /usr/local

# select compilers and optimized compiler-options, add -fopenmp or similar to use OpenMP
CC = gcc
CXX = g++ # -stdlib=libc++
FC = gfortran
FCEXTRAFLAGS = -fno-second-underscore
# FLAGS= -mfpmath=sse -msse2
FLAGS= -O3 -funroll-loops -pipe -Wall -Wextra \
  -Wsign-compare -Wundef -Wimplicit-fallthrough -fdiagnostics-color
LD = $(CXX)
INSTALL = install
RSYNC = rsync
SED = sed
GENGETOPT = gengetopt

# The path to the external libraries: metis, taucs, lapack, blas and zlib, often the same as $(PREFIX)
LEXTPREFIX = $(PREFIX)

# The external libraries to link with
#EXTLIBS1 = -L$(PREFIX)/lib -L$(LEXTPREFIX)/lib -lgsl -lgslcblas
EXTLIBS1 = -L$(PREFIX)/lib -L$(LEXTPREFIX)/lib
EXTLIBS2 =
EXTLIBS3 =

EXTLIBS = $(EXTLIBS1) $(EXTLIBS2) $(EXTLIBS3)



OS_TYPE:=$(shell uname -s)
ifeq ($(OS_TYPE), Linux)
  OS_FLAGS=-DLINUX
# Only activate X if specifically asked for:
ifeq ($(FMESHER_WITH_X), 1)
  OS_LDFLAGS=-lX11
  OS_FLAGS:=$(OS_FLAGS) -DFMESHER_WITH_X
  ADD_X_O = xtmpl.o x11utils.o
else
  OS_LDFLAGS=
  ADD_X_O =
endif
endif
# What does uname -s say on Mac/OSX: Darwin.
ifeq ($(OS_TYPE), Darwin)
  OS_FLAGS=-DMAC
# Activate X if specifically asked for:
ifeq ($(FMESHER_WITH_X), 1)
  OS_LDFLAGS=-L/usr/X11R6/lib -L/usr/X11/lib -lX11
  OS_FLAGS:=$(OS_FLAGS) -DFMESHER_WITH_X
  ADD_X_O = xtmpl.o x11utils.o
else
  OS_LDFLAGS=
  ADD_X_O =
endif
endif
ifeq ($(OS_TYPE), Windows)
  OS_FLAGS=-DWINDOWS
  OS_LDFLAGS=
  ADD_X_O =
endif


# Activate GSL, only needed for spherical harmonics, which can be computed with the R gsl package instead
ifeq ($(FMESHER_WITH_GSL), 1)
  OS_LDFLAGS:=$(OS_LDFLAGS) -lgsl -lgslcblas
  OS_FLAGS:=$(OS_FLAGS) -DFMESHER_WITH_GSL -DFMESHER_WITH_SPHERICAL_HARMONICS
endif

################################################################################
#
# no changes usually needed below this point
#
################################################################################

RM    = rm -f
RMR   = rm -rf
RMDIR = rmdir

GITCOMMIT:= $(shell git log --pretty=format:"%H - %cd" | head -n1)
INCL    = -I$(PREFIX)/include
##LDFLAGS = $(FLAGS) $(OS_LDFLAGS) -lgsl -lgslcblas
LDFLAGSXX = $(FLAGS) $(OS_LDFLAGS)
FFLAGS  = $(FLAGS) $(FCEXTRAFLAGS) $(INCL)
CFLAGS  = -std=gnu99 $(FLAGS) $(OS_FLAGS) $(INCL) -I$(LEXTPREFIX)/include -DCMDLINE_VERSION="\"${GITCOMMIT}\""
CXXFLAGS  = -std=c++17 $(FLAGS) $(OS_FLAGS) $(INCL) -I$(LEXTPREFIX)/include -DCMDLINE_VERSION="\"${GITCOMMIT}\""

FMESHER = fmesher
FMESHER_TEST = fmesher_test
FMESHER_TEST2 = treetest

OBJ = $(ADD_X_O) predicates.o mesh.o meshc.o vector.o ioutils.o cmdline.o \
  basis.o locator.o trees.o fmesher_helpers.o
CMDLINE_CH = cmdline.h cmdline.c
all: $(FMESHER)


# !!! Do not change! /FL
DEPDIR = .deps/.dummy
$(DEPDIR):
	-mkdir $(dir $(DEPDIR))
	touch $(DEPDIR)

%.o: %.c $(DEPDIR)
	$(COMPILE.c) -MMD -MF $(dir $(DEPDIR))$(notdir $*.d.tmp) -o $@ $<
	-@$(SED) -e ':more;/\\/N;s/[ ]*\\\r*\n[ ]*/ /;tmore' \
	      -e 's,\($*\.o\)[ :]*\(.*\),$@ : $$\(wildcard \2\),g' \
	  < $(dir $(DEPDIR))$(notdir $*.d.tmp) > $(dir $(DEPDIR))$(notdir $*.d)
	-$(RM) $(dir $(DEPDIR))$(notdir $*.d.tmp)

%.o: %.cc $(DEPDIR)
	$(COMPILE.C) -MMD -MF $(dir $(DEPDIR))$(notdir $*.d.tmp) -o $@ $<
	-@$(SED) -e ':more;/\\/N;s/[ ]*\\\r*\n[ ]*/ /;tmore' \
	      -e 's,\($*\.o\)[ :]*\(.*\),$@ : $$\(wildcard \2\),g' \
	  < $(dir $(DEPDIR))$(notdir $*.d.tmp) > $(dir $(DEPDIR))$(notdir $*.d)
	-$(RM) $(dir $(DEPDIR))$(notdir $*.d.tmp)



$(FMESHER) : $(OBJ) fmesher.o
	$(LD) -o $@ $^ $(LDFLAGSXX) -L$(PREFIX)/lib  $(EXTLIBS)

$(FMESHER_TEST) : $(OBJ) fmesher_test.o
	$(LD) -o $@ $^ $(LDFLAGSXX) -L$(PREFIX)/lib  $(EXTLIBS)

$(FMESHER_TEST2) : $(OBJ) treetest.o
	$(LD) -o $@ $^ $(LDFLAGSXX) -L$(PREFIX)/lib  $(EXTLIBS)

cmdline.o: $(CMDLINE_CH)
$(OBJ) : $(CMDLINE_CH)
$(CMDLINE_CH) : cmdline.ggo
	$(GENGETOPT) --file-name=cmdline --conf-parser --unamed-opts=PREFIX \
		--long-help --func-name=cmdline < cmdline.ggo

tags :; etags --output=TAGS *.c *.h *.cc

clean:; -$(RM) $(OBJ) $(FMESHER).o $(FMESHER_TEST).o $(FMESHER_TEST2).o $(FMESHER) $(FMESHER_TEST) $(FMESHER_TEST2) $(CMDLINE_CH) Rcpp*.o *.so *.o
clean-deps:
	-$(RMR) $(dir $(DEPDIR))

depend: $(dir $(DEPDIR))dependencies.d
$(dir $(DEPDIR))dependencies.d: Makefile $(DEPDIR) $(wildcard *.c) $(wildcard *.cc) $(wildcard *.h)
	-$(RM) $(dir $(DEPDIR))dependencies.d
	touch $(dir $(DEPDIR))dependencies.d
	for file in *.c *.h *.cc ; do \
	  if [ -n "`grep CMDLINE_VERSION $$file`" ]; then \
	    echo $$file | \
	      $(SED) 's/\(.*\)\(\.[hc]\)$$/\1\.o: gitstamp/' \
	      >> $(dir $(DEPDIR))dependencies.d ; \
	  fi ; \
	done

gitstamp: dummytarget
	[ -f $@ ] || touch $@
	echo '$(GITCOMMIT)' | cmp -s $@ - || echo '$(GITCOMMIT)' > $@

dummytarget: ;

.PHONY: depend clean clean-deps uninstall install tags dummytarget

include $(dir $(DEPDIR))dependencies.d $(wildcard $(dir $(DEPDIR))*.d)
