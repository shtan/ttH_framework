# System type
UNAME := $(shell uname)

LIBDIR = lib

CXX = $(shell root-config --cxx)
CXXFLAGS = $(shell root-config --cflags) -std=c++11 -Wall -Wextra -pedantic -O2 \
-Wshadow -Wunused-variable -Werror=sign-compare \
-Wunused-but-set-variable -Werror=return-type -Werror=missing-braces \
-Werror=delete-non-virtual-dtor  -fPIC \
$(INCLUDE_RULES)
# needs to be added: -Werror=maybe-uninitialized
INCLUDE_RULES = -I $(shell root-config --incdir) -I . \
-I include/ -I include/ttHReco/ttbar_reco_C/
LD = $(shell root-config --ld)
LDFLAGS = $(shell root-config --ldflags)
LDLIBS =  $(shell root-config --glibs) -lMinuit2 -lMathMore -l GenVector

OBJECTS = src/WDaughterEllipseCalculator.o src/lightJetChiSquareMinimumSolver.o \
	  src/topEventMinimizer.o src/topSystemChiSquare.o \
	  src/ttbar_reco_C.o

LINK = $(LD) $(LDFLAGS)
LINKEND = $(OBJECTS) $(LDLIBS)




###############
### Targets ###
###############

libname = ttbar_reco_C
targetlib = $(LIBDIR)/lib$(libname).so

all: $(targetlib)




############################
### Linux building block ###
############################

ifeq ($(UNAME), Linux)
# do something Linux-y

$(targetlib): $(OBJECTS)
	$(LINK) -fPIC -shared -Wl,-soname,lib$(libname).so \
-o $(targetlib) $(LINKEND)
endif




#############################
### Mac OS building block ###
#############################

ifeq ($(UNAME), $(filter $(UNAME),OSX Darwin))
# do something OSX-y

$(targetlib): $(OBJECTS)
	$(LINK) -shared -Wl,-install_name,lib$(libname).so \
-o $(targetlib) $(LINKEND)
endif




clean:
	rm -f $(targetlib); \
	rm -f $(OBJECTS)

.PHONY: clean

.SUFFIXES: .cpp .cc .cxx .c

.cxx.o:
	$(CXX) -c $(CXXFLAGS) -o $@ $<
