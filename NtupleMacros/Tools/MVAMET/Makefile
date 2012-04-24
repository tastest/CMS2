# All targets with # symbol are self-documenting, i.e. make help or simply make will
# show the targets among available options
#
# User targets are at the bottom
#
ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else

all: 	help
help:
	@echo "Available Parameters:";\
	echo -e "\tVERSION - skim name for processing data.";\
	echo -e "\tDataDir - alternative location for samples. Default is data/";\
	echo -e "\tSkimDir - alternative location for output skims. Default is data/";\
	echo -e "\tSkimName - name of the skim to produce when you skim data";\
	echo -e "\tMerge - merge output skims during skimming by datasets/version";\
	echo -e "\tSync - show detail information about cut application for synchronization needs";\
	echo "Available Targets:";\
	cat Makefile | perl -ne 'printf("\t%-15s %s\n",$$1,$$2) if(/^(\S+):[^#]+(#.*)$$/)';\
	echo "Example:";\
	echo -e "\tmake DataDir=skims/ VERSION=wwskim table";\
	echo -e "\tmake DataDir=data/ SkimDir=skims/ SkimName=ww202020 skim"
ifndef VERBOSE
  QUIET := @
endif

CC = g++
CMSROOT = ../../
ROOFITINCLUDE = 
# ifdef CMSSW_VERSION
# 	ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
# endif
INCLUDE = -I$(CMSROOT) $(ROOFITINCLUDE) -I$(CMSROOT)/CORE
CFLAGS = -Wall -Wno-unused-function -g -O2 -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS)

DICTINCLUDE = $(ROOTSYS)/include/Math/QuantFuncMathCore.h $(ROOTSYS)/include/TLorentzVector.h $(ROOTSYS)/include/Math/Vector4D.h GBRForest.h GBRTree.h

LIBFLAGS =  -Wl,-rpath,.

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags)
ifdef CMSSW_RELEASE_BASE
	LINKERLIBS = -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH)/ -L$(CMSSW_RELEASE_BASE)/external/$(SCRAM_ARCH)/lib $(shell root-config --libs) -lEG -lGenVector -lTMVA
else
	LINKERLIBS = $(shell root-config --libs) -lEG -lGenVector -lTMVA
endif

ifeq ($(shell root-config --platform),macosx)
	MACOSXFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace
endif

#
SOURCES = $(wildcard *.cc) 
OBJECTS = $(SOURCES:.cc=.o) LinkDef_out.o
LIB = libCMS2NtupleMacrosMVAMet.so

FWLIB = ../MiniFWLite/libMiniFWLite.so

.PHONY: all help compile clean cms2env

libs:	$(LIB)

$(LIB):	$(SOURCES) $(OBJECTS) 
	$(QUIET) echo "Linking $(LIB)"; \
	$(LINKER) $(LINKERFLAGS) $(MACOSXFLAGS) $(LINKERLIBS) -shared $(OBJECTS) -o $@ 2>&1|perl -ne 'print if(!/skipping incompatible/)'

$(FWLIB):
	$(QUIET) echo "making MiniFWLite"; \
	cd ../MiniFWLite; \
	$(MAKE) -f Makefile;

LinkDef_out.cxx: PHMetAnalysisLinkDef.h 
	$(QUIET) echo "Making CINT dictionaries"; \
	rootcint -f LinkDef_out.cc -c -p $(DICTINCLUDE) PHMetAnalysisLinkDef.h; \
	cat PHMetAnalysisLinkDef.h LinkDef_out.cc > LinkDef_out.cxx; rm LinkDef_out.cc


# General rule for making object files
%.d:	%.cc
	$(QUIET) echo "Checking dependencies for $<"; \
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@

%.o: 	%.cc 
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@ -DPFISOFROMNTUPLE

%.o: 	%.cc %.h
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@ -DPFISOFROMNTUPLE

%.o: 	%.cxx 
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

data:
	@echo 'Make a directory or a symbolic link called "data" pointing to sample location'; echo; exit 1

#################################################################
#                       User Targets 
#################################################################
build: $(LIB) # compile code

b: build

clean: # clean up 
	$(QUIET) rm -v -f \
	*.o *.d *.so *.exe LinkDef_out*; echo "Done"

endif
