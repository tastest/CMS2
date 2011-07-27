# All targets with # symbol are self-documenting, i.e. make help or simply make will
# show the targets among available options
#
# User targets are at the bottom
#

#
ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else

CC = g++
CMSROOT = ../
ROOFITINCLUDE = 
ifdef CMSSW_VERSION
	ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
endif
INCLUDE = -I$(CMSROOT) $(ROOFITINCLUDE) -I$(CMSROOT)/CORE
CFLAGS = -Wall -Wno-unused-function -g -O2 -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS)

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags)

ifeq ($(shell root-config --platform),macosx)
ifdef CMSSW_RELEASE_BASE
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH)/ -L$(CMSSW_RELEASE_BASE)/external/$(SCRAM_ARCH)/lib $(shell root-config --libs) -lEG -lGenVector
else
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
endif
endif


CORESOURCES = $(CMSROOT)/CORE/CMS2.cc $(CMSROOT)/CORE/muonSelections.cc $(CMSROOT)/CORE/electronSelections.cc $(CMSROOT)/CORE/electronSelectionsParameters.cc $(CMSROOT)/CORE/jetSelections.cc $(CMSROOT)/CORE/MITConversionUtilities.cc $(CMSROOT)/CORE/trackSelections.cc $(CMSROOT)/CORE/eventSelections.cc $(CMSROOT)/CORE/metSelections.cc $(CMSROOT)/CORE/mcSelections.cc  
COREOBJECTS = $(CORESOURCES:.cc=.o) 
CORELIB = libCMS2NtupleMacrosCORE.so

#SOURCES = $(wildcard *.cc) 
SOURCES = doAnalysis.cc 
OBJECTS = $(SOURCES:.cc=.o) LinkDef_out.o
LIB = libCMS2NtupleMacrosLooper.so

FWLIB = ../Tools/MiniFWLite/libMiniFWLite.so

LIBS = $(CORELIB) $(LIB)

.PHONY: all help compile clean cms2env

libs:	$(LIBS)

$(CORELIB):	$(CORESOURCES) $(COREOBJECTS) 
	$(QUIET) echo "Linking $(CORELIB)"; \
	$(LINKER) $(LINKERFLAGS) $(MACOSXFLAGS) $(LINKERLIBS) -shared $(COREOBJECTS) -o $@ 2>&1|perl -ne 'print if(!/skipping incompatible/)'

$(LIB):	$(OBJECTS) 
	$(QUIET) echo "Linking $(LIB)"; \
	$(LINKER) $(LINKERFLAGS) $(MACOSXFLAGS) $(LINKERLIBS) -shared $(OBJECTS) -o $@ 2>&1|perl -ne 'print if(!/skipping incompatible/)'

$(FWLIB):
	$(QUIET) echo "making MiniFWLite"; \
	cd ../Tools/MiniFWLite; \
	$(MAKE) -f Makefile;

LinkDef_out.cxx: LinkDef.h
	$(QUIET) echo "Making CINT dictionaries"; \
	rootcint -f LinkDef_out.cc -c -p $(INCLUDE) LinkDef.h; \
	cat LinkDef.h LinkDef_out.cc > LinkDef_out.cxx; rm LinkDef_out.cc

# General rule for making object files
%.d:	%.cc
	$(QUIET) echo "Checking dependencies for $<"; \
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@
%.o: 	%.cc 
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.cxx 
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

Util_HWW.root: doAllHWW.C $(CORELIB) $(LIB) $(FWLIB)
	$(QUIET) root -b -q doAllHWW.C

#################################################################
#                       User Targets 
#################################################################

build: $(LIBS) # compile code
	@
b: build

clean:   
	$(QUIET) rm -v -f \
	$(CMSROOT)/CORE/*.o \
	$(CMSROOT)/CORE/*.d \
	$(CMSROOT)/CORE/*.so \
	$(CMSROOT)/Tools/*.so \
	$(CORELIB) $(LOOPERLIB) $(FWLIB) \
	*.so *.o *.d LinkDef_out*; echo "Done"

ww:  Util_HWW.root
zz:  Util_HZZ.root


endif
