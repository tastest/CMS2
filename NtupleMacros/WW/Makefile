# All targets with # symbol are self-documenting, i.e. make help or simply make will
# show the targets among available options
#
# User targets are at the bottom
#
all: 	help
help:
	@echo "Available Targets:";\
	cat Makefile | perl -ne 'printf("\t%-15s %s\n",$$1,$$2) if(/^(\S+):[^#]+(#.*)$$/)'
ifndef VERBOSE
  QUIET := @
endif

CC = g++
CMSROOT = ../
ROOFITINCLUDE = 
ifdef CMSSW_VERSION
	ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
endif
INCLUDE = -I$(CMSROOT) $(ROOFITINCLUDE)
CFLAGS = -Wall -Wno-unused-function -g -O2 -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS)

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags)

ifeq ($(shell root-config --platform),macosx)
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
endif

CORESOURCES = $(wildcard $(CMSROOT)/CORE/*.cc) 
COREOBJECTS = $(CORESOURCES:.cc=.o) 
CORELIB = libCMS2NtupleMacrosCORE.so

SOURCES = $(wildcard *.cc) 
OBJECTS = $(SOURCES:.cc=.o) LinkDef_out.o
LIB = libCMS2NtupleMacrosLooper.so

FWLIB = ../Tools/MiniFWLite/libMiniFWLite.so

LIBS = $(CORELIB) $(LIB)


.PHONY: all help compile clean cms2env rootenv

libs:	$(LIBS)

$(CORELIB):	$(COREOBJECTS) 
	$(QUIET) echo "Linking $(CORELIB)"; \
	$(LINKER) $(LINKERFLAGS) -shared $(COREOBJECTS) -o $@ 2>&1|perl -ne 'print if(!/skipping incompatible/)'

$(LIB):	$(OBJECTS) 
	$(QUIET) echo "Linking $(LIB)"; \
	$(LINKER) $(LINKERFLAGS) -shared $(OBJECTS) -o $@ 2>&1|perl -ne 'print if(!/skipping incompatible/)'

$(FWLIB):
	$(QUIET) echo "making MiniFWLite"; \
	cd ../Tools/MiniFWLite; \
	$(MAKE) -f Makefile;

LinkDef_out.cxx: LinkDef.h
	$(QUIET) echo "Making CINT dictionaries"; \
	rootcint -f LinkDef_out.cc -c -p $(INCLUDE) LinkDef.h; \
	cat LinkDef.h LinkDef_out.cc > LinkDef_out.cxx; rm LinkDef_out.cc
# [ -e LinkDef_out.cc ] || [ -e LinkDef_out.h ] || [ -e LinkDef_out.o ] && rm LinkDef_out.* \

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

rootenv:
ifndef ROOTSYS
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo; exit 1
endif

processed_data.root: processData.C $(CORELIB) $(LIB) $(FWLIB)
ifndef ROOTSYS
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo; exit 1
else
	$(QUIET) root -b -q processData.C
endif

fakeIsoControlSamples.root: makeFakeIsoControlSamples.C
ifndef ROOTSYS
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo; exit 1
else
	$(QUIET) root -b -q makeFakeIsoControlSamples.C+
endif

#################################################################
#                       User Targets 
#################################################################
build: rootenv $(LIBS) # compile code
	@#root -b -q ".L doAnalysis.C++"
b: build

clean: # clean up 
	$(QUIET) rm -f \
	$(CMSROOT)/CORE/*.o \
	$(CMSROOT)/CORE/*.d \
	$(CMSROOT)/CORE/*.so \
	$(CORELIB) $(LOOPERLIB) $(FWLIB) \
	*.o *.d *.so LinkDef_out*; echo "Done"

table: processed_data.root # show results
ifndef ROOTSYS
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo; exit 1
else
	$(QUIET) root -q -b showResults.C+
endif
t: table

jetbkg: fakeIsoControlSamples.root processed_data.root  # estimate jet induced background
	$(QUIET) echo "ok"

skim: processData.C $(CORELIB) $(LIB) $(FWLIB) # skim data
ifndef ROOTSYS
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo; exit 1
else
	$(QUIET) export SkimSamples=1; root -b -q processData.C
endif
