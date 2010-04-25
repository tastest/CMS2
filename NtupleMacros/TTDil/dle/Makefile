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
	@echo "Available Targets:";\
	cat Makefile | perl -ne 'printf("\t%-15s %s\n",$$1,$$2) if(/^(\S+):[^#]+(#.*)$$/)'
ifndef VERBOSE
  QUIET := @
endif

CC = g++
CMSROOT = ../../
ROOFITINCLUDE = 
ifdef CMSSW_VERSION
	ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
endif
INCLUDE = -I$(CMSROOT) $(ROOFITINCLUDE) -I$(CMSROOT)/CORE -I$(CMSROOT)/Tools
CFLAGS = -Wall -Wno-unused-function -g -O2 -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS)

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags)

ifeq ($(shell root-config --platform),macosx)
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
endif

#TOOLSSOURCES = $(wildcard $(CMSROOT)/Tools/*.cc)
TOOLSSOURCES = $(CMSROOT)/Tools/tools.cc
TOOLSOBJECTS = $(TOOLSSOURCES:.cc=.o)
TOOLSLIB = libCMS2NtupleMacrosTOOLS.so

CORESOURCES = \
    $(CMSROOT)/CORE/CMS2.cc \
    $(CMSROOT)/CORE/electronSelections.cc \
	$(CMSROOT)/CORE/electronSelectionsCICParameters.cc \
    $(CMSROOT)/CORE/muonSelections.cc \
    $(CMSROOT)/CORE/trackSelections.cc \
    $(CMSROOT)/CORE/utilities.cc \
    $(CMSROOT)/CORE/eventSelections.cc \
    $(CMSROOT)/CORE/mcSelections.cc \
	$(CMSROOT)/CORE/ttbarSelections.cc \
	$(CMSROOT)/CORE/jetSelections.cc

COREOBJECTS = $(CORESOURCES:.cc=.o) 
CORELIB = libCMS2NtupleMacrosCORE.so

SOURCES = $(wildcard *.cc) 
OBJECTS = $(SOURCES:.cc=.o) LinkDef_out.o
LIB = libCMS2NtupleMacrosLooper.so

FWLIB = ../../Tools/MiniFWLite/libMiniFWLite.so

LIBS = $(CORELIB) $(TOOLSLIB) $(LIB) $(FWLIB)

.PHONY: all help compile clean cms2env

libs:	$(LIBS)

$(CORELIB):	$(COREOBJECTS) 
	$(QUIET) echo "Linking $(CORELIB)"; \
	$(LINKER) $(LINKERFLAGS) -shared $(COREOBJECTS) -o $@ 2>&1|perl -ne 'print if(!/skipping incompatible/)'

$(TOOLSLIB): $(TOOLSOBJECTS) 
	$(QUIET) echo "Linking $(TOOLSLIB)"; \
	$(LINKER) $(LINKERFLAGS) -shared $(TOOLSOBJECTS) -o $@ 2>&1|perl -ne 'print if(!/skipping incompatible/)'

$(LIB):	$(OBJECTS) 
	$(QUIET) echo "Linking $(LIB)"; \
	$(LINKER) $(LINKERFLAGS) -shared $(OBJECTS) -o $@ 2>&1|perl -ne 'print if(!/skipping incompatible/)'

$(FWLIB) :
	$(QUITE) echo "making MiniFWLite"; \
	cd ../../Tools/MiniFWLite; \
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


#################################################################
#                       User Targets 
#################################################################
build: $(LIBS) # compile code
	@#root -b -q ".L MyScanChain.C++"
b: build

clean: # clean up 
	$(QUIET) rm -f \
	$(CMSROOT)/CORE/*.o \
	$(CMSROOT)/CORE/*.d \
	$(CMSROOT)/CORE/*.so \
	$(CMSROOT)/Tools/*.o \
	$(CMSROOT)/Tools/*.d \
	$(CMSROOT)/Tools/*.so \
	$(CORELIB) $(TOOLSLIB) $(LOOPERLIB) $(FWLIB) \
	*.o *.d *.so LinkDef_out*; echo "Done"

endif

