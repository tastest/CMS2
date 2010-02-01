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
INCLUDE = -I${CMS2_LOCATION}/NtupleMacros
# CFLAGS = -Wall -Wno-unused-function -g -O2 -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS) -DTOOLSLIB
CFLAGS = -Wall -Wno-unused-function -g -O2 -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS)

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags)
# ifeq ($(shell root-config --platform),macosx)
# 	LINKERFLAGS = -bundle -undefined suppress -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
# endif

ifeq ($(shell root-config --platform),macosx)
#	LINKERFLAGS = -dynamiclib -undefined suppress -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
endif

CORESOURCES = $(wildcard ${CMS2_LOCATION}/NtupleMacros/CORE/*.cc) 
COREOBJECTS = $(CORESOURCES:.cc=.o) 
CORELIB = libCMS2NtupleMacrosCORE.so

TOOLSSOURCES = $(wildcard ${CMS2_LOCATION}/NtupleMacros/Tools/*.cc) 
TOOLSOBJECTS = $(TOOLSSOURCES:.cc=.o) 
TOOLSLIB = libCMS2NtupleMacrosTools.so

SOURCES = $(wildcard *.cc) 
OBJECTS = $(SOURCES:.cc=.o) LinkDef_out.o
LIB = libCMS2NtupleMacrosLooper.so

# LIBS = $(CORELIB) $(TOOLSLIB) $(LIB)
LIBS = $(CORELIB) $(LIB)


.PHONY: all help compile clean cms2env rootenv

libs:	$(LIBS)

$(TOOLSLIB):	$(TOOLSOBJECTS) 
	$(QUIET) echo "Linking $(TOOLSLIB)"; \
	$(LINKER) $(LINKERFLAGS) -shared $(TOOLSOBJECTS) -o $@ 2>&1|perl -ne 'print if(!/skipping incompatible/)'

$(CORELIB):	$(COREOBJECTS) 
	$(QUIET) echo "Linking $(CORELIB)"; \
	$(LINKER) $(LINKERFLAGS) -shared $(COREOBJECTS) -o $@ 2>&1|perl -ne 'print if(!/skipping incompatible/)'

$(LIB):	$(OBJECTS) 
	$(QUIET) echo "Linking $(LIB)"; \
	$(LINKER) $(LINKERFLAGS) -shared $(OBJECTS) -o $@ 2>&1|perl -ne 'print if(!/skipping incompatible/)'

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

cms2env:
ifndef CMS2_LOCATION
	@echo "Please set CMS2_LOCATION"; echo; exit 1
endif

rootenv:
ifndef ROOTSYS
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo; exit 1
endif

processed_data.root: processData.C $(LIB) $(CORELIB)
ifndef CMS2_NTUPLE_LOCATION
	@echo "Please set CMS2_NTUPLE_LOCATION"; echo; exit 1
endif
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
build: rootenv cms2env $(LIBS) # compile code
	@#root -b -q ".L doAnalysis.C++"

clean: cms2env # clean up 
	$(QUIET) rm -f \
	${CMS2_LOCATION}/NtupleMacros/Tools/*.o \
	${CMS2_LOCATION}/NtupleMacros/Tools/*.d \
	${CMS2_LOCATION}/NtupleMacros/Tools/*.so \
	${CMS2_LOCATION}/NtupleMacros/CORE/*.o \
	${CMS2_LOCATION}/NtupleMacros/CORE/*.d \
	${CMS2_LOCATION}/NtupleMacros/CORE/*.so \
	$(TOOLSLIB) $(CORELIB) $(LOOPERLIB) \
	*.o *.d *.so LinkDef_out*; echo "Done"

table: processed_data.root # show results
ifndef ROOTSYS
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo; exit 1
else
	$(QUIET) root -q -b showResults.C+
endif

jetbkg: fakeIsoControlSamples.root processed_data.root  # estimate jet induced background
	$(QUIET) echo "ok"