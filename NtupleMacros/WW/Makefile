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

CORESOURCES = $(CMSROOT)/CORE/CMS2.cc $(CMSROOT)/CORE/muonSelections.cc $(CMSROOT)/CORE/electronSelections.cc $(CMSROOT)/CORE/electronSelectionsParameters.cc $(CMSROOT)/CORE/jetSelections.cc 
#$(CMSROOT)/CORE/fakerates.cc
# CORESOURCES = $(wildcard $(CMSROOT)/CORE/*.cc) 
COREOBJECTS = $(CORESOURCES:.cc=.o) 
CORELIB = libCMS2NtupleMacrosCORE.so

SOURCES = $(wildcard *.cc) 
OBJECTS = $(SOURCES:.cc=.o) LinkDef_out.o
LIB = libCMS2NtupleMacrosLooper.so

FWLIB = ../Tools/MiniFWLite/libMiniFWLite.so

LIBS = $(CORELIB) $(LIB)

.PHONY: all help compile clean cms2env

libs:	$(LIBS)

$(CORELIB):	$(CORESOURCES) $(COREOBJECTS) 
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

data:
	@echo 'Make a directory or a symbolic link called "data" pointing to sample location'; echo; exit 1

processed_data.root: data processData.C $(CORELIB) $(LIB) $(FWLIB)
	$(QUIET) root -b -q processData.C

fakeIsoControlSamples.root: makeFakeIsoControlSamples.C
	$(QUIET) root -b -q makeFakeIsoControlSamples.C+

processed_data.eps: processed_data.root makePSFile.C
	$(QUIET) root -q -b 'makePSFile.C+("processed_data.root")'

#################################################################
#                       User Targets 
#################################################################
build: $(LIBS) # compile code
	@#root -b -q ".L doAnalysis.C++"
b: build

clean: # clean up 
	$(QUIET) rm -v -f \
	$(CMSROOT)/CORE/*.o \
	$(CMSROOT)/CORE/*.d \
	$(CMSROOT)/CORE/*.so \
	$(CMSROOT)/Tools/*.so \
	$(CORELIB) $(LOOPERLIB) $(FWLIB) \
	processed_data.* *.list \
	*.o *.d *.so LinkDef_out*; echo "Done"

table: processed_data.root # show results
	$(QUIET) root -q -l showResults.C+
t: table

jetbkg: processed_data.root  # estimate jet induced background
	$(QUIET) root -q -l estimateJetBkg_fakes.C+; \
	root -q -l estimateJetBkg_sideband.C+

ElFrVer = el_fr_v2_wwV1
MuFrVer = mu_fr_fo_wwV1_10
ElFrFile = ww_el_fr_EGandEGMon.root
MuFrFile = ww_mu_fr_Mu.root

jetbkg2: processed_data.root # fake rate method
	$(QUIET) root -q -l 'estimateJetBkg_fakes.C+("processed_data.root","files/$(ElFrFile)","$(ElFrVer)","files/$(MuFrFile)","$(MuFrVer)")'

topbkg: processed_data.root  # estimate top background
	$(QUIET) root -q -l estimateTopBkg.C+

plots: processed_data.eps # make plots based on processed data
	$(QUIET) epstopdf processed_data.eps

skim: data processData.C $(CORELIB) $(LIB) $(FWLIB)# skim data
	$(QUIET) if [ ! -n "$(SkimName)" ]; then echo "Please set SkimName variable in order to proceed"; exit 1; fi; \
	export SkimSamples=1; root -b -q processData.C
#	export SkimSamples=1; valgrind --track-origins=yes root.exe -b -q processData.C
lumi: processed_data.root processed.json # estimated integrated luminosity looking at processed runs
	$(QUIET) lumiCalc.py -c frontier://LumiProd/CMS_LUMI_PROD -i processed.json overview
lumi2: processed_data.root processed_detailed.json # estimated integrated luminosity looking at processed runs and lumis
	$(QUIET) lumiCalc.py -c frontier://LumiProd/CMS_LUMI_PROD -i processed_detailed.json overview
endif
