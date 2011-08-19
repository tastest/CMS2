# $Id: Makefile,v 1.5 2009/07/25 17:21:52 ibloch Exp $

CC = g++
INCLUDE = -I${CMS2_LOCATION}/NtupleMacros
CFLAGS = -Wall -Wno-unused-function -g -O0 -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS) -DTOOLSLIB

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags)
ifeq ($(shell root-config --platform),macosx)
	# linker option for Mac OS X, got it from old project, seems to work ;-)
	LINKERFLAGS = -bundle -undefined suppress -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
endif

CORESOURCES = $(wildcard ${CMS2_LOCATION}/NtupleMacros/CORE/*.cc) 
COREOBJECTS = $(CORESOURCES:.cc=.o) 
CORELIB = libCMS2NtupleMacrosCORE.so

TOOLSSOURCES = $(wildcard ${CMS2_LOCATION}/NtupleMacros/Tools/*.cc) 
TOOLSOBJECTS = $(TOOLSSOURCES:.cc=.o) 
TOOLSLIB = libCMS2NtupleMacrosTools.so

LOOPERSOURCES = $(wildcard *.cc) 
LOOPEROBJECTS = $(LOOPERSOURCES:.cc=.o) linkdef_out.o
LOOPERLIB = libCMS2NtupleMacrosLooper.so

LIBS = $(CORELIB) $(TOOLSLIB) $(LOOPERLIB)

RESULTTABLES = \
SSSUSY.tbl\
OSSUSY.tbl\
TtbarBaseline_SS_FlipPredicted.tbl\
TtbarBaseline_SS_FlipObserved.tbl\
FRPredOSSUSY.tbl\
FRPredSSSUSY.tbl

FLIPTABLES = \
TtbarLoose_SS_FlipPredicted.tbl\
TtbarLoose_SS_FlipObserved.tbl\
TtbarBaseline_SS_FlipPredicted.tbl\
TtbarBaseline_SS_FlipObserved.tbl

FAKETABLES = \
Wjets_Numerator.tbl\
Wjets_FOs_Not_Numerator.tbl\
Wjets_Fakerate.tbl\
Wjets_SS_Numerator.tbl\
Wjets_SS_FOs_Not_Numerator.tbl\
Wjets_SS_Fakerate.tbl\
Ttbar_SS_Numerator.tbl\
Ttbar_SS_FOs_Not_Numerator.tbl\
Ttbar_SS_Fakerate.tbl\
Wjets_MuonNumerator.tbl\
Wjets_FOs_Not_MuonNumerator.tbl\
Wjets_MuonFakerate.tbl\
Wjets_SS_MuonNumerator.tbl\
Wjets_SS_FOs_Not_MuonNumerator.tbl\
Wjets_SS_MuonFakerate.tbl\
Ttbar_SS_MuonNumerator.tbl\
Ttbar_SS_FOs_Not_MuonNumerator.tbl\
Ttbar_SS_MuonFakerate.tbl

.PHONY: all
all:	$(LIBS) 

.PHONY: libs
libs:	$(LIBS)

.PHONY: tables
tables:	resulttables fliptables faketables

.PHONY: resulttables
resulttables: $(RESULTTABLES)

.PHONY: fliptables
fliptables: $(FLIPTABLES)

.PHONY: faketables
faketables: $(FAKETABLES)


$(TOOLSLIB):	$(TOOLSOBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(TOOLSOBJECTS) -o $@ 

$(CORELIB):	$(COREOBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(COREOBJECTS) -o $@ 

$(LOOPERLIB):	$(LOOPEROBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(LOOPEROBJECTS) -o $@ 

linkdef_out.cxx: ${CMS2_LOCATION}/NtupleMacros/linkdef.h ConfigAndRun.h
	rootcint -f $@ -c $(INCLUDE) ConfigAndRun.h $<

# General rule for making object files
%.d:	%.cc
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@
%.d:	%.cxx
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@

%.o: 	%.cc 
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.cxx
	$(CC) $(CFLAGS) $< -c -o $@

.PHONY: clean
clean:  
	rm -f ${CMS2_LOCATION}/NtupleMacros/Tools/*.o ${CMS2_LOCATION}/NtupleMacros/Tools/*.d $(TOOLSLIB) \
		${CMS2_LOCATION}/NtupleMacros/CORE/*.o ${CMS2_LOCATION}/NtupleMacros/CORE/*.d $(CORELIB) \
		*.o *.d linkdef_out.* $(LOOPERLIB) \
		$(TABLES) 

# how to make tables
$(TABLES): $(LIBS)

%.tbl %.root %.log:	$(LIBS)
	echo $(@:.tbl=)"()" | root.exe -b setup.C

-include $(TOOLSSOURCES:.cc=.d)
-include $(CORESOURCES:.cc=.d)
-include $(LOOPERSOURCES:.cc=.d)
-include linkdef_out.d

.PRECIOUS: $(wildcard *.tbl) $(wildcard *.log) $(wildcard *.root)
