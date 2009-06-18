# $Id: Makefile,v 1.2 2009/06/18 17:08:12 fgolf Exp $

CC = g++
INCLUDE = -I${CMS2_LOCATION}/NtupleMacros -I$(CLHEP_PARAM_PATH)/include
CFLAGS = -Wall -Wno-unused-function -g -O2 -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS) -DTOOLSLIB

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

TABLES = Results.tbl Feb_Results.tbl Feb_Results_With_Ntrks.tbl		\
Feb_Results_With_Ntrks.tbl Feb_Results_With_TcMET.tbl			\
Feb_Results_With_TrackJets.tbl Feb_Results_With_Btags.tbl		\
Feb_Results_NoPass4MET.tbl Feb_Results_With_CaloIso.tbl			\
Oct_Results.tbl Results_NoTcMET.tbl Results_NoTrackJets.tbl		\
Results_NoBtags.tbl Results_NoCaloIso.tbl Results_NoPass4MET.tbl	\
Results_NoNtrks.tbl Wjets_Numerator.tbl Wjets_FOs_Not_Numerator.tbl	\
Wjets_Fakerate.tbl

.PHONY: all
all:	$(LIBS) 

.PHONY: libs
libs:	$(LIBS)

.PHONY: tables
tables:	$(TABLES)

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