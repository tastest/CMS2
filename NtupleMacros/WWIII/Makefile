# $Id: Makefile,v 1.11 2009/06/22 13:18:58 jmuelmen Exp $

CC = g++
INCLUDE = -I${CMS2_LOCATION}/NtupleMacros
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

FAKETABLES = Wjets_Numerator.tbl Wjets_FOs_Not_Numerator.tbl	\
Wjets_Fakerate.tbl Wjets_Oingo.tbl \
Wjets_SS_Fakerate.tbl Wjets_SS_FOs_Not_Numerator.tbl Wjets_SS_Numerator.tbl \
Wjets_Histat.tbl Wjets_Histat_Fakerate.tbl Wjets_Histat_SS.tbl Wjets_Histat_SS_Fakerate.tbl

NM1TABLES = Results_NoBtags.tbl Results_NoCaloIso.tbl Results_CaloIso_1_6.tbl	\
Calojet.tbl Calojet_Trkjet.tbl Results_NoJetVeto.tbl Results_NoPass4MET.tbl\
Results_NoPass2MET.tbl Results_NoZVeto.tbl

RESULTTABLES = Results.tbl Results_W.tbl Results_DY.tbl Results_Vgamma.tbl

EFFTABLES = Efficiency_base.tbl Efficiency_trigger.tbl Efficiency_id.tbl Efficiency_iso.tbl \
Efficiency_tcmet.tbl Efficiency_jet.tbl Efficiency_zveto.tbl Efficiency_muveto.tbl

.PHONY: all
all:	$(LIBS) 

.PHONY: libs
libs:	$(LIBS)

.PHONY: tables
tables:	$(RESULTTABLES) $(NM1TABLES)  $(FAKETABLES) 

.PHONY: faketables
faketables: $(FAKETABLES)

.PHONY: efftables
efftables: $(EFFTABLES)

.PHONY: nm1tables
nm1tables: $(NM1TABLES)

.PHONY: resulttables
resulttables: $(RESULTTABLES)

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
