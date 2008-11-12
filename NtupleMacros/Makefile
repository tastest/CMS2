# $Id: Makefile,v 1.5 2008/11/12 23:52:50 jmuelmen Exp $

CC = g++
INCLUDE = -ITools -IWW -IDSG
CFLAGS = -Wall -g -O2 -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS) -DTOOLSLIB

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags)

TOOLSSOURCES = $(wildcard Tools/*.cc) 
TOOLSOBJECTS = $(TOOLSSOURCES:.cc=.o) # Tools/linkdef_out.o
TOOLSLIB = libCMS2NtupleMacrosTools.so

WWSOURCES = $(wildcard WW/*.cc) 
WWOBJECTS = $(WWSOURCES:.cc=.o) 
WWLIB = libCMS2NtupleMacrosWW.so

DSGSOURCES = $(wildcard DSG/*.cc) 
DSGOBJECTS = $(DSGSOURCES:.cc=.o) 
DSGLIB = libCMS2NtupleMacrosDSG.so

DICTOBJECTS = linkdef_out.o
DICTLIB = libCMS2Dict.so

LIBS = $(TOOLSLIB) $(WWLIB) $(DSGLIB) $(DICTLIB)

WWTABLES = WW_Results.tbl WW_Results_softmu.tbl WW_NoMuTag.tbl WW_NoTrackJets.tbl WW_NoCaloIso.tbl \
	   WW_OldResults.tbl WW_OldResults_trkjets.tbl WW_OldResults_caloiso.tbl WW_OldResults_muveto.tbl \
	   WW_OldResults_caloiso_trkjets.tbl WW_OldResults_caloiso_muveto.tbl \
	   WW_TopEstimate.tbl WW_TopEstimate_noptcut.tbl 
DSGTABLES = DSG_BaseLine.tbl DSG_MET_10.tbl DSG_MET_1.tbl

.PHONY: all
all:	$(TOOLSLIB) $(WWLIB) $(WWTABLES) $(DSGLIB) $(DSGTABLES)

.PHONY: libs
libs:	$(LIBS)

$(TOOLSLIB):	$(TOOLSOBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(TOOLSOBJECTS) -o $@ 

$(WWLIB):	$(WWOBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(WWOBJECTS) -o $@ 

$(DSGLIB):	$(DSGOBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(DSGOBJECTS) -o $@ 

$(DICTLIB):	$(DICTOBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(DICTOBJECTS) -o $@ 

linkdef_out.cxx: linkdef.h Tools/Sample.h Tools/LooperBase.h WW/WWDriver.h DSG/DSGDriver.h
	rootcint -f $@ -c $(INCLUDE) Sample.h LooperBase.h WWDriver.h DSGDriver.h $<

Tools/linkdef_out.cxx: Tools/linkdef.h Tools/Sample.h Tools/LooperBase.h
	rootcint -f $@ -c $(INCLUDE) Sample.h LooperBase.h $<

WW/linkdef_out.cxx: WW/linkdef.h WW/WWDriver.h
	rootcint -f $@ -c $(INCLUDE) WWDriver.h $<

DSG/linkdef_out.cxx: DSG/linkdef.h DSG/DSGDriver.h
	rootcint -f $@ -c $(INCLUDE) DSGDriver.h $<

# General rule for making object files
%.d:	%.cc
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@

%.o: 	%.cc 
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.cxx
	$(CC) $(CFLAGS) $< -c -o $@

.PHONY: clean
clean:  
	rm -f Tools/*.o Tools/*.d Tools/linkdef_out.* $(TOOLSLIB) \
		WW/*.o WW/*.d WW/linkdef_out.* $(WWLIB) \
		DSG/*.o DSG/*.d DSG/linkdef_out.* $(DSGLIB) \
		$(WWTABLES) $(DSGTABLES)

# how to make tables
$(WWTABLES): $(LIBS)

# how to make tables
$(DSGTABLES): $(LIBS)

%.tbl:	$(LIBS)
	echo $(@:.tbl=)"()" | root.exe -b


-include $(TOOLSSOURCES:.cc=.d)
-include $(WWSOURCES:.cc=.d)
-include $(DSGSOURCES:.cc=.d)

