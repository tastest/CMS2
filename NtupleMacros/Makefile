# $Id: Makefile,v 1.2 2008/10/15 08:41:28 jmuelmen Exp $

CC = g++
INCLUDE = -ITools -IWW
CFLAGS = -Wall -g -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS) -DTOOLSLIB

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags)

TOOLSSOURCES = $(wildcard Tools/*.cc) 
TOOLSOBJECTS = $(TOOLSSOURCES:.cc=.o) # Tools/linkdef_out.o
TOOLSLIB = libCMS2NtupleMacrosTools.so

WWSOURCES = $(wildcard WW/*.cc) 
WWOBJECTS = $(WWSOURCES:.cc=.o) WW/linkdef_out.o
WWLIB = libCMS2NtupleMacrosWW.so

LIBS = $(TOOLSLIB) $(WWLIB)

TABLES = WW_Results.tbl WW_Results_softmu.tbl WW_NoMuTag.tbl WW_NoTrackJets.tbl WW_NoCaloIso.tbl \
	WW_OldResults.tbl WW_OldResults_trkjets.tbl WW_OldResults_caloiso.tbl WW_OldResults_muveto.tbl \
	WW_OldResults_caloiso_trkjets.tbl WW_OldResults_caloiso_muveto.tbl \
	WW_TopEstimate.tbl WW_TopEstimate_noptcut.tbl 

.PHONY: all
all:	$(TOOLSLIB) $(WWLIB) $(TABLES)

.PHONY: libs
libs:	$(LIBS)

$(TOOLSLIB):	$(TOOLSOBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(TOOLSOBJECTS) -o $@ 

$(WWLIB):	$(WWOBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(WWOBJECTS) -o $@ 

Tools/linkdef_out.cxx: Tools/linkdef.h Tools/Sample.h Tools/LooperBase.h
	rootcint -f $@ -c $(INCLUDE) Sample.h LooperBase.h $<

WW/linkdef_out.cxx: WW/linkdef.h WW/WWDriver.h WW/WWLooper.h
	rootcint -f $@ -c $(INCLUDE) WWDriver.h WWLooper.h $<

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
		$(TABLES)

# how to make tables
$(TABLES): $(LIBS)

%.tbl:	$(LIBS)
	echo $(@:.tbl=)"()" | root.exe -b

-include $(TOOLSSOURCES:.cc=.d)
-include $(WWSOURCES:.cc=.d)

