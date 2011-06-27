# $Id: Makefile,v 1.5 2010/08/15 18:44:30 jmuelmen Exp $

CC = g++
INCLUDE = -I${CMS2_LOCATION}/NtupleMacros
CFLAGS = -Wall -g -O0 -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS) -DTOOLSLIB
TXT2JSON = ./convertGoodRunsToJSON.py
JSON2LUMI = ${CMSSW_BASE}/src/RecoLuminosity/LumiDB/scripts/lumiCalc.py

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

TABLES = Results.tbl

.PHONY: all
all:	check_env $(LIBS) 

.PHONY: check_env
check_env:
	if [ -z $(CMS2_LOCATION) ] ; then echo "first set CMS2_LOCATION" ; exit 1 ; fi

.PHONY: check_env_data
check_env_data:
	if [ -z $(DATA) ] ; then echo "first set DATA" ; exit 1 ; fi
	for i in `ls $(DATA)` ; do if ! [ -e $$i ] ; then echo "data file \"$$i\" does not exist" ; exit 1 ; fi ; done
	if [ -z $(GOODRUN) ] ; then echo "first set GOODRUN" ; exit 1 ; fi

.PHONY: libs
libs:	$(LIBS)

$(TOOLSLIB):	$(TOOLSOBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(TOOLSOBJECTS) -o $@ 

$(CORELIB):	$(COREOBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(COREOBJECTS) -o $@ 

$(LOOPERLIB):	$(LOOPEROBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(LOOPEROBJECTS) -o $@ -lncurses

linkdef_out.cxx: linkdef.h ConfigAndRun.h DSGTable.h DSGSearchWindow.h
	rootcint -f $@ -c $(INCLUDE) ConfigAndRun.h DSGTable.h DSGSearchWindow.h $<

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

# General rule for making JSONs from good-run lists
%.json:	%.txt
	$(TXT2JSON) $< > $@ || rm -f $@

# General rule for calculating lumi from json
%.lumi:	%.json
	$(JSON2LUMI) --nowarning -c frontier://LumiProd/CMS_LUMI_PROD -i $< overview 2>&1 > $@ || rm -f $@

.PHONY: clean
clean:  
	rm -f ${CMS2_LOCATION}/NtupleMacros/Tools/*.o ${CMS2_LOCATION}/NtupleMacros/Tools/*.d $(TOOLSLIB) \
		${CMS2_LOCATION}/NtupleMacros/CORE/*.o ${CMS2_LOCATION}/NtupleMacros/CORE/*.d $(CORELIB) \
		*.o *.d linkdef_out.* $(LOOPERLIB) Results*

# how to make tables
$(TABLES): $(LIBS)

%.tbl %.root %.log:	$(LIBS)
	echo $(@:.tbl=)"()" | root.exe -b setup.C

.PHONY: data
data:	check_env_data $(LIBS) $(GOODRUN) $(GOODRUN:.txt=.lumi)
	export LUMI=`grep -A 5 Total $(GOODRUN:.txt=.lumi) | tail -2 | awk -F \| '{print $$5}' `; \
	echo got $(LIBS), $(DATA), $(GOODRUN), $(GOODRUN:.txt=.lumi), lumi is $${LUMI}; \
	echo "Results_data()" | root.exe -b setup.C

-include $(TOOLSSOURCES:.cc=.d)
-include $(CORESOURCES:.cc=.d)
-include $(LOOPERSOURCES:.cc=.d)
-include linkdef_out.d

.PRECIOUS: $(wildcard *.tbl) $(wildcard *.log) $(wildcard *.root)
