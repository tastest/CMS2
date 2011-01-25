
CC = g++
INCLUDE = -I../
CFLAGS = -Wall -g -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS) -DTOOLSLIB
LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags)
CMSROOT=../

CORESOURCES = \
    $(CMSROOT)/CORE/CMS2.cc \
    $(CMSROOT)/CORE/electronSelections.cc \
	$(CMSROOT)/CORE/electronSelectionsParameters.cc \
    $(CMSROOT)/CORE/muonSelections.cc \
    $(CMSROOT)/CORE/utilities.cc \
    $(CMSROOT)/CORE/eventSelections.cc \
    $(CMSROOT)/CORE/mcSelections.cc \
	$(CMSROOT)/CORE/trackSelections.cc \
	$(CMSROOT)/CORE/jetSelections.cc \
	$(CMSROOT)/CORE/triggerUtils.cc \
    $(CMSROOT)/CORE/metSelections.cc \
    $(CMSROOT)/CORE/ssSelections.cc
COREOBJECTS = $(CORESOURCES:.cc=.o) 
CORELIB = libCMS2NtupleMacrosCORE.so

TOOLSSOURCES = $(CMSROOT)/Tools/tools.cc \
			   $(CMSROOT)/Tools/goodrun.cc 
TOOLSOBJECTS = $(TOOLSSOURCES:.cc=.o) 
TOOLSLIB = libCMS2NtupleMacrosTools.so

MT2SOURCES = $(CMSROOT)/CORE/MT2/MT2.cc
MT2OBJECTS = $(MT2SOURCES:.cc=.o)
MT2LIB = libCMS2NtupleMacrosMT2.so

LOOPERSOURCES = $(wildcard *.cc) 
LOOPEROBJECTS = $(LOOPERSOURCES:.cc=.o) LinkDef_out.o
LOOPERLIB = libCMS2NtupleMacrosLooper.so

FWLIB = ../Tools/MiniFWLite/libMiniFWLite.so

$(TOOLSLIB):	$(TOOLSOBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(TOOLSOBJECTS) -o $@ 

$(CORELIB):	$(COREOBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(COREOBJECTS) -o $@ 

$(LOOPERLIB):	$(LOOPEROBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(LOOPEROBJECTS) -o $@ 

$(MT2LIB):   $(MT2OBJECTS)
	$(LINKER) $(LINKERFLAGS) -shared $(MT2OBJECTS) -o $@ 

$(FWLIB) :
	$(QUITE) echo "making MiniFWLite"; \
	cd ../Tools/MiniFWLite; \
	$(MAKE) -f Makefile;

LinkDef_out.cxx: LinkDef.h emubabymaker.h dilepbabymaker.h
	rootcint -f $@ -c $(INCLUDE) emubabymaker.h dilepbabymaker.h $<

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

LIBS = $(CORELIB) $(TOOLSLIB) $(LOOPERLIB) $(FWLIB) $(MT2LIB)

.PHONY: all
all:	$(LIBS)  

.PHONY: clean
clean:  
	rm -f $(CMSROOT)/Tools/*.o $(CMSROOT)/Tools/*.d $(TOOLSLIB) \
		$(CMSROOT)/CORE/*.o $(CMSROOT)/NtupleMacros/CORE/*.d $(CORELIB) \
		$(CMSROOT)/CORE/MT2/*.o $(CMSROOT)/NtupleMacros/CORE/MT2/*.d $(MT2LIB) \
		*.o *.d LinkDef_out.* $(LOOPERLIB)

-include $(TOOLSSOURCES:.cc=.d)
-include $(CORESOURCES:.cc=.d)
-include $(LOOPERSOURCES:.cc=.d)
-include $(MT2SOURCES:.cc=.d)
-include LinkDef_out.d

.PRECIOUS: $(wildcard *.tbl) $(wildcard *.log) $(wildcard *.root)

