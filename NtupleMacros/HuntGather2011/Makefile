
CC = g++
INCLUDE = -I../ -I./
CFLAGS = -Wall -g -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS) -DTOOLSLIB
LINKER = g++

LINKERFLAGS = $(shell root-config --ldflags)
ifeq ($(shell root-config --platform),macosx)
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
endif

CMSROOT=../
PLOTTERSRC=./plotter/
BABYMAKERSRC=./babymaker/
LIBDIR=libs

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
CORELIB = $(LIBDIR)/libCMS2NtupleMacrosCORE.so

TOOLSSOURCES = $(CMSROOT)/Tools/tools.cc \
			   $(CMSROOT)/Tools/goodrun.cc 
TOOLSOBJECTS = $(TOOLSSOURCES:.cc=.o)
TOOLSLIB = $(LIBDIR)/libCMS2NtupleMacrosTools.so

MT2SOURCES = $(CMSROOT)/CORE/MT2/MT2.cc
MT2OBJECTS = $(MT2SOURCES:.cc=.o)
MT2LIB = $(LIBDIR)/libCMS2NtupleMacrosCOREMT2.so

BABYMAKERSOURCES = $(wildcard $(BABYMAKERSRC)/*.cc) 
BABYMAKEROBJECTS = $(BABYMAKERSOURCES:.cc=.o) $(LIBDIR)/LinkDef_out.o
BABYMAKERLIB = $(LIBDIR)/libHuntGather2011Babymaker.so

PLOTTERSOURCES = $(wildcard $(PLOTTERSRC)/*.cc)
PLOTTEROBJECTS = $(PLOTTERSOURCES:.cc=.o) 
PLOTTERLIB = $(LIBDIR)/libHuntGather2011Plotter.so

FWLIB = ../Tools/MiniFWLite/libMiniFWLite.so

$(TOOLSLIB):	$(TOOLSOBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(TOOLSOBJECTS) -o $@ 

$(CORELIB):	$(COREOBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(COREOBJECTS) -o $@ 

$(BABYMAKERLIB):	$(BABYMAKEROBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(BABYMAKEROBJECTS) -o $@ 

$(MT2LIB):   $(MT2OBJECTS)
	$(LINKER) $(LINKERFLAGS) -shared $(MT2OBJECTS) -o $@ 

$(PLOTTERLIB):   $(PLOTTEROBJECTS)
	$(LINKER) $(LINKERFLAGS) -shared $(PLOTTEROBJECTS) -o $@ 


$(FWLIB) :
	$(QUITE) echo "making MiniFWLite"; \
	cd ../Tools/MiniFWLite; \
	$(MAKE) -f Makefile;

$(LIBDIR)/LinkDef_out.cxx: $(LIBDIR)/LinkDef.h  $(CMSROOT)/Tools/goodrun.h $(PLOTTERSRC)/BabyDorkIdentifier.h $(PLOTTERSRC)/BabySample.h $(PLOTTERSRC)/gather.h $(BABYMAKERSRC)/emubabymaker.h $(BABYMAKERSRC)/dilepbabymaker.h 
	rootcint -f $@ -c $(INCLUDE)  $(CMSROOT)/Tools/goodrun.h $(PLOTTERSRC)/BabyDorkIdentifier.h $(PLOTTERSRC)/BabySample.h $(PLOTTERSRC)/gather.h $(BABYMAKERSRC)/emubabymaker.h $(BABYMAKERSRC)/dilepbabymaker.h $<

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

LIBS = $(CORELIB) $(TOOLSLIB) $(BABYMAKERLIB) $(FWLIB) $(MT2LIB) $(PLOTTERLIB)

.PHONY: all
all:	$(LIBS)  

.PHONY: plotter
plotter: $(PLOTTERLIB)

.PHONY: clean
clean:  
	rm -f $(CMSROOT)/Tools/*.o $(CMSROOT)/Tools/*.d $(TOOLSLIB) \
		$(CMSROOT)/CORE/*.o $(CMSROOT)/NtupleMacros/CORE/*.d $(CORELIB) \
		$(CMSROOT)/CORE/MT2/*.o $(CMSROOT)/NtupleMacros/CORE/MT2/*.d $(MT2LIB) \
		plotter/*.o plotter/*.d $(PLOTTERLIB) \
		babymaker/*.o babymaker/*.d $(LIBDIR)/LinkDef_out.* $(BABYMAKERLIB)

-include $(TOOLSSOURCES:.cc=.d)
-include $(CORESOURCES:.cc=.d)
-include $(BABYMAKERSOURCES:.cc=.d)
-include $(MT2SOURCES:.cc=.d)
-include $(PLOTTERSOURCES:.cc=.d)
-include $(LIBDIR)/LinkDef_out.d

