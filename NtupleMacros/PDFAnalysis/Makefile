CC = g++
LHAPDF = /tas/dlevans/HWW2012/CMSSW_5_2_3/src/LHAPDF-5.8.92b/
INCLUDE = -I../ -I./ -I$(LHAPDF)/include
CFLAGS = -Wall -g -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS) -DTOOLSLIB
LINKER = g++

LINKERFLAGS = $(shell root-config --ldflags)
ifeq ($(shell root-config --platform),macosx)
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
endif

SOURCES = MyScanChain.cc PlotUtilities.cc ../CORE/CMS2.cc
OBJECTS = $(SOURCES:.cc=.o) LinkDef_out.o
LIB = libCMS2NtupleMacrosLooper.so

$(LIB):	$(OBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(OBJECTS) -o $@ 

LIBS = $(LIB) 

LinkDef_out.cxx: LinkDef.h MyScanChain.h PlotUtilities.h
	rootcint -f $@ -c $(INCLUDE) MyScanChain.h PlotUtilities.h $<

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

.PHONY: clean all
all:  $(LIBS)

clean:  
	rm -f *.d \
	rm -f *.o \
	rm -f *.so \
	rm -f *.cxx
	rm -f ../CORE/CMS2.o \
	rm -f ../CORE/CMS2.d 

-include $(SOURCES:.cc=.d)
-include $(LIBDIR)/LinkDef_out.d

