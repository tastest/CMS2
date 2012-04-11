CC = g++
ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
INCLUDE = -I../ -I./ $(ROOFITINCLUDE)
CFLAGS = -Wall -g -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS) -DTOOLSLIB
LINKER = g++

LINKERFLAGS = $(shell root-config --ldflags)
ifeq ($(shell root-config --platform),macosx)
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
endif

SOURCES = LeptonTreeLooper.cc PlotUtilities.cc ../Tools/goodrun.cc
OBJECTS = $(SOURCES:.cc=.o) LinkDef_out.o
LIB = libLeptonTreeLooper.so

$(LIB):	$(OBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(OBJECTS) -o $@ 

LinkDef_out.cxx: LinkDef.h LeptonTreeLooper.h PlotUtilities.h
	rootcint -f $@ -c $(INCLUDE) LeptonTreeLooper.h PlotUtilities.h $<

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

LIBS = $(LIB) 

.PHONY: all
all:	$(LIBS)  

.PHONY: clean
clean:  
	rm -f *.d \
	rm -f *.o \
	rm -f *.so

-include $(SOURCES:.cc=.d)
-include $(LIBDIR)/LinkDef_out.d

