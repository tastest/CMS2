# $Id: Makefile,v 1.1 2008/11/20 19:46:18 jmuelmen Exp $

CC = g++
INCLUDE = -I../Tools -I../CORE
CFLAGS = -Wall -g -O2 -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS) -DTOOLSLIB

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags)

CORESOURCES = $(wildcard ../CORE/*.cc) 
COREOBJECTS = $(CORESOURCES:.cc=.o) 
CORELIB = libCMS2NtupleMacrosCORE.so

TOOLSSOURCES = $(wildcard ../Tools/*.cc) 
TOOLSOBJECTS = $(TOOLSSOURCES:.cc=.o) 
TOOLSLIB = libCMS2NtupleMacrosTools.so

LOOPERSOURCES = $(wildcard *.cc) 
LOOPEROBJECTS = $(LOOPERSOURCES:.cc=.o) linkdef_out.o
LOOPERLIB = libCMS2NtupleMacrosLooper.so

LIBS = $(CORELIB) $(TOOLSLIB) $(LOOPERLIB)

TABLES = Results.tbl

.PHONY: all
all:	$(LIBS) 

.PHONY: libs
libs:	$(LIBS)

$(TOOLSLIB):	$(TOOLSOBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(TOOLSOBJECTS) -o $@ 

$(CORELIB):	$(COREOBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(COREOBJECTS) -o $@ 

$(LOOPERLIB):	$(LOOPEROBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(LOOPEROBJECTS) -o $@ 

linkdef_out.cxx: ../linkdef.h ConfigAndRun.h
	rootcint -f $@ -c $(INCLUDE) ConfigAndRun.h $<

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
	rm -f ../Tools/*.o ../Tools/*.d $(TOOLSLIB) \
		../CORE/*.o ../CORE/*.d $(CORELIB) \
		*.o *.d linkdef_out.* $(LOOPERLIB) \
		$(TABLES) 

# how to make tables
$(TABLES): $(LIBS)

%.tbl %.root %.log:	$(LIBS)
	echo $(@:.tbl=)"()" | root.exe -b

-include $(TOOLSSOURCES:.cc=.d)
-include $(CORESOURCES:.cc=.d)
-include $(LOOPERSOURCES:.cc=.d)

