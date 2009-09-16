include Makefile.arch

ROOTINCDIR   := $(shell $(ROOTCONFIG) --incdir)
LDFLAGS      += -fPIC

FORWRAPPERH   = vector.h Wrapper.h Math/LorentzVector.h Math/PositionVector3D.h
MINIFWLITEH   = $(FORWRAPPERH) DataFormats.h

all:	libMiniFWLite.so

libMiniFWLite.so: MiniFWLiteLinkDef.h
		@rootcint -v3 -f MiniFWLiteDict.cc -c -p $(MINIFWLITEH) MiniFWLiteLinkDef.h
		$(LD) $(LDFLAGS) -Wno-deprecated $(SOFLAGS) -I$(ROOTINCDIR) MiniFWLiteDict.cc TaSBanner.cc -o $@

clean:
		@rm -f core *.o *.so MiniFWLiteDict.*
