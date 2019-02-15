include ../../config.mk

OUTPUTDIR = ../../libc/

OBJECTC = wrappers.o ../waveletcdf97_3d/waveletcdf97_3d.o ../rangecod/rangecod.o
CXXSOURCES = wrappers.cpp

all: wrappers

wrappers: 
	$(CXX) $(CXXSOURCES) $(CDNLFLAG) $(CXXFLAGS) 
	$(AR) rcs $(OUTPUTDIR)libwaverange.a $(OBJECTC)
.PHONY: clean
clean:
	$(RM) ./*.gc??
	$(RM) ./*.o
