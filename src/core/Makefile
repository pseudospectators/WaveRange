include ../../config.mk

OUTPUTDIR = ../../libc/

OBJECTC = wrappers.o ../waveletcdf97_3d/waveletcdf97_3d.o ../rangecod/rangecod.o
CXXSOURCES = wrappers.cpp

ifeq ($(CC),gcc)
  CPICFLAG = -fPIC
  SHARELIB = gcc -shared -o $(OUTPUTDIR)libwaverange.so $(OBJECTC)
else
  CPICFLAG = 
  SHARELIB =
endif

all: wrappers

wrappers: 
	$(CXX) $(CXXSOURCES) $(CDNLFLAG) $(CPICFLAG) $(CXXFLAGS) 
	$(AR) rcs $(OUTPUTDIR)libwaverange.a $(OBJECTC)
	$(SHARELIB)
.PHONY: clean
clean:
	$(RM) ./*.gc??
	$(RM) ./*.o
