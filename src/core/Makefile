include ../../config.mk

OUTPUTDIR = ../../bin/
AUXDIR = ../../libc/

OBJECTC = wrappers.o ../waveletcdf97_3d/waveletcdf97_3d.o ../rangecod/rangecod.o
CXXSOURCES = wrappers.cpp

ifeq ($(CC),gcc)
  CPICFLAG = -fPIC
  SODIR = $(OUTPUTDIR)lib/
  SHARELIB = gcc -shared -o $(SODIR)libwaverange.so $(OBJECTC)
else
  CPICFLAG = 
  SHARELIB =
endif

all: wrappers

wrappers: 
	$(CXX) $(CXXSOURCES) $(CDNLFLAG) $(CPICFLAG) $(CXXFLAGS) 
	$(AR) rcs $(AUXDIR)libwaverange.a $(OBJECTC)
	$(SHARELIB)
.PHONY: clean
clean:
	$(RM) ./*.gc??
	$(RM) ./*.o
