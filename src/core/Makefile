include ../../config.mk

OUTPUTDIR = ../../libc/

OBJECTC = wrappers.o
CXXSOURCES = wrappers.cpp

all: wrappers

wrappers: 
	$(MKDIR)
	$(CXX) $(CXXSOURCES) $(CXXFLAGS) 
	$(AR) rcs $(OUTPUTDIR)libwaverange.a $(OBJECTC)
.PHONY: clean
clean:
	$(RM) ./*.gc??
	$(RM) ./*.o
