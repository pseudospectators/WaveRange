FC = f95
OUTPUTDIR = ./
# LIBDIR is the path to the WaveRange library file libwaverange.a or libwaverange.so
LIBDIR = ../../bin/lib/
MKDIR = mkdir -p  
#FFLAGS = -Wpadded -Wpacked -malign-double -mpreferred-stack-boundary=8 -fno-underscoring
FFLAGS = -fno-underscoring

all:
	$(MKDIR) $(OUTPUTDIR)
	$(FC) example_fort.f90 $(FFLAGS) -L$(LIBDIR) -lwaverange -lstdc++ -o $(OUTPUTDIR)example_fort
.PHONY: clean
clean:
	$(RM) ./*.gc??
	$(RM) ./*.o
