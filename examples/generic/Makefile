FC = f95
OUTPUTDIR = ./
MKDIR = mkdir -p  
#FFLAGS = -Wpadded -Wpacked -malign-double -mpreferred-stack-boundary=8 -fno-underscoring
FFLAGS = -fno-underscoring

all:
	$(MKDIR) $(OUTPUTDIR)
	$(FC) create_in_field.f90 $(FFLAGS) -o $(OUTPUTDIR)create_in_field
.PHONY: clean
clean:
	$(RM) ./*.o
	$(RM) create_in_field
