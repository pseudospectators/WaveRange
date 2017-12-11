include ../../config.mk

OUTPUTDIR = ../../libc/

OBJECTC = waveletcdf97_3d.o
CSOURCES = waveletcdf97_3d.c

all: waveletcdf

waveletcdf: 
	$(MKDIR)
	$(CC) $(CSOURCES) $(CFLAGS) 
	$(AR) rcs $(OUTPUTDIR)libwcdf.a $(OBJECTC)
.PHONY: clean
clean:
	$(RM) ./*.gc??
	$(RM) ./*.o
