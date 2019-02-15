include ../../config.mk

OBJECTC = waveletcdf97_3d.o
CSOURCES = waveletcdf97_3d.c

all: waveletcdf

waveletcdf: 
	$(CC) $(CSOURCES) $(CDNLFLAG) $(CFLAGS) 
.PHONY: clean
clean:
	$(RM) ./*.gc??
	$(RM) ./*.o
