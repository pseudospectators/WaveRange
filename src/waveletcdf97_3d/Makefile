include ../../config.mk

OBJECTC = waveletcdf97_3d.o
CSOURCES = waveletcdf97_3d.c

ifeq ($(CC),gcc)
  CPICFLAG = -fPIC
else
  CPICFLAG =
endif

all: waveletcdf

waveletcdf: 
	$(CC) $(CSOURCES) $(CDNLFLAG) $(CPICFLAG) $(CFLAGS) 
.PHONY: clean
clean:
	$(RM) ./*.gc??
	$(RM) ./*.o
