include ../../config.mk

OBJECTC = rangecod.o 
CSOURCES = rangecod.c

ifeq ($(CC),gcc)
  CPICFLAG = -fPIC
else
  CPICFLAG = 
endif

all: rangecoder

rangecoder: 
	$(CC) $(CSOURCES) $(CDNLFLAG) $(CPICFLAG) $(CFLAGS) 
.PHONY: clean
clean:
	$(RM) ./*.gc??
	$(RM) ./*.o
