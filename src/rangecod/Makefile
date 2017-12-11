include ../../config.mk

OUTPUTDIR = ../../libc/

OBJECTC = rangecod.o 
CSOURCES = rangecod.c

all: rangecoder

rangecoder: 
	$(CC) $(CSOURCES) $(CFLAGS) 
	$(AR) rcs $(OUTPUTDIR)libcoders.a $(OBJECTC)
.PHONY: clean
clean:
	$(RM) ./*.gc??
	$(RM) ./*.o
