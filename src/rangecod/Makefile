include ../../config.mk

OBJECTC = rangecod.o 
CSOURCES = rangecod.c

all: rangecoder

rangecoder: 
	$(CC) $(CSOURCES) $(CDNLFLAG) $(CFLAGS) 
.PHONY: clean
clean:
	$(RM) ./*.gc??
	$(RM) ./*.o
