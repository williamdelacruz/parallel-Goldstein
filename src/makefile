IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR)

ODIR=obj
LDIR =../lib

_DEPS = pi.h util.h file.h extract.h grad.h

DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o util.o extract.o grad.o

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -fopenmp -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	gcc -fopenmp -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
