CC = mpicc
CFLAGS = 
LFLAGS = -lrt -lm

all: blurc thresc

clean:
	-$(RM) *.o blurc thresc

blurc: ppmio.o gaussw.o blurfilter.o blurmain.o
	$(CC) -o $@ ppmio.o gaussw.o blurfilter.o blurmain.o $(LFLAGS) $(CFLAGS)

thresc:: thresmain.o ppmio.o thresfilter.o
	$(CC) -o $@ thresmain.o ppmio.o thresfilter.o $(LFLAGS) $(CFLAGS)

arc:
	tar cf - *.c *.h Makefile data/* | gzip - > filters.tar.gz
