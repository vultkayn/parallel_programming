CC = gcc
CFLAGS = 
LFLAGS = -lpthread -lrt -lm

all: blurc thresc

clean:
	-$(RM) *.o blurc thresc

blurc: ppmio.c gaussw.c blurfilter.c blurmain.c
	$(CC) -o $@ ppmio.c gaussw.c blurfilter.c blurmain.c $(LFLAGS) $(CFLAGS)

thresc: thresmain.c ppmio.c thresfilter.c
	$(CC) -o $@ thresmain.c ppmio.c thresfilter.c $(LFLAGS) $(CFLAGS)

arc:
	tar cf - *.c *.h Makefile data/* | gzip - > filters.tar.gz
