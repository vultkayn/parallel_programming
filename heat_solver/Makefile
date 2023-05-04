CC = gcc
LFLAGS = -fopenmp -lm

laplsolv-f: laplsolv.f90
	ifort -fopenmp $^ -o $@

laplsolv: laplsolv.c
	$(CC) -std=c11 $^ -o $@ $(LFLAGS) 
	
laplsolv-back: laplsolv.back.c
	$(CC) -std=c11 $^ -o $@ $(LFLAGS) 

clean:
	-$(RM) laplsolv laplsolv-f laplsolv-back