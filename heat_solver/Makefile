laplsolv-f: laplsolv.f90
	ifort $(CFLAGS) -fopenmp $^ -o $@

laplsolv: laplsolv.c
	$(CC) $(CFLAGS) -std=c11 -fopenmp $^ -o $@ -lm
	
laplsolv-back: laplsolv.back.c
	$(CC) $(CFLAGS) -std=c11 -fopenmp $^ -o $@ -lm

clean:
	-$(RM) laplsolv laplsolv-f