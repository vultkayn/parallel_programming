laplsolv-f: laplsolv.f90
	ifort -fopenmp $^ -o $@

laplsolv: laplsolv.c
	$(CC) -std=c11 -fopenmp $^ -o $@ -lm
	
laplsolv-back: laplsolv.back.c
	$(CC) -std=c11 -fopenmp $^ -o $@ -lm

clean:
	-$(RM) laplsolv laplsolv-f