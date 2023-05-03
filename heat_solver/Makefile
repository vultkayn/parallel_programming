laplsolv-f: laplsolv.f90
	ifort -fopenmp $^ -o $@

laplsolv: laplsolv.c
	$(CC) -std=c11 -fopenmp $^ -o $@
	
clean:
	-$(RM) laplsolv laplsolv-f