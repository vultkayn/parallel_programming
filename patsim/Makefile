parallel:main_parallel.c  physics.c
	mpicc -std=c99 -g -o $@ $^ -lrt -lm

sequential:main_sequential.c  physics.c
	gcc -std=c99 -g -o $@ $^ -lrt -lm

clean:
	-$(RM) parallel sequential
