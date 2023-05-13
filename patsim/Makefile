LFLAGS = -lrt -lm
CPPFLAGS := -O3

DEPS := $(wildcard physics.c)

parallel-c:main_parallel.c  physics.c
	mpicc -std=c99 -g -o $@ $^ $(LFLAGS)

sequential-c:main_sequential.c  physics.c
	$(CC) -std=c99 -g -o $@ $^ $(LFLAGS)

parallel: main_parallel.cpp
	mpic++ -std=c++14 -Wall -o $@ $^ $(DEPS) $(CPPFLAGS) $(LFLAGS)

sequential:main_sequential.cpp 
	$(CXX) -std=c++14  -Wall -o $@ $^ $(DEPS) $(CPPFLAGS) $(LFLAGS)

clean:
	-$(RM) parallel sequential parallel-c sequential-c
