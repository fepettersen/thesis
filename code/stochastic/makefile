OBJS = Random.o Walk.o walksolver.o
CC = g++
DEBUG = -g
CFLAGS = -O3 $(DEBUG)

walk_solver: alles
	$(CC) $(OBJS) -o walk_solver

alles:
	$(CC) -c *.cpp $(CFLAGS)

Walk.o: Walk.h Walk.cpp
	$(CC) -c Walk.cpp $(CFLAGS)

Random.o: Random.h Random.cpp
	$(CC) -c Random.cpp $(CFLAGS)

walksolver.o: walksolver.h walksolver.cpp
	$(CC) -c walksolver.cpp $(CFLAGS)

clean:
	rm *.o