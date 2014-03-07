# makefile for project

SRCS = Combine.cpp Diffusion.cpp Tridiag.cpp Random.cpp Dendrite.cpp main_walk.cpp
OBJS = Combine.o Diffusion.o Tridiag.o Random.o Dendrite.o main_walk.o
CC = g++
DEBUG = -g
CFLAGS = -O2 
LFLAGS = -larmadillo -llapack -lblas $(DEBUG)

main_walk: libs
	$(CC) $(OBJS) $(LFLAGS) -o main_walk

libs:
	$(CC) -c $(SRCS) $(CFLAGS) $(LFLAGS) 

main_walk.o: main_walk.h main_walk.cpp
	$(CC) -c $(CFLAGS) $(LFLAGS) main_walk.cpp 

Combine.o: Combine.h Combine.cpp
	$(CC) $(CFLAGS) $(LFLAGS) Combine.cpp

Diffusion.o: Diffusion.h Diffusion.cpp
	$(CC) $(CFLAGS) $(LFLAGS) Diffusion.cpp

Tridiag.o: Tridiag.h Tridiag.cpp
	$(CC) $(CFLAGS) $(LFLAGS) Tridiag.cpp

Random.o: Random.h Random.cpp
	$(CC) $(CFLAGS) Random.cpp

Dendrite.o: main_walk.h Combine.h Dendrite.h Dendrite.cpp
	$(CC) $(CFLAGS) $(LFLAGS) Dendrite.cpp

clean:
	rm $(OBJS) main_walk