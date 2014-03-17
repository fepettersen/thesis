# makefile for project

#SRCS = Combine.cpp Diffusion.cpp Tridiag.cpp Random.cpp Dendrite.cpp Spine.cpp Walker.h 
SRCS = *.cpp Walker.h
OBJS = Combine.o Diffusion.o Tridiag.o Random.o Dendrite.o main_walk.o
CC = g++
DEBUG = -g
CFLAGS = -O2 
LFLAGS = -larmadillo -llapack -lblas $(DEBUG)

main_walk: libraries
	$(CC) main_walk.cpp $(wildcard libraries/*) $(LFLAGS) -o main_walk

libraries:
	mkdir libraries/
	for file in $(SRCS); do \
		$(CC) -c $$file; \
		done
	mv *.o libraries/; \


clean:
	rm -rf libraries/; \
	rm main_walk

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
