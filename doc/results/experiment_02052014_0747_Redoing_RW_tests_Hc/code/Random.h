#ifndef RANDOM_H
#define RANDOM_H

// #include <climits>
// #include <stdlib.h>
// #include <time.h>

class Random
{
public:
	Random(int seed = 1);
	double uniform(void);
	double ran0(void);
	/* data */
	unsigned long x,y,z,w,v;
	long idum;
};


#endif // RANDOM_H