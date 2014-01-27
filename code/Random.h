#ifndef RANDOM_H
#define RANDOM_H

#include <climits>
#include <stdlib.h>
#include <time.h>

class Random
{
public:
	Random(void);
	double uniform(void);
	/* data */
	unsigned long x,y,z,w,v;
};
#endif // RANDOM_H