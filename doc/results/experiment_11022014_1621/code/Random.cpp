#include "main_walk.h"
using namespace std;

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

Random::Random(void){
	srand (time(NULL));
	x = rand(); y = rand(); z = rand();
	w = rand(); v = rand();
	idum = -1*rand();
}

double Random::uniform(void){
	/*This is a modified xorshift prng from George Marsaglia
	with 5 random seeds. The period of a general xorshift algorithm 
	is approximatiely 2^(32k) where k is the number of seeds. This one 
	should have a period of around 2^160*/
	unsigned long t;
 	t=(x^(x>>7)); 
 	x=y; y=z; z=w; w=v;
 	v=(v^(v<<6))^(t^(t<<13)); 
 	return double((y+y+1)*v)/ULLONG_MAX;
}



double Random::ran0(void)
{
   long     k;
   double   ans;

   idum ^= MASK;
   k = (idum)/IQ;
   idum = IA*(idum - k*IQ) - IR*k;
   if(idum < 0) idum += IM;
   ans=AM*(idum);
   idum ^= MASK;
   return ans;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK