#include "Random.h"
using namespace std;

Random::Random(void){
	srand (time(NULL));
	x = rand(); y = rand(); z = rand();
	w = rand(); v = rand();
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