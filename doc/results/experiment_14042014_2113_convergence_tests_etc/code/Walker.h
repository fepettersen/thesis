#ifndef WALKER_H
#define WALKER_H
class Walker
{
	/*Walker object containing ony a 2d position so it can be removed easier*/
public:
	Walker(double x,double y){
		r = new double[2];
		r[0] = x; r[1] = y;};
	~Walker();

	/* data */
	double* r;
};
#endif // WALKER_H