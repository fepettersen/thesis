#include <iostream>
#include <cmath>
#include <stdlib.h>

class Walk
{
	public:
		Walk(int d);

		int nwalkers;
		double **area;
		double **walkers;	//Array of the positions of all the walkers [[0.1,0.2],[0.2,0.3],...]
		int d,m,n;
		double factor,dt;
		double x0, x1, y0, y1, z0, z1;
		double x0_, x1_, y0_, y1_; 		/*x0_ = x0 - (dx/2.0); etc*/
		double *x;
		double *y;
		double dx,dy;
		long idum;
		

		bool HasLeftArea(double *);
		int **advance(int **);
		double *Step(double*, double*);
		int InitializeTimestep(int **);
		void PutWalkers(int, int, int);
		double **ReturnBoundary();	//should have a better name
		int *FindPosition(double *);
		double **CalculateGradient();
		double *checkpos(double*,double*);
		void SetInitialCondition(int **, int, int);

};

// #endif // WALK_H