#ifndef WALK_H
#define WALK_H

class Walk
{
	public:
		Walk(int d,double _dt=0.01);

		int nwalkers, steps;
		double **area, **aD;
		double **walkers;				/*Array of the positions of all the walkers*/
		// std::vector<double*> walkers;
		std::vector<double**> aDx;
		int d,m,n;
		double factor,dt,D;
		double x0, x1, y0, y1, z0, z1;
		double x0_, x1_, y0_, y1_; 		/*x0_ = x0 - (dx/2.0); etc*/
		double *x;
		double *y;
		double dx,dy,drift;
		bool inhomogenous;
		Random *rng;					/*Random Number Generator*/
		

		
		void FindPosition(double *, int*);
		
		void checkpos(double*);

		void SetDiffusionTensor(double **, int, int);
		void SetDiffusionConstant(double);

		void Advance(void);
		void Load(std::string,int,int);


};

#endif // WALK_H