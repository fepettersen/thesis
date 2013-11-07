#ifndef WALK_H
#define WALK_H

class Walk
{
	public:
		Walk(int d);

		int nwalkers, steps;
		double **area;
		double **walkers;	//Array of the positions of all the walkers [[0.1,0.2],[0.2,0.3],...]
		int d,m,n;
		double factor,dt;
		double x0, x1, y0, y1, z0, z1;
		double x0_, x1_, y0_, y1_; 		/*x0_ = x0 - (dx/2.0); etc*/
		double *x;
		double *y;
		double dx,dy,drift;
		// long Idum;
		

		bool HasLeftArea(double *);
		void advance(int **);
		void InhomogenousAdvance(int **C, double **D, double _dt);
		double *Step(double*, double*);
		double *InhomogenousStep(double*, double*);
		int InitializeTimestep(int **);
		void PutWalkers(int, int, int);
		double **ReturnBoundary();			//should have a better name
		int *FindPosition(double *);
		double **CalculateGradient();
		double *checkpos(double*,double*);
		void SetInitialCondition(int **, int, int);
		void ResetInitialCondition(int **);
		void ResetWalkers(){
			for(int k=0;k<nwalkers;k++)
				for(int l=0; l<d; l++)
					walkers[k][l] = 0;
		};
		// double ran0(long*);

};

#endif // WALK_H