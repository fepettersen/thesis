#ifndef WALK_H
#define WALK_H

class Walk
{
	public:
		Walk(int d,double _dt=0.01);

		int nwalkers, steps;
		double **area, **aD;
		// double **walkers;	//Array of the positions of all the walkers
		std::vector<double*> walkers;
		int d,m,n;
		double factor,dt,D;
		double x0, x1, y0, y1, z0, z1;
		double x0_, x1_, y0_, y1_; 		/*x0_ = x0 - (dx/2.0); etc*/
		double *x;
		double *y;
		double dx,dy,drift;
		bool inhomogenous;
		

		bool HasLeftArea(double *);
		void advance(int **);
		void InhomogenousAdvance(int **C, double _dt);
		double *Step(double*, double*);
		double *InhomogenousStep(double*, double*);
		int InitializeTimestep(int **);
		void PutWalkers(int, int, int);
		void FindPosition(double *, int*);
		double *checkpos(double*,double*);
		void SetInitialCondition(int **, int, int);
		void ResetInitialCondition(int **);
		void SetDiffusionTensor(double **, int, int);
		void SetDiffusionConstant(double);

		void ResetWalkers(){
			walkers.clear();
			// for(int k=0;k<nwalkers;k++)
			// 	delete [] walkers[k];
				// for(int l=0; l<d; l++)
				// 	walkers[k][l] = 0;
			// delete [] walkers;
		};
		// double ran0(long*);

};

#endif // WALK_H