class Walk
{
	public:
		Walk(int n, int d);
		// double Hc;
		int nwalkers;
		double **area;
		double **walkers;	//Array of the positions of all the walkers [[0.1,0.2],[0.2,0.3],...]
		int d;
		// double M;
		double factor;
		double x0, x1, y0, y1, z0, z1;
		
		bool HasLeftArea(double *);
		double **advance(double **);
		int InitializeTimestep(int **);
		void PutWalkers();
		double **ReturnBoundary();	//should have a better name
		double *FindPosition();
		double **CalculateGradient();
		double *checkpos();
		void SetInitialCondition(double **);

};

// #endif // WALK_H