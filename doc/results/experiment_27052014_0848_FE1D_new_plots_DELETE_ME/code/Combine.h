#ifndef COMBINE_H
#define COMBINE_H

class Combine
{
	public:
		Combine(int, int, double, double, double, double, double, double, double);
		Combine(int, int, double, double, double, double, double**, double, double);

		void advance();
		void AddWalkArea(double*, double*);
		void Solve();
		void ConvertToWalkers(double **, std::string, int **);
		void ConvertFromWalkers(double **, std::string, int **);
		void MapAreaToIndex(double *, double *, int **);
		void SetInitialCondition(double**,int,int);
		void polyreg(double *x,double *y,int m);
		void CubicSpline(double*,double*,double*,int,double,double,int);
		double abs_max(double **array,int m, int n);
		double norm(double **U,double **Up,int m, int n);
		void TestRWConvergence(int steps,std::string path);
		void SaveDiffusionTensor(double**,int, int, int);

		double dx, dy, D;
		double *X,*Y;
		double **aD;			/*anisotropic diffusion*/
		double x0, x1, y0, y1;
		int m,n,d;				/*resolution in x and y dir. and dim. n=0 => d=1*/
		double **U, **Up;		/*Actual solution*/
		double Hc;				/*conversion factors*/
		int walk_areas, walk_steps;
		int **signmap;
		bool inhomogenous;
		std::string prgm;

		std::vector<int**> indeces;				/*Indeces of walk-area i*/
		std::vector<std::string> inifilenames;
		std::vector<int*> parameters;

		Diffusion *pde_solver;
		Random* rng;
};


#endif // COMBINE_H