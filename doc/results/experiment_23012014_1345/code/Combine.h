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
		void ConvertToWalkers(double **, int**, int **);
		void ConvertFromWalkers(double **, int**, int **);
		void MapAreaToIndex(double *, double *, int **);
		void SetInitialCondition(double**,int,int);
		void CubicSpline(double*,double*,int,double,double);
		double abs_max(double **array,int m, int n);
		double norm(double **U,double **Up,int m, int n);
		void TestRWConvergence(int steps,std::string path);

		double dx, dy, D;
		double *X,*Y;
		double **aD;			/*anisotropic diffusion*/
		double x0, x1, y0, y1;
		int m,n,d;				/*resolution in x and y dir. and dim. n=0 => d=1*/
		double **U, **Up;		/*Actual solution*/
		int **C;				/*Distribution of walkers -- replaced?*/
		double Hc;			/*conversion factors*/
		int walk_areas;
		int **signmap;
		bool inhomogenous;

		std::vector<Walk*> walk_solvers; 	/*A linked list of the walk-solvers*/
		std::vector<int**> c;				/*Linked list of walker-distr. for area i*/
		std::vector<int**> indeces;			/*Indeces of walk-area i*/

		Diffusion *pde_solver;
};


#endif // COMBINE_H