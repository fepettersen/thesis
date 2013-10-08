#ifndef COMBINE_H
#define COMBINE_H

class Combine
{
	public:
		Combine(int, int, double, double, double, double, double);

		void advance();
		void AddWalkArea(int );
		void Solve();
		void ConvertToWalkers(int**C,int i);
		void ConvertFromWalkers(int**C,int i);
		void MapAreaToIndex();

		double dx, dy, D;
		double x0, x1, y0, y1;
		int m,n,d;
		double **U, **Up;		/*Actual solution*/
		int **C;				/*Distribution of walkers -- Should 
									support more than one walk-solver*/
		int walk_areas;
		std::vector<Walk*> walk_solvers; 	/*A linked list of the walk-solvers*/
		Diffusion pde_solver(double, double, double);
};


#endif // COMBINE_H