#include "main_walk.h"
using namespace std;


Combine::Combine(int M, int N, double X0, double X1, double Y0, double Y1,double DiffusionConstant)
{
	m = M; n = N;
	x0 = X0; x1 = X1;
	y0 = Y0; y1 = Y1;
	D = DiffusionConstant;
	walk_areas = 0;
	dx = (x1-x0)/(m-1);
	dy = (y1-y0)/(n-1);
	pde_solver = new Diffusion(dx,dy,D);
	int **tmp = new int*[n];
	for(int i=0; i<n; i++){
		tmp[i] = new int[n];
	}
	C = tmp;

};

void Combine::Solve(){
	int counter = 0;
	pde_solver.advance();
	for(vector<Walk*>::iterator it1 = walk_solvers.begin(); it1 != walk_solvers.end(); it1++){
		ConvertToWalkers(C,counter);
		it1->advance();
		ConvertFromWalkers(C,counter);
	}
}

void Combine::AddWalkArea(int d){
	walk_areas++;
	Walk tmp(d);
	tmp.SetInitialCondition(C,M,N);		/*The solution in the relevant area converted to walkers*/
	walk_solvers.push_back(*tmp);
}

void Combine::ConvertToWalkers(int**C,int i){

}

void Combine::ConvertFromWalkers(int**C,int i){

	for(int k=0; k<len(C); k++){
		for(int j=0; j<len(C[i]); j++){
			/*This is where we insert least squares or similar*/
			U[k+indeces[i]][j+indeces[i]] += C[i][j]/(Hc*M);
		}
	}
}

void Combine::MapAreaToIndex(){
	// indeces = [[],[]]
	// if self.d ==1:
	// 	xcoor = [area[0][0],area[0][1]]
	// 	for i in xcoor:
	// 		for j in xrange(len(self.mesh[0])):
	// 			if i-self.mesh[0][j]<eps:
	// 				break
	// 		indeces[0].append(j)
	// 	indeces[-1] = 0	
	// elif self.d==2:
	// 	xcoor = [area[0][0],area[1][0]]
	// 	ycoor = [area[0][1],area[1][1]]
	// 	for i in xcoor:
	// 		for j in xrange(len(self.mesh[0])):
	// 			if i-self.mesh[0][j]<eps:
	// 				break
	// 		indeces[0].append(j)
	// 	for i in ycoor:
	// 		for j in xrange(len(self.mesh[-1])):
	// 			if i-self.mesh[-1][j]<eps:
	// 				break
	// 		indeces[-1].append(j)
}