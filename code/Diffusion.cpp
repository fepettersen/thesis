#include "Diffusion.h"
using namespace std;

Diffusion::Diffusion(double dx, double dy, double D){
	dx = dx;
	dy = dy;
	D = D;
	dt = dx*dy/5.0;
	_D = D*dt/(dx*dy);
};

double **Diffusion::advance(double **U,double **Up,int m, int n){
	if (solver==0){
		/*Forward Euler*/
		for(int i=1; i<m-1;i++){
			for(int j=1; j<n-1; j++){
				U[i][j] = _D*(Up[i+1][j]-2*Up[i][j] +Up[i-1][j]) +
				_D*(Up[i][j+1]-2*Up[i][j]+Up[i][j-1]) + Up[i][j];
			}
			U = boundary(U,Up,m,n);
		}
	}
}

double **Diffusion::boundary(double **U,double **Up,int m, int n){
	for(int i=1; i<(m-1); i++){
		U[0][i] = 2*_D*(Up[0][i]-Up[1][i]) + 
		_D*(Up[0][i-1]-2*Up[0][i] + Up[0][i+1]) + Up[0][i];
		U[m][i] = 2*_D*(Up[m][i]-Up[m-1][i]) + 
		_D*(Up[m][i-1]-2*Up[m][i] + Up[m][i+1]) + Up[m][i];
	}
	for(int j=0; j<(n-1); j++){
		U[j][0] = _D*(Up[j-1][0]-2*Up[j][0]+Up[j+1][0]) + 
		2*_D*(Up[j][0]-Up[j][1]) +Up[j][0];
		U[j][n] = _D*(Up[j-1][n]-2*Up[j][n]+Up[j+1][n]) + 
		2*_D*(Up[j][n]-Up[j][n-1]) +Up[j][n];
	}
	U[0][0] = 2*_D*(Up[1][0]-Up[0][0]) +Up[0][0] + 
	2*_D*(Up[0][1]-Up[0][0]);

	U[m][0] = 2*_D*(Up[m-1][0]-Up[m][0]) +Up[m][0] + 
	2*_D*(Up[m][1]-Up[m][0]);
	
	U[0][n] = 2*_D*(Up[1][n]-Up[0][n]) + Up[0][n] + 
	2*_D*(Up[0][n-1]-Up[0][n]);
	
	U[m][n] = 2*_D*(Up[m-1][n]-Up[m][n]) +Up[m][n] +
	 2*_D*(Up[m][n-1]-Up[m][n]);
}