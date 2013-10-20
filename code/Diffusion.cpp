#include "main_walk.h"
using namespace std;



Diffusion::Diffusion(double dx, double dy, double D){
	dx = dx;
	dy = dy;
	D = D;
	d = (dy>0)?2:1;
	// cout<<"dx,dy = "<<dx<<","<<dy<<endl;
	if(d==2){
		dt = dx*dy/5.0;
		_D = D*dt/(dx*dy);
	}
	else if(d==1){
		dt = dx*dx/3.0;
		_D = D*dt/(dx*dx);
	}
	solver = 0;
};

void Diffusion::advance(double **U,double **Up, int m, int n){
	if (solver==0){
		/*Forward Euler*/
		if(d==2){
			for(int i=1; i<(m-1);i++){
				for(int j=1; j<(n-1); j++){
					U[i][j] = _D*(Up[i+1][j]-2*Up[i][j] +Up[i-1][j]) +
					_D*(Up[i][j+1]-2*Up[i][j]+Up[i][j-1]) + Up[i][j];
				}
				boundary(U,Up,m,n);
			}
		}
		else if(d==1){
			for(int i=1; i<(m-1);i++){
				for(int j=0; j<1; j++){
					U[i][j] = _D*(Up[i+1][j]-2*Up[i][j] +Up[i-1][j]) + Up[i][j];
				}
				boundary(U,Up,m,n);
			}
		}
	}
}

void Diffusion::boundary(double **U,double **Up,int m, int n){
	if(d==2){
		for(int i=1; i<(m-1); i++){
			U[0][i] = 2*_D*(Up[1][i]-Up[0][i]) + 
			_D*(Up[0][i-1]-2*Up[0][i] + Up[0][i+1]) + Up[0][i];
			U[m-1][i] = 2*_D*(Up[m-2][i]-Up[m-1][i]) + 
			_D*(Up[m-1][i-1]-2*Up[m-1][i] + Up[m-1][i+1]) + Up[m-1][i];
		}
		for(int j=1; j<(n-1); j++){
			U[j][0] = _D*(Up[j-1][0]-2*Up[j][0]+Up[j+1][0]) + 
			2*_D*(Up[j][1]-Up[j][0]) +Up[j][0];
			U[j][n-1] = _D*(Up[j-1][n-1]-2*Up[j][n-1]+Up[j+1][n-1]) + 
			2*_D*(Up[j][n-2]-Up[j][n-1]) +Up[j][n-1];
		}
		U[0][0] = 2*_D*(Up[1][0]-Up[0][0]) +Up[0][0] + 
		2*_D*(Up[0][1]-Up[0][0]);

		U[m-1][0] = 2*_D*(Up[m-2][0]-Up[m-1][0]) +Up[m-1][0] + 
		2*_D*(Up[m-1][1]-Up[m-1][0]);
		
		U[0][n-1] = 2*_D*(Up[1][n-1]-Up[0][n-1]) + Up[0][n-1] + 
		2*_D*(Up[0][n-2]-Up[0][n-1]);
		
		U[m-1][n-1] = 2*_D*(Up[m-2][n-1]-Up[m-1][n-1]) +Up[m-1][n-1] +
		 2*_D*(Up[m-1][n-2]-Up[m-1][n-1]);
	}
	else if(d==1){
		U[0][0] =  2*_D*(Up[1][0]-Up[0][0]) +Up[0][0];
		U[m-1][0] = 2*_D*(Up[m-2][0]-Up[m-1][0]) +Up[m-1][0];
	}
}