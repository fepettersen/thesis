#include "main_walk.h"
#include <algorithm>
using namespace std;



Diffusion::Diffusion(double _dx, double _dy, double _D, double Dt){
	dx = _dx;
	dy = _dy;
	D = _D;
	d = (dy>0)?2:1;
	solver = 0;
	if(d==2 && solver==0){
		dt = (Dt>(dx*dy/4.0))? (dx*dy/(5.0)):Dt;
		_Dx = D*dt/(dx*dx);
		_Dy = D*dt/(dy*dy);
	}
	else if(d==1 && solver==0){
		dt = (Dt>(dx*dx/2.0))? (dx*dx/(3.0)):Dt;
		_Dx = D*dt/(dx*dx);
	}
	t=0;
};
Diffusion::Diffusion(double _dx, double _dy, double **D, double Dt){
	dx = _dx;
	dy = _dy;
	dt = Dt;
	_Dx = dt/(dx*dx);
	aD = D;
	d = (dy>0)?2:1;
	solver = 2;
	if(d==2 && solver==1){
		dt = (Dt>(dx*dy/4.0))? (dx*dy/(5.0)):Dt;
		// _Dx = D*dt/(dx*dx);
		// _Dy = D*dt/(dy*dy);
	}
	else if(d==1 && solver==1){
		double Dtmp = 3.1415926535897932;			/*Should be max(aD)*/
		dt = (Dt>(dx*dx/(2.0*Dtmp)))? (dx*dx/(10.0*Dtmp)):Dt;
		// dt = (dx*dx/(10.0*Dtmp));
		cout<<dt<<setprecision(16)<<endl;
		_Dx = dt/(dx*dx);
	}
	t=0;
}

void Diffusion::advance(double **U,double **Up, int m, int n){
	if (solver==0){
		/*Forward Euler*/
		if(d==2){
			for(int i=1; i<(m-1);i++){
				for(int j=1; j<(n-1); j++){
					U[i][j] = _Dx*(Up[i+1][j]-2*Up[i][j] +Up[i-1][j]) +
					_Dy*(Up[i][j+1]-2*Up[i][j]+Up[i][j-1]) + Up[i][j];
				}
			}
			boundary(U,Up,m,n);
		}
		else if(d==1){
			for(int i=1; i<(m-1);i++){
				for(int j=0; j<1; j++){
					U[i][j] = _Dx*(Up[i+1][j]-2*Up[i][j] +Up[i-1][j]) + Up[i][j];
				}
			}
			boundary(U,Up,m,n);
		}
	}
	else if(solver==1){
		/*Anisotropy - Forward Euler*/
		if(d==1){
			// cout<<aD[m-1][0]<<endl;
			for(int i=1;i<(m-1);i++)
				for(int j=0;j<1;j++){
					U[i][j] = (_Dx/2.0)*((aD[i+1][j]+aD[i][j])*(Up[i+1][j]-Up[i][j]) -
						(aD[i][j]+aD[i-1][j])*(Up[i][j]-Up[i-1][j])) + Up[i][j] + dt*f(i*dx,0,t*dt);
				}
			boundary(U,Up,m,n);
		}
	}
	else{
		/*Anisotropic convection diffusion - Forward Euler*/
		if(d==1){
			double v = 1;
			for(int i=1;i<(m-1);i++){
				for(int j=0;j<1;j++){
					U[i][j] = (_Dx/2.0)*((aD[i+1][j]+aD[i][j])*(Up[i+1][j]-Up[i][j]) -
						(aD[i][j]+aD[i-1][j])*(Up[i][j]-Up[i-1][j])) +Up[i][j] - 
					((dt*v)/(2.0*dx))*(Up[i+1][j]-Up[i-1][j]) + dt*f(i*dx,0,t*dt);
					// U[i][j] = _Dx*(Up[i+1][j]-2*Up[i][j]+Up[i-1][j]) +Up[i][j] - 
					// (dt*v)/(2.0*dx)*(Up[i+1][j]-Up[i-1][j]) + dt*f(i*dx,0,t*dt);
				}
			}
			boundary(U,Up,m,n);
		}
	}
	t++;
}
// U[1:-1] = factor*(Up[2:]-2*Up[1:-1]+Up[:-2]) + Up[1:-1] - (dt*v)/(2*dx)*(Up[2:]-Up[:-2]) +dt*F
void Diffusion::boundary(double **U,double **Up,int m, int n){
	if(solver==0){
		if(d==2){
			for(int j=1; j<(m-1); j++){
				U[j][0] = _Dx*(Up[j+1][0]-2*Up[j][0]+Up[j-1][0]) + 
				2*_Dy*(Up[j][1]-Up[j][0]) +Up[j][0];
				U[j][n-1] = _Dx*(Up[j+1][n-1]-2*Up[j][n-1]+Up[j-1][n-1]) + 
				2*_Dy*(Up[j][n-2]-Up[j][n-1]) +Up[j][n-1];
			}
			for(int i=1; i<(n-1); i++){
				U[0][i] = 2*_Dx*(Up[1][i]-Up[0][i]) + 
				_Dy*(Up[0][i+1]-2*Up[0][i] + Up[0][i-1]) + Up[0][i];
				U[m-1][i] = 2*_Dx*(Up[m-2][i]-Up[m-1][i]) + 
				_Dy*(Up[m-1][i+1]-2*Up[m-1][i] + Up[m-1][i-1]) + Up[m-1][i];
			}
			U[0][0] = 2*_Dx*(Up[1][0]-Up[0][0]) +Up[0][0] + 
			2*_Dy*(Up[0][1]-Up[0][0]);

			U[m-1][0] = 2*_Dx*(Up[m-2][0]-Up[m-1][0]) +Up[m-1][0] + 
			2*_Dy*(Up[m-1][1]-Up[m-1][0]);
			
			U[0][n-1] = 2*_Dx*(Up[1][n-1]-Up[0][n-1]) + Up[0][n-1] + 
			2*_Dy*(Up[0][n-2]-Up[0][n-1]);
			
			U[m-1][n-1] = 2*_Dx*(Up[m-2][n-1]-Up[m-1][n-1]) +Up[m-1][n-1] +
			 2*_Dy*(Up[m-1][n-2]-Up[m-1][n-1]);
		}
		else if(d==1){
			U[0][0] =  2*_Dx*(Up[1][0]-Up[0][0]) +Up[0][0];
			U[m-1][0] = 2*_Dx*(Up[m-2][0]-Up[m-1][0]) +Up[m-1][0];
		}
	}
	else if(solver==1 || solver==2){
		if(d==1){
			U[0][0] =  _Dx*((aD[1][0]+aD[0][0])*(Up[1][0]-Up[0][0])) +Up[0][0]+dt*f(0,0,t*dt);
			U[m-1][0] = _Dx*((aD[m-2][0]+aD[m-1][0])*(Up[m-2][0]-Up[m-1][0])) +Up[m-1][0] +dt*f((m-1)*dx,0,t*dt);
		}
	}
}

double Diffusion::f(double x,double y, double t){
	double pi = 3.1415926535897932;
	// return 2*pi*sin(pi*x)*exp(-t*pi*pi);
	return -pi*sin(pi*x)*exp(-t*pi*pi);
}