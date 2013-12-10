#include "main_walk.h"
#include <algorithm>
using namespace std;



Diffusion::Diffusion(double _dx, double _dy, double _D, double Dt, double _v){
	dx = _dx;
	dy = _dy;
	D = _D;
	d = (dy>0)?2:1;
	solver = 0;		/*solver=3 ==> Backward Euler; solver=0 ==> Forward Euler*/
	if(d==2 && solver==0){
		dt = (Dt>(dx*dy/4.0))? (dx*dy/(5.0)):Dt;
		_Dx = D*dt/(dx*dx);
		_Dy = D*dt/(dy*dy);
	}
	else if(d==1 && solver==0){
		dt = (Dt>(dx*dx/2.0))? (dx*dx/(3.0)):Dt;
		_Dx = D*dt/(dx*dx);
	}
	else{
		dt = Dt;
	}
	t=1;
	v = _v;
	cout<<"dt,dx,dy = "<<dt<<","<<dx<<","<<dy<<	endl;
};

Diffusion::Diffusion(double _dx, double _dy, double **D, double Dt, double _v){
	dx = _dx;
	dy = _dy;
	dt = Dt;
	_Dx = dt/(dx*dx);
	aD = D;
	d = (dy>0)?2:1;
	solver = 2;
	if(d==2){
		dt = (Dt>(dx*dy/4.0))? (dx*dy/(5.0)):Dt;
		_Dx = dt/(dx*dx);
		_Dy = dt/(dy*dy);
	}
	else if(d==1 && solver==1){
		double Dtmp = 3.1415926535897932;			/*Should be max(aD)*/
		dt = (Dt>(dx*dx/(2.0*Dtmp)))? (dx*dx/(10.0*Dtmp)):Dt;
		// dt = (dx*dx/(10.0*Dtmp));
		cout<<dt<<setprecision(16)<<endl;
		_Dx = dt/(dx*dx);
	}
	t=1;
	v = _v;
	vdtdx2 = (v*dt)/(2*dx);
	vdtdy2 = (v*dt)/(2*dy);
	cout<<"dt,vdtdx2,vdtdy2 = "<<dt<<","<<vdtdx2<<","<<vdtdy2<<endl;
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
			for(int i=1;i<(m-1);i++)
				for(int j=0;j<1;j++){
					U[i][j] = (_Dx/2.0)*((aD[i+1][j]+aD[i][j])*(Up[i+1][j]-Up[i][j]) -
						(aD[i][j]+aD[i-1][j])*(Up[i][j]-Up[i-1][j])) + Up[i][j] + dt*f(i*dx,0,t*dt);
				}
			boundary(U,Up,m,n);
		}
	}
	else if(solver==3){
		/*Isotropic Bacward Euler discretization (1d)*/
		if(d==1){
			double a = -dt/(dx*dx);
			double c = a;
			double b = 1.0-2.0*a;
		    double *bv = new double[m];
		    for(int i=0;i<m;i++){
		    	bv[i] = 0;
		    }
		    double temp = 0;

		    bv[0]=2*c/b;
		    U[0][0] = Up[0][0]/b;
		    for(int i=1; i<m-1; i++){
		        //forward substitution
		       	temp = 1.0/(b-a*bv[i-1]);
		        bv[i] = c*temp;
		        U[i][0] = (Up[i][0] -a*U[i-1][0])*temp;
		    }
		    /*Boundary condition (Neumann)*/
	       	temp = 1.0/(1-c-a*bv[m-2]);
	        bv[m-1] = c*temp;
	        U[m-1][0] = (Up[m-1][0] -a*U[m-2][0])*temp;
		    for(int i=(m-2);i>0;i--){
		        //Backward substitution 
		        U[i][0] -= bv[i]*U[i+1][0];
		    }
		    boundary(U,Up,m,n);
		    delete [] bv;
		}
		else if(d==2){
			BE2D(U,Up,m,n);
			// boundary(U,Up,m,n);
		}
	}
	else{
		/*Anisotropic convection diffusion - Forward Euler*/
		if(d==1){
			// double v = 1;
			for(int i=1;i<(m-1);i++){
				for(int j=0;j<1;j++){
					U[i][j] = (_Dx/2.0)*((aD[i+1][j]+aD[i][j])*(Up[i+1][j]-Up[i][j]) -
						(aD[i][j]+aD[i-1][j])*(Up[i][j]-Up[i-1][j])) +Up[i][j] - 
					((dt*v)/(2.0*dx))*(Up[i+1][j]-Up[i-1][j]) + dt*f(i*dx,0,t*dt);
				}
			}
			boundary(U,Up,m,n);
		}
		else if(d==2){
			for(int i=1;i<(m-1);i++){
				for(int j=1;j<(n-1);j++){
					U[i][j] = (_Dx/2.0)*((aD[i+1][j]+aD[i][j])*(Up[i+1][j]-Up[i][j]) -(aD[i][j]+aD[i-1][j])*(Up[i][j]-Up[i-1][j])) +
					 (_Dy/2.0)*((aD[i][j+1]+aD[i][j])*(Up[i][j+1]-Up[i][j]) - (aD[i][j]+aD[i][j-1])*(Up[i][j]-Up[i][j-1])) + 
					 vdtdx2*(Up[i+1][j]-Up[i-1][j]) + vdtdy2*(Up[i][j+1]-Up[i][j-1])+Up[i][j] + dt*f(i*dx,j*dy,t*dt);
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
		else if(d==2){
			for(int j=1; j<(m-1); j++){
				U[j][0] = (_Dx/2.0)*((aD[j+1][0]+aD[j][0])*(Up[j+1][0]-Up[j][0])-(aD[j][0]+aD[j-1][0])*(Up[j][0]-Up[j-1][0])) + 
				_Dy*((aD[j][1]+aD[j][0])*(Up[j][1]-Up[j][0])) +Up[j][0] +vdtdx2*(Up[j+1][0]-Up[j-1][0])+dt*f(dx*j,0,t*dt);
				
				U[j][n-1] = (_Dx/2.0)*((aD[j+1][n-1]+aD[j][n-1])*(Up[j+1][n-1]-Up[j][n-1]) -(aD[j][n-1]+aD[j-1][n-1])*(Up[j][n-1]-Up[j-1][n-1])) + 
				_Dy*(aD[j][n-2]+aD[j][n-1])*(Up[j][n-2]-Up[j][n-1]) +Up[j][n-1]+vdtdx2*(Up[j+1][n-1]-Up[j-1][n-1])+dt*f(dx*j,(n-1)*dy,t*dt);
			}
			for(int i=1; i<(n-1); i++){
				U[0][i] = _Dx*((aD[1][i]+aD[0][i])*(Up[1][i]-Up[0][i])) + 
				(_Dy/2.0)*((aD[0][i+1]+aD[0][i])*(Up[0][i+1]-Up[0][i])-(aD[0][i]+aD[0][i-1])*(Up[0][i]-Up[0][i-1])) + 
				Up[0][i]+vdtdy2*(Up[0][i+1]-Up[0][i-1])+dt*f(0,i*dy,t*dt);
				
				U[m-1][i] = _Dx*((aD[m-2][i]+aD[m-1][i])*(Up[m-2][i]-Up[m-1][i])) + 
				(_Dy/2.0)*((aD[m-1][i+1]+aD[m-1][i])*(Up[m-1][i+1]-Up[m-1][i])-(aD[m-1][i]+aD[m-1][i-1])*(Up[m-1][i]-Up[m-1][i-1])) + 
				Up[m-1][i]+vdtdy2*(Up[m-1][i+1]-Up[m-1][i-1])+dt*f((m-1)*dx,i*dy,t*dt);
			}
			U[0][0] = _Dx*(aD[1][0]+aD[0][0])*(Up[1][0]-Up[0][0]) +Up[0][0] + 
			_Dy*(aD[0][1]+aD[0][0])*(Up[0][1]-Up[0][0])+dt*f(0,0,t*dt);

			U[m-1][0] = _Dx*(aD[m-2][0]+aD[m-1][0])*(Up[m-2][0]-Up[m-1][0]) + Up[m-1][0] + 
			_Dy*(aD[m-1][1]+aD[m-1][0])*(Up[m-1][1]-Up[m-1][0])+dt*f((m-1)*dx,0,t*dt);
			
			U[0][n-1] = _Dx*(aD[1][n-1]+aD[0][n-1])*(Up[1][n-1]-Up[0][n-1]) + Up[0][n-1] + 
			_Dy*(aD[0][n-2]+aD[0][n-1])*(Up[0][n-2]-Up[0][n-1])+dt*f(0,(n-1)*dy,t*dt);
			
			U[m-1][n-1] = _Dx*(aD[m-2][n-1]+aD[m-1][n-1])*(Up[m-2][n-1]-Up[m-1][n-1]) +Up[m-1][n-1] +
			 _Dy*(aD[m-1][n-2]+aD[m-1][n-1])*(Up[m-1][n-2]-Up[m-1][n-1])+dt*f((m-1)*dx,(n-1)*dy,t*dt);
			// cout<<aD[0][0]<<","<<aD[1][0]<<","<<aD[0][1]<<endl;
		}
	}
	else if(solver==3){
		/*Backward Euler*/
		if(d==1){
			double Q = (2.0*D*dt/(dx*dx));
			U[0][0] = (Q*U[1][0]+Up[0][0])/(1.0+Q);
			U[m-1][0] = (Q*U[m-2][0] +Up[m-1][0])/(1.0+Q);
			// U[0][0] = 1.0;
			// U[m-1][0] = -1.0;
		}
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
	}
}



double Diffusion::f(double x,double y, double t){
	double pi = 3.1415926535897932;
	// return exp(-t*pi*pi)*pi*pi*(sin(pi*x)+cos(pi*x)*(pi*x-1));
	// return -pi*sin(pi*x)*exp(-t*pi*pi);
	// double py = pi*y;
	// double px = pi*x;
	// return pi*exp(-t*pi*pi)*(2*pi*cos(px)*cos(py)*(x+y-0.5) + cos(px)*sin(py) +sin(px)*cos(py));
	return 0.0;
}

void Diffusion::tridiag(double *u, double *up, int N, double *di, double *ab, double *bel){
    double *temp = new double[N];
    for(int i=0;i<N;i++){
    	temp[i] = 0;
    }
    double btemp = di[0];

    // bv[0]=2*ab[1]/btemp;
    u[0] = up[0]/btemp;
    for(int i=1; i<N-1; i++){
        //forward substitution
       	// btemp = 1.0/(di[i]-bel[i]*bv[i-1]);
       	temp[i] = ab[i-1]/btemp;
       	btemp = di[i]-bel[i]*temp[i];
        // bv[i] = ab[i]*btemp;
        // u[i] = (up[i] -bel[i]*u[i-1])*btemp;
        u[i] = (up[i] -bel[i]*u[i-1])/btemp;
    }
    /*Boundary condition (Neumann)*/
   	// btemp = 1.0/(1-ab[N-2]-bel[N-2]*bv[N-2]);
    // bv[N-1] = ab[N-2]*btemp;
    // u[N-1] = (up[N-1] -bel[N-2]*u[N-2])*btemp;
    for(int i=(N-2);i>0;i--){
        //Backward substitution 
        // u[i] -= bv[i]*u[i+1];
        u[i] -= temp[i+1]*u[i+1];
    }
}

void Diffusion::BE2D(double **U, double **Up, int m, int n){
	double Utmp[m*n];
	double Uptmp[m*n];
	double diag[m*n];
	double lower[m*n-2];
	double upper[m*n-2];
	double alpha = D*dt/(2*dx*dx);
	double beta = D*dt/(2*dy*dy);
	int k=0;
	double b = (1+alpha);
	double b1 = (1-beta);

	for(int ii=0; ii<n;ii++){
		if(ii==0)
			Uptmp[k] = 2*beta*Up[0][1] + b1*Up[0][0];
		else if(ii=n-1)
			Uptmp[k] = 2*beta*Up[0][n-2] + b1*Up[0][n-1];
		else
			Uptmp[k] = beta*Up[0][ii+1] + b1*Up[0][ii] + beta*Up[0][ii-1];
		upper[k] = -2*alpha; 	
		lower[k] = 0;
		diag[k] = b;
		k++;
	}
	for(int i=1;i<(m-1);i++){
		Uptmp[k] = 2*beta*Up[i][1] + b1*Up[i][0];
		upper[k] = -2*alpha; 	
		lower[k] = 0;
		diag[k] = b;
		k++;
		for(int j=1;j<(n-1);j++){
			Uptmp[k] = beta*Up[i][j+1] + b1*Up[i][j] + beta*Up[i][j-1];
			upper[k] = -alpha; 	
			lower[k] = -alpha;
			diag[k] = b;
			k++;
		}
		Uptmp[k] = 2*beta*Up[i][n-2] + b1*Up[i][n-1];
		upper[k] = 0; 
		lower[k] = -2*alpha;
		diag[k] = b;
		k++;
	}
	for(int ii=0; ii<n;ii++){
		if(ii==0)
			Uptmp[k] = 2*beta*Up[m-1][1] + b1*Up[m-1][0];
		else if(ii=n-1)
			Uptmp[k] = 2*beta*Up[m-1][n-2] + b1*Up[m-1][n-1];
		else
			Uptmp[k] = beta*Up[m-1][ii+1] + b1*Up[m-1][ii] + beta*Up[m-1][ii-1];
		upper[k] = 0; 
		lower[k] = -2*alpha;
		diag[k] = b;
		k++;
	}
	tridiag(Utmp,Uptmp,m*n,diag,upper,lower);
	b = 2*(1+beta);
	b1 = 2*(1-alpha);
	k=0;
	int kk=0;
	for(int i=1; i<(m-1);i++){
		for(int l=1; l<(n-1); l++){
			Up[i][l] = Utmp[kk];
			kk++;
		}
	}
	for(int ii=0; ii<n;ii++){
		if(ii==0)
			Uptmp[k] = 2*alpha*Up[1][0] + b1*Up[0][0];
		else if(ii=n-1)
			Uptmp[k] = 2*alpha*Up[m-2][0] + b1*Up[m-1][0];
		else
			Uptmp[k] = alpha*Up[ii+1][0] + b1*Up[ii][0] + alpha*Up[ii-1][0];
		upper[k] = -2*beta; 
		lower[k] = 0;	
		diag[k] = b;
		k++;
	}
	for(int i=1;i<(m-1);i++){
		Uptmp[k] = 2*alpha*Up[1][i] + b1*Up[0][i];
		upper[k] = -2*beta;
		lower[k] = 0;
		diag[k] = b;
		k++;
		for(int j=1;j<(n-1);j++){
			Uptmp[k] = alpha*Up[i+1][j] + b1*Up[i][j] + alpha*Up[i-1][j];
			upper[k] = -beta;
			lower[k] = -beta;
			diag[k] = b;
			k++;
		}
		Uptmp[k] = 2*alpha*Up[m-2][i] + b1*Up[m-1][i];
		upper[k] = 0;
		lower[k] = -2*beta;
		diag[k] = b;
		k++;
	}
	tridiag(Utmp,Uptmp,m*n,diag,upper,lower);
	k=0;
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			Up[i][j]=Uptmp[k];
			U[i][j] = Utmp[k];
		}
	}
}