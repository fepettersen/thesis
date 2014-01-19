
#include "main_walk.h"
using namespace std;
using namespace arma;



Diffusion::Diffusion(double _dx, double _dy, double _D, double Dt, double _v){
	/*Initializer for Isotropic diffusion
	_Dx and _Dy are collections of variables to save some FLOPs
	v is for the advection term vC*du/dx*/
	dx = _dx;
	dy = _dy;
	D = _D;
	d = (dy>0)?2:1;
	solver = 3;		/*solver=3 ==> Backward Euler; solver=0 ==> Forward Euler*/
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
	isotropic = 1;
};

Diffusion::Diffusion(double _dx, double _dy, double **D, double Dt, double _v){
	/*Initializer for anisotropic diffusion. 
	aD denotes the diffusion "tensor" 
	other variables are as before or collections of variables to save FLOPs*/
	dx = _dx;
	dy = _dy;
	dt = Dt;
	_Dx = dt/(dx*dx);
	aD = D;
	d = (dy>0)?2:1;
	solver = 3;		/*solver=2 ==> Forward Euler; solver=3 ==> BE*/
	if(d==2 && solver==2){
		dt = (Dt>(dx*dy/4.0))? (dx*dy/(5.0)):Dt;
		_Dx = dt/(dx*dx);
		_Dy = dt/(dy*dy);
	}
	else if(d==1 && solver==1){
		double Dtmp = 3.1415926535897932;			/*Should be max(aD)*/
		dt = (Dt>(dx*dx/(2.0*Dtmp)))? (dx*dx/(10.0*Dtmp)):Dt;
		_Dx = dt/(dx*dx);
	}
	t=1;
	v = _v;
	vdtdx2 = (v*dt)/(2*dx);
	vdtdy2 = (v*dt)/(2*dy);
	isotropic = 0;
}

void Diffusion::advance(double **U,double **Up, int m, int n){
	/*Advaces the previous solution one timestep using the solvers 
	specified earlier. 
	TODO: It should be possible to change solver from python script...*/
	if (solver==0){
		/*Isotropic Forward Euler*/
		_Dx = dt/(dx*dx);
		if(d==2){
			_Dy = dt/(dy*dy);
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
	else if(solver==3){
		/*anisotropic BE in 2d is under development and may not work
		1D version seems to work perfectly*/
		if(isotropic && d==1){
			/*should call tridiag in stead of this mess!*/
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
		else{
			BE2D(U,Up,m,n);
		}
	}
	else{
		/*Anisotropic convection diffusion - Forward Euler*/
		_Dx = dt/(dx*dx);
		if(d==1){
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
			_Dy = dt/(dy*dy);
			vdtdy2 = v*dt/(dy*dy);
			vdtdx2 = v*dt/(dx*dx);
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

void Diffusion::boundary(double **U,double **Up,int m, int n){
	/*Boundary conditions for the explicit schemes.*/
	if(solver==0){
		/*Isotropic FE*/
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
			/*corners*/
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
		/*Anisotropic FE (with advection)*/
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
			/*Corners*/
			U[0][0] = _Dx*(aD[1][0]+aD[0][0])*(Up[1][0]-Up[0][0]) +Up[0][0] + 
			_Dy*(aD[0][1]+aD[0][0])*(Up[0][1]-Up[0][0])+dt*f(0,0,t*dt);

			U[m-1][0] = _Dx*(aD[m-2][0]+aD[m-1][0])*(Up[m-2][0]-Up[m-1][0]) + Up[m-1][0] + 
			_Dy*(aD[m-1][1]+aD[m-1][0])*(Up[m-1][1]-Up[m-1][0])+dt*f((m-1)*dx,0,t*dt);
			
			U[0][n-1] = _Dx*(aD[1][n-1]+aD[0][n-1])*(Up[1][n-1]-Up[0][n-1]) + Up[0][n-1] + 
			_Dy*(aD[0][n-2]+aD[0][n-1])*(Up[0][n-2]-Up[0][n-1])+dt*f(0,(n-1)*dy,t*dt);
			
			U[m-1][n-1] = _Dx*(aD[m-2][n-1]+aD[m-1][n-1])*(Up[m-2][n-1]-Up[m-1][n-1]) +Up[m-1][n-1] +
			 _Dy*(aD[m-1][n-2]+aD[m-1][n-1])*(Up[m-1][n-2]-Up[m-1][n-1])+dt*f((m-1)*dx,(n-1)*dy,t*dt);
		}
	}
}



double Diffusion::f(double x,double y, double t){
	/*Source term. This should be possible to specify from somewhere. 
	Perhaps by inheritence*/
	double pi = 3.1415926535897932;
	// return exp(-t*pi*pi)*pi*pi*(sin(pi*x)+cos(pi*x)*(pi*x-1));
	// return -pi*sin(pi*x)*exp(-t*pi*pi);
	double py = pi*y;
	double px = pi*x;
	// return pi*exp(-t*pi*pi)*(2*pi*cos(px)*cos(py)*(x+y-0.5) + cos(px)*sin(py) +sin(px)*cos(py));
	return 0.0;
}

void Diffusion::tridiag(double *u, double *f, int N, double *a, double *b, double *c){
	/*Specialized gaussian elimination for tridiagonal matrices (BE1D)*/
    double *temp = new double[N];
    for(int i=0;i<N;i++){
    	temp[i] = 0;
    }
    double btemp = b[0];
    u[0] = f[0]/btemp;
    for(int i=1; i<N; i++){
        //forward substitution
       	temp[i] = c[i-1]/btemp;
       	btemp = b[i]-a[i]*temp[i];
        u[i] = (f[i] -a[i]*u[i-1])/btemp;
    }
    for(int i=(N-2);i>=0;i--){
        //Backward substitution
        u[i] -= temp[i+1]*u[i+1];
    }
    delete[] temp;
}
vec Diffusion::BlockTridiag(mat A, vec Up, int m, int n){
	/*Specialized Gaussian elimination for block-tridiagonal systems A*x = Up. 
	Input is the block-tridiagonal matrix A, the solution in the vector Up 
	and the unknown x which will be filled with the solution*/
	std::vector<mat> tmp;
	std::vector<mat> H;
	std::vector<vec> g;
	mat a = zeros(m,n);
	mat b = zeros(m,n);
	mat c = zeros(m,n);
	vec K = zeros(m);
	vec gtmp;
	int i = 0;
	for(int j=0;j<m;j++){
		for(int k=0;k<n;k++){
			b(j,k) = A(i*n+j,i*n+k);
			c(j,k) = A(i*n+j,(i+1)*n+k);
		}
		K(j) = Up[i*n+j];
	}
	tmp.push_back(inv(b));
	H.push_back(-1*tmp[i]*c);
	g.push_back(tmp[i]*K);
	for(i=1;i<m-1;i++){
		for(int j=0;j<m;j++){
			for(int k=0;k<n;k++){
				a(j,k) = A(i*n+j,(i-1)*n+k);
				b(j,k) = A(i*n+j,i*n+k);
				c(j,k) = A(i*n+j,(i+1)*n+k);
			}
			K(j) = Up[i*n+j];
		}
		tmp.push_back(inv(b+a*H[i-1]));
		H.push_back(-1*tmp[i]*c);
		gtmp = g[i-1];
		g.push_back(tmp[i]*(K-a*gtmp));
	}
	i = m-1;
	for(int j=0;j<m;j++){
		for(int k=0;k<n;k++){
			a(j,k) = A(i*n+j,(i-1)*n+k);
			b(j,k) = A(i*n+j,i*n+k);
		}
		K(j) = Up[i*n+j];
	}
	tmp.push_back(inv(b+a*H[i-1]));
	gtmp = g[i-1];
	g.push_back(tmp[i]*(K-a*gtmp));
	gtmp = g[i];

	vec x = zeros(m*n);
	for(int j=0;j<m;j++){
		x(i*n+j) = gtmp(j);
	}
	vec xtmp;
	for(i=n-2;i>-1;i--){
		xtmp = g[i]+H[i]*gtmp;
		for(int j=0;j<m;j++){
			x(i*n+j) = xtmp(j);
			gtmp(j) = xtmp(j);
		}
	}
	return x;
}

void Diffusion::BE2D(double **U, double **Up, int m, int n){
	/*Backward Euler scheme in 1d & 2d. Assembles a matrix if necessary and 
	saves a LU decomposition for further use. This will reduce the number of 
	FLOPs needed by m*m (m must be equal to n) or one order.
	Because we can reuse a lot of code, the tridiag solver is not used in 1d, 
	but this could save another order of FLOPs*/
	int N = m*n;
	if(t==1 || beta != D*dt/(dy*dy) ){
		beta = D*dt/(dy*dy);
		alpha = D*dt/(dx*dx);
		mat A = zeros(N,N);
		if(not isotropic){
			Lower = AssembleAnisotropic(dt/(2*dx*dx),dt/(2*dx*dx),m,n);
		}
		else{
			Lower = Assemble(alpha,beta,m,n);
		}
		// lu(Lower, Upper, Permutation, A);
		// Lower = A.copy();
	}
	vec Uptmp = zeros(N);
	vec Utmp = zeros(N);
    /*Backward substitution for the equations 
    Lower*y = Uptmp(+f)
    Upper*Utmp = y*/
    
    double y[N];
    double LUtemp;
	int k=0;
	for(int i=0; i<m;i++){
		for(int j=0; j<n; j++){
			Uptmp[k] = Up[i][j] +dt*f(i*dx,j*dy,t*dt);
			k++;
		}
	}
	Utmp = BlockTridiag(Lower,Uptmp,m,n);
    // y[0] = Uptmp[0];
    // for(int i=1;i<N;i++){
    //     /* This loop calculates y from L*y = f. It does not consider calculations
    //     with 0. Results are only valid for L matrix in LU decomposition. */
    //     LUtemp=0;
    //     for(int j=0;j<i;j++){
    //         LUtemp += Lower(i,j)*y[j];
    //     }
    //     y[i] = Uptmp[i] -LUtemp;
    // }
    
    // double factor = 0;
    // Utmp[N-1] = y[N-1]/Upper(N-1,N-1);
   
    // for(int i=N-2;i>=0;i--){
    //     /*This loop continues the calculations from the previous one.
    //      Calculates U*z = y, where z = v from earlier in the program.
    //      Resluts are only valid for U matrix in LU decomposition. */
    //     LUtemp=0;
    //     for(int j=(N-1);j>=0;j--){
    //         if(j==i){
    //             factor = 1.0/Upper(i,j);
    //             break;
    //         }
    //         LUtemp += Upper(i,j)*Utmp[j];
    //     }
    //     Utmp[i] =(y[i] -LUtemp)*factor;
    // }
	// Utmp = solve(A,Uptmp);
	k=0;
	for(int i=0; i<m;i++){
		for(int j=0; j<n; j++){
			U[i][j] = Utmp[k];
			k++;
		}
	}
}

mat Diffusion::Assemble(double alpha, double beta, int m, int n){
    /*Assemble the matrix A which will be constant as long as dt is constant
    Note that this matrix will implement Neumann boundary conditions*/
    double gamma = 1+2*alpha+2*beta;
    int N = m*n;
	mat A = eye(N,N);
	A *= gamma;
	for(int i=0;i<m;i++){
		A(i,i+n) = -2*alpha;
		A(N-i-1,N-i-n-1) = -2*alpha;
	}
	int k = 0;
	for (int i=0;i<m;i++){
		A(k,k+1) = -2*beta;
		if (k>(m-1) && k<(N-m)){
			A(k,k+m) = -alpha;
			A(k,k-m) = -alpha;
		}
		k++;
		for (int j=1;j<(n-1);j++){
			A(k,k+1) = -beta;
			A(k,k-1) = -beta;
			if (k>(m-1) && k<(N-m)){
				A(k,k+m) = -alpha;
				A(k,k-m) = -alpha;
			}
			k++;
		}
		A(k,k-1) = -2*beta;
		if(k>(m-1) && k<(N-m)){
			A(k,k+m) = -alpha;
			A(k,k-m) = -alpha;
		}
		k++;
	}
	return A;
}

mat Diffusion::AssembleAnisotropic(double a, double b, int m, int n){
	/*Assembles the matrix for further use in the anisotropic BE 1D & 2D solver.
	Note that this assembler implements Neumann boundary conditions.*/
	
	int N = m*n;
	mat A = zeros(N,N);
	if(d==1){
		int k=0;
		A(k,k) = 1+2*a*(aD[1][0]+aD[0][0]);
		A(k,k+1) = -2*a*(aD[0][0]+aD[1][0]);
		for(k=1;k<(m-1);k++){
			A(k,k) = 1+a*aD[k+1][0]+2*a*aD[k][0] +a*aD[k-1][0];
			A(k,k+1) = -a*(aD[k+1][0]+aD[k][0]);
			A(k,k-1) = -a*(aD[k-1][0]+aD[k][0]);
		}
		A(k,k) = 1+2*a*(aD[m-2][0]+aD[m-1][0]);
		A(k,k-1) = -2*a*(aD[m-1][0]+aD[m-2][0]);
	}
	else if(d==2){
		int k=0;
		for (int i=0;i<m;i++){
			A(i,i+n) = -2*a*(aD[0][i]+aD[1][i]);
			A(N-i-1,N-i-n-1) = -2*a*(aD[m-1][n-1-i]+aD[m-2][n-1-i]);
			for (int j=0;j<n;j++){
				if (j==0){
					A(k,k+1) = -2*b*(aD[i][j+1]+aD[i][j]);
					if (i==0){
						A(k,k) = 1+2*a*aD[1][j] + aD[i][j]*(2*a+2*b) + 2*b*aD[i][1];
					}
					else if(i==m-1){
						A(k,k-1) = 0;
						A(k,k) = 1+2*a*aD[m-2][j] + aD[i][j]*(2*a+2*b) + 2*b*aD[i][1];
					}
					else{
						A(k,k-1) = 0;
						A(k,k) = 1+a*aD[i+1][j]+a*aD[i-1][j] + aD[i][j]*(2*a+2*b) + 2*b*aD[i][1];
					}
				}
				else if (j==n-1){
					if(k != (m*n-1)){
						A(k,k+1) = 0;
					}
					A(k,k-1) = -2*b*(aD[i][j-1]+aD[i][j]);
					if (i==0){
						A(k,k) = 1+2*a*aD[1][j] + aD[i][j]*(2*a+2*b) + 2*b*aD[i][j-1];
					}
					else if (i==m-1){
						A(k,k) = 1+2*a*aD[m-2][j] + aD[i][j]*(2*a+2*b) +2*b*aD[i][j-1];
					}
					else{
						A(k,k) = 1+a*aD[i+1][j]+a*aD[i-1][j] + aD[i][j]*(2*a+2*b) +2*b*aD[i][j-1];
					}
				}
				else if(j!=0){
					A(k,k+1) = -b*(aD[i][j+1]+aD[i][j]);
					A(k,k-1) = -b*(aD[i][j-1]+aD[i][j]);
					if (i==0){
						A(k,k) = 1+2*a*aD[1][j] + aD[i][j]*(2*a+2*b) + b*aD[i][j+1]+b*aD[i][j-1];
					}
					else if( i==m-1){
						A(k,k) = 1+2*a*aD[m-2][j] + aD[i][j]*(2*a+2*b) + b*aD[i][j+1]+b*aD[i][j-1];
					}
					else{
						A(k,k) = 1+a*aD[i+1][j]+a*aD[i-1][j] + aD[i][j]*(2*a+2*b) + b*aD[i][j+1]+b*aD[i][j-1];
					}
				}
				if (k>(m-1) && k<(N-m)){
					A(k,k+m) = -a*(aD[i+1][j]+aD[i][j]);
					A(k,k-m) = -a*(aD[i-1][j]+aD[i][j]);
				}
				k+=1;
			}
		}
	}
	return A;
}