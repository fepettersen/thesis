#include "main_walk.h"
using namespace std;
using namespace arma;

Tridiag::Tridiag(void){
	
}

void Tridiag::precondition(mat A, int m, int n){
	/*Preconditioner for the efficient tridiag method.*/
	D.clear();
	H.clear();
	aa.clear();
	mat a = zeros(n,n);
	mat b = zeros(n,n);
	mat c = zeros(n,n);
	
	int i = 0;
	for(int j=0;j<n;j++){
		for(int k=0;k<n;k++){
			b(j,k) = A(i*n+j,i*n+k);
			c(j,k) = A(i*n+j,(i+1)*n+k);
		}
	}

	D.push_back(inv(b));
	H.push_back(-1*D[i]*c);
	
	for(i=1;i<m-1;i++){
		for(int j=0;j<n;j++){
			for(int k=0;k<n;k++){
				a(j,k) = A(i*n+j,(i-1)*n+k);
				b(j,k) = A(i*n+j,i*n+k);
				c(j,k) = A(i*n+j,(i+1)*n+k);
			}
		}
		aa.push_back(a);
		D.push_back(inv(b+a*H[i-1]));
		H.push_back(-1*D[i]*c);

	}
	i = m-1;
	for(int j=0;j<n;j++){
		for(int k=0;k<n;k++){
			a(j,k) = A(i*n+j,(i-1)*n+k);
			b(j,k) = A(i*n+j,i*n+k);
		}
	}
	aa.push_back(a);
	D.push_back(inv(b+a*H[i-1]));
}



vec Tridiag::blockTridiag(mat A, vec Up, int m, int n){
	std::vector<mat> tmp;
	std::vector<mat> HH;
	std::vector<vec> g;
	mat a = zeros(n,n);
	mat b = zeros(n,n);
	mat c = zeros(n,n);
	vec K = zeros(n);
	vec gtmp;
	int i = 0;
	for(int j=0;j<n;j++){
		for(int k=0;k<n;k++){
			b(j,k) = A(i*n+j,i*n+k);
			c(j,k) = A(i*n+j,(i+1)*n+k);
		}
		K(j) = Up[i*n+j];
	}
	tmp.push_back(inv(b));
	HH.push_back(-1*tmp[i]*c);
	g.push_back(tmp[i]*K);
	for(i=1;i<m-1;i++){
		for(int j=0;j<n;j++){
			for(int k=0;k<n;k++){
				a(j,k) = A(i*n+j,(i-1)*n+k);
				b(j,k) = A(i*n+j,i*n+k);
				c(j,k) = A(i*n+j,(i+1)*n+k);
			}
			K(j) = Up[i*n+j];
		}
		tmp.push_back(inv(b+a*HH[i-1]));
		HH.push_back(-1*tmp[i]*c);
		gtmp = g[i-1];
		g.push_back(tmp[i]*(K-a*gtmp));
	}
	i = m-1;
	for(int j=0;j<n;j++){
		for(int k=0;k<n;k++){
			a(j,k) = A(i*n+j,(i-1)*n+k);
			b(j,k) = A(i*n+j,i*n+k);
		}
		K(j) = Up[i*n+j];
	}
	tmp.push_back(inv(b+a*HH[i-1]));
	gtmp = g[i-1];
	g.push_back(tmp[i]*(K-a*gtmp));
	gtmp = g[i];

	vec x = zeros(m*n);
	for(int j=0;j<n;j++){
		x(i*n+j) = gtmp(j);
	}
	vec xtmp;
	for(i=m-2;i>-1;i--){
		xtmp = g[i]+HH[i]*gtmp;
		for(int j=0;j<n;j++){
			x(i*n+j) = xtmp(j);
			gtmp(j) = xtmp(j);
		}
	}
	return x;
}

vec Tridiag::efficient_blockTridiag(vec Up, int m, int n){
	/*efficient version of the blockTridiag version taking advantage 
	of the fact that the matrix A  stays the same as long as the 
	parameters are unchanged. 
	Requires that the preconditioner has ben called!*/
	std::vector<vec> g;
	vec K = zeros(n);
	vec gtmp;
	int i = 0;
	for(int j=0;j<n;j++){
		K(j) = Up[i*n+j];
	}

	g.push_back(D[i]*K);
	for(i=1;i<m-1;i++){
		for(int j=0;j<n;j++){
			K(j) = Up[i*n+j];
		}
		gtmp = g[i-1];
		g.push_back(D[i]*(K-aa[i-1]*gtmp));
	}
	i = m-1;
	for(int j=0;j<n;j++){
		K(j) = Up[i*n+j];
	}

	gtmp = g[i-1];
	g.push_back(D[i]*(K-aa[i-1]*gtmp));
	gtmp = g[i];

	vec x = zeros(m*n);
	for(int j=0;j<n;j++){
		x(i*n+j) = gtmp(j);
	}
	vec xtmp;
	for(i=m-2;i>-1;i--){
		xtmp = g[i]+H[i]*gtmp;
		for(int j=0;j<n;j++){
			x(i*n+j) = xtmp(j);
			gtmp(j) = xtmp(j);
		}
	}
	return x;
}

void Tridiag::tridiag(double *u, double *f, int N, double *a, double *b, double *c){
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