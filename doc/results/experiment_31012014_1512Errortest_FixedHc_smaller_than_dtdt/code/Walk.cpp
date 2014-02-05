#include "main_walk.h"
using namespace std;

bool debug_walk = false;

Walk::Walk(int dimension, double _dt)
{
	if(debug_walk){cout<<"Walk::Walk"<<endl;}
	d = dimension;
	steps = 100;		/*Should not be specified here!*/
	x0 = 0; x1 = 1;
	y0 = 0; y1 = 1;
	z0 = 0; z1 = 1;
	D = 1;
	dt = _dt/steps;
	factor = sqrt(2*D*dt);
	rng = new Random();
	drift = factor/steps;
	inhomogenous = false;
};

void Walk::SetInitialCondition(int **C, int M, int N){
	/*Takes the initial condition as an array C[m,n] of ints(?) 
	describing how many walkers in each entry. */
	if(debug_walk){cout<<"Walk::SetInitialCondition"<<endl;}
	nwalkers = 0;
	m = M;
	n = N;
	for(int i=0; i<m;i++){
		for(int j=0; j<n;j++){
			nwalkers += C[i][j];
		}
	}
	walkers.resize(nwalkers);
	for(int i=0; i<nwalkers;i++){
		walkers[i] = new double[d];
		for(int l=0; l<d;l++){
			walkers[i][l] = 0;
		}
	}
	x = new double[m];
	y = new double[n];
	dx = (x1-x0)/(m-1);
	dy = (n>1)?((y1-y0)/(n-1)):0;
	for(int k=0;k<m;k++){
		x[k] = k*dx;
	}
	for(int k=0;k<n;k++){
		y[k] = k*dy;
	}
	int counter = 0;

	for(int i=0; i<m;i++){
		for(int j=0; j<n;j++){
			for(int l=0; l<C[i][j]; l++){
				PutWalkers(i,j,counter);
				counter ++;
			}
		}
	}
	x0_ = x0 - (dx/2.0);
	x1_ = x1 + (dx/2.0);
	y0_ = y0 - (dx/2.0);
	y1_ = y1 + (dx/2.0);
}

void Walk::ResetInitialCondition(int **C){
	if(debug_walk){cout<<"Walk:ResetInitialCondition"<<endl;}
	int nwalkers_prev = nwalkers;
	nwalkers = 0;
	for(int i=0; i<m;i++){
		for(int j=0; j<n;j++){
			nwalkers += C[i][j];
		}
	}
	if (nwalkers>nwalkers_prev){
		/*Do something smart*/
		walkers.resize(nwalkers);
	}
	// ResetWalkers();
	// cout<<"nwalkers = "<<nwalkers<<endl;
	for(int i=nwalkers_prev;i<nwalkers;i++)
		walkers[i] = new double[d];

	int counter = 0;
	for(int i=0; i<m;i++){
		for(int j=0; j<n;j++){
			for(int l=0; l<C[i][j]; l++){
				PutWalkers(i,j,counter);
				counter ++;
			}
		}
	}
}
bool Walk::HasLeftArea(double *pos){
	// if(debug_walk){cout<<"Walk::HasLeftArea"<<endl;}
	/*pos = [x] in 1d;
	pos = [x,y] in 2d etc*/
	if(d==1){
		if(pos[0]<x0 || pos[0]>x1){
			return true;
		}
	}
	else if(d==2){
		if(pos[0]>x1 || pos[0]<x0){
			return true;
		}
		if(pos[1]<y0 || pos[1]>y1){
			return true;
		}
	}
	return false;
}

void Walk::InhomogenousAdvance(int **C, double _dt){
	/*This should now implement the normal advance-function as well*/
	if(debug_walk){cout<<"Walk::InhomogenousAdvance"<<endl;}
	dt = _dt/steps;

	for(int i=0;i<m; i++){
		//Empty the array to conserve energy
		for(int j=0; j<n; j++){
			C[i][j] = 0;
		}
	}
	double DX[d];
	DX[0] = dx;
	if(d==2){DX[1] = dx;}
	// #pragma omp parallel 
	// {
	double *newPos, **s;
	newPos = new double[d];
	int *index = new int[d];
	double L = 0;
	double L_deriv = 0;
	double L0 = sqrt(2*dt);
	double L_deriv0 = (dt/d)/(2*sqrt(2*(dt/d)));
	double Tr,Tl,Delta_m,Delta_p,r,stepvector[d];
	if(not inhomogenous){
		L = L0*sqrt(D);
		Tr = 0.5;
		Delta_p = L;
		Delta_m = L;
	}
	int p=0;
	// #pragma omp for
	for(int i=0; i<nwalkers; i++){
		/*For every walker: */
		for(int j=0;j<steps;j++){
			/*Advance n steps*/
			FindPosition(walkers[i],index);
			p = int(round(rng->uniform()*(d-1)));	/*pick a spatial dimension to advance the walker*/
			// for(p=0; p<d; p++){
				if(inhomogenous){
					L = (d>1)?(L0*sqrt(aD[index[0]][index[1]])):(L0*sqrt(aD[index[0]][0]));
					/*This might work and is slightly better*/
					L_deriv = (d>1)?(L_deriv0/(DX[p]*sqrt(aD[index[0]][index[1]]))*(aDx[p][index[0]][index[1]])):
					(L_deriv0/(DX[p]*sqrt(aD[index[0]][0])))*(aDx[p][index[0]][0]);
					Tr = (1+0.5*L_deriv);
					Tl = (1-0.5*L_deriv);
					Delta_p = L*Tr;
					Delta_m = L*Tl;
					Tr /= 2.0;
					// cout<<"Delta_p,Delta_m = "<<Delta_p<<","<<Delta_m<<endl;
				}
				// r = ran0(&Idum);
				r = rng->uniform();
				if(r>Tr){
					stepvector[p] = Delta_p;
				}
				else{
					stepvector[p] = -Delta_m;
				}
			// }
			newPos = InhomogenousStep(walkers[i],stepvector);
			walkers[i] = checkpos(newPos,stepvector);
			stepvector[p] = 0;
		}
		FindPosition(walkers[i],index);
		if(d==1){
			C[index[0]][0] += 1;
		}
		else if(d==2){
			C[index[0]][index[1]] += 1;
		}
	}
// }
	// delete [] newPos,index;
}


double *Walk::InhomogenousStep(double *r, double *s){
	for(int l=0;l<d;l++){
		r[l] += s[l];
	}
	return r;
}

void Walk::SetDiffusionTensor(double **Diff, int M, int N){
	if(debug_walk){cout<<"Walk::SetDiffusionTensor"<<endl;}
	aD = new double*[M];
	dx = (x1-x0)/(M-1);
	dy = (N>1)?((y1-y0)/(N-1)):0;
	aDx.resize(d);
	for(int p=0;p<d;p++){
		aDx[p] = new double*[M];
		for(int i=0;i<M;i++){
			aDx[p][i] = new double[N];
		}
	}
	for(int i=0;i<M;i++){
		aD[i] = new double[N];
	}
	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			aD[i][j] = Diff[i][j];
		}
	}
	if(d==1){
		for(int i=1;i<M-1;i++){
			for(int j=0;j<1;j++){
				aDx[0][i][j] = (1.0/(2*dx))*(aD[i+1][j]-aD[i-1][j]);
			}
		}
		if(M!=0){
			aDx[0][0][0] = (1.0/dx)*(aD[1][0]-aD[0][0]);
			aDx[0][M-1][0] = (1.0/dx)*(aD[M-1][0]-aD[M-2][0]);
		}
	}
	if(d==2){
		for(int i=1;i<M-1;i++){
			for(int j=1;j<N-1;j++){
				aDx[0][i][j] = (1.0/(2*dx))*(aD[i+1][j]-aD[i-1][j]);
				aDx[1][i][j] = (1.0/(2*dy))*(aD[i][j+1]-aD[i][j-1]);
			}
		}
		for(int i=0;i<N;i++){
			aDx[0][0][i] = (1.0/dx)*(aD[1][i]-aD[0][i]);
			aDx[0][M-1][i] = (1.0/dx)*(aD[M-1][i]-aD[M-2][i]);
			if(i>0 && i<N-1){
				aDx[1][0][i] = (1.0/(2*dy))*(aD[0][i+1]-aD[0][i-1]);
				aDx[1][M-1][i] = (1.0/(2*dy))*(aD[M-1][i+1]-aD[M-1][i-1]);
			}
		}
		for(int i=0;i<M;i++){
			aDx[1][i][0] = (1.0/dy)*(aD[i][1]-aD[i][0]);
			aDx[1][i][N-1] = (1.0/dy)*(aD[i][N-2]-aD[i][N-1]);
			if(i>0 && i<M-1){
				aDx[0][i][0] = (1.0/(2*dx))*(aD[i+1][0]-aD[i-1][0]);
				aDx[0][i][N-1] = (1.0/(2*dx))*(aD[i+1][N-1]-aD[i-1][N-1]);
			}
		}
	}
	inhomogenous = true;
}

void Walk::SetDiffusionConstant(double Diff){
	if(debug_walk){cout<<"Walk::SetDiffusionConstant"<<endl;}
	D = Diff;
	inhomogenous = false;
}

void Walk::PutWalkers(int i, int j, int counter){
	// if(debug_walk){cout<<"Walk::PutWalkers"<<endl;}
	// walkers[counter][0] = x[i]+dx*(0.5-ran0(&Idum));
	walkers[counter][0] = x[i]+dx*(0.5-rng->uniform());
	if(d==2){
		// walkers[counter][1] = y[j]+dy*(0.5-ran0(&Idum));
		walkers[counter][1] = y[j]+dy*(0.5-rng->uniform());
	}
}


void Walk::FindPosition(double *pos, int *indx){
	/*Maps the walkers position to its index*/
	indx[0] = int(round(pos[0]/dx));
	if(d==2){
		indx[1] = int(round(pos[1]/dy));
	}
}

double *Walk::checkpos(double *r,double *s){
	/*Implements reflecting boundaries -- Need to adjust the boundaries by dx/2 so each thing has 
	as much space*/
	// tmp = r+s 	/*This is r in this case*/
	if(not HasLeftArea(r)){
		return r;
	}
	else{
		if(d ==1){
			if (r[0]<x0_){
				r[0] = x0_ - (r[0]-x0_);
				}
			else if(r[0]>x1_){
				r[0] = x1_ - (r[0]-x1_);
				}
			}
		else if(d == 2){
			if(r[0]<x0_){
				r[0] = x0_ - (r[0]-x0_);
				}
			else if( r[0]>x1_){
				r[0] = x1_ - (r[0]-x1_);
				}
			if (r[1]<y0_){
				r[1] = y0_ - (r[1]-y0_);
				}
			else if (r[1]>y1_){
				r[1] = y1_ - (r[1]-y1_);
				}
			}
		return r;
	}
}
