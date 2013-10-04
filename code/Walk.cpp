#include "Walk.h"
using namespace std;

Walk::Walk(int dimension)
{
	d = dimension;
	x0 = 0; x1 = 1;
	y0 = 0; y1 = 1;
	z0 = 0; z1 = 1;
	dt = 0.001;
};

void Walk::SetInitialCondition(int **C, int m, int n){
	/*Takes the initial condition as an array C[m,n] of ints(?) 
	describing how many walkers in each entry. */
	int nwalkers = 0;
	for(int i=0; i<m;i++){
		for(int j=0; j<n;j++){
			nwalkers += C[i][j];
		}
	}
	double **tmp = new double*[nwalkers];
	
	for(int i=0; i<nwalkers;i++){
		tmp[i] = new double[d];
	}
	walkers = tmp;
	x = new double[m];
	y = new double[n];
	double dx = (x1-x0)/(m-1);
	double dy = (y1-y0)/(n-1);
	for(int k=0;k<m;k++){
		x[k] = k*dx;
	}
	for(int k=0;k<m;k++){
		y[k] = k*dy;
	}
	int counter = 0;

	for(int i=0; i<m;i++){
		for(int j=0; j<n;j++){
			for(int n=0;n<C[i][j];n++){
				PutWalkers(i,j,counter);
				counter ++;
			}
		}
	}
	double x0_, x1_, y0_, y1_;
	x0_ = x0 - (dx/2.0);
	x1_ = x1 + (dx/2.0);
	y0_ = y0 - (dx/2.0);
	y1_ = y1 + (dx/2.0);
}

bool Walk::HasLeftArea(double *pos){
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

int **Walk::advance(int **C){
	int steps = 100;
	double *newPos, **s;
	newPos = new double[d];
	s = new double*[steps];
	for(int k=0; k<d; k++){
		s[k] = new double[d];
	}
	for(int i=0; i<nwalkers; i++){
		for(int k=0;k<steps; k++){
			for(int n=0; n<d; n++){
				s[k][n] = factor*(0.5-rand());
			}
		}
		for(int j=0;j<steps;j++){
			newPos = Step(walkers[i],s[j]);
			walkers[i] = checkpos(newPos,s[j]);
		}
	}
	return C; 
}

double *Walk::Step(double *r,double *s){
	/*Possibility of changing the algorithm*/
	if(true){
		for(int i=0; i<d; i++){
			r[i] += s[i];
		}
	}
	else if(false){
		for(int i=0; i<d; i++){
			r[i] += dt*dt*s[i];		//Need velovity as well
		}
	}
	return r;
}

int Walk::InitializeTimestep(int **C){
	return 0;
}

void Walk::PutWalkers(int i, int j, int counter){
	for(int k=0; k<d; k++){
		if(k==0){
			walkers[counter][k] = x[i]+factor*(0.5-rand());
		}
		else if(k==1){
			walkers[counter][k] = y[i]+factor*(0.5-rand());
		}
	}
}

double **Walk::ReturnBoundary(){
	//should have a better name
	double **tmp;
	return tmp;
}

int *Walk::FindPosition(double *pos){
	/*Remember to make dx and dy*/
	int *indx;
	if(d==1){
		indx = new int[1];
		for(int i=0; i<m; i++){
			if(fabs(pos[0]-x[i])<dx/2.0){
				indx[0] = i;
				return indx;
			}
		}
	}
	else if(d==2){
		indx = new int[2];
		indx[0] = indx[1] = -1;
		for(int i=0; i<m; i++){
			if(fabs(pos[0]-x[i])<dx/2.0){
				indx[0] = i;
				break;
			}
		}
		for(int j=0; j<n; j++){
			if(fabs(pos[1]-y[j])<dy/2.0){
				indx[1] = j;
				break;
			}
		}
		if(indx[0]==-1||indx[1]==-1){
			cout<<"panic!!"<<endl;
			exit(1);
		}
	}
	return indx;
}
double **Walk::CalculateGradient(){
	double **tmp;
	return tmp;
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

// int len(double*){
// 	/*Length function in python*/
// 	double*::iterator it;
// 	while(it != null){
// 		it++;
// 	}
// }