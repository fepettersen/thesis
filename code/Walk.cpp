#include "main_walk.h"
using namespace std;

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double ran0(long *idum)
{
   long     k;
   double   ans;

   *idum ^= MASK;
   k = (*idum)/IQ;
   *idum = IA*(*idum - k*IQ) - IR*k;
   if(*idum < 0) *idum += IM;
   ans=AM*(*idum);
   *idum ^= MASK;
   return ans;
}

// long Idum = -1*rand();
long Idum = -1*time(0);
bool debug_walk = false;

Walk::Walk(int dimension)
{
	if(debug_walk){cout<<"Walk::Walk"<<endl;}
	d = dimension;
	steps = 100;
	x0 = 0; x1 = 1;
	y0 = 0; y1 = 1;
	z0 = 0; z1 = 1;
	int D = 1;
	dt = 1.0/steps;
	factor = sqrt(2*D*dt);
	drift = factor/steps;
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
	walkers = new double*[nwalkers];
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
	ResetWalkers();
	nwalkers = 0;
	for(int i=0; i<m;i++){
		for(int j=0; j<n;j++){
			nwalkers += C[i][j];
		}
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
}
bool Walk::HasLeftArea(double *pos){
	if(debug_walk){cout<<"Walk::HasLeftArea"<<endl;}
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

void Walk::advance(int **C){
	if(debug_walk){cout<<"Walk::advance"<<endl;}
	double *newPos, **s;
	newPos = new double[d];			//delete at the end!
	s = new double*[steps];			//delete at the end!
	int *index = new int[d];		//delete at the end!
	for(int k=0; k<steps; k++){
		s[k] = new double[d];
		for(int l=0; l<d; l++){
			s[k][l] = 0;
		}
	}
	for(int i=0; i<nwalkers; i++){
		/*For every walker: */
		for(int k=0; k<steps; k++){
			/*Fill array s with steps*d random numbers*/
			for(int l=0; l<d; l++){
				s[k][l] = factor*(0.5-ran0(&Idum));
			}
		}
		for(int j=0;j<steps;j++){
			/*Advance n steps*/
			newPos = Step(walkers[i],s[j]);
			walkers[i] = checkpos(newPos,s[j]);
		}
	}

	if(d==1){
		for(int i=0;i<nwalkers;i++){
			index = FindPosition(walkers[i]);
			C[index[0]] += 1;
		}
	}
	else if(d==2){
		for(int i=0; i<nwalkers; i++){
			index = FindPosition(walkers[i]);
			C[index[0]][index[1]] += 1;
		}
	}
}

double *Walk::Step(double *r,double *s){
	/*Possibility of changing the algorithm*/
	if(debug_walk){cout<<"Walk::Step"<<endl;}
	if(false){
		for(int i=0; i<d; i++){
			r[i] += s[i];
		}
	}
	else if(false){
		for(int i=0; i<d; i++){
			r[i] += dt*dt*s[i];		//Need velocity as well
		}
	}
	else if(true){
		if(d==1){
			r[0] = (s[0]>0) ?(r[0]+factor):(r[0]-factor);
		}
		else if(d==2){
			if(s[0]>s[1]){
				r[0] = (s[0]>0) ?(r[0]+factor):(r[0]-factor);
			}
			else{
				r[1] = (s[1]>0) ?(r[1]+factor):(r[1]-factor);
			}
		}
	}
	r[0] += drift;		//slight drift in x direction
	return r;
}

int Walk::InitializeTimestep(int **C){
	return 0;
}

void Walk::PutWalkers(int i, int j, int counter){
	if(debug_walk){cout<<"Walk::PutWalkers"<<endl;}
	for(int k=0; k<d; k++){
		if(k==0){
			/*Should this be something with the steplength? - cange factor in that case*/
			walkers[counter][k] = x[i]+factor*(0.5-ran0(&Idum));
		}
		else if(k==1){
			walkers[counter][k] = y[j]+factor*(0.5-ran0(&Idum));
		}
	}
}

double **Walk::ReturnBoundary(){
	//should have a better name
	double **tmp;
	return tmp;
}

int *Walk::FindPosition(double *pos){
	/*Maps the walkers position to its index*/
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
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)