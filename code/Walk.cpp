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
	ResetWalkers();
	nwalkers = 0;
	for(int i=0; i<m;i++){
		for(int j=0; j<n;j++){
			nwalkers += C[i][j];
		}
	}
	// cout<<"nwalkers = "<<nwalkers<<endl;
	walkers.resize(nwalkers);
	for(int i=0;i<nwalkers;i++)
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
	double *newPos, **s;
	newPos = new double[d];			//delete at the end!
	int *index = new int[d];		//delete at the end!

	for(int i=0;i<m; i++){
		//Empty the array to conserve energy
		for(int j=0; j<n; j++){
			C[i][j] = 0;
		}
	}
	double L = 0;
	double L_deriv = 0;
	double L0 = sqrt(2*dt);
	double L_deriv0 = (dt/d)/(2*dx*sqrt(2*(dt/d)));
	double Tr,Tl,Delta_m,Delta_p,r,stepvector[d];
	if(not inhomogenous){
		L = L0*sqrt(D);
		Tr = 0.5;
		Delta_p = L;
		Delta_m = L;
	}

	for(int i=0; i<nwalkers; i++){
		/*For every walker: */
		if(inhomogenous){
			FindPosition(walkers[i],index);
			L = (d>1)?(L0*sqrt(aD[index[0]][index[1]])):(L0*sqrt(aD[index[0]][0]));
			L_deriv = (d>1)?(L_deriv0/sqrt(aD[index[0]][index[1]])*(aD[index[0]+1][index[1]]-aD[index[0]+1][index[1]])):(L_deriv0/sqrt(aD[index[0]][0])*(aD[index[0]+1][0]-aD[index[0]+1][0])); 	/*This should not work and is horrible programming*/
			Tr = (1+0.5*L_deriv);
			Tl = (1-0.5*L_deriv);
			Delta_p = L*Tr;
			Delta_m = L*Tl;
			Tr /= 2.0;
		}
		for(int j=0;j<steps;j++){
			/*Advance n steps*/
			for(int p=0; p<d; p++){
				r = ran0(&Idum);
				if(r>Tr){
					stepvector[p] = Delta_p;
				}
				else{
					stepvector[p] = -Delta_m;
				}
			}
			newPos = InhomogenousStep(walkers[i],stepvector);
			walkers[i] = checkpos(newPos,stepvector);
		}
		FindPosition(walkers[i],index);
		if(d==1){
			C[index[0]][0] += 1;
		}
		else if(d==2){
			C[index[0]][index[1]] += 1;
		}
	}
	// delete [] newPos,index;
}


double *Walk::InhomogenousStep(double *r, double *s){
	for(int l=0;l<d;l++){
		r[l] += s[l];
	}
	return r;
}

void Walk::SetDiffusionTensor(double **Diff, int m, int n){
	aD = new double*[m];
	for(int i=0;i<m;i++){
		aD[i] = new double[n];
	}
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++)
		aD[i][j] = Diff[i][j];
	}
	inhomogenous = true;
}
void Walk::SetDiffusionConstant(double Diff){
	D = Diff;
	inhomogenous = false;
}

void Walk::PutWalkers(int i, int j, int counter){
	// if(debug_walk){cout<<"Walk::PutWalkers"<<endl;}
	for(int k=0; k<d; k++){
		if(k==0){
			walkers[counter][k] = x[i]+dx*(0.5-ran0(&Idum));
		}
		else if(k==1){
			walkers[counter][k] = y[j]+dy*(0.5-ran0(&Idum));
		}
	}
}


void Walk::FindPosition(double *pos, int *indx){
	/*Maps the walkers position to its index*/
	// int indx[d];
	if(d==1){
		indx[0] = -1;
		for(int q=0; q<m; q++){
			if(fabs(pos[0]-x[q])<dx/2.0){
				indx[0] = q;
				break;
			}
		}
		if(indx[0]==-1){
			cout<<"panic!!"<<endl;
			exit(1);
		}
	}
	else if(d==2){
		indx[0] = indx[1] = -1;
		for(int q=0; q<m; q++){
			if(fabs(pos[0]-x[q])<dx/2.0){
				indx[0] = q;
				break;
			}
		}
		for(int h=0; h<n; h++){
			if(fabs(pos[1]-y[h])<dy/2.0){
				indx[1] = h;
				break;
			}
		}
		if(indx[0]==-1||indx[1]==-1){
			cout<<"panic!!"<<endl;
			cout<<pos[0]<<","<<pos[1]<<endl;
			exit(1);
		}
	}
	// return indx;
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