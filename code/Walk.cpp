#include "Walk.h"
#include <iostream>
using namespace std;

Walk::Walk(int n, int dimension)
{
	d = dimension;
	x0 = 0; x1 = 1;
	y0 = 0; y1 = 1;
	z0 = 0; z1 = 1;
	nwalkers = n;
	double *tmp[nwalkers];
	for(int i=0; i<nwalkers;i++){
		tmp[i] = new double[d];
	}
	walkers = tmp;
	cout<<walkers<<endl;
};

void Walk::SetInitialCondition(double **C){
	int a = 0;
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

double **Walk::advance(double **C){
	double **i;

	return i; 
}


int Walk::InitializeTimestep(int **C){
	

	return 0;

}
void Walk::PutWalkers(){

}
double **Walk::ReturnBoundary(){
	//should have a better name
	double **tmp;
	return tmp;
}	
double *Walk::FindPosition(){
	double *tmp;
	return tmp;
}
double **Walk::CalculateGradient(){
	double **tmp;
	return tmp;
}
double *Walk::checkpos(){
	double *tmp;
	return tmp;
}

int len(double*){
	/*Length function in python*/
	double*::iterator it;
	while(it != null){
		it++;
	}
}