#include "main_walk.h"
using namespace std;

int debug = 0;

Combine::Combine(int M, int N, double X0, double X1, double Y0, double Y1,double DiffusionConstant, double factor, double Dt)
{
	if(debug){cout<<"Combine::Combine"<<endl;}
	m = M; n = N;
	d = (n>1) ? 2:1; // x = (condition) ? (value_if_true) : (value_if_false);
	x0 = X0; x1 = X1;
	y0 = Y0; y1 = Y1;
	D = DiffusionConstant;
	walk_areas = 0;
	dx = (x1-x0)/(m-1);
	dy = (n>1)?(y1-y0)/(n-1):0;

	X = new double[m];
	Y = new double[n];
	C = new int*[m];
	U = new double*[m];
	Up = new double*[m];
	for(int i=0; i<m; i++){
		C[i] = new int[n];
		U[i] = new double[n];
		Up[i] = new double[n];
		X[i] = i*dx;
	}
	for(int j=0; j<n;j++){
		Y[j] = j*dy;
	}

	pde_solver = new Diffusion(dx,dy,D,Dt);
	Hc = factor;
};
Combine::Combine(int M, int N, double X0, double X1, double Y0, double Y1,double **DiffTensor, double factor, double Dt)
{
	if(debug){cout<<"Combine::Combine"<<endl;}
	m = M; n = N;
	d = (n>1) ? 2:1; // x = (condition) ? (value_if_true) : (value_if_false);
	x0 = X0; x1 = X1;
	y0 = Y0; y1 = Y1;
	aD = DiffTensor;
	walk_areas = 0;
	dx = (x1-x0)/(m-1);
	dy = (n>1)?(y1-y0)/(n-1):0;

	X = new double[m];
	Y = new double[n];
	C = new int*[m];
	U = new double*[m];
	Up = new double*[m];
	for(int i=0; i<m; i++){
		C[i] = new int[n];
		U[i] = new double[n];
		Up[i] = new double[n];
		X[i] = i*dx;
	}
	for(int j=0; j<n;j++){
		Y[j] = j*dy;
	}

	double stabilitycriteria = dx*dx/(4*abs_max(aD,m,n));
	if(Dt>stabilitycriteria){
		Dt = stabilitycriteria;
	}
	pde_solver = new Diffusion(dx,dy,aD,Dt);
	Hc = factor;
};


void Combine::Solve(){
	if(debug){cout<<"Combine::Solve"<<endl;}
	int counter = 0;
	pde_solver->advance(U,Up,m,n);
	for(vector<Walk*>::iterator it1 = walk_solvers.begin(); it1 != walk_solvers.end(); it1++){
		(*it1)->drift = 0;
		ConvertToWalkers(U,c[counter],indeces[counter]);
		(*it1)->ResetInitialCondition(C);
		(*it1)->advance(c[counter]);
		// (*it1)->InhomogenousAdvance(c[counter],ad[indeces[counter]]);	/*Maybe we should have a walk.SetDiffusionConstant(D) and let the advance function figure out what to do*/
		ConvertFromWalkers(U,c[counter],indeces[counter]);
		counter ++;
	}

	for(int k=0; k<m; k++){
		for(int l=0; l<n; l++){
			Up[k][l] = U[k][l];
		}
	}
}

void Combine::AddWalkArea(double *x, double *y){
	if(debug){cout<<"Combine::AddWalkArea"<<endl;}
	walk_areas++;
	int M = 0; int N = 0;
	int **index = new int*[2];
	for(int k=0; k<2;k++){
		index[k] = new int[d];
	}
	Walk *tmp = new Walk(d);
	// tmp->SetDiffusionConstant();
	// tmp->SetDiffusionTensor();
	MapAreaToIndex(x,y,index);
	if(d==1){
		M = index[1][0]-index[0][0];
		N = 1;
	}
	else if(d==2){
		M = index[0][1]-index[0][0];
		N = index[1][1]-index[1][0];
	}
	int **Ctmp = new int*[M];
	for(int k=0; k<M;k++){
		Ctmp[k] = new int[N];
	}
	ConvertToWalkers(Up,Ctmp,index);
	tmp->SetInitialCondition(Ctmp,M,N);		/*The solution in the relevant area converted to walkers*/
	walk_solvers.push_back(tmp);		/*Throws std::bad_alloc*/
	indeces.push_back(index);
	c.push_back(Ctmp);
}

void Combine::ConvertToWalkers(double **u, int **Conc, int **index){
	if(debug){cout<<"Combine::ConvertToWalkers"<<endl;}
	int M,N,m0,n0;
	m0 = index[0][0];
	if(d==1){
		M = index[1][0]-index[0][0];
		N = 1;
		n0 = 0;
	}
	else if(d==2){
		M = index[0][1]-index[0][0];
		N = index[1][1]-index[1][0];
		n0 = index[1][0];
	}
	signmap = new int*[M];
	for(int e=0;e<M;e++){
		signmap[e] = new int[N];
	}
	for(int k=0; k<M; k++){
		for(int l=0; l<N; l++){
			Conc[k][l] = (int) (fabs(u[k+m0][l+n0]*Hc));
			(u[k+m0][l+n0]>0)?(signmap[k][l]=1):(signmap[k][l]=-1);
		}
	}
}

void Combine::ConvertFromWalkers(double **u, int**Conc, int **index){
	if(debug){cout<<"Combine::ConvertFromWalkers"<<endl;}
	int M,N,m0,n0;
	m0 = index[0][0];
	if(d==1){
		M = index[1][0]-index[0][0];
		N = 1;
		n0 = 0;
	}
	else if(d==2){
		M = index[0][1]-index[0][0];
		N = index[1][1]-index[1][0];
		n0 = index[1][0];
	}

	for(int k=0; k<M; k++){
		for(int j=0; j<N; j++){
			/*This is where we insert least squares or similar*/
			u[k+m0][j+n0] = 0.5*((signmap[k][j]*Conc[k][j]/Hc)+u[k+m0][j+n0]);
			// u[k+m0][j+n0] = Conc[k][j]/Hc;
		}
	}
	delete[] signmap;
}

void Combine::MapAreaToIndex(double *x,double *y, int **index){
	if(debug){cout<<"Combine::MapAreaToIndex"<<endl;}
	/*x and y holds the limits of the area:
	x[0],y[0]; x[0],y[1]; x[1],y[0]; x[1],y[1]*/
	double eps = 0.0001;
	if(d ==1){
		for(int k=0; k<m;k++){ 		//loop over the x-mesh
			if(fabs(x[0]-X[k])<eps){
				index[0][0] = k;
			}
			if(fabs(x[1]-X[k])<eps)
				index[1][0] = k;
		}
		/*Test that this worked*/
	}
	else if(d==2){
		for(int k=0; k<m;k++){
			if(fabs(x[0]-X[k])<eps){
				index[0][0] = k;
			}
			if(fabs(x[1]-X[k])<eps)
				index[0][1] = k;
		}
		for(int l=0; l<n; l++){
			if(fabs(y[0]-Y[l])<eps){
				index[1][0] = l;
			}
			if(fabs(y[1]-Y[l])<eps){
				index[1][1] = l;
			}
		}
	}
	// cout<<"MapAreaToIndex:"<<endl;
	// cout<<index[0][0]<<","<<index[0][1]<<" , "<<index[1][0]<<","<<index[1][1]<<endl;
}

void Combine::SetInitialCondition(double** U0,int x,int y){	
	if(debug){cout<<"Combine::SetInitialCondition"<<endl;}
	if(x !=m || y!=n){
		cout<<"Initial condition must have same dimensions as rest of problem"<<endl;
		exit(1);
	}
	else{
		for(int i=0; i<x; i++){
			for(int j=0; j<y; j++){
				Up[i][j] = U0[i][j];
			}
		}
	}
}

double Combine::abs_max(double **array,int m, int n){
	double max = array[0][0];
	double tmp = 0;
	for(int i=0; i<m; i++){
		for(int j=0; j<n; j++){
			tmp = fabs(array[i][j]);
			if(tmp>max){
				max = tmp;
			}
		}
	}
	return max;
}