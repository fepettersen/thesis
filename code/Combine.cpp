#include "main_walk.h"
using namespace std;

bool DEBUG = false;

Combine::Combine(int M, int N, double X0, double X1, double Y0, double Y1,double DiffusionConstant)
{
	if(DEBUG){cout<<"Combine::Combine"<<endl;}
	m = M; n = N;
	d = (n>0) ? 2:1; // x = (condition) ? (value_if_true) : (value_if_false);
	x0 = X0; x1 = X1;
	y0 = Y0; y1 = Y1;
	D = DiffusionConstant;
	walk_areas = 0;
	dx = (x1-x0)/(m-1);
	dy = (y1-y0)/(n-1);
	Hc = 100;

	X = new double[m];
	Y = new double[n];
	int **tmp = new int*[m];
	U = new double*[m];
	Up = new double*[m];
	for(int i=0; i<m; i++){
		tmp[i] = new int[n];
		U[i] = new double[n];
		Up[i] = new double[n];
		X[i] = i*dx;
	}
	for(int j=0; j<n;j++){
		Y[j] = j*dy;
	}
	C = tmp;
	pde_solver = new Diffusion(dx,dy,D);
};

void Combine::Solve(){
	if(DEBUG){cout<<"Combine::Solve"<<endl;}
	int counter = 0;
	pde_solver->advance(U,Up,m,n);
	for(vector<Walk*>::iterator it1 = walk_solvers.begin(); it1 != walk_solvers.end(); it1++){
		ConvertToWalkers(U,c[counter],indeces[counter]);
		(*it1)->ResetInitialCondition(C);
		(*it1)->advance(c[counter]);
		ConvertFromWalkers(U,c[counter],indeces[counter]);
		counter ++;
	}
	for(int k=0; k<m; k++){
		for(int l=0; l<m; l++){
			Up[k][l] = U[k][l];
		}
	}
}

void Combine::AddWalkArea(double *x, double *y){
	if(DEBUG){cout<<"Combine::AddWalkArea"<<endl;}
	walk_areas++;
	int M = m; int N = n;
	int **index = new int*[2];
	for(int k=0; k<2;k++){
		index[k] = new int[d];
	}
	Walk *tmp = new Walk(d);
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
	// for(int k=0; k<M; k++){
	// 	for(int l=0; l<N; l++){
	// 		cout<<C[k][l]<<"  ";
	// 	}
	// 	cout<<endl;
	// }
	ConvertToWalkers(Up,Ctmp,index);
	tmp->SetInitialCondition(Ctmp,M,N);		/*The solution in the relevant area converted to walkers*/
	walk_solvers.push_back(tmp);		/*Throws std::bad_alloc*/
	indeces.push_back(index);
	c.push_back(Ctmp);
}

void Combine::ConvertToWalkers(double **u, int **Conc, int **index){
	if(DEBUG){cout<<"Combine::ConvertToWalkers"<<endl;}
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
		for(int l=0; l<N; l++){
			Conc[k][l] = (int) (u[k+m0][l+n0]*Hc);
		}
	}
}

void Combine::ConvertFromWalkers(double **u, int**Conc, int **index){
	if(DEBUG){cout<<"Combine::ConvertFromWalkers"<<endl;}
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
			u[k+m0][j+n0] = 0.5*((Conc[k][j]/Hc)+u[k+m0][j+n0]);
			// u[k+m0][j+n0] = Conc[k][j]/Hc;
		}
	}
}

void Combine::MapAreaToIndex(double *x,double *y, int **index){
	if(DEBUG){cout<<"Combine::MapAreaToIndex"<<endl;}
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
	if(DEBUG){cout<<"Combine::SetInitialCondition"<<endl;}
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