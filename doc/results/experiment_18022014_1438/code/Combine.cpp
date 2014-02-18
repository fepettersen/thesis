#include "main_walk.h"
using namespace std;
using namespace arma;

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
	U = new double*[m];
	Up = new double*[m];
	for(int i=0; i<m; i++){
		U[i] = new double[n];
		Up[i] = new double[n];
		X[i] = i*dx;
	}
	for(int j=0; j<n;j++){
		Y[j] = j*dy;
	}

	pde_solver = new Diffusion(dx,dy,D,Dt);
	rng = new Random();
	Hc = factor;
	inhomogenous = false;
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
	U = new double*[m];
	Up = new double*[m];
	for(int i=0; i<m; i++){
		U[i] = new double[n];
		Up[i] = new double[n];
		X[i] = i*dx;
	}
	for(int j=0; j<n;j++){
		Y[j] = j*dy;
	}

	// double stabilitycriteria = dx*dx/(4*abs_max(aD,m,n));
	// if(Dt>stabilitycriteria){
	// 	Dt = stabilitycriteria;
	// }
	pde_solver = new Diffusion(dx,dy,aD,Dt);
	rng = new Random();
	Hc = factor;
	inhomogenous = true;
};


void Combine::Solve(){
	if(debug){cout<<"Combine::Solve"<<endl;}
	int num_procs = 2;
	pde_solver->advance(U,Up,m,n);
	double diffnorm = norm(U,Up,m,n);
	char cmd[100];
	char diffT[40];
	for(int i=0; i<walk_areas;i++){
		// (*it1)->drift = 0;
		ConvertToWalkers(Up,inifilenames[i],indeces[i]);
		// (*it1)->ResetInitialCondition(c[counter]);
		// (*it1)->InhomogenousAdvance(c[counter],pde_solver->dt);
		sprintf(diffT,"stochastic/DiffusionTensor_%d.txt",i+1);
		// sprintf(cmd,"mpirun -np %d stochastic/%s %d %d %d %g",num_procs,prgm.c_str(),parameters[i][0],parameters[i][1],walk_steps,pde_solver->dt);
		sprintf(cmd,"./stochastic/%s %d %d %d %g %s %s",prgm.c_str(),parameters[i][0],parameters[i][1],walk_steps,pde_solver->dt,inifilenames[i].c_str(),diffT);
		int failure = system(cmd);
		if(failure){
			cout<<endl;
			cout<<"command \""<<cmd<<"\" gave return value: "<<failure<<endl;
			cout<<endl<<endl;
			exit(1);
		}
		ConvertFromWalkers(U,inifilenames[i],indeces[i]);
	}
	for(int k=0; k<m; k++){
		for(int l=0; l<n; l++){
			Up[k][l] = U[k][l];
		}
	}
	if (diffnorm<0.1*pde_solver->dt && pde_solver->dt <0.01)
	{
		int lkjh=0;
		// pde_solver->dt *=10;
		// cout<<"changing dt"<<endl;
	}
}

void Combine::AddWalkArea(double *x, double *y){
	/*sets up a new walk-solver for the area defined by the rectangle x,y. 
	Based on the start and end points in x and y we find the relevant indeces in 
	the solution array, U, and calculate a distribution of walkers C = Hc*U. 
	The isotropic/anisotropic diffusion constant is set for the area.*/
	if(debug){cout<<"Combine::AddWalkArea"<<endl;}
	walk_areas++;
	int M = 0; int N = 0;
	int **index = new int*[2];
	for(int k=0; k<2;k++){
		index[k] = new int[2];
	}
	//Walk *tmp = new Walk(d,pde_solver->dt);
	MapAreaToIndex(x,y,index);
	if(d==1){
		M = index[0][1]-index[0][0];
		N = 1;
	}
	else if(d==2){
		M = index[0][1]-index[0][0];
		N = index[1][1]-index[1][0];
	}

	if(d==2){
		double **temp = new double*[M];
		for(int k=0; k<(M);k++){
			temp[k] = new double[N];
		}
		for(int k=0; k<M;k++){
			for(int l=0; l<N;l++){
				temp[k][l] = aD[index[0][0]+k][index[1][0]+l];
			}
		}
		SaveDiffusionTensor(temp,M,N,walk_areas);
	}
	else if(d==1){
		double **temp = new double*[M];
		for(int k=0; k<M;k++){
			temp[k] = new double[1];
		}
		for(int k=0; k<M;k++){
			for(int l=0; l<1;l++){
				temp[k][l] = aD[index[0][0]+k][0];
			}
		}
		SaveDiffusionTensor(temp,M,N,walk_areas);
	}
	char tmp[100];
	sprintf(tmp,"stochastic/inifile_%d.bin",walk_areas);
	string name = tmp;
	inifilenames.push_back(name);
	int* parameter_tmp = new int[2];
	parameter_tmp[0] = M; parameter_tmp[1] = N;
	parameters.push_back(parameter_tmp);
	walk_steps = 10;
	prgm = "walk_solver";
	indeces.push_back(index);
	// int **Ctmp = new int*[M];
	// signmap = new int*[M];
	// for(int k=0; k<M;k++){
	// 	Ctmp[k] = new int[N];
	// 	signmap[k] = new int[N];
	// }
	// ConvertToWalkers(Up,Ctmp,index);
	// tmp->SetInitialCondition(Ctmp,M,N);		/*The solution in the relevant area converted to walkers*/
	// walk_solvers.push_back(tmp);		/*Throws std::bad_alloc*/
	// c.push_back(Ctmp);
	/*deallocate index?*/

}
void Combine::SaveDiffusionTensor(double** tmp,int M, int N, int no){
	if(debug){cout<<"Combine::SaveDiffusionTensor"<<endl;}
	/*Save the relevant part of the diffusion tensor to 
	a file.*/
	ofstream outfile;
	char neame[100];
	sprintf(neame,"stochastic/DiffusionTensor_%d.txt",no);
	outfile.open(neame);
	for(int i =0;i<M;i++){
		for(int j=0;j<N; j++){
			outfile<<tmp[i][j]<<" ";
		}
		outfile<<endl;
	}
	outfile.close();
}
void Combine::ConvertToWalkers(double **u, string filename, int **index){
	/*Converts the solution, U, to a distribution of walkers for these particular 
	indeces.*/
	if(debug){cout<<"Combine::ConvertToWalkers"<<endl;}
	int M,N,m0,n0;
	m0 = index[0][0];
	double DX=0,DY=0;
	if(d==1){
		M = index[0][1]-index[0][0];
		N = 1;
		n0 = 0;
		DX = 1.0/(M-1);
	}
	else if(d==2){
		M = index[0][1]-index[0][0];
		N = index[1][1]-index[1][0];
		n0 = index[1][0];
		DX = 1.0/(M-1);
		DY = 1.0/(N-1);
	}
	int** Conc = new int*[M];
	int nwalkers = 0;
	for(int k=0; k<M; k++){
		Conc[k] = new int[N];
		for(int l=0; l<N; l++){
			Conc[k][l] = (int) (round(fabs(u[k+m0][l+n0]*Hc)));
			nwalkers += Conc[k][l];
		}
	}	
	char* asdf = new char[100];
	sprintf(asdf,"stochastic/%s",filename.c_str());
	ofstream inifile (filename.c_str(), ios::out | ios::binary);
	

	/* Need to make new x & y to make sure they are mapped to unit square*/
	vec x = linspace(0,1,M);
	vec y = linspace(0,1,N);
	inifile<<nwalkers<<endl<<endl;		/*Write header*/
	double thingy = 0;

	for(int i=0; i<M;i++){
		for(int j=0; j<N;j++){
			nwalkers = 0;
			for(int l=0; l<Conc[i][j]; l++){
				thingy = x[i]+0.99*DX*(0.5-rng->uniform());
				inifile<<0<<" "<<thingy<<" ";
				if(d==2){
					thingy = y[j]+0.99*DY*(0.5-rng->uniform());
					inifile<<thingy<<" "<<0;
				}
				else{
					inifile<<0<<" "<<0;	
				}
				inifile<<endl;
				nwalkers++;
			}
		}
	}
	inifile.close();
	// delete [] asdf;
}

void Combine::ConvertFromWalkers(double **u, string filename, int **index){
	if(debug){cout<<"Combine::ConvertFromWalkers"<<endl;}
	int M,N,m0,n0; 
	double DX=0,DY=0;
	
	ifstream infile (filename.c_str(), ios::in | ios::binary);
	string line, sub;

	M = index[0][1]-index[0][0];
	N = (d>1)?(index[1][1]-index[1][0]):(1);
	m0 = index[0][0];
	n0 = (d>1)?(index[1][0]):(0);
	DX = 1.0/(M-1);
	DY = 1.0/(N-1);
	double* pos = new double[d];
	int indx[3]; indx[0] = indx[1] = indx[2] = 0;
	// int** C = new int*[M];
	// for(int i=0; i<M;i++){
	// 	C[i] = new int[N];
	// 	for(int j=0;j<N;j++){
	// 		C[i][j] = 0;
	// 	}
	// }
	mat C = zeros(M,N);
	/*Read header which contains the number of walkers*/
	int j = 0;
	getline(infile,line);
	int nwalkers = atoi(line.c_str());
	getline(infile,line);		
	for(int i=0;i<nwalkers;i++){
		getline(infile,line);
		istringstream iss(line);
		iss >> sub;
		for (j = 0; j < 3; j++){
			iss >> sub;
			pos[j] = atof(sub.c_str());
		}
		// cout<<pos[0]<<","<<pos[1]<<endl;
		indx[0] = int(round(pos[0]/DX));
		if(d==2){
			indx[1] = int(round(pos[1]/DY));
		}
		try
		{
		// C[indx[0]][indx[1]] += 1;		/*segfault possibilities*/
		C(indx[0],indx[1]) +=1;
		// throw 20;
		}
		catch(logic_error){
			cout<<"Segmentation fault. indx = "<<indx[0]<<","<<indx[1]<<" :: "<<pos[0]<<","<<pos[1]<<endl;
			exit(1);
		}
	}
	cout<<"n0,m0 = "<<n0<<","<<m0<<endl;
	for(int k=0; k<M; k++){
		for(int j=0; j<N; j++){
			/*This is where we insert least squares or similar*/
			u[k+m0][j+n0] = 0.5*((C(k,j)/Hc)+u[k+m0][j+n0]);
			// u[k+m0][j+n0] = (C(k,j)/Hc);
			// u[k+m0][j+n0] = 0.5*((C[k][j]/Hc)+u[k+m0][j+n0]);
		}
	}
}

void Combine::MapAreaToIndex(double *x,double *y, int **index){
	if(debug){cout<<"Combine::MapAreaToIndex"<<endl;}
	/*x and y holds the limits of the area:
	x[0],y[0]; x[0],y[1]; x[1],y[0]; x[1],y[1]*/
	if(d ==1){
		index[0][0] = int(round(x[0]/dx));
		index[0][1] = int(round(x[1]/dx));

		if(index[0][0]<0 || index[0][1]>m-1){
			cout<<"Invalid index, exiting"<<endl;
			exit(1);
		}
	}
	else if(d==2){
		index[0][0] = int(round(x[0]/dx));
		index[0][1] = int(round(x[1]/dx));

		if(index[0][0]<0 || index[0][1]>m-1){
			cout<<"Invalid index, exiting"<<endl;
			exit(1);
		}
		index[1][0] = int(round(y[0]/dy));
		index[1][1] = int(round(y[1]/dy));

		if(index[1][0]<0 || index[1][1]>m-1){
			cout<<"Invalid index, exiting"<<endl;
			exit(1);
		}
	}
}

void Combine::polyreg(double *x,double *y,int m){
	/*Polynomial regression of degree < degmax*/
	int deg = m;
	int degmax = 4;
	if(m>degmax){
		deg=degmax;
	}
	mat A = zeros(deg,deg);
	vec b = zeros(deg);
	double tmp;
	for (int i=0;i<deg;i++){
		for (int l=0; l<m;l++){
			b[i] += pow(x[l],i)*y[l];
		}
		for(int j=i;j<deg;j++){
			tmp = 0;
			for(int l=0; l<m;l++){
		 		tmp += pow(x[l],(i+j));
		 	}
			A(i,j) = A(j,i) = tmp;
		}
	}
	vec a = arma::solve(A,b);
	for(int l=0; l<m;l++){
		y[l] = 0;
		for(int i=0;i<deg;i++){
			y[l] += a[i]*pow(x[l],i);
		}
	}
}

void Combine::CubicSpline(double *x, double *y, double *v, int M, double yp1, double yp2, int m0){
	/*This function is only implemented in 1d. It does a Cubic Spline interpolation of 
	y to the equally spaced mesh with spacing dx. The results are returned in v. 
	yp1 & yp2 are the first derivatives of y at the start- and end- points, or even better, 
	the first derivatives of the deterministic solution in the same points.
	Based on the function spline from lib.cpp written by Morten Hjorth Jensen*/
	int i,k;
	double inf = 1e23;
	double p, qn, sig, un, *bu;
	m0 = 0;
	bu = new double[M];
	// if(!bu) {
	// 	cout<<"Error in function CubicSpline():"<<endl;
	// 	cout<<"Not enough memory for u["<<M<<"]. Exiting"<<endl;
	// 	exit(1);
	// }
	cout<<"yp1,yp2 = "<<yp1<<","<<yp2<<endl;
	if(yp1>inf) {v[0] = bu[0] = 0.0;}
	else {
		v[0] = -0.5;
		bu[0] = (3.0/(x[m0+1]-x[m0+0]))*((y[1]-y[0])/(x[m0+1]-x[m0+0]) - yp1);
	}
	for(i=0; i<M-1; i++){
		sig = (x[m0+i] - x[m0+i-1])/(x[m0+i+1] - x[m0+i-1]);
		p = sig*v[i-1] +2.0;
		v[i] = (sig-1)/p;
		bu[i] = (y[i+1]-y[i])/(x[m0+i+1]-x[m0+i]) - (y[i]-y[i-1])/(x[m0+i]-x[m0+i-1]);
		bu[i]  = (6.0*bu[i]/(x[m0+i+1]-x[m0+i-1]) - sig*bu[i-1])/p;
	}
	if(yp2 > inf) { qn = un = 0.0;}
    else {
      qn = 0.5;
      un = (3.0/(x[m0+M-1]-x[m0+M-2]))*(yp2 - (y[M-1] - y[M-2])/(x[m0+M-1]-x[m0+M-2]));
   }
   for(k=(M-2); k>=0; k--){
   	v[k] = v[k]*v[k+1] + bu[k];
   }
   delete [] bu;
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

// void Combine::TestRWConvergence(int steps,string path){
// 	Walk walks(d,pde_solver->dt);
// 	int **distr = new int*[m];
// 	for(int i=0;i<m;i++){
// 		distr[i] = new int[n];
// 		for(int j=0;j<n;j++){
// 			distr[i][j] = Hc*Up[i][j];
// 		}
// 	}
// 	ofstream ofile;
// 	char *name = new char[120];
// 	walks.SetInitialCondition(distr,m,n);
// 	walks.SetDiffusionTensor(aD,m,n);
// 	// walks.drift = 0.05;
// 	for(int t=0;t<steps;t++){
// 		walks.InhomogenousAdvance(distr,pde_solver->dt);
// 		sprintf(name,"%s/results_FE_Hc%d_n%04d.txt",path.c_str(),(int) Hc,t);
// 		ofile.open(name);
// 		for(int i=0;i<m;i++){
// 			for(int j=0;j<n;j++){
// 				U[i][j] = distr[i][j]/Hc;
// 				ofile<<U[i][j]<<" ";
// 			}
// 			ofile<<endl;
// 		}
// 		ofile.close();
// 		walks.ResetInitialCondition(distr);
// 		cout<<"t = "<<t<<endl;
// 	}
// }


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

double Combine::norm(double **U,double **Up,int m, int n){
	double norm = 0;
	for (int i = 0; i < m; ++i){
		for (int j = 0; j < n; ++j){
			norm += (U[i][j]-Up[i][j])*(U[i][j]-Up[i][j]);
		}
	}
	return sqrt(norm);
}