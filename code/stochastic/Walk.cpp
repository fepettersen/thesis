#include "walksolver.h"

using namespace std;

bool debug_walk = false;

Walk::Walk(int dimension, double _dt)
{
	if(debug_walk){cout<<"Walk::Walk"<<endl;}
	d = dimension;
	steps = 1;		/*Should not be specified here!*/
	x0 = 0; x1 = 1;
	y0 = 0; y1 = 1;
	z0 = 0; z1 = 1;
	D = 1;
	dt = _dt/steps;
	rng = new Random();
	drift = 0; //factor/steps;
	inhomogenous = false;
};


void Walk::Advance(void){
/*This should now implement the normal advance-function as well*/
	if(debug_walk){cout<<"Walk::Advance"<<endl;}
	double DX[d];
	DX[0] = dx;
	if(d==2){DX[1] = dy;}

	int *index = new int[2];
	double L = 0;
	double L_deriv = 0;
	double L0 = sqrt(2*dt);						/*should it be sqrt(2*d*dt); ???*/
	double L_deriv0 = (dt/d)/(2*sqrt(2*(dt/d)));
	double Tr,Tl,Delta_m,Delta_p,r,stepvector[d];

	if(not inhomogenous){
		L = L0*sqrt(D);
		Tr = 0.5;
		Delta_p = L;
		Delta_m = L;
	}
	int p = 0;
	double temp= 0;
	index[1] = 0;
	for(int i=0; i<nwalkers; i++){
		/*Loop over all walkers: */
		FindPosition(walkers[i],index);
		p = int(rng->ran0()*d);		/*pick a spatial dimension to advance the walker*/
		if(inhomogenous){
			/*Calculate step-length and step-probabilities based D(x,y)*/
			temp = sqrt(aD[index[0]][index[1]]);
			L = L0*temp;
			L_deriv = L_deriv0/(DX[p]*temp)*(aDx[p][index[0]][index[1]]);
			Tr = (1+0.5*L_deriv);
			Tl = (1-0.5*L_deriv);
			Delta_p = L*Tr;
			Delta_m = L*Tl;
			Tr *= 0.5;
		}
		walkers[i][p] += (rng->ran0()>Tr)?(Delta_p):(-Delta_m);
		checkpos(walkers[i]);		/*Make sure walker does not leave area*/
	}	
}

void Walk::Load(std::string filename,int M, int N){
	if (debug_walk)	{cout<<"Walk::Load"<<endl;}
	/*Loads the file "filename" which is a (binary) .xyz file 
	with the positions of all walkers and saves them in the 
	vector walkers. More or less a replacement for the 
	SetInitialCondition(C,m,n) - function*/
	m = M;
	n = N;
	o = 0;
	x = new double[m];
	y = new double[n];
	z = new double[o];
	dx = (x1-x0)/(m-1);
	dy = (n>1)?((y1-y0)/(n-1)):0;
	dz = (o>1)?((z1-z0)/(o-1)):0;

	x0_ = x0 - 0.99*(dx/2.0);
	x1_ = x1 + 0.99*(dx/2.0);
	y0_ = y0 - 0.99*(dy/2.0);
	y1_ = y1 + 0.99*(dy/2.0);
	z0_ = z0 + 0.99*(dz/2.0);
	z1_ = z1 + 0.99*(dz/2.0);
	
	for(int k=0;k<m;k++){
		x[k] = k*dx;
	}
	for(int k=0;k<n;k++){
		y[k] = k*dy;
	}
	ifstream infile (filename.c_str(), ios::in | ios::binary);

	/*Could read the file like this. Might be faster*/	
	// file.read(reinterpret_cast<char*>(&N),sizeof(int));
 //    data_type *tmp_data = new data_type[6*N];

 //    file.read(reinterpret_cast<char*>(tmp_data), 6*N*sizeof(data_type));
 //    file.close();
 //    vector<data_type> &r = system->r;
 //    vector<data_type> &v = system->v;

 //    for(int n=0;n<N;n++) {
 //        r.at(3*n+0) = tmp_data[6*n+0];
 //        r.at(3*n+1) = tmp_data[6*n+1];
 //        r.at(3*n+2) = tmp_data[6*n+2];
 //        v.at(3*n+0) = tmp_data[6*n+3];
 //        v.at(3*n+1) = tmp_data[6*n+4];
 //        v.at(3*n+2) = tmp_data[6*n+5];
 //        int mapped_cell_index = system->cell_index_map[system->cell_index_from_position(n)];
 //        Cell *cell = system->active_cells.at(mapped_cell_index);
 //        cell->add_molecule(n,system->molecule_index_in_cell,system->molecule_cell_index);
 //    }

 //    delete filename;
 //    delete tmp_data;
	
	string line,sub;
	getline(infile,line);
	nwalkers = atoi(line.c_str());
	walkers = new double*[nwalkers];
	// walkers.resize(nwalkers);

	getline(infile,line);
	int j;
	for (int i = 0; i < nwalkers; ++i){
		walkers[i] = new double[d];
		getline(infile,line);

		istringstream iss(line);
		iss >> sub;
		for(j = 0; j<d; j++){
			iss >> sub;
			walkers[i][j] = atof(sub.c_str());
		}
	}
	infile.close();
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



void Walk::FindPosition(double *pos, int *indx){
	/*Maps the walkers position to its index*/
	indx[0] = int(round(pos[0]/dx));
	// indx[1] = (d>1)?(int(round(pos[1]/dy))):(0);
	if(d>1){
		indx[1] = int(round(pos[1]/dy));
	}
	if(d>2){
		indx[2] = int(round(pos[2]/dz));
	}
}
void Walk::checkpos(double *r){
	/*Implements reflecting boundaries.
	r[0] = x-coordinate
	r[1] = y-coordinate	...*/
	if (r[0]<x0_){
		r[0] = x0_ - (r[0]-x0_);
		}
	else if(r[0]>x1_){
		r[0] = x1_ - (r[0]-x1_);
	}
	if(d>1){
		if (r[1]<y0_){
			r[1] = y0_ - (r[1]-y0_);
		}
		else if (r[1]>y1_){
			r[1] = y1_ - (r[1]-y1_);
		}
	}
	if(d>2){
		if (r[2]<z0_){
			r[2] = z0_ - (r[2]-z0_);
		}
		else if (r[2]>z1_){
			r[2] = z1_ - (r[2]-z1_);
		}
	}
}
