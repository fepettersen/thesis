#include "main_walk.h"
using namespace std;
using namespace arma;

Dendrite::Dendrite(int M, int N, double X0, double X1, double Y0, double Y1,double **DiffTensor, double factor, double Dt):Combine(M,N,X0,X1,Y0,Y1,DiffTensor,factor,Dt){
	min_spine_contact_point = 0.09/dx;		/*spine neck diameter min = 0.09 micrometer (Arellano et.al. 2007)*/
	max_spine_contact_point = 0.51/dx;		/*spine neck diameter max = 0.51 micrometer (Arellano et.al. 2007)*/

	left_spine_pos_limit = 2;
	right_spine_pos_limit = m - 2;
	dt = Dt;
	diffuse_into_spine_probability = 0.2*dt*dx;
	RW_timesteps = 100;
	num_spines = t = 0;
	ammount_last_step = 0;

	ofstream spine_info;
	spine_info.open("spine_info.txt");
	spine_info.close();
	cout<<"diffuse_into_spine_probability = "<<diffuse_into_spine_probability<<endl;
	}

void Dendrite::ConvertToWalkers(double **u, string filename, int **index){
	/*Converts the solution, U, to a distribution of walkers for these particular 
	indeces.*/
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
	double test = 0;
	ofstream inifile (filename.c_str());
	// ofstream inifile (filename.c_str(), ios::out | ios::binary);

	
	vec x = linspace(minimum,maximum,M);
	vec y = linspace(minimum,maximum,N);
	double thingy = 0;
	for(int k=0; k<M; k++){
		Conc[k] = new int[N];
		for(int l=0; l<N; l++){
			Conc[k][l] = (int) (round(fabs(u[k+m0][l+n0]*Hc)));
			nwalkers += Conc[k][l];
		}
	}	

	inifile<<nwalkers<<endl<<endl;		/*Write xyz-header*/
	
	for(int i=0; i<M; i++){
		for(int j=0; j<N; j++){
			for(int l=0; l<Conc[i][j]; l++){
				thingy = x[i]+0.99*DX*(0.5-rng->uniform());
				Walker* tmp = new Walker(thingy,0);
				dendrite_walkers.push_back(tmp);
				inifile<<"Ar "<<thingy<<" ";
				if(d==2){
					thingy = y[i]+0.99*DY*(0.5-rng->uniform());
					inifile<<thingy<<" "<<minimum;
				}
				else{
					inifile<<minimum<<" "<<minimum;
				}
				inifile<<endl;
			}
		}
	}
	inifile.close();
}

void Dendrite::AddSpine(double drift){
	int spine_position = left_spine_pos_limit+ right_spine_pos_limit*rng->uniform();
	int spine_length_in_gridpoints = m + 1;
	int max_spine_length_this_spine = min(max_spine_contact_point,m-spine_position);

	while(spine_length_in_gridpoints >= max_spine_length_this_spine){
		spine_length_in_gridpoints = min_spine_contact_point +rng->uniform()*max_spine_contact_point;
	}
	// spine_placements.push_back(spine_position);
	Spine* tmp = new Spine(spine_position,spine_length_in_gridpoints + spine_position,0.33,dt/RW_timesteps, dx);
	spines.push_back(tmp);
	spines[num_spines]->SetDrift(drift);
	num_spines ++;
	cout<<"Added spine at "<<spine_position<<" to "<<spine_length_in_gridpoints + spine_position<<endl;
}

void Dendrite::AddWalkArea(double* x, double* y, string cmd, string path, string excecutable_name){
	/*make sure that cmd is */
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
	char tmp[100],balle[100];
	sprintf(tmp,"stochastic/inifile_%d.bin",walk_areas);
	sprintf(balle,"%s%s%s",cmd.c_str(),path.c_str(),excecutable_name.c_str());
	// sprintf(tmp,"state.xyz");
	string name = tmp;
	inifilenames.push_back(name);
	int* parameter_tmp = new int[2];
	parameter_tmp[0] = M; parameter_tmp[1] = N;
	parameters.push_back(parameter_tmp);
	command.push_back(balle);
	walk_steps = 250;
	prgm = "walk_solver";
	indeces.push_back(index);
	ConvertToWalkers(Up,name,index);
}

void Dendrite::Solve(void){
	pde_solver->advance(U,Up,m,n);
	char diffT[40];
	// for(int i=0; i<walk_areas;i++){
	// 	ConvertToWalkers(U,inifilenames[i],indeces[i]);
	// 	int failure = system(command[i].c_str());
	// 	if(failure){
	// 		cout<<endl;
	// 		cout<<"command \""<<command[i]<<"\" gave return value: "<<failure<<endl;
	// 		cout<<endl<<endl;
	// 		exit(1);
	// 	}
	// 	ConvertFromWalkers(U,inifilenames[i],indeces[i]);
	// }
		// for(vector<Walker*>::iterator Ion = dendrite_walkers.begin(); Ion != dendrite_walkers.end(); ++Ion){
		// 	/*This loop is currently not used, and is ment for only using walkers on the dendrite*/
		// 	(*Ion)->r[0] += pow(-1,int(2*rng->uniform()))*(step_length);
		// 	if ((*Ion)->r[0] < x0 - dx/2.0){
		// 		(*Ion)->r[0] += (x0-dx/2.0) - (*Ion)->r[0];
		// 	}
		// 	if ((*Ion)->r[0] < x1 + dx/2.0){
		// 		(*Ion)->r[0] += (x1 + dx/2.0) - (*Ion)->r[0];
		// 	}
		// }
	t += 1;
	double step_length = sqrt(2*d*aD[0][0]*dt);
	for(int i=0;i<RW_timesteps;i++){
		for (vector<Spine*>::iterator spine = spines.begin(); spine != spines.end(); ++spine){

			if(rng->uniform()<(*spine)->probability_factor*diffuse_into_spine_probability){
				/*Particle diffusing into spine. Probabiliy increases with spine neck width*/
				// int position = (*spine)->pos + int((*spine)->dendrite_gridpoints*rng->uniform());
				double integral = 0;
				for(int j=(*spine)->pos; j < ((*spine)->pos + (*spine)->dendrite_gridpoints);j++){
					integral += U[j][0];
				}
				if (integral*dx>=1.0/Hc){
					fprintf(stderr,"Particle diffusing into spine\n");

					(*spine)->AddIonFromDendrite();
					for(int j=(*spine)->pos; j < ((*spine)->pos + (*spine)->dendrite_gridpoints);j++){
						U[j][0] -= 1.0/((*spine)->dendrite_gridpoints*Hc);
					}
				}
				/*Pick a random random walker at relevant position and delete it*/
				// dendrite_walkers.delete()
			}
			(*spine)->Solve();
			/*Particles diffusing from spine*/
			SpineBoundary((*spine));
		}
	}
	/*Write information about spine heads to a separate file*/

	int counter = 0;
	bool flag = false;
	for(vector<Spine*>::iterator spine = spines.begin(); spine != spines.end(); ++spine){
		counter += (*spine)->ions_in_spine_head;
		double integral = 0;
		if(t>1){	
			for(int j=(*spine)->pos; j < ((*spine)->pos + (*spine)->dendrite_gridpoints);j++){
				integral += U[j][0];
			}
			if(integral*dx>=1.0/Hc && not (*spine)->is_reported){
				(*spine)->to_report = true;
				flag = true;
			}
		}
	}
	if(counter != ammount_last_step || t==1 || flag){
		ofstream spine_info;
		spine_info.open("spine_info.txt",ios_base::app);
		sprintf(diffT,"t = %04d  %d",t,int(spines.size()));
		spine_info<<diffT<<endl;
		if (t==1){
			sprintf(diffT,"dx = %g  dt = %g",dx,dt);
			spine_info<<diffT<<endl;
		}
		for(vector<Spine*>::iterator spine = spines.begin(); spine != spines.end(); ++spine){
			// cout<<"Hei, jeg har ID "<<(*spine)->pos<<endl;
			if ((*spine)->ions_in_spine_head > 0 && t>1){
				spine_info<<(*spine)->ions_in_spine_head<<"  "<<(*spine)->neck_length<<"  "<<(*spine)->pos<<"  "<<(*spine)->dendrite_gridpoints<<endl;
			}
			else if (t==1){
				spine_info<<(*spine)->ions_in_spine_head<<"  "<<(*spine)->neck_length<<"  "<<(*spine)->pos<<"  "<<(*spine)->dendrite_gridpoints<<endl;
			}
			else if((*spine)->to_report){
				spine_info<<"spine at "<<(*spine)->pos<<" has enough available"<<endl;
				(*spine)->is_reported = true;
				(*spine)->to_report = false;
			}
		}		
		spine_info.close();
	}
	ammount_last_step = counter;
	for(int k=0; k<m; k++){
		for(int l=0; l<n; l++){
			Up[k][l] = U[k][l];
		}
	}
}

void Dendrite::SpineBoundary(Spine* spine){
	int ions = 0;
	for (vector<Walker*>::iterator i = spine->dendrite_boundary.begin(); i != spine->dendrite_boundary.end(); ++i){
		U[spine->pos + int(round((*i)->r[0]/spine->dx))][0] += 1.0/Hc;
		(*i)->r[1] = 0;
		dendrite_walkers.push_back((*i));
	}
	// cout<<"Dendrite::SpineBoundary:"<<endl<<"boundary ions: "<<spine->dendrite_boundary.size()<<endl;
	spine->dendrite_boundary.clear();
}