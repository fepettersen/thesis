#include "main_walk.h"
using namespace std;
using namespace arma;

Spine::Spine(int spine_start, int spine_stop, double D, double dt){
	rng = new Random();
	double some_factor = 0.1;
	max_spike_size = 15;
	d = 2;
	dendrite_gridpoints = spine_stop - spine_start;
	pos = spine_start;
	spike_probability = 0.0;

	neck_width = some_factor*dendrite_gridpoints;
	neck_length = (1-0.5*rng->uniform())*neck_width;
	head_height = 1-neck_length;
	head_width = 0.5*(1-neck_width);
	a = head_height/head_width;
	left_neck_limit = 0.5*(1 - neck_width);
	right_neck_limit = 0.5*(1 + neck_width);
	drift = 0;
	step_length = sqrt(2*D*dt);
	
	_x1 = _y1 = 1.0;
	_x0  =_y0 = 0.0;
	dx = (dendrite_gridpoints>1)?(1.0/(dendrite_gridpoints-1)):(1.0);
	cout<<"dx = "<<dx<<endl;
	// cout<<endl<<"Values: "<<endl<<"neck_width = "<<neck_width<<endl;
	// cout<<"neck_length = "<<neck_length<<endl<<"a = "<<a<<endl;
}

void Spine::AddSpike(void){
	int amplitude = max_spike_size*rng->uniform();
	for (int i = 0; i < amplitude; ++i){
		/*Add a random walker in the proximity of the PSD*/
		Walker* tmp = new Walker(rng->uniform(),1.0-rng->uniform());
		Ions.push_back(tmp);
	}
	cout<<" of size "<<amplitude<<endl;
}

void Spine::Solve(void){
	bool delete_me;
	ions_in_spine_head = 0;
	if (rng->uniform()<spike_probability){
		cout<<"Spike";
		AddSpike();
	}
	// cout<<Ions.size()<<"  ";
	// for(auto Ion:Ions){
	for(vector<Walker*>::iterator Ion = Ions.begin(); Ion != Ions.end(); ++Ion){
		int index = int(d*rng->uniform());
		(*Ion)->r[index] += pow(-1,int(d*rng->uniform()))*(step_length);
		AddDrift((*Ion));
		delete_me = checkpos((*Ion));
		if (delete_me){
			Ions.erase(Ion);
			Ion--;
		}
	}
}

void Spine::AddIonFromDendrite(void){
	Walker* tmp = new Walker(left_neck_limit+neck_width*rng->uniform(),0.1);
	Ions.push_back(tmp);
}

bool Spine::checkpos(Walker* ion){
	if (ion->r[1]>=neck_length){
		/*The "Ion" is located in the spine head*/
		ions_in_spine_head += 1;
		double minimum_y_value_right = right_limit(ion->r[0]);
		double minimum_y_value_left = left_limit(ion->r[0]);
		if (ion->r[1]<minimum_y_value_right){
			ion->r[0] = (ion->r[1]-neck_length)/a + 0.5*(1+neck_width);
			ion->r[1] += (minimum_y_value_right - ion->r[1]);
		}
		else if (ion->r[1]<minimum_y_value_left){
			ion->r[0] = (ion->r[1]-1)/(-1*a);
			ion->r[1] += (minimum_y_value_left - ion->r[1]);
		}
		if (ion->r[1]>_y1){
			ion->r[1] = _y1 - (ion->r[1] - _y1);
		}
		if (ion->r[0]<_x0){
			ion->r[0] = _x0 + (ion->r[0] - _x0);
		}
		else if (ion->r[0]>_x1){
			ion->r[0] = _x1 - (ion->r[0] - _x1);
		}
	}
	else{
		/*The "Ion" is located in the spine neck*/
		if (ion->r[0]<left_neck_limit){
			ion->r[0] = left_neck_limit + (left_neck_limit - ion->r[0]);
		}
		else if (ion->r[0]>right_neck_limit){
			ion->r[0] = right_neck_limit - (ion->r[0] - right_neck_limit);
		}
	}
	if (ion->r[1]<=_y0){
		dendrite_boundary.push_back(ion);
		// Ions.erase(ion);
		return true;
	}
	return false;
}

void Spine::AddDrift(Walker* ion){
	ion->r[1] -= drift;
}

void Spine::Write(string filename){
	ofstream outfile (filename.c_str());
	for (int i = 0; i < Ions.size(); ++i){
		outfile<<Ions.at(i)->r[0]<<" "<<Ions.at(i)->r[1]<<endl;
	}
	outfile.close();
}