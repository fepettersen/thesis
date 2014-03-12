#include "main_walk.h"
using namespace std;
using namespace arma;

Spine::Spine(int spine_start, int spine_stop, double D, double dt){
	rng = new Random();
	double some_factor = 0.1;
	max_spike_size = 15;
	d = 2;
	dendrite_gridpoints = spine_stop - spine_start;

	neck_width = some_factor*dendrite_gridpoints;
	neck_length = (1-rng->uniform())*neck_width;
	head_height = 1-neck_length;
	head_width = 0.5*(1-neck_width);
	spike_probability = 1e-5;
	a = head_height/head_width;
	left_neck_limit = 0.5*(1 - neck_width);
	right_neck_limit = 0.5*(1 + neck_width);
	drift = 0;
	step_length = sqrt(2*D*dt);
	
	_x1 = _y1 = 1.0;
	_x0  =_y0 = 0.0;
	dx = 1.0/(dendrite_gridpoints-1);
}

void Spine::AddSpike(void){
	int amplitude = max_spike_size*rng->uniform();
	for (int i = 0; i < amplitude; ++i){
		/*Add a random walker in the proximity of the PSD*/
		Walker* tmp = new Walker(rng->uniform(),1-rng->uniform());
		Ions.push_back(tmp);
	}
}

void Spine::Solve(void){
	if (rng->uniform()<spike_probability){
		AddSpike();
	}
	for(vector<Walker*>::iterator Ion = Ions.begin(); Ion != Ions.end(); ++Ion){
		(*Ion)->r[int((d+1)*rng->uniform())] += pow(-1,int((d+1)*rng->uniform()))*(step_length);
		AddDrift((*Ion));
		checkpos((*Ion));
	}
}

void Spine::AddIonFromDendrite(void){
	Walker* tmp = new Walker(left_neck_limit+neck_width*rng->uniform(),0.1);
	Ions.push_back(tmp);
}

void Spine::checkpos(Walker* ion){
	if (ion->r[1]>=neck_length){
		double tmpr = right_limit(ion->r[0]);
		double tmpl = left_limit(ion->r[0]);
		if (ion->r[1]<tmpr){
			ion->r[0] -= (ion->r[1]-tmpr)/a;
			ion->r[1] += (ion->r[1]-tmpr);
		}
		else if (ion->r[1]<tmpl){
			ion->r[0] += (ion->r[1]-tmpr)/a;
			ion->r[1] += (ion->r[1]-tmpr);
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
		if (ion->r[0]<left_neck_limit){
			ion->r[0] = left_neck_limit + (left_neck_limit - ion->r[0]);
		}
		else if (ion->r[0]>right_neck_limit){
			ion->r[0] = right_neck_limit - (ion->r[0] - right_neck_limit);
		}
	}
	if (ion->r[1]<=_y0){
		dendrite_boundary.push_back(ion);
		// Ions.remove(ion);
	}
}

void Spine::AddDrift(Walker* ion){
	ion->r[1] -= drift;
}