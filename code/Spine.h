#ifndef SPINE_H
#define SPINE_H
class Spine
{
public:
	Spine(int,int,double,double,double);
	~Spine();

	void SetDrift(double val){
		drift = val;
	};
	double right_limit(double x){
		// return a*(x-0.5*(1+neck_width))+neck_length;
		return b*(x - _x1) + _y1;
	};
	double left_limit(double x){
		return _y1-a*x;
	};
	bool checkpos(Walker*);
	void Solve(void);
	void AddSpike(void);
	void AddDrift(Walker*);
	void AddIonFromDendrite(void);
	void Write(std::string);

	std::vector<Walker*> Ions;
	std::vector<Walker*> dendrite_boundary;

	Random* rng;
	int max_spike_size, d, pos, dendrite_gridpoints, ions_in_spine_head;
	double dx,neck_length, probability_factor;
	bool is_reported, to_report;
	/* data */
private:
	/* data */
	double drift, neck_width, head_height, head_width, left_neck_limit, right_neck_limit, spike_probability;
	double step_length, a, b,_x0, _x1, _y0, _y1;

};
#endif // SPINE_H