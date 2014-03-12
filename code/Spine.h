#ifndef SPINE_H
#define SPINE_H
class Spine
{
public:
	Spine(int,int,double,double);
	~Spine();

	void SetDrift(double val){
		drift = val;
	};
	double right_limit(double x){
		return a*(x-0.5*(1+neck_width))+neck_length;
	};
	double left_limit(double x){
		return 1.0-a*x;
	};
	void checkpos(Walker*);
	void Solve(void);
	void AddSpike(void);
	void AddDrift(Walker*);
	void AddIonFromDendrite(void);

	std::vector<Walker*> Ions;
	std::vector<Walker*> dendrite_boundary;

	Random* rng;
	int max_spike_size, d, pos, dendrite_gridpoints;
	double dx;
	/* data */
private:
	/* data */
	double drift, neck_width, neck_length, head_height, head_width, left_neck_limit, right_neck_limit, spike_probability;
	double step_length, a, _x0, _x1, _y0, _y1;

};
#endif // SPINE_H