#ifndef DENDRITE_H
#define DENDRITE_H

class Dendrite: public Combine
{
public:
	Dendrite(int, int, double, double, double, double, double**, double, double);
	void ConvertToWalkers(double**, std::string, int**);
	void SetMin(double value){
		minimum = value;
	}
	void SetMax(double value){
		maximum = value;
	}
	void AddWalkArea(double*, double*, std::string, std::string, std::string);
	void Solve(void);
	void AddSpine(double drift=0);
	void SpineBoundary(Spine*);
	double minimum,maximum,dt;
	
	int max_spine_contact_point, num_spines, left_spine_pos_limit, right_spine_pos_limit;
	double diffusie_into_spine_probability;

	std::vector<Spine*> spines;
	std::vector<int*> spine_placements;
	std::vector<std::string> command;
	std::vector<Walker*> dendrite_walkers;
};
#endif // DENDRITE_H