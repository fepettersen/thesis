class Dendrite: public Combine
{
public:
	Dendrite(int, int, double, double, double, double, double**, double, double);
	void ConvertToWalkers(double**, std::string, int**);
	void setmin(double value){
		minimum = value;
	}
	void setmax(double value){
		maximum = value;
	}
	void AddWalkArea(double*, double*, std::string, std::string, std::string);
	void Solve(void);
	double minimum,maximum;
	std::vector<std::string> command;
};