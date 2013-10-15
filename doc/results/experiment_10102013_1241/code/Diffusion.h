class Diffusion
{
	public:
		Diffusion(double dx, double dy, double D);

		int d;
		double dt,dx,dy;
		double D,_D;	/*should be double **D ?*/
		int solver; /*Which solver to use*/
		int t; 		/*Timestep number*/

		void advance(double **U,double **Up, int, int);
		void boundary(double **,double **, int, int);
};