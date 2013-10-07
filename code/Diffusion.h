class Diffusion
{
	public:
		int d;
		double dt,dx,dy;
		double D,_D;	/*should be double **D ?*/
		int solver; /*Which solver to use*/
		int t; 		/*Timestep number*/

		Diffusion();
		double **advance(double **U,double **Up, int, int);
		double **boundary(double **,int, int);
};