#ifndef DIFFUSION_H
#define DIFFUSION_H

class Diffusion
{
	public:
		Diffusion(double dx, double dy, double D, double Dt, double _v=0);
		Diffusion(double dx, double dy, double **D, double Dt, double _v=0);

		double **aD;							/*anisotropic D*/
		int d;									/*spatial dimension*/
		double dt,dx,dy, v;
		double D,_Dx,_Dy,vdtdx2,vdtdy2;
		int solver; 							/*Which solver to use*/
		int t; 									/*Timestep number*/
		int isotropic;							/*set to zero if not*/
		arma::mat Lower, Upper, Permutation;	/*used for the LU decomposition*/
		double alpha,beta;						/*Saved in case dt is changed*/

		void advance(double **U,double **Up, int, int);
		void boundary(double **,double **, int, int);
		double f(double x,double y, double t);
		void BE2D(double **U, double **Up, int m, int n);
		void tridiag(double *u, double *up, int N, double *di, double *ab, double *bel);
		arma::mat Assemble(double, double, int, int);
		arma::mat AssembleAnisotropic(double a, double b, int m, int n);
};

#endif // DIFFUSION_H