#ifndef TRIDIAG_H
#define TRIDIAG_H

class Tridiag
{
public:
	Tridiag(void);
	void precondition(arma::mat, int, int);
	arma::vec blockTridiag(arma::mat, arma::vec, int, int);
	arma::vec efficient_blockTridiag(arma::vec,int, int);
	void tridiag(double *u, double *f, int N, double *a, double *b, double *c);
	/* data */
	std::vector<arma::mat> D,H,aa;
};


#endif // TRIDIAG_H