#ifndef TRIDIAG_H
#define TRIDIAG_H

class Tridiag
{
public:
	Tridiag(void);
	// void precondition(arma::mat, int, int);
	void precondition(std::vector<arma::mat> A, std::vector<arma::mat> B,std::vector<arma::mat> C, int m, int n);
	arma::vec blockTridiag(arma::mat, arma::vec, int, int);
	arma::vec efficient_blockTridiag(arma::vec,int, int);
	void tridiag(double *u, double *f, int N, double *a, double *b, double *c);
	/* data */
	std::vector<arma::mat> D,H,aa;
};


#endif // TRIDIAG_H