#include "main_walk.h"
using namespace std;


char *make_filename(char* buffer,int n){
    //Returns a filename saying something about the particular run.
	sprintf(buffer,"results_FE_n%03d.txt",n);
    return buffer;
}

void output(ofstream* outfile, double **u, char* buffer,int n, int N){
    /*outfile is an ofstram-object letting us open a file
    **u is an armadillo-object containing the solution at time n
    **n is the timestep number
    **scheme is an integer telling what scheme is used to obtain the solution
    **N is the size of the array (in one direction)*/
    outfile->open(make_filename(buffer,n));
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            *outfile <<u[i][j]<<setprecision(16)<<"  ";
            }
        if(i<N){*outfile <<endl;}
    }
    outfile->close();
}
#define PI 3,1415926535897932; //38462643383279502884197169399375105820974944592307816406;
int main()
{
    char* buffer = new char[60];
    ofstream outfile;
	int T = 20;
	int n = 11;
	int **C = new int*[n];
	double **Up = new double*[n];
	double **U = new double*[n];

	for(int i=0; i<n; i++){
		C[i] = new int[n];
		Up[i] = new double[n];
		U[i] = new double[n];
	}
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(i<n/2 && j<n/2){
				C[i][j] = 10;
				Up[i][j] = PI;
			}
			else{
				C[i][j] = 0;
				Up[i][j] = 0;
			}
			U[i][j] = 0;
		}
	}
	Walk balle(2);
	Diffusion solver(0.1,0.1,1.0);
	balle.SetInitialCondition(C,n,n);

	for(int t=0; t<T; t++){
		C = balle.advance(C);
		solver.advance(U,Up,n,n);
		// output(&outfile,Up,buffer,t,n);
		for(int i=0; i<n; i++){
			for(int j=0; j<n; j++){
				Up[i][j] = U[i][j];
			}
		}
	}
	cout<<"Hello, world! "<<endl;
	return 0;
}