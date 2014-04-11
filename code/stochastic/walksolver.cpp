#include "walksolver.h"
using namespace std;

void FromFile(double **array, string filename,int m, int n){
	/*Reads a .txt file containing a m*n matrix and saves the 
	matrix to array.*/
	ifstream infile (filename.c_str());
	string sub;
	string line;
	for(int i=0;i<m;i++){
		getline(infile,line);
		istringstream iss(line);
		for(int j=0;j<n;j++){
			iss >>sub;
			array[i][j] = atof(sub.c_str());
		}
	}
	infile.close();
}

int main(int argc, char** argv)
{

    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int T = atoi(argv[3]);
    double Dt = strtod(argv[4],NULL);

    Dt /=T;
    string infile = argv[5];
    string DiffusionTensor = argv[6];

    ofstream outfile;
    int d;
    (n>1)?(d=2):(d=1);
	double dx = 1.0/(m-1);
	double dy = 1.0/(n-1);

	double **aD = new double*[m];

	for(int i=0; i<m; i++){
		aD[i] = new double[n];
	}

	FromFile(aD,DiffusionTensor.c_str(),m,n);

	Walk solver(d,Dt);

	solver.Load(infile,m,n);
	solver.SetDiffusionTensor(aD,m,n);
	for (int i = 0; i < m; ++i){
		for (int j = 0; j < n; ++j)	{
		}
	}
	int t=0;
    // cout<<"RW-solver; dt = "<<Dt<<", running "<<T<<" steps...";
	for(t=0; t<T; t++){
		solver.Advance();
	}
	// cout<<"..done! "<<t<<" steps completed"<<endl;
	outfile.open(infile.c_str(),ios::out | ios::binary);
	outfile<<solver.nwalkers<<endl<<endl;
	for(int i=0; i<solver.nwalkers; i++){
		outfile<<0<<" "<<solver.walkers[i][0]<<" ";
		if(d==2){outfile<<solver.walkers[i][1]<<" "<<0;}
		else {outfile<<0<<" "<<0;}
		outfile<<endl;
	}
	outfile.close();
	return 0;
}