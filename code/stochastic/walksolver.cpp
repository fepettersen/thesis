#include "walksolver.h"
using namespace std;

void FromFile(double **array, string filename,int m, int n){
	/*Reads a .txt file containing a m*n matrix and saves the 
	matrix to array.*/
	ifstream infile (filename.c_str());
	string sub;
	string line;
	int k = 0;
	for(int i=0;i<m;i++){
		getline(infile,line);
		istringstream iss(line);
		k =0;
		for(int j=0;j<n;j++){
			sub = "";
			iss >>sub;
			array[i][k] = atof(sub.c_str());
			k++;
		}
	}
	infile.close();
}

int main(int argc, char** argv)
{
	/*These should be changed*/

    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int T = atoi(argv[3]);
    double Dt = strtod(argv[4],NULL);

    string buffer = argv[5];
    ofstream outfile;
    int d;
    (n>1)?(d=2):(d=1);
	double dx = 1.0/(m-1);
	double dy = 1.0/(n-1);

	double **aD = new double*[m];

	for(int i=0; i<m; i++){
		aD[i] = new double[n];
	}
	// cout<<"balle"<<endl;
	// FromFile(Up,"InitialCondition.txt",m,n);
	FromFile(aD,"DiffusionTensor_1.txt",m,n);
	// cout<<"balle2"<<endl;
	Walk solver(d,Dt);

	solver.Load(buffer,m,n);
	solver.SetDiffusionTensor(aD,m,n);
	
	for(int t=0; t<T; t++){
		solver.Advance();
	}

	outfile.open(buffer.c_str(),ios::out | ios::binary);
	outfile<<solver.nwalkers<<endl<<endl;
	for(int i=0; i<solver.nwalkers; i++){
		outfile<<solver.walkers[i][0]<<" ";
		if(d==2){outfile<<solver.walkers[i][1]<<" "<<0;}
		else {outfile<<0<<" "<<0;}
		outfile<<endl;
	}
	outfile.close();
	return 0;
}