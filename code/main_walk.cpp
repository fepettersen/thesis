#include "main_walk.h"
using namespace std;


string make_filename(string buffer,string filename,int n){
    //Returns a filename saying something about the particular run.
	char buff[100];
	if(filename ==""){
		sprintf(buff,"/results_FE_n%03d.txt",n);
	}
	else{
		sprintf(buff,"/%s_n%03d.txt",filename.c_str(),n);
	}
  	buffer = buff;
  	// delete(buff);
  	return buffer;
}

void output(ofstream* outfile, double **u, string buffer, string path,string filename,int n, int N){
    /*outfile is an ofstram-object letting us open a file
    **u is an armadillo-object containing the solution at time n
    **n is the timestep number
    **scheme is an integer telling what scheme is used to obtain the solution
    **N is the size of the array (in one direction)*/
    string tmp = path;
    tmp.append(make_filename(buffer,filename,n));
    outfile->open(tmp.c_str());
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            *outfile <<u[i][j]<<setprecision(16)<<"  ";
            }
        if(i<N){*outfile <<endl;}
    }
    outfile->close();
    // delete(&tmp);
}
#define PI 3,1415926535897932;



int main(int argc, char** argv)
{
    int tofile = atoi(argv[1]);
    double x0 = atof(argv[2]);
    double x1 = atof(argv[3]);
    double y0 = atof(argv[4]);
    double y1 = atof(argv[5]);

    int n = atoi(argv[6]);
    int T = atoi(argv[7]);

    string path = argv[8];
    string filename = argv[9];
    string buffer;

    ofstream outfile;
	// int T = 20;
	// int n = 11;
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
			if(i>n/2 && j<n/2){
				C[i][j] = 1000;
				Up[i][j] = PI;
			}
			else{
				C[i][j] = 0;
				Up[i][j] = 0;
			}
			U[i][j] = 0;
		}
	}
	double *x = new double[2];
	double *y = new double[2];
	x[0] = x0; x[1] = x1;
	y[0] = y0; y[1] = y1;

	Walk snow(2);
	snow.SetInitialCondition(C,n,n);

	Combine gremlin(n,n,0,1,0,1,1);
	gremlin.SetInitialCondition(Up,n,n);
	gremlin.AddWalkArea(x,y);
	for(int t=0; t<T; t++){
		cout<<"Step "<<t<<" of "<<T-1<<endl;
		gremlin.Solve();
		// snow.advance(C);
		// for(int g=0; g<n; g++)
		// 	for(int h=0; h<n; h++)
		// 		U[g][h] = C[g][h]/1000.0;
		if(tofile){
			// output(&outfile,U,buffer,path,filename,t,n);
			output(&outfile,gremlin.Up,buffer,path,filename,t,n);
		}
	}
	cout<<"Hello, world! "<<endl;
	return 0;
}