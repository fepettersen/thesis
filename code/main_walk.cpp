#include "main_walk.h"
using namespace std;


string make_filename(string buffer,string filename,int step_no){
    //Returns a filename saying something about the particular run.
	char buff[100];
	if(filename ==""){
		sprintf(buff,"/results_FE_n%03d.bin",step_no);
	}
	else{
		sprintf(buff,"/%s_n%03d.bin",filename.c_str(),step_no);
	}
  	buffer = buff;
  	// delete(buff);
  	return buffer;
}

void output(ofstream* outfile, double **u, string buffer, string path,string filename,int m,int n, int N){
    /*outfile is an ofstram-object letting us open a file
    **u is an armadillo-object containing the solution at time n
    **n is the timestep number
    **scheme is an integer telling what scheme is used to obtain the solution
    **N is the size of the array (in one direction)*/
    string tmp = path;
    tmp.append(make_filename(buffer,filename,N));
    outfile->open(tmp.c_str(),ios::binary);
    // outfile->open(tmp.c_str());
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            *outfile <<u[i][j]<<setprecision(16)<<" ";
            }
        if(i<m){*outfile <<endl;}
    }
    outfile->close();
    // delete(&tmp);
}
#define PI 3.1415926535897932;



int main(int argc, char** argv)
{
    int tofile = atoi(argv[1]);
    double x0 = atof(argv[2]);
    double x1 = atof(argv[3]);
    double y0 = atof(argv[4]);
    double y1 = atof(argv[5]);

    int m = atoi(argv[6]);
    int n = atoi(argv[7]);
    int T = atoi(argv[8]);

    string path = argv[9];
    string filename = argv[10];
    double factor =atof(argv[11]);
    string buffer;

    ofstream outfile;

	double *x = new double[2];
	double *y = new double[2];
	x[0] = x0; x[1] = x1;
	y[0] = y0; y[1] = y1;
	double dx = 1.0/(m-1);
	double wth = 0;

	int **C = new int*[m];
	double **Up = new double*[m];
	double **U = new double*[m];
	double *X = new double[m];

	for(int i=0; i<m; i++){
		C[i] = new int[n];
		Up[i] = new double[n];
		U[i] = new double[n];
		X[i] = i*dx;
	}

	for(int i=0; i<m; i++){
		for(int j=0; j<n; j++){
			if(i>=m/2 && j<=n/2){
				C[i][j] = (int) (factor);
			}
			else{
				C[i][j] = 0;
				// Up[i][j] = 0;
			}
			wth = X[i]*PI;
			Up[i][j] = cos(wth) +1;
			U[i][j] = 0;
		}
	}

	// for(int k=0;k<m;k++){
	// 	cout<<Up[k][0]<<"  ";
	// }
	// cout<<endl;

	// Walk snow(2);
	// snow.SetInitialCondition(C,m,n);
	// Diffusion giraffe(dx,0,1);

	Combine gremlin(m,n,0,1,0,1,1,factor);
	gremlin.SetInitialCondition(Up,m,n);
	gremlin.AddWalkArea(x,y);

	for(int t=0; t<T; t++){
		gremlin.Solve();

		// giraffe.advance(U,Up,m,n);
		// snow.advance(C);
		// for(int g=0; g<m; g++)
		// 	for(int h=0; h<n; h++)
		// 		Up[g][h] = U[g][h];
		if(tofile){
			// output(&outfile,U,buffer,path,filename,m,n,t);
			output(&outfile,gremlin.Up,buffer,path,filename,m,n,t);
		}
	}
	return 0;
}