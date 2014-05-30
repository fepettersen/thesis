#include "main_walk.h"
using namespace std;
using namespace arma;

string make_filename(string buffer,string filename,int conversion_factor,int step_no){
    //Returns a filename saying something about the particular run.
	char buff[120];
	if(filename =="a"){
		sprintf(buff,"/results_FE_Hc%d_n%04d.txt",conversion_factor,step_no);
	}
	else{
		sprintf(buff,"/%s_n%04d.txt",filename.c_str(),step_no);
	}
  	buffer = buff;
  	// delete[] buff;
  	return buffer;
}

void output(ofstream* outfile, double **u, string buffer, string path,string filename,int m,int n,int conversion_factor, int N){
    /*outfile is an ofstram-object letting us open a file
    **u is a double** containing the solution at time n
    **n is the timestep number
    **N is the size of the array (in one direction)*/
    string tmp = path;
    tmp.append(make_filename(buffer,filename,conversion_factor,N));
    // outfile->open(tmp.c_str(),ios::binary);
    outfile->open(tmp.c_str());
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            *outfile <<u[i][j]<<setprecision(16)<<" ";
            }
        if(i<m){*outfile <<endl;}
    }
    outfile->close();
    // delete(&tmp);
}

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
			// array[i][k] = strtod(sub.c_str(),NULL);
			array[i][k] = atof(sub.c_str());
			k++;
		}
	}
	infile.close();
}
#define PI 3.1415926535897932;



int main(int argc, char** argv)
{
	cout<<"initializing... ";
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
    double conversion_factor =atof(argv[11]);
    double Dt = strtod(argv[12],NULL);

    string buffer;
    ofstream outfile;
	double *x = new double[2];
	double *y = new double[2];
	x[0] = x0; x[1] = x1;
	y[0] = y0; y[1] = y1;

	double x_start = 0;
	double x_end = 1.0;

	double dx = (x_end-x_start)/(m-1);
	double dy = 1.0/(n-1);

	double **Up = new double*[m];
	double **aD = new double*[m];

	for(int i=0; i<m; i++){
		Up[i] = new double[n];
		aD[i] = new double[n];
	}

	FromFile(Up,"InitialCondition.txt",m,n);
	FromFile(aD,"DiffusionTensor.txt",m,n);

	string RWname = "spine";
	bool test_convergence = false;

	// Spine* s = new Spine(5,8,1.0,0.01);
	// s->AddSpike();

	// for (int i = 0; i < T; ++i){
	// 	s->Solve();
	// 	s->Write(make_filename(buffer,RWname,conversion_factor,i));
	// }
	Combine BlackBox(m,n,x_start,x_end,0,1,aD,conversion_factor,Dt);
	// Dendrite BlackBox(m,n,x_start,x_end,0,1,aD,conversion_factor,Dt);
	BlackBox.SetInitialCondition(Up,m,n);

	BlackBox.AddWalkArea(x,y);
	
/*	for (int i = 0; i < 20; ++i){
		// Adds some number of spines to the dendrite
		BlackBox.AddSpine();
	}*/
	
	cout<<"done!"<<endl;
	for(int t=0; t<T; t++){
		BlackBox.Solve();
		if(tofile){
			// output(&outfile,BlackBox.U,buffer,path,filename,m,n,int(1/Dt),t);
			output(&outfile,BlackBox.U,buffer,path,filename,m,n,conversion_factor,t);
		}
		cout<<"t = "<<t+1<<" of "<<T<<endl;
	}

	return 0;
}