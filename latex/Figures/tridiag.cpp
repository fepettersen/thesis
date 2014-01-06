void tridiag(double *u, double *f, int N, double *a, double *b, double *c){
    double *temp = new double[N];
    for(int i=0;i<N;i++){
    	temp[i] = 0;
    }
    double btemp = b[0];
    u[0] = f[0]/btemp;
    for(int i=1; i<N; i++){
        //forward substitution
       	temp[i] = c[i-1]/btemp;
       	btemp = b[i]-a[i]*temp[i];
        u[i] = (f[i] -a[i]*u[i-1])/btemp;
    }
    for(int i=(N-2);i>=0;i--){
        //Backward substitution
        u[i] -= temp[i+1]*u[i+1];
    }
    delete[] temp;
}
int main()
{
    return 0;
}