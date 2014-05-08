void tridiag(double *u, double *up, int N, double *a, double *b, double *c){
    double *H = new double[N];
    double *g = new double[N];
    for(int i=0;i<N;i++){
    	H[i] = 0;
        g[i] = 0;
    }
    g[0] = up[0]/b[0];
    H[0] = c[0]/b[0];
    for(int i=1; i<N; i++){
        //forward substitution
        H[i] = -c[i]/(b[i] + a[i]*H[i-1]);
        g[i] = (up[i] - a[i]*g[i-1])/(b[i] + a[i]*H[i-1]);
    }
    u[N-1] = g[N-1];
    for(int i=(N-2);i>=0;i--){
        //Backward substitution
        u[i] = g[i] - H[i]*u[i+1];
    }
}