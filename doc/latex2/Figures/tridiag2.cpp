g[0] = up[0]/b[0];
H[0] = c[0]/b[0];
for(int i=1; i<n; i++){
    //forward substitution
    H[i] = -c[i]/(b[i] + a[i]*H[i-1]);
    g[i] = (up[i] - a[i]*g[i-1])/(b[i] + a[i]*H[i-1]);
}
u[n-1] = g[n-1];
for(int i=(n-2);i>=0;i--){
    //Backward substitution
    u[i] = g[i] - H[i]*u[i+1];
}