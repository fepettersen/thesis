#include "Walk.h"
#include <iostream>
using namespace std;

int main()
{
	int T = 20;
	int n = 11;
	int **C = new int*[n];

	for(int i=0; i<n; i++){
		C[i] = new int[n];
	}
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(i<n/2 && j<n/2){
				C[i][j] = 1;
			}
			else{
				C[i][j] = 0;
			}
		}
	}
	Walk balle(2);
	balle.SetInitialCondition(C,n,n);
	for(int t=0; t<T; t++){
		C = balle.advance(C);
	}
	cout<<"Hello, world!"<<endl;
	return 0;
}