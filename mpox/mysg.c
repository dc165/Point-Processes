#include <R.h>#include <Rmath.h>void sgc (double *lon, double *lat, double *t, int *n,	  double *T, double *theta, int *grid, int *m, double *result){       int i,j,w;    double sum1, r2, mu, K, a, b, lam, area;    mu = theta[0]; K = theta[1]; a = theta[2]; b = theta[3];    *result = 0.0;    area = *T / *m;    

    for(i=0; i < *m; i++){
       sum1 = 0.0;  
       for(j=0; j < *n; j++){
           if(grid[j] == i){ 
		lam = mu; 
		for(w = 0; w < j; w++){
			r2 = (lon[j]-lon[w])*(lon[j]-lon[w]) + (lat[j]-lat[w])*(lat[j]-lat[w]); 
			lam += K * b * exp(-1.0 * b * (t[j]-t[w])) * a / 3.141593 * exp(-1.0 * a * r2);
			}
		sum1 += 1/lam; 
	        } 
           }
       }
    *result += (sum1 - area)*(sum1-area); 
}
