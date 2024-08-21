#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>

int largerornot(double x,double y, int tiesmethod) // Include tiesmethod
{
double a;
if (x < y) {return(0);}
if (x > y) {return(1);}
if (x == y) {
    if (tiesmethod == 1) { // New method "first"; in favor of increasing patterns
        return(0);
    }
    else if (tiesmethod == 0) { // Original method "random"
        GetRNGstate();
        a=unif_rand();
        PutRNGstate();
        if (a> 0.5) {return(1);}
        return(0);
    }
	}
return(0);
}




static void verg(int n, int p, double*data, double*aus, int tiesmethod) {
int speicherve[p*(p+1)/2];
int w;
int i;
int j;
int k;

for (j=0;j<p-1;j++) {
	for (i=j+1;i<p;i++) {
		w = largerornot(data[j],data[i],tiesmethod);
		speicherve[j*p-j*(j-1)/2+i-j] = 1-w;
		speicherve[j*p-j*(j+1)/2+j] = w+speicherve[j*p-j*(j+1)/2+j];
		aus[i*(n-p+1)] = aus[i*(n-p+1)]+1-w;
		aus[j*(n-p+1)] = aus[j*(n-p+1)]+w;
		}
	}

for (i=1;i<n-p+1;i++){
	for(j=0;j<p-1;j++) {
		aus[i+j*(n-p+1)] = aus[(i-1)+(j+1)*(n-p+1)]-speicherve[j+1];
		}
	for(j=0;j<p-1;j++) {
		for(k=0;k<p-j-1;k++) {
			speicherve[j*p-j*(j+1)/2+k+j] = speicherve[(j+1)*p-(j+1)*(j+2)/2+k+(j+1)];
			}
		}
    for (j=0;j<p-1;j++) {
	w = largerornot(data[i+p-1],data[i+j],tiesmethod);
	speicherve[(j+1)*p-(j+1)*(j+2)/2+j] = w;
	aus[i+(p-1)*(n-p+1)] = w+aus[i+(p-1)*(n-p+1)];
	speicherve[j*p-j*(j-1)/2] = speicherve[j*p-j*(j-1)/2]+1-w;
        aus[i+j*(n-p+1)] = aus[i+j*(n-p+1)]+1-w;
	}	
    }
}


SEXP vergleich(SEXP data, SEXP ergebnis, SEXP vektor, SEXP tiesmethod)
{
    int n = LENGTH(data);
    int p = LENGTH(vektor);
    SEXP ans = duplicate(ergebnis);
    PROTECT(ans);
    verg(n,p,REAL(data),REAL(ans),REAL(tiesmethod)[0]);
    UNPROTECT(1);
    return ans;
}

static const R_CallMethodDef CallEntries[] = {
    {"vergleich", (DL_FUNC) &vergleich, 4},
    {NULL, NULL, 0}
};


void R_init_ordinalpattern(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}



