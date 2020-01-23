// spread_cur_test.c
// example usage of spread_cur.c
//
// under linux use Makefile or 
// compile by:  gcc -o spread_cur_test -lm spread_cur_test.c spread_cur.c
// execute:     ./spread_cur_test

#include <stdio.h>
#include "spread_cur.h"

int main()
{
#define NPTS	11		// # of data points

	// output variables
	double	Ispread[NPTS], jv[NPTS], Jv[NPTS];
	double	Voc;
	int		i0;
	
	// input variables
	// [V] applied  voltage (vector) :
	double V[NPTS]		= {-1.0, -0.8, -0.6, -0.4, -0.2,  0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
	// [A] measured current (vector) :
	double Imeas[NPTS]	= {-0.001, -0.0008, -0.0006, -0.0004, -0.0002, 0.0, 0.002, 0.004, 0.006, 0.008, 0.01};
	double A			= 5e-3*5e-3; // [m2] active area (electrode overlap area)
	double C			= 4*5e-3;    // [m] circumference of active area
	double Rsq			= 100e3;     // [Ohm] sheet resistance of device outside active area

	
	// ==========================================================================
	//            function call for spreading calculation 
	// ==========================================================================
	int retc = spread_cur(NPTS, V, Imeas, A, C, Rsq, Ispread, jv, &Voc, &i0, Jv);
	// ==========================================================================
	//
	// ==========================================================================
	
	if(retc<0) {
		fprintf(stderr, "error exit (retc = %d)!\n", retc);
		return retc;
	}

	printf("\n spread_cur.c test :\n\n");
	printf("A   = %.2le m2\n", A);
	printf("C   = %.2le m\n", C);
	printf("Rsq = %.2le Ohm\n", Rsq);
	printf("\n");	
	printf("Voc = %.6lf V\n", Voc);
	printf("voltage    measured     true value (= w/o spreading)\n");
	for(int i=0; i<NPTS; i++) printf("%5.2lf V  %7.3lf mA -> %11.7lf mA\n", V[i], Imeas[i]*1e3, (Imeas[i] - Ispread[i])*1e3);

	printf("\n<press return>");
	getc(stdin);
	
	return 0;
}
