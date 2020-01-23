// spread_cur.c, version 08.07.2019
//
//  written bei Veit Wagner, Jacobs University Bremen gGmbH (https://www.jacobs-university.de/directory/vwagner)
//
//  This software is an C-program code (https://en.wikipedia.org/wiki/ANSI_C) and is free to use.
//  If you publish data analyzed using this code or approach cite our paper below.
//  This file has to remain unchanged (code, author and paper reference information must be present).
//  The routine below calculate the spreading current analytically as described in the paper.
//
//    "Modeling of photoactive area spreading in unstructured photovoltaic cells"
//    M. Gruber, V. Jovanov, V. Wagner
//    Solar Energy Materials and Solar Cells v200 (2019) 110011
//    https://doi.org/10.1016/j.solmat.2019.110011
//
// 
//  usage:
//    Ispread = spread_cur(V, Imeas, A, C, Rsq)
//
//      V, Imeas  experimental data given as vectors for voltage [V] and current [A], respectively
//      A, C      area [m^2] and circumference [m] of device given as scalars 
//      Rsq       sheet resistance [Ohm.sq] of lateral conducting layers given as scalar
//
//  corrected experimental current values are calculated by:
//    I = Imeas - spread_cur(V, Imeas, A, C, Rsq)
//
// example data:
//    V     = [-1:0.2:1];  // [V] applied  voltage (row-vector)
//    Imeas = [[-1e-3:2e-4:0], [2e-3:2e-3:1e-2]];  // [A] measured current (row-vector)
//    A     = 5e-3*5e-3; // [m2] active area (electrode overlap area)
//    C     = 4*5e-3;    // [m] circumference of active area
//    Rsq   = 100e3;     // [Ohm] sheet resistance of device outside active area

// tested with SciLab v6.0.2 (64-bit)  [obtainable via https://www.scilab.org/download/6.0.2 ]

// =================================
// --- calculate spreading current (see doi.org/10.1016/j.solmat.2019.110011) ---
// input:	n		# of points (must be >=2) in V, Imeas, and available output space in Ispread, jv, Jv
//			V		voltage values (vector with strictly increasing values) [V]
//			Imeas	measured curent values (vector, must contain a zero crossing) [A]
//			A		device area (scalar) [m2]
//			C		device circumference (scalar) [m]
//			Rsq		device sheet resistance (scalar) [Ohm.sq] 
// output:	Ispread	spreading current component in Imeas (vector) [A]
//			jv		vertical diode curent values (vector) [A/m2]
//			*pVoc	voltage corresponding to Imeas=0 (scalar) [V]      , is [] if not found
//			*pi0	index, voltage Voc is in interval [ V(i0)..V(i0+1) ), is [] if not found
//			Jv		voltage integrated jv (=power) (vector) [AV/m2]
// return: 0=ok, <0=error  (values of output variables or undefined if error is indicated)
#include <stdio.h>
#include <math.h>

#include "spread_cur.h"

//function [Ispread, jv, Voc, i0, Jv] = spread_cur(V, Imeas, A, C, Rsq)
int spread_cur(int n, const double *V, const double *Imeas, double A, double C, double Rsq,
			   double *Ispread, double *jv, double *pVoc, int *pi0, double *Jv)
{
int	i;

    // --- input data check ---
    if(n < 2) {fprintf(stderr, "spread_cur(): Error: wrong format of input data\n"); return -1;}
	for(i=0; i<n-1; i++)
		if( V[i+1] <= V[i] ) {fprintf(stderr, "spread_cur(): Error: voltage values V not strictly increasing.\n"); return -2;}

	// --- calculation ---
	// -- find (first) Imeas=0 position -> Voc, i0 --
	for(i=0; i<n-1; i++)
		if( Imeas[i+1] * Imeas[i] <= 0) break;
	*pi0 = i;
	if(*pi0 >= n-2) {
		fprintf(stderr, "spread_cur(): Error: can't find zero crossing before last value of Imeas .\n");
		return -3;	 // error case "not found"
	}
	if(Imeas[*pi0+1]==0) {
		(*pi0)++;      // special case, have exact zero-point @ i0
		*pVoc = V[*pi0];
	} else {           // general case, zero-point in interval [i0,i0+1)
		*pVoc = V[*pi0] - Imeas[*pi0] * (V[*pi0+1] - V[*pi0]) / (Imeas[*pi0+1] - Imeas[*pi0]); // <=> Voc = interp1(Imeas(i0:i0+1), V(i0:i0+1), 0)
	}
    // -- deconvolute Imeas -> jv --
    //x jv      = zeros(V);     // init
    //x IvA     = Imeas/A;      // precalc
	double CA2_Rsq = (C*C/(A*A))/Rsq;	// precalc
	for(int istep = -1; istep <=1; istep+=2) {	// start with downwards (istep=-1), thereafter do upwards (istep=+1) from Imeas=0-position
		int		inxt = *pi0 + ((istep>0) ? 1 : 0); 
		double	jv_  = 0;
		double	IvA_ = 0;
		double	dV   = V[inxt] - *pVoc;
		double	dIvA = Imeas[inxt]/A - IvA_;
		int		sign_V_Voc = istep; // <=> sign(dV), but for dV=0 case problematic;
		while(1) {
			// double p        = Imeas[inxt]/A + 0.5*CA2_Rsq * dV;                        // eqn 14
			// double q        = Imeas[inxt]*Imeas[inxt]/(A*A) - (IvA_ - jv_)*(IvA_ - jv_) - CA2_Rsq * jv_ * dV;   // eqn 15
			// jv[inxt] = p - sign_V_Voc * sqrt(p*p - q);                      // eqn 13
			// identical but numerically more stable version: -> solve for diff. to best guess: jv(i) + ( Imeas(i+1) - Imeas(i) )/A
			double	p        = IvA_ - jv_ + 0.5*CA2_Rsq * dV;
			double	q        = -2*CA2_Rsq * (jv_ + 0.5*dIvA) * dV;
			jv[inxt] = (jv_ + dIvA) + p - sign_V_Voc * sqrt(p*p - q);
			// next step
			i        = inxt;
			inxt     = i + istep;
			if( inxt < 0 || inxt >= n) break;
			jv_      = jv[i];
			IvA_     = Imeas[i]/A;
			dV       = V[inxt] - V[i];
			dIvA     = Imeas[inxt]/A - IvA_;
		}
	}
	
    // -- numerical integration of jv (starting at Voc) -> Jv --
	double	Vlast, jvlast, Jvlast;
	// [Vstart..Voc] :
	Vlast  = *pVoc;
	jvlast = 0;
	Jvlast = 0;
	for(i=*pi0; i>=0; i--) {
		Jv[i]  = Jvlast + .5 * (jv[i] + jvlast) * (V[i] - Vlast);	// eqn 12
		Vlast  = V[i];
		jvlast = jv[i];
		Jvlast = Jv[i];
	}
	// [Voc..Vend] :
	Vlast  = *pVoc;
	jvlast = 0;
	Jvlast = 0;
	for(i=*pi0+1; i<n; i++) {
		Jv[i]  = Jvlast + .5 * (jv[i] + jvlast) * (V[i] - Vlast);	// eqn 12
		Vlast  = V[i];
		jvlast = jv[i];
		Jvlast = Jv[i];
	}

	// -- final formula for spreading current -> Ispread --
	for(i=0; i<n; i++)
		Ispread[i]  = C * ((V[i] >= *pVoc) ? 1 : -1) * sqrt( (2/Rsq) * Jv[i] );	// eqn 7
	
	return 0;
}

