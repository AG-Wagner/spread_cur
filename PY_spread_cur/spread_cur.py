# -*- coding: utf-8 -*-
# spread_cur.py, version 08.07.2019
#
#  written bei Veit Wagner, Jacobs University Bremen gGmbH (https://www.jacobs-university.de/directory/vwagner)
#
#  This software is a phyton-code (https://www.python.org/) and is free to use.
#  If you publish data analyzed using this code or approach cite our paper below.
#  This file has to remain unchanged (code, author and paper reference information must be present).
#  The routine below calculate the spreading current analytically as described in the paper.
#
#    "Modeling of photoactive area spreading in unstructured photovoltaic cells"
#    M. Gruber, V. Jovanov, V. Wagner
#    Solar Energy Materials and Solar Cells v200 (2019) 110011
#    https://doi.org/10.1016/j.solmat.2019.110011
#
# 
#  usage:
#    (Ispread, jv, Voc, i0, Jv) = spread_cur(V, Imeas, A, C, Rsq)
#
#      V, Imeas  experimental data given as vectors for voltage [V] and current [A], respectively
#      A, C      area [m^2] and circumference [m] of device given as scalars 
#      Rsq       sheet resistance [Ohm.sq] of lateral conducting layers given as scalar
#
#  corrected experimental current values are calculated by:
#    I = Imeas - spread_cur(V, Imeas, A, C, Rsq)[0]
#
# example data:
#    V     = np.linspace(-1.0, 1.0, num=11)  # [V] applied  voltage (row-vector)
#    Imeas = np.append(np.linspace(-1e-3,0.0, num=5, endpoint=False), np.linspace(0.0,1e-2, num=6))  # [A] measured current (row-vector)
#    A     = 5e-3*5e-3  # [m2] active area (electrode overlap area)
#    C     = 4*5e-3     # [m] circumference of active area
#    Rsq   = 100e3      # [Ohm] sheet resistance of device outside active area

# tested with SciLab v3.7.0 (64-bit)  [obtainable via https://www.python.org/downloads/ ]
import sys
import numpy as np

# =================================
# --- calculate spreading current (see doi.org/10.1016/j.solmat.2019.110011) ---
# input:  V      voltage values (at least two) (np.array-vector with strictly increasing values) [V]
#         Imeas  measured curent values (np.array-vector, must contain a zero crossing) [A]
#         A      device area (scalar) [m2]
#         C      device circumference (scalar) [m]
#         Rsq    device sheet resistance (scalar) [Ohm.sq] 
# return: Tupel with data below   or  None if error occurred
#         Ispread spreading current component in Imeas (np.array vector) [A]
#         jv     vertical diode curent values (np.array-vector) [A/m2]
#         Voc    voltage corresponding to Imeas=0 (scalar) [V]
#         i0     index, voltage Voc is in interval [ V[i0]..V[i0+1] )
#         Jv     voltage integrated jv (=power) (np.array-vector) [AV/m2]
def spread_cur(V, Imeas, A, C, Rsq):
    """Calculate spreading current (see doi.org/10.1016/j.solmat.2019.110011)"""
    # --- input data check ---
    err = np.ndim(V)!=1 or np.ndim(Imeas)!=1 or np.ndim(A)!=0 or np.ndim(C)!=0 or np.ndim(Rsq)!=0 or V.size != Imeas.size or V.size < 2
    if(err):
        print('spread_cur(): Error: wrong format of input data\n', file=sys.stderr)
    else:
        err = np.all(V[1:] <= V[:-1])
        if(err):
            print('spread_cur(): Error: voltage values V not strictly increasing.', file=sys.stderr)
    if(err):
        return   # input data error

    # --- calculation ---
    # -- find (first) Imeas=0 position -> Voc, i0 --
    i0 = np.argmax(Imeas[1:] * Imeas[:-1] <= 0)
    if((i0 == 0 and Imeas[0] * Imeas[1] > 0) or (i0==Imeas.size-2)):
        print('spread_cur(): Error: cant find zero crossing before last value of Imeas .\n', file=sys.stderr)
        return    # error case "not found"
    if(Imeas[i0+1]==0):
        i0 = i0+1       # special case, have exact zero-point @ i0
        Voc = V[i0]
    else:               # general case, zero-point in interval [i0,i0+1)
        Voc = V[i0] - Imeas[i0] *(V[i0+1] - V[i0]) / (Imeas[i0+1] - Imeas[i0])  # <=> Voc = interp1(Imeas(i0:i0+1), V(i0:i0+1), 0)
    # -- deconvolute Imeas -> jv --
    jv      = np.zeros(V.size)  # init
    IvA     = Imeas/A           # precalc
    CA2_Rsq = (C/A)**2/Rsq      # precalc
    for istep in (-1,1):        # start with downwards (istep=-1), thereafter do upwards (istep=+1) from Imeas=0-position
        inxt = i0 + (istep>0)*1
        jv_  = 0
        IvA_ = 0
        dV   = V[inxt] - Voc
        dIvA = IvA[inxt] - IvA_
        sign_V_Voc = istep  # <=> sign(dV), but for dV=0 case problematic;
        while True:
            # p        = IvA[inxt] + 0.5*CA2_Rsq * dV                          # eqn 14
            # q        = IvA[inxt]**2 - (IvA_ - jv_)**2 - CA2_Rsq * jv_ * dV   # eqn 15
            # jv[inxt] = p - sign_V_Voc * np.sqrt(p*p - q)                     # eqn 13
            # identical but numerically more stable version: -> solve for diff. to best guess: jv[i] + ( Imeas[i+1] - Imeas[i] )/A
            p        = IvA_ - jv_ + 0.5*CA2_Rsq * dV
            q        = -2*CA2_Rsq * (jv_ + 0.5*dIvA) * dV
            jv[inxt] = (jv_ + dIvA) + p - sign_V_Voc * np.sqrt(p*p - q)
            # next step
            i        = inxt
            inxt     = i + istep
            if( inxt < 0 or inxt >= V.size):
                break
            jv_      = jv[i]
            IvA_     = IvA[i]
            dV       = V[inxt] - V[i]
            dIvA     = IvA[inxt] - IvA_
    # -- numerical integration of jv (starting at Voc) -> Jv --
    VV  = np.flip(np.append(V[0:i0+1],Voc))    # <- [Vstart..Voc]
    jvx = np.flip(np.append(jv[0:i0+1], 0.0))
    Jv  = np.cumsum(.5*(jvx[1:] + jvx[:-1]) * (VV[1:] - VV[:-1]))                  # eqn 12
    Jv  = np.flip(Jv)
    VV  = np.append(Voc, V[i0+1:])           # <- [Voc..Vend]
    jvx = np.append(0.0, jv[i0+1:])
    Jv  = np.append(Jv , np.cumsum(.5*(jvx[1:] + jvx[:-1]) * (VV[1:] - VV[:-1])))  # eqn 12
    # -- final formula for spreading current -> Ispread --
    Ispread  = C * np.sign(V-Voc) * np.sqrt( (2/Rsq) * Jv )                        # eqn 7
    # --- finished
    return (Ispread, jv, Voc, i0, Jv)
