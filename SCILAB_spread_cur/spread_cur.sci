// spread_cur.sci, version 08.07.2019
//
//  written bei Veit Wagner, Jacobs University Bremen gGmbH (https://www.jacobs-university.de/directory/vwagner)
//
//  This software is a Scilab-code (https://www.scilab.org) and is free to use.
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
// input:  V      voltage values (at least two) (row-vector with strictly increasing values) [V]
//         Imeas  measured curent values (row-vector, must contain a zero crossing) [A]
//         A      device area (scalar) [m2]
//         C      device circumference (scalar) [m]
//         Rsq    device sheet resistance (scalar) [Ohm.sq] 
// output: Ispread spreading current component in Imeas (row vector) [A]
//         jv     vertical diode curent values (row-vector) [A/m2]
//         Voc    voltage corresponding to Imeas=0 (scalar) [V]      , is [] if not found
//         i0     index, voltage Voc is in interval [ V(i0)..V(i0+1) ), is [] if not found
//         Jv     voltage integrated jv (=power) (row-vector) [AV/m2]
function [Ispread, jv, Voc, i0, Jv] = spread_cur(V, Imeas, A, C, Rsq)

    // --- input data check ---
    err = or([size(V,1),size(Imeas,1),size(A),size(C),size(Rsq)]~=1) | size(V,2) ~= size(Imeas,2) | size(V,2) < 2;
    if(err) then disp("spread_cur(): Error: wrong format of input data");
    else
        err = or(V(2:$) <= V(1:$-1));
        if(err) then disp("spread_cur(): Error: voltage values V not strictly increasing."); end
    end
    if(err) then Ispread=[]; jv=[]; Voc=[]; i0=[]; Jv=[]; return; end

    // --- calculation ---
    // -- find (first) Imeas=0 position -> Voc, i0 --
    i0 = find(Imeas(2:$).*Imeas(1:$-1) <= 0, 1);
    if(i0 == [] | i0==size(Imeas,2)-1) then disp("spread_cur(): Error: can''t find zero crossing before last value of Imeas ."); Ispread=[]; jv=[]; Voc=[]; Jv=[]; return; end // error case "not found"
    if(Imeas(i0+1)==0) then
        i0 = i0+1;      // special case, have exact zero-point @ i0
        Voc = V(i0);
    else                // general case, zero-point in interval [i0,i0+1)
        Voc = V(i0) - Imeas(i0) *(V(i0+1) - V(i0)) / (Imeas(i0+1) - Imeas(i0)); // <=> Voc = interp1(Imeas(i0:i0+1), V(i0:i0+1), 0)
    end
    // -- deconvolute Imeas -> jv --
    jv      = zeros(V);     // init
    IvA     = Imeas/A;      // precalc
    CA2_Rsq = (C/A)^2/Rsq;  // precalc
    for istep = -1:2:1     // start with downwards (istep=-1), thereafter do upwards (istep=+1) from Imeas=0-position
        inxt = i0 + (istep>0)*1; 
        jv_  = 0;
        IvA_ = 0;
        dV   = V(inxt) - Voc;
        dIvA = IvA(inxt) - IvA_;
        sign_V_Voc = istep; // <=> sign(dV), but for dV=0 case problematic;
        while(%T)
            // p        = IvA(inxt) + 0.5*CA2_Rsq * dV;                        // eqn 14
            // q        = IvA(inxt)^2 - (IvA_ - jv_)^2 - CA2_Rsq * jv_ * dV;   // eqn 15
            // jv(inxt) = p - sign_V_Voc * sqrt(p*p - q);                      // eqn 13
            // identical but numerically more stable version: -> solve for diff. to best guess: jv(i) + ( Imeas(i+1) - Imeas(i) )/A
            p        = IvA_ - jv_ + 0.5*CA2_Rsq * dV;
            q        = -2*CA2_Rsq * (jv_ + 0.5*dIvA) * dV;
            jv(inxt) = (jv_ + dIvA) + p - sign_V_Voc * sqrt(p*p - q);
            // next step
            i        = inxt;
            inxt     = i + istep;
            if( inxt < 1 | inxt > size(V,'*')) then break; end
            jv_      = jv(i);
            IvA_     = IvA(i);
            dV       = V(inxt) - V(i);
            dIvA     = IvA(inxt) - IvA_;
        end
    end
    // -- numerical integration of jv (starting at Voc) -> Jv --
    VV  = [Voc, V(i0:-1:1)];    // <- [Vstart..Voc]
    jvx = [0, jv(i0:-1:1)];
    Jv  = cumsum(.5*(jvx(2:$) + jvx(1:$-1)) .* (VV(2:$) - VV(1:$-1)));         // eqn 12
    Jv  = Jv($:-1:1);
    VV  = [Voc, V(i0+1:$)];     // <- [Voc..Vend]
    jvx = [0, jv(i0+1:$)];
    Jv  = [Jv , cumsum(.5*(jvx(2:$) + jvx(1:$-1)) .* (VV(2:$) - VV(1:$-1)))];  // eqn 12
    // -- final formula for spreading current -> Ispread --
    Ispread  = C * sign(V-Voc) .* sqrt( (2/Rsq) * Jv );                        // eqn 7

endfunction
