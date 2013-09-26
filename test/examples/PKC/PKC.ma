(* initial values *)
iPIP2 = 10;
iIP3 = 0;
iDAG = 0;
iCa_ER = 1100;
iCa = 0.1;
irop = 0;
iPKC = 1;
iPKCt = 0;

(* constants *)
C_GtPLC = 0.07;
k6_s = 0.4;
k6_m = 0.4;
k7 = 2;
F = 0.02;
k12 = 0.7;
k11 = 0.1;
b11 = 70;

(* expressions *)
k6_p = 10-PIP2;
k6 = k6_s*C_GtPLC*Ca/(k6_m+Ca*C_GtPLC);
v8 = (20*(IP3/(0.4+IP3))^4+0.01)*Ca_ER*(1-rop);
v10 = 40*Ca^2/(0.15^2+Ca^2);
kg = Ca^4*(1-rop);

(* ode for independent species *)
dPIP2dt = + k6_p - k6 PIP2 ;
dIP3dt = + k6 PIP2 - k7 IP3 ;
dDAGdt = + k6 PIP2 + b11 PKCt - k12 DAG - k11 DAG PKC ;
dCa_ERdt = + v10 - v8 ;
dCadt = + v8 - v10 ;
dropdt = + kg - F rop ;
dPKCdt = + b11 PKCt - k11 DAG PKC ;
dPKCtdt = + k11 DAG PKC - b11 PKCt ;
