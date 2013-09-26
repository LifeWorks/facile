(* initial values *)
ig = 2;
ia = 0;
ih = 2;
ib = 2;
ic = 0;
id = 0;
ik = 0;
if = 0;
ij = 0;
ie = 0;
ii = 0;

(* constants *)
v1 = 1;
v2 = 1;
v3 = 1;
v4 = 1;
v7 = 1;
v5 = 1;
v6 = 1;
v8 = 1;
v9 = 1;

(* ode for independent species *)
dgdt = + v6 + v8 + v9 - v1 - v2 ;
dadt = + v1 - v2 ;
dhdt = + v1 + v2 - v6 - v8 - v9 - v9 ;
dbdt = + v2 - v3 ;
dcdt = + v3 - v4 - v7 ;
dddt = + v3 + v4 - v5 ;
dkdt = + v5 - v7 ;
dfdt = + v7 - v8 ;
djdt = + v7 - v5 ;
dedt = + v5 - v6 ;
didt = + v9 ;
