(* initial values *)
ib = 100;
ic = 0;
id = 0;
ie = 0;

(* constants *)
f1 = 1;

(* ode for independent species *)
dbdt = - f1 b - f1 b ;
dcdt = + f1 b ;
dddt = + f1 b ;
dedt = + f1 b ;
