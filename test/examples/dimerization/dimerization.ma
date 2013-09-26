(* initial values *)
iA = 10;
iB = 1;

(* constants *)
f1 = 1.0;
b1 = 2.0;

(* ode for independent species *)
dAdt = + b1 B + b1 B - f1 A A - f1 A A ;
dBdt = + f1 A A - b1 B ;
