(* initial values *)
iE = 1e-06;
iS = 0.001;
iC = 0;
iP = 0;

(* constants *)
f = 1.5e6;
b = 1e3;
k = 1e3;

(* ode for independent species *)
dEdt = + b C + k C - f E S ;
dSdt = + b C - f E S ;
dCdt = + f E S - b C - k C ;
dPdt = + k C ;
