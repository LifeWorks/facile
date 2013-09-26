(* initial values *)
iS0 = 0;
iS1 = 0;
iA = 0;
iB = 0;
iC = 0;
iD = 0;
iX = 0;
iY = 0.44e+1;
iW = 0;
iZ = 5.5;
iE = 10.2e+1;
iF = 0;
iG = 0;
iM1 = 0;
iM2 = 0;

(* constants *)
sourceS0f = 1;
fa0 = 10;
fb0 = 8;
b1 = 0.1;
f1 = 5;
f2 = 60e-1;
b2 = 0.02e+1;
sinkCf = 10;
sourceXb = 99;
sourceXf = 1;
sinkYb = 11;
sinkYf = 9;
sinkZb = 1.3;
sinkZf = 1;
gg = 0;
hh = 10;
ii = 5;
m1 = 0;
m2 = 0;
v2 = 55.5;

(* expressions *)
func1 = 0.5*(square(2*pi*t/0.2) + 1);
func2 = 0.5*(square(2*pi*t/0.2) + 1);
sourceWb = cos(t);
sourceWf = sin(2*pi*t/0.1);
ff = fa0;
v1 = E + G;
v3 = E + E;

(* ode for independent species *)
dS0dt = + sourceS0f - fa0 S0 - fa0 S0 ;
dS1dt = + func1 - fb0 S1 ;
dAdt = + fa0 S0 + fa0 S0 + b2 D - f1 A B - f2 A B ;
dBdt = + fb0 S1 + b2 D - f1 A B - f2 A B ;
dCdt = + f1 A B - sinkCf C ;
dDdt = + f2 A B + func2 - b2 D ;
dXdt = + sourceXf ;
dYdt = - sinkYf Y ;
dWdt = + sourceWf ;
dZdt = - sinkZf Z ;
dEdt = - ff E - gg E ;
dFdt = + ff E + ii G - hh F ;
dGdt = + gg E + hh F - ii G ;
dM1dt = + m2 M2 - m1 M1 ;
dM2dt = + m1 M1 - m2 M2 ;
