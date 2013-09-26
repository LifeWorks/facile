(* initial values *)
iW1 = 1.0;
iX1 = 0;
iY1 = 0;
iZ1 = 0;
iW2 = 2.0;
iX2 = 0;
iY2 = 0;
iZ2 = 0;
iE0 = 0.10;
iS0 = 10.0;
iC0 = 0;
iS1 = 0;
iE1 = 0.11;
iC1 = 0;
iS2 = 0;
iE2 = 0.12;
iC2 = 0;
iS3 = 0;

(* constants *)
f1 = 0.10;
kf0 = 1.0;
kb0 = 0.1;
kp0 = 10.0;
kf1 = 1.1;
kb1 = 0.11;
kp1 = 11.0;
kf2 = 1.2;
kb2 = 0.12;

(* expressions *)
f2 = f1;
f3 = f1;
f4 = f1;
kp2 = 12.0*1/X1^2;

(* ode for independent species *)
dW1dt = + f4 Z1 - f1 W1 ;
dX1dt = + f1 W1 - f2 X1 ;
dY1dt = + f2 X1 - f3 Y1 ;
dZ1dt = + f3 Y1 - f4 Z1 ;
dW2dt = + f4 Z2 - f1 W2 ;
dX2dt = + f1 W2 - f2 X2 ;
dY2dt = + f2 X2 - f3 Y2 ;
dZ2dt = + f3 Y2 - f4 Z2 ;
dE0dt = + kb0 C0 + kp0 C0 - kf0 E0 S0 ;
dS0dt = + kb0 C0 - kf0 E0 S0 ;
dC0dt = + kf0 E0 S0 - kb0 C0 - kp0 C0 ;
dS1dt = + kp0 C0 + kb1 C1 - kf1 E1 S1 ;
dE1dt = + kb1 C1 + kp1 C1 - kf1 E1 S1 ;
dC1dt = + kf1 E1 S1 - kb1 C1 - kp1 C1 ;
dS2dt = + kp1 C1 + kb2 C2 - kf2 E2 S2 ;
dE2dt = + kb2 C2 + kp2 C2 - kf2 E2 S2 ;
dC2dt = + kf2 E2 S2 - kb2 C2 - kp2 C2 ;
dS3dt = + kp2 C2 ;
