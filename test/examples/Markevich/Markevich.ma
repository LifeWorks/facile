(* initial values *)
iM = 0;
iK = 0;
iMK = 0;
iMp = 0;
iMpK = 0;
iMpp = 0;
iP3 = 0;
iMppP3 = 0;
iMpP3 = 0;

(* constants *)
k1 = 1;
k3 = 1;
h1 = 1;
h_1 = 1;
my_h3 = 1;
h4 = 1;
my_h6 = 1;

(* expressions *)
k_1 = 3;
k2 = k_1;
k_3 = k1;
k4 = sin(k1^2.3*3)^k2;
h_4 = 1;

(* ode for independent species *)
dMdt = + k_1 MK + my_h6 MpP3 - k1 M K ;
dKdt = + k_1 MK + k2 MK + k_3 MpK + k4 MpK - k1 M K - k3 Mp K ;
dMKdt = + k1 M K - k_1 MK - k2 MK ;
dMpdt = + k2 MK + k_3 MpK + my_h3 MppP3 + h_4 MpP3 - k3 Mp K - h4 Mp P3 ;
dMpKdt = + k3 Mp K - k_3 MpK - k4 MpK ;
dMppdt = + k4 MpK + h_1 MppP3 - h1 Mpp P3 ;
dP3dt = + h_1 MppP3 + my_h3 MppP3 + h_4 MpP3 + my_h6 MpP3 - h1 Mpp P3 - h4 Mp P3 ;
dMppP3dt = + h1 Mpp P3 - h_1 MppP3 - my_h3 MppP3 ;
dMpP3dt = + h4 Mp P3 - h_4 MpP3 - my_h6 MpP3 ;
