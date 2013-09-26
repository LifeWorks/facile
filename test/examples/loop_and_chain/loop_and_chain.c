
#include "auto_f2c.h"
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*   ab :            The A --> B reaction */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int func (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *f, doublereal *dfdu, doublereal *dfdp) {
doublereal X2,Y2,Z2,X1,Y1,Z1,C0,S1,C1,S2,C2,S3,f1,kf0,kb0,kp0,kf1,kb1,kp1,kf2,kb2,W1_tot,E1_tot,E2_tot,S0_tot,E0_tot,W2_tot,W1,E1,E2,S0,E0,W2,f2,f3,f4,kp2;
  /* Evaluates the algebraic equations or ODE right hand side */

  /* Input arguments : */
  /*      ndim   :   Dimension of the ODE system */
  /*      u      :   State variables */
  /*      icp    :   Array indicating the free parameter(s) */
  /*      par    :   Equation parameters */

  /* Values to be returned : */
  /*      f      :   ODE right hand side values */

  /* Normally unused Jacobian arguments : IJAC, DFDU, DFDP (see manual) */

/* free species */
X2 = u[0];
Y2 = u[1];
Z2 = u[2];
X1 = u[3];
Y1 = u[4];
Z1 = u[5];
C0 = u[6];
S1 = u[7];
C1 = u[8];
S2 = u[9];
C2 = u[10];
S3 = u[11];

/* constants */
f1=par[0];
kf0=par[1];
kb0=par[2];
kp0=par[3];
kf1=par[4];
kb1=par[5];
kp1=par[6];
kf2=par[7];
kb2=par[8];

/* moiety totals */
W1_tot=par[9];
E1_tot=par[10];
E2_tot=par[11];
S0_tot=par[12];
E0_tot=par[13];
W2_tot=par[14];

/* dependent species */
W1=  - Z1 - Y1 - X1 + W1_tot;
E1=  - C1 + E1_tot;
E2=  - C2 + E2_tot;
S0=  - S3 - C2 - S2 - C1 - S1 - C0 + S0_tot;
E0=  - C0 + E0_tot;
W2=  - Z2 - Y2 - X2 + W2_tot;

/* expressions */
f2= f1;
f3= f1;
f4= f1;
kp2= 12.0*1/pow(X1,2);

/* ode for independent species */
f[0] =  f1*W2 - f2*X2                   ;   /* d/dt (X2) */
f[1] =  f2*X2 - f3*Y2                   ;   /* d/dt (Y2) */
f[2] =  f3*Y2 - f4*Z2                   ;   /* d/dt (Z2) */
f[3] =  f1*W1 - f2*X1                   ;   /* d/dt (X1) */
f[4] =  f2*X1 - f3*Y1                   ;   /* d/dt (Y1) */
f[5] =  f3*Y1 - f4*Z1                   ;   /* d/dt (Z1) */
f[6] =  kf0*E0*S0 - kb0*C0 - kp0*C0     ;   /* d/dt (C0) */
f[7] =  kp0*C0 + kb1*C1 - kf1*E1*S1     ;   /* d/dt (S1) */
f[8] =  kf1*E1*S1 - kb1*C1 - kp1*C1     ;   /* d/dt (C1) */
f[9] =  kp1*C1 + kb2*C2 - kf2*E2*S2     ;   /* d/dt (S2) */
f[10] =  kf2*E2*S2 - kb2*C2 - kp2*C2    ;   /* d/dt (C2) */
f[11] =  kp2*C2                         ;   /* d/dt (S3) */

return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int stpnt (integer ndim, doublereal t,
           doublereal *u, doublereal *par) {
  /* Input arguments : */
  /*      ndim   :   Dimension of the ODE system */

  /* Values to be returned : */
  /*      u      :   A starting solution vector */
  /*      par    :   The corresponding equation-parameter values */


  /* Initialize the equation parameters */

/* free species */
u[0] = (doublereal) 0; /* X2 */
u[1] = (doublereal) 0; /* Y2 */
u[2] = (doublereal) 0; /* Z2 */
u[3] = (doublereal) 0; /* X1 */
u[4] = (doublereal) 0; /* Y1 */
u[5] = (doublereal) 0; /* Z1 */
u[6] = (doublereal) 0; /* C0 */
u[7] = (doublereal) 0; /* S1 */
u[8] = (doublereal) 0; /* C1 */
u[9] = (doublereal) 0; /* S2 */
u[10] = (doublereal) 0; /* C2 */
u[11] = (doublereal) 0; /* S3 */

/* constants */
par[0]=0.10;         /* f1 */
par[1]=1.0;          /* kf0 */
par[2]=0.1;          /* kb0 */
par[3]=10.0;         /* kp0 */
par[4]=1.1;          /* kf1 */
par[5]=0.11;         /* kb1 */
par[6]=11.0;         /* kp1 */
par[7]=1.2;          /* kf2 */
par[8]=0.12;         /* kb2 */

/* moiety totals */
par[9]=1;            /* W1_tot */
par[10]=0.11;        /* E1_tot */
par[11]=0.12;        /* E2_tot */
par[12]=10;          /* S0_tot */
par[13]=0.1;         /* E0_tot */
par[14]=2;           /* W2_tot */

 return 0;
}
/* The following subroutines are not used here, */
/* but they must be supplied as dummy routines */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int bcnd (integer ndim, const doublereal *par, const integer *icp,
          integer nbc, const doublereal *u0, const doublereal *u1, integer ijac,
          doublereal *fb, doublereal *dbc) {
  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int icnd (integer ndim, const doublereal *par, const integer *icp,
          integer nint, const doublereal *u, const doublereal *uold,
          const doublereal *udot, const doublereal *upold, integer ijac,
          doublereal *fi, doublereal *dint) {
    return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int fopt (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *fs, doublereal *dfdu, doublereal *dfdp) {
    return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int pvls (integer ndim, const doublereal *u,
          doublereal *par) {
    return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
