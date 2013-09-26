
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
doublereal MK,Mp,MpK,Mpp,MppP3,MpP3,k1,k3,h1,h_1,my_h3,h4,my_h6,K_tot,M_tot,P3_tot,K,M,P3,k_1,k2,k_3,k4,h_4;
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
MK = u[0];
Mp = u[1];
MpK = u[2];
Mpp = u[3];
MppP3 = u[4];
MpP3 = u[5];

/* constants */
k1=(doublereal)1;
k3=(doublereal)1;
h1=(doublereal)1;
h_1=(doublereal)1;
my_h3=(doublereal)1;
h4=(doublereal)1;
my_h6=(doublereal)1;

/* moiety totals */
K_tot=par[0];
M_tot= 0;
P3_tot= 0;

/* dependent species */
K=  - MK - MpK + K_tot;
M=  - MK - MpP3 - MppP3 - Mpp - MpK - Mp + M_tot;
P3=  - MpP3 - MppP3 + P3_tot;

/* expressions */
k_1= 3;
k2= k_1;
k_3= k1;
k4= pow(sin(pow(k1,2.3)*3),k2);
h_4= 1;

/* ode for independent species */
f[0] =  k1*M*K - k_1*MK - k2*MK         ;   /* d/dt (MK) */
f[1] =  k2*MK + k_3*MpK + my_h3*MppP3 + h_4*MpP3 - k3*Mp*K - h4*Mp*P3 ;   /* d/dt (Mp) */
f[2] =  k3*Mp*K - k_3*MpK - k4*MpK      ;   /* d/dt (MpK) */
f[3] =  k4*MpK + h_1*MppP3 - h1*Mpp*P3  ;   /* d/dt (Mpp) */
f[4] =  h1*Mpp*P3 - h_1*MppP3 - my_h3*MppP3 ;   /* d/dt (MppP3) */
f[5] =  h4*Mp*P3 - h_4*MpP3 - my_h6*MpP3 ;   /* d/dt (MpP3) */

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
u[0] = (doublereal) 0; /* MK */
u[1] = (doublereal) 0; /* Mp */
u[2] = (doublereal) 0; /* MpK */
u[3] = (doublereal) 0; /* Mpp */
u[4] = (doublereal) 0; /* MppP3 */
u[5] = (doublereal) 0; /* MpP3 */

/* constants */

/* moiety totals */
par[0]=0;            /* K_tot */

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
