
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
doublereal Ca,PKCt,PIP2,IP3,DAG,rop,C_GtPLC,k6_s,k6_m,k7,F,k12,k11,b11,PKC_tot,Ca_ER_tot,PKC,Ca_ER,k6_p,k6,v8,v10,kg;
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
Ca = u[0];
PKCt = u[1];
PIP2 = u[2];
IP3 = u[3];
DAG = u[4];
rop = u[5];

/* constants */
C_GtPLC=par[0];
k6_s=par[1];
k6_m=par[2];
k7=(doublereal)2;
F=(doublereal)0.02;
k12=(doublereal)0.7;
k11=(doublereal)0.1;
b11=(doublereal)70;

/* moiety totals */
PKC_tot= 1;
Ca_ER_tot= 1100.1;

/* dependent species */
PKC=  - PKCt + PKC_tot;
Ca_ER=  - Ca + Ca_ER_tot;

/* expressions */
k6_p= 10-PIP2;
k6= k6_s*C_GtPLC*Ca/(k6_m+Ca*C_GtPLC);
v8= (20*pow((IP3/(0.4+IP3)),4)+0.01)*Ca_ER*(1-rop);
v10= 40*pow(Ca,2)/(pow(0.15,2)+pow(Ca,2));
kg= pow(Ca,4)*(1-rop);

/* ode for independent species */
f[0] =  v8 - v10                        ;   /* d/dt (Ca) */
f[1] =  k11*DAG*PKC - b11*PKCt          ;   /* d/dt (PKCt) */
f[2] =  k6_p - k6*PIP2                  ;   /* d/dt (PIP2) */
f[3] =  k6*PIP2 - k7*IP3                ;   /* d/dt (IP3) */
f[4] =  k6*PIP2 + b11*PKCt - k12*DAG - k11*DAG*PKC ;   /* d/dt (DAG) */
f[5] =  kg - F*rop                      ;   /* d/dt (rop) */

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
u[0] = (doublereal) 0.1; /* Ca */
u[1] = (doublereal) 0; /* PKCt */
u[2] = (doublereal) 10; /* PIP2 */
u[3] = (doublereal) 0; /* IP3 */
u[4] = (doublereal) 0; /* DAG */
u[5] = (doublereal) 0; /* rop */

/* constants */
par[0]=0.07;         /* C_GtPLC */
par[1]=0.4;          /* k6_s */
par[2]=0.4;          /* k6_m */

/* moiety totals */

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
