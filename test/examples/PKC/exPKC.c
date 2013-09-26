
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
doublereal Ca_ER,PKC,rop,DAG,PIP2,IP3,C_GtPLC,k6_p,k6_s,k6_m,k7,v8,v10,kg,F,k12,k11,b11,Ca_tot,Ca,PKCt_tot,PKCt,k6;
  /* Evaluates the algebraic equations or ODE right hand side */

  /* Input arguments : */
  /*      ndim   :   Dimension of the ODE system */
  /*      u      :   State variables */
  /*      icp    :   Array indicating the free parameter(s) */
  /*      par    :   Equation parameters */

  /* Values to be returned : */
  /*      f      :   ODE right hand side values */

  /* Normally unused Jacobian arguments : IJAC, DFDU, DFDP (see manual) */

Ca_ER = u[0];
PKC = u[1];
rop = u[2];
DAG = u[3];
PIP2 = u[4];
IP3 = u[5];
C_GtPLC=par[0];
k6_p= 10-PIP2;
k6_s=par[1];
k6_m=par[2];
k7=(doublereal)2;
v8= (20*pow((IP3/(0.4+IP3)),4)+0.01)*Ca_ER*(1-rop);

/* v10 should be defined after ca below !!!!!!!!!!!!!!*/

v10= 40*pow(Ca,2)/(pow(0.15,2)+pow(Ca,2));
kg= pow(Ca,4)*(1-rop);
F=(doublereal)0.02;
k12=(doublereal)0.7;
k11=(doublereal)0.1;
b11=(doublereal)70;
Ca_tot=(doublereal)1100.1;
Ca=  - Ca_ER + Ca_tot;
PKCt_tot=(doublereal)1;
PKCt=  - PKC + PKCt_tot;
k6= k6_s*C_GtPLC*Ca/(k6_m+Ca*C_GtPLC);
f[0] =  v10 - v8 ;
f[1] =  b11*PKCt - k11*DAG*PKC ;
f[2] =  kg - F*rop ;
f[3] =  k6*PIP2 + b11*PKCt - k12*DAG - k11*DAG*PKC ;
f[4] =  k6_p - k6*PIP2 ;
f[5] =  k6*PIP2 - k7*IP3 ;

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

u[0] = (doublereal)1100;
u[1] = (doublereal)1;
u[2] = (doublereal)0;
u[3] = (doublereal)0;
u[4] = (doublereal)10;
u[5] = (doublereal)0;
par[0]=0.07;
par[1]=0.4;
par[2]=0.4;

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
