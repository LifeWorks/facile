
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
doublereal g,a,h,b,c,d,k,f,v1,v2,v3,v4,v7,v5,v6,v8,v9,i_tot,e_tot,j_tot,i,e,j;
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
g = u[0];
a = u[1];
h = u[2];
b = u[3];
c = u[4];
d = u[5];
k = u[6];
f = u[7];

/* constants */
v1=par[0];
v2=par[1];
v3=par[2];
v4=par[3];
v7=par[4];
v5=par[5];
v6=par[6];
v8=par[7];
v9=par[8];

/* moiety totals */
i_tot=par[9];
e_tot=par[10];
j_tot=par[11];

/* dependent species */
i=  - h - g + i_tot;
e=  - f - d - c - (2)*b - h - a - (2)*g + e_tot;
j=  - k + j_tot;

/* expressions */

/* ode for independent species */
f[0] =  v6 + v8 + v9 - v1 - v2          ;   /* d/dt (g) */
f[1] =  v1 - v2                         ;   /* d/dt (a) */
f[2] =  v1 + v2 - v6 - v8 - v9 - v9     ;   /* d/dt (h) */
f[3] =  v2 - v3                         ;   /* d/dt (b) */
f[4] =  v3 - v4 - v7                    ;   /* d/dt (c) */
f[5] =  v3 + v4 - v5                    ;   /* d/dt (d) */
f[6] =  v5 - v7                         ;   /* d/dt (k) */
f[7] =  v7 - v8                         ;   /* d/dt (f) */

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
u[0] = (doublereal) 2; /* g */
u[1] = (doublereal) 0; /* a */
u[2] = (doublereal) 2; /* h */
u[3] = (doublereal) 2; /* b */
u[4] = (doublereal) 0; /* c */
u[5] = (doublereal) 0; /* d */
u[6] = (doublereal) 0; /* k */
u[7] = (doublereal) 0; /* f */

/* constants */
par[0]=1;            /* v1 */
par[1]=1;            /* v2 */
par[2]=1;            /* v3 */
par[3]=1;            /* v4 */
par[4]=1;            /* v7 */
par[5]=1;            /* v5 */
par[6]=1;            /* v6 */
par[7]=1;            /* v8 */
par[8]=1;            /* v9 */

/* moiety totals */
par[9]=4;            /* i_tot */
par[10]=10;          /* e_tot */
par[11]=0;           /* j_tot */

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
