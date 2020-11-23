#include "mshdist.h"

extern char  ddb;

/* Initialize the signed distance function from entities existing in mesh */
int iniencdomain_s(Info info,pMesh mesh, pSol sol){
  int k,np0;
  
  /* Starting point (for now) */
  np0 = 1901;
  
  /* Large initial value */
  for (k=1; k<=sol->np; k++)
    sol->val[k] = INIVAL_3d;
  
  sol->val[np0] = 0.0;
  
  return(1);
}

/* Propagation of the signed distance function by the Fast Marching Method */
int ppgdistfmm_s(pMesh mesh,pSol sol) {
  
  return(1);
}