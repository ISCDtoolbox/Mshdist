#include "mshdist.h"

extern char  ddb;

/* Initialize the signed distance function from entities existing in mesh */
int iniencdomain_s(Info info,pMesh mesh, pSol sol){
  pTria       pt;
  pPoint      p0,p1;
  double      dd;
  int         k,ip0,ip,iel,ilist,list[LONMAX],*adja;
  char        i0,i;
  
  /* Starting point (for now) */
  ip0 = 1901;
  
  /* Large initial value */
  for (k=1; k<=sol->np; k++)
    sol->val[k] = INIVAL_3d;
  
  /* Travel the ball of the starting point */
  p0 = &mesh->point[ip0];
  pt = &mesh->tria[p0->s];
  for (i0=0; i0<3; i0++)
    if ( pt->v[i0] == ip0 ) break;
  
  ilist = boulet_2d(mesh,p0->s,i0,list);
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    pt = &mesh->tria[iel];

    for(i=0; i<3; i++) {
      ip = pt->v[i];
      p1 = &mesh->point[ip];
      dd = (p1->c[0]-p0->c[0])* (p1->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])* (p1->c[1]-p0->c[1]) + (p1->c[2]-p0->c[2])* (p1->c[2]-p0->c[2]);
      sol->val[ip] = D_MIN(sol->val[ip],dd);
    }
  }
  
  /* Take square root */
  for (k=1; k<=sol->np; k++)
    sol->val[k] = sqrt(sol->val[k]);
  
  return(1);
}

/* Calculate a (positive) active value at vertex i in triangle k based on the values in the other two vertices */
double actival_s(pMesh mesh,pSol sol,int k,int i) {
  double dist;
  
  
  return(dist);
}

/* Propagation of the signed distance function by the Fast Marching Method */
int ppgdistfmm_s(pMesh mesh,pSol sol) {
  
  return(1);
}