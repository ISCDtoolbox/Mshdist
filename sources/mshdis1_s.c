#include "mshdist.h"

extern char  ddb;
extern unsigned char inxt2[5];

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
  pTria         pt;
  pPoint        p0,p1,p2;
  double        dist,d0,d1,d2,ll,defval1,defval2,alpha,det,ialpha,idet,m[2][2],im[2][2],gre[2][3],Gr[3][3],g[3],a[3],r[2];
  double        ps1,ps2,rmin,rmax,g1g2,ng0,ng1,ng2;
  int           ip,ip1,ip2,nr;
  char          i1,i2;
  
  i1 = inxt2[i];
  i2 = inxt2[i1];
  
  pt = &mesh->tria[k];
  ip  = pt->v[i];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];
  p0  = &mesh->point[ip];
  p1  = &mesh->point[ip1];
  p2  = &mesh->point[ip2];
  
  d0  = sol->val[ip];
  d1  = sol->val[ip1];
  d2  = sol->val[ip2];
  
  dist = INIVAL_3d;
  
  /* If ip1 is not accepted, calculate a trial value based on ip2 */
  if ( p1->tag != 1 ) {
    ll = (p2->c[0]-p0->c[0])*(p2->c[0]-p0->c[0]) + (p2->c[1]-p0->c[1])*(p2->c[1]-p0->c[1]) + (p2->c[2]-p0->c[2])*(p2->c[2]-p0->c[2]);
    ll = sqrt(ll);
    dist = (d0 < 0.0) ? d2 -ll : d2 + ll;
    return( fabs(dist) );
  }
  /* Else if ip2 is not accepted, calculate a trial value based on ip1 */
  else if ( p2->tag != 1 ) {
    ll = (p1->c[0]-p0->c[0])*(p1->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])*(p1->c[1]-p0->c[1]) + (p1->c[2]-p0->c[2])*(p1->c[2]-p0->c[2]);
    ll = sqrt(ll);
    dist = (d0 < 0.0) ? d1 -ll : d1 + ll;
    return( fabs(dist) );
  }
  
  /* At this point, both values d1 and d2 are accepted; calculation of the gradient of basis functions */
  m[0][0]            = (p1->c[0]-p0->c[0])*(p1->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])*(p1->c[1]-p0->c[1]) + (p1->c[2]-p0->c[2])*(p1->c[2]-p0->c[2]);
  m[0][1] = m[1][0]  = (p1->c[0]-p0->c[0])*(p2->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])*(p2->c[1]-p0->c[1]) + (p1->c[2]-p0->c[2])*(p2->c[2]-p0->c[2]);
  m[1][1]            = (p2->c[0]-p0->c[0])*(p2->c[0]-p0->c[0]) + (p2->c[1]-p0->c[1])*(p2->c[1]-p0->c[1]) + (p2->c[2]-p0->c[2])*(p2->c[2]-p0->c[2]);
  
  /* Divide matrix by the mean value of the coefficients to avoid degenerate values */
  alpha = fabs(m[0][0]) + 2.0*fabs(m[0][1]) + fabs(m[1][1]);
  if ( alpha < EPS1 ) return(INIVAL_3d);
  ialpha = 4.0 / alpha;
  
  m[0][0] *= ialpha ; m[0][1] *= ialpha;
  m[1][0] *= ialpha ; m[1][1] *= ialpha;
  
  /* Inverse matrix */
  det = m[0][0]*m[1][1] - m[1][0]*m[0][1];
  if ( det < EPS1 ) return(INIVAL_3d);
  idet = 1.0 / det;
  
  im[0][0] = idet*m[1][1];     im[0][1] = -idet*m[0][1];
  im[1][0] = -idet*m[1][0];    im[1][1] = idet*m[0][0];
  
  /* Reduced gradients = coefficients in the basis (p1-p0,p2-p0) */
  gre[0][0] = ialpha*(-im[0][0] -im[0][1])  ;  gre[0][1] = ialpha*im[0][0] ; gre[0][2] = ialpha*im[0][1] ;
  gre[1][0] = -ialpha*(-im[1][0] -im[1][1]) ;  gre[1][1] = ialpha*im[1][0] ; gre[1][2] = ialpha*im[1][1] ;
  
  /* Gr[0,1,2][i] = coordinates 0,1,2 of ith basis function */
  Gr[0][0] = gre[0][0]*(p1->c[0]-p0->c[0]) + gre[1][0]*(p2->c[0]-p0->c[0]);
  Gr[1][0] = gre[0][0]*(p1->c[1]-p0->c[1]) + gre[1][0]*(p2->c[1]-p0->c[1]);
  Gr[2][0] = gre[0][0]*(p1->c[2]-p0->c[2]) + gre[1][0]*(p2->c[2]-p0->c[2]);
  
  Gr[0][1] = gre[0][1]*(p1->c[0]-p0->c[0]) + gre[1][1]*(p2->c[0]-p0->c[0]);
  Gr[1][1] = gre[0][1]*(p1->c[1]-p0->c[1]) + gre[1][1]*(p2->c[1]-p0->c[1]);
  Gr[2][1] = gre[0][1]*(p1->c[2]-p0->c[2]) + gre[1][1]*(p2->c[2]-p0->c[2]);
  
  Gr[0][2] = gre[0][2]*(p1->c[0]-p0->c[0]) + gre[1][2]*(p2->c[0]-p0->c[0]);
  Gr[1][2] = gre[0][2]*(p1->c[1]-p0->c[1]) + gre[1][2]*(p2->c[1]-p0->c[1]);
  Gr[2][2] = gre[0][2]*(p1->c[2]-p0->c[2]) + gre[1][2]*(p2->c[2]-p0->c[2]);
  
  /* Local solver for Eikonal equation */
  ng0  = Gr[0][0]*Gr[0][0] + Gr[1][0]*Gr[1][0] + Gr[2][0]*Gr[2][0];
  ng1  = Gr[0][1]*Gr[0][1] + Gr[1][1]*Gr[1][1] + Gr[2][1]*Gr[2][1];
  ng2  = Gr[0][2]*Gr[0][2] + Gr[1][2]*Gr[1][2] + Gr[2][2]*Gr[2][2];
  g1g2 = Gr[0][1]*Gr[0][2] + Gr[1][1]*Gr[1][2] + Gr[2][1]*Gr[2][2];
  
  a[2] = ng0;
  a[1] = -2.0*(d1*ng1 + d2*ng2 + (d1+d2)*g1g2);
  a[0] = d1*d1*ng1 + d2*d2*ng2 + 2.0*d1*d2*g1g2 - 1.0;
  
  nr = eqquad(a,r);
  
  /* Only one real root */
  if ( nr == 1 ) {
    if ( (d0 >= 0.0 && r[0] >= d1 && r[0] >= d2) || (d0 <= 0.0 && r[0] <= d1 && r[0] <= d2) ) {
      /* Check the direction of the gradient in this case */
      g[0] = (d1-r[0])*Gr[0][1] + (d2-r[0])*Gr[0][2];
      g[1] = (d1-r[0])*Gr[1][1] + (d2-r[0])*Gr[1][2];
      g[2] = (d1-r[0])*Gr[2][1] + (d2-r[0])*Gr[2][2];
      ps1  = g[0]*Gr[0][1] + g[1]*Gr[1][1] + g[2]*Gr[2][1];
      ps2  = g[0]*Gr[0][2] + g[1]*Gr[1][2] + g[2]*Gr[2][2];
      
      if ( ( d0 >= 0.0 && ps1 < EPS1 && ps2 < EPS1 ) || ( d0 <= 0.0 && ps1 >- EPS1 && ps2 >- EPS1 )) dist = fabs(r[0]);
    }
  }
  
  /* Two real roots */
  else if ( nr == 2 ) {
    /* Sort the roots */
    if ( d0 >= 0.0 ) {
      if ( r[0] < r[1] ) {
        rmin = r[0];
        rmax = r[1];
      }
      else {
        rmin = r[1];
        rmax = r[0];
      }
    }
    else {
      if ( r[0] < r[1] ) {
        rmin = r[1];
        rmax = r[0];
      }
      else {
        rmin = r[1];
        rmax = r[0];
      }
    }
    
    /* Try first rmin */
    if ( (d0 >= 0.0 && rmin >= d1 && rmin >= d2) || (d0 <= 0.0 && rmin <= d1 && rmin <= d2) ) {
      g[0] = (d1-rmin)*Gr[0][1] + (d2-rmin)*Gr[0][2];
      g[1] = (d1-rmin)*Gr[1][1] + (d2-rmin)*Gr[1][2];
      g[2] = (d1-rmin)*Gr[2][1] + (d2-rmin)*Gr[2][2];
      ps1  = g[0]*Gr[0][1] + g[1]*Gr[1][1] + g[2]*Gr[2][1];
      ps2  = g[0]*Gr[0][2] + g[1]*Gr[1][2] + g[2]*Gr[2][2];
      
      if ( ( d0 >= 0.0 && ps1 < EPS1 && ps2 < EPS1 ) || ( d0 <= 0.0 && ps1 >- EPS1 && ps2 >- EPS1 )) dist = fabs(rmin);
      
    }
    /* Then try rmax if no update from rmin */
    if ( (fabs( dist - INIVAL_2d ) < EPS2) && ((d0 >= 0.0 && rmax >= d1 && rmax >= d2) || (d0 <= 0.0 && rmax <= d1 && rmax <= d2)) ) {
      g[0] = (d1-rmax)*Gr[0][1] + (d2-rmax)*Gr[0][2];
      g[1] = (d1-rmax)*Gr[1][1] + (d2-rmax)*Gr[1][2];
      g[2] = (d1-rmax)*Gr[2][1] + (d2-rmax)*Gr[2][2];
      ps1  = g[0]*Gr[0][1] + g[1]*Gr[1][1] + g[2]*Gr[2][1];
      ps2  = g[0]*Gr[0][2] + g[1]*Gr[1][2] + g[2]*Gr[2][2];
      
      if ( ( d0 >= 0.0 && ps1 < EPS1 && ps2 < EPS1 ) || ( d0 <= 0.0 && ps1 >- EPS1 && ps2 >- EPS1 )) dist = fabs(rmax);
    }
    
    /* If no other value has been assigned to dist, calculate a trial value based on both triangle edges */
    if ( fabs( dist - INIVAL_2d ) < EPS2 ) {
      ll = (p1->c[0]-p0->c[0])*(p1->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])*(p1->c[1]-p0->c[1]) + (p1->c[2]-p0->c[2])*(p1->c[2]-p0->c[2]);
      ll = sqrt(ll);
      defval1 = (d0 < 0.0) ? d1 -ll : d1 + ll;
      
      ll = (p2->c[0]-p0->c[0])*(p2->c[0]-p0->c[0]) + (p2->c[1]-p0->c[1])*(p2->c[1]-p0->c[1]) + (p2->c[2]-p0->c[2])*(p2->c[2]-p0->c[2]);
      ll = sqrt(ll);
      defval2 = (d0 < 0.0) ? d2 -ll : d2 + ll;
      
      dist = D_MIN(fabs(defval1),fabs(defval2));
    }
  }
  
  return(dist);
}

/* Propagation of the signed distance function by the Fast Marching Method */
int ppgdistfmm_s(pMesh mesh,pSol sol) {

  
  return(1);
}