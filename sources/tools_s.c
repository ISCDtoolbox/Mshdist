#include "mshdist.h"

extern unsigned char inxt2[5];

/* Return (squared) distance from point pa to segment (p1,p2);
   proj = 2 if distance is realized by p1 or p2,
          1 if it is realized by the orthogonal projection of pa on (p1p2) */
double distpt_s(pPoint p0,pPoint p1,pPoint pa,int *proj) {
  double   ux,uy,uz,p1p0,pap0,ps,lambda;
  
  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];
  p1p0 = ux*ux + uy*uy + uz*uz;
  pap0 = (pa->c[0]-p0->c[0])*(pa->c[0]-p0->c[0]) + (pa->c[1]-p0->c[1])*(pa->c[1]-p0->c[1]) + (pa->c[2]-p0->c[2])*(pa->c[2]-p0->c[2]);

  /* If p0p1 is too short, return squared distance to p0 */
  if ( p1p0 < EPS1 )
    return ( pap0 );
  
  ps = (pa->c[0]-p0->c[0])*ux + (pa->c[1]-p0->c[1])*uy + (pa->c[2]-p0->c[2])*uz;
  lambda = ps / p1p0;
  
  /* Closest point is p0 */
  if ( lambda < 0.0 ) {
    *proj = 2;
    return (pap0);
  }
  /* Closest point is p1 */
  else if ( lambda > 1.0 ) {
    *proj = 2;
    return ( (pa->c[0]-p1->c[0])*(pa->c[0]-p1->c[0]) + (pa->c[1]-p1->c[1])*(pa->c[1]-p1->c[1]) + (pa->c[2]-p1->c[2])*(pa->c[2]-p1->c[2]) );
  }
  /* Closest point is q = p0 + lambda*(p1-p0) \in [p0,p1] */
  else {
    *proj = 1;
    return ( fabs(pap0 - lambda*lambda*p1p0) );
  }
}

/* Compute (squared) distance from pa to the 0 level set in triangle k */
double distnv0_s(pMesh mesh, pSol sol, int k, pPoint pa, int *proj) {
  pTria   pt;
  pPoint  p0,p1,p2;
  Point   q,r;
  double  aux,lambda,mu,v0,v1,v2;
  int     i0,i1,i2;
  
  pt = &mesh->tria[k];
  i0 = pt->v[0];
  i1 = pt->v[1];
  i2 = pt->v[2];
  p0 = &mesh->point[i0];
  p1 = &mesh->point[i1];
  p2 = &mesh->point[i2];
  v0 = sol->val[i0];
  v1 = sol->val[i1];
  v2 = sol->val[i2];
  
  /* Wrong configuration : flat level set function */
  if ( !((v0 != 0.0) || (v1 != 0.0) || (v2 != 0.0)) )
    printf("k = %d: %E %E %E\n",k,v0,v1,v2);
  assert ( (v0 != 0.0) || (v1 != 0.0) || (v2 != 0.0) );
  
  /* Case where 2 of the vertices are 0.*/
  if ( (v0 == 0.0) && (v1 == 0.0) ){
    return(distpt_2d(p0, p1, pa, proj));
  }
  
  if ( (v0 == 0.0) && (v2 == 0.0) ){
    return(distpt_2d(p0, p2, pa, proj));
  }
  
  if ( (v1 == 0.0) && (v2 == 0.0) ){
    return(distpt_2d(p1, p2, pa, proj));
  }
      
  /* Case where exactly 1 of the vertices is 0. */
  if ( v0 == 0.0 ){
    if ( v1*v2 < 0.0 ){
      lambda = v1 / (v1-v2);
      q.c[0] = p1->c[0] + lambda*(p2->c[0] - p1->c[0]);
      q.c[1] = p1->c[1] + lambda*(p2->c[1] - p1->c[1]);
      return (distpt_s(p0, &q, pa, proj));
    }
    else {
      *proj = 2;
      return((pa->c[0]-p0->c[0])*(pa->c[0]-p0->c[0]) + (pa->c[1]-p0->c[1])*(pa->c[1]-p0->c[1]));
    }
  }
  
  if ( v1 == 0.0 ){
    if ( v0*v2 < 0.0 ){
      lambda = v0/(v0-v2);
      q.c[0] = p0->c[0] + lambda*(p2->c[0] - p0->c[0]);
      q.c[1] = p0->c[1] + lambda*(p2->c[1] - p0->c[1]);
      return (distpt_s(p1, &q, pa, proj));
    }
    else {
      *proj = 2;
      return((pa->c[0]-p1->c[0])*(pa->c[0]-p1->c[0]) + (pa->c[1]-p1->c[1])*(pa->c[1]-p1->c[1]));
    }
  }
  
  if ( v2 == 0.0 ){
    if ( v0*v1 < 0.0 ){
      lambda = v0/(v0-v1);
      q.c[0] = p0->c[0] + lambda*(p1->c[0] - p0->c[0]);
      q.c[1] = p0->c[1] + lambda*(p1->c[1] - p0->c[1]);
      return (distpt_s(p2, &q, pa, proj));
    }
    else {
      *proj = 2;
      return((pa->c[0]-p2->c[0])*(pa->c[0]-p2->c[0]) + (pa->c[1]-p2->c[1])*(pa->c[1]-p2->c[1]));
    }
  }
      
  /* General case; make a permutation of indices so that p1 and p2 have same sign*/
  if ( v0*v1 >= 0.0 ){
    p0  = &mesh->point[i2];
    aux = v0;
    v0  = v2;
    
    p2 = &mesh->point[i0];
    v2 = aux;
  }
  else if ( v0*v2 >= 0.0 ){
    p0  = &mesh->point[i1];
    aux = v0;
    v0  = v1;
    
    p1 = &mesh->point[i0];
    v1 = aux;
  }
  
  assert( (v0 != v1) && (v0 != v2) );
  lambda = v0/(v0-v1);
  mu     = v0/(v0-v2);
  
  q.c[0] = p0->c[0] + lambda*(p1->c[0] - p0->c[0]);
  q.c[1] = p0->c[1] + lambda*(p1->c[1] - p0->c[1]);
  q.c[2] = p0->c[2] + lambda*(p1->c[2] - p0->c[2]);

  r.c[0] = p0->c[0] + mu*(p2->c[0] - p0->c[0]);
  r.c[1] = p0->c[1] + mu*(p2->c[1] - p0->c[1]);
  r.c[2] = p0->c[2] + mu*(p2->c[2] - p0->c[2]);

  return( distpt_s(&q, &r, pa, proj) );
}
