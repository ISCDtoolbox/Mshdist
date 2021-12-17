#include "mshdist.h"

extern unsigned char inxt2[5];
extern char  ddb;

/* Calculate the roots of the second order polynomial stored in a, and put the result in r; return number of real roots (-1 if polynomial is null) */
int eqquad(double *a,double *r) {
  double b,c,d,D;
  
  /* Polynomial is at most of degree 1*/
  if ( fabs(a[2]) < EPS1 ) {
    if ( fabs(a[1]) < EPS1 ) {
      if ( fabs(a[0]) < EPS1 ) return(-1);
      else                    return(0);
    }
    else {
      r[0] = -a[0] / a[1];
      return(1);
    }
  }
  
  /* Normalize coefficients */
  b = a[1] / a[2];
  c = a[0] / a[2];
  D = b*b - 4.0*c;
  
  if ( D < -EPS1 )
    return(0);
  else if ( fabs(D) < EPS1 ) {
    r[0] = -0.5*b;
    return(1);
  }
  else {
    d    = sqrt(D);
    r[0] = 0.5*(-b-d);
    r[1] = 0.5*(-b+d);
    return(2);
  }
}

/* Check whether edge (pa,pb) intersects (p1,p2) */
int intersec_2d(pPoint p1,pPoint p2,pPoint pa,pPoint pb) {
  double  det1,det2,det3,det4,ux,uy,vx,vy,wx,wy,zx,zy;
  
  ux = p2->c[0] - p1->c[0];
  uy = p2->c[1] - p1->c[1];
  vx = pa->c[0] - p1->c[0];
  vy = pa->c[1] - p1->c[1];
  det1 = ux*vy - vx*uy;
  
  wx = pb->c[0] - p1->c[0];
  wy = pb->c[1] - p1->c[1];
  det2 = ux*wy - wx*uy;
  
  /* coplanarity */
  if ( fabs(det1) < EPS1 && fabs(det2) < EPS1 ) {
    if ( fabs(pa->c[0]-pb->c[0]) < EPS1 ) {
      if ( wy*vy <= 0.0 )  return(1);
      if ( (uy-vy)*(uy-wy) <= 0.0 )  return(1);
      if ( wy*(uy-wy) >= 0.0 )  return(1);
    }
    else if ( fabs(pa->c[1]-pb->c[1]) < EPS1 ) {
      if ( wx*vx <= 0.0 )  return(1);
      if ( (ux-vx)*(ux-wx) <= 0.0 )  return(1);
      if ( wx*(ux-wx) >= 0.0 )  return(1);
    }
    
    return(0);
  }
  if ( det1*det2 > 0.0 )  return(0);
  
  zx = pb->c[0] - pa->c[0];
  zy = pb->c[1] - pa->c[1];
  det3 = vx*zy - zx*vy;
  det4 = zx*(uy-vy) - (ux-vx)*zy;
  
  return( det3*det4 <= 0.0 );
}

/* Return (squared) distance from point pa to segment (p1,p2);
   proj = 2 if distance is realized by p1 or p2,
          1 if it is realized by the orthogonal projection of pa on (p1p2) */
double distpt_2d(pPoint p1,pPoint p2,pPoint pa,int *proj) {
  double   a,b,c,d,dd,ux,uy,vx,vy,wx,wy,xp,yp,lambda;
  
  *proj = 1;
  a = p1->c[1] - p2->c[1];
  b = p2->c[0] - p1->c[0];
  c = -b*p1->c[1] - a*p1->c[0];
  d = INIVAL_2d;
  
  dd = a*a + b*b;
  if ( dd < EPS1 ) {
    d = (pa->c[0]-p1->c[0])*(pa->c[0]-p1->c[0]) + (pa->c[1]-p1->c[1])*(pa->c[1]-p1->c[1]);
    *proj =2;
    return(d);
  }
  
  lambda = (pa->c[0] - p2->c[0])*(-b) + (pa->c[1] - p2->c[1])*a;	
  dd = 1.0 / dd;
  lambda *=dd;
  xp = p2->c[0] + lambda*(-b);
  yp = p2->c[1] + lambda*a;
  
  ux = xp - p1->c[0];
  uy = yp - p1->c[1];
  vx = xp - p2->c[0];
  vy = yp - p2->c[1];
  wx = p2->c[0] - p1->c[0];
  wy = p2->c[1] - p1->c[1];
  
  if ( fabs(b) < EPS1 ) {
    if ( uy*wy <= 0.0 ) {
      d = (pa->c[0]-p1->c[0])*(pa->c[0]-p1->c[0]) + (pa->c[1]-p1->c[1])*(pa->c[1]-p1->c[1]);
      *proj = 2;
    }
    else if ( vy*wy <= 0.0 )
      d = (pa->c[0]-xp)*(pa->c[0]-xp) + (pa->c[1]-yp)*(pa->c[1]-yp);
    else {
      d = (pa->c[0]-p2->c[0])*(pa->c[0]-p2->c[0]) + (pa->c[1]-p2->c[1])*(pa->c[1]-p2->c[1]);
      *proj = 2;
    }
  }
  else {
    if ( ux*wx <= 0.0 ) {
      d = (pa->c[0]-p1->c[0])*(pa->c[0]-p1->c[0]) + (pa->c[1]-p1->c[1])*(pa->c[1]-p1->c[1]);
      *proj = 2;
    }
    else if ( vx*wx <= 0.0 )
      d = (pa->c[0]-xp)*(pa->c[0]-xp) + (pa->c[1]-yp)*(pa->c[1]-yp);
    else {
      d = (pa->c[0]-p2->c[0])*(pa->c[0]-p2->c[0]) + (pa->c[1]-p2->c[1])*(pa->c[1]-p2->c[1]);
      *proj = 2;
    }
  }
  
  return(d);
}

/* compute (squared) distance from pa to the 0 level set in triangle ntria ; same use of proj as before */
double distnv0_2d(pMesh mesh, pSol sol, int ntria, pPoint pa, int *proj) {
  pTria   pt; 
  pPoint  p0,p1,p2;
  Point   q,r;
  double  aux,lambda,mu,v0,v1,v2;
  int     i0,i1,i2;
  
  pt = &mesh->tria[ntria];
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
  if ( !((v0 != 0.0)||(v1 != 0.0)||(v2 != 0.0)) )
    printf("ntria = %d: %E %E %E\n",ntria,v0,v1,v2);
  assert( (v0 != 0.0)||(v1 != 0.0)||(v2 != 0.0) );
	
  /* 2 of the vertices are 0.*/
  if ( (v0 == 0.0) && (v1 == 0.0) ){
    return(distpt_2d(p0, p1, pa, proj));
  }
	
  if ( (v0 == 0.0) && (v2 == 0.0) ){
    return(distpt_2d(p0, p2, pa, proj));
  }
	
  if ( (v1 == 0.0) && (v2 == 0.0) ){
    return(distpt_2d(p1, p2, pa, proj));
  }
	    
  /* Only 1 of the vertices is 0. */
  if ( v0 == 0.0 ){
    if ( v1*v2 < 0.0 ){
      lambda = v1 / (v1-v2);
      q.c[0] = p1->c[0] + lambda*(p2->c[0] - p1->c[0]);
      q.c[1] = p1->c[1] + lambda*(p2->c[1] - p1->c[1]);
      return (distpt_2d(p0, &q, pa, proj));
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
      return (distpt_2d(p1, &q, pa, proj));
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
      return (distpt_2d(p2, &q, pa, proj));
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
	
  r.c[0] = p0->c[0] + mu*(p2->c[0] - p0->c[0]);
  r.c[1] = p0->c[1] + mu*(p2->c[1] - p0->c[1]);
	
  return( distpt_2d(&q, &r, pa, proj) );
}

/* Compute Hausdorff distance between two set of edges already appearing in mesh */
double hausdorff(pMesh mesh1, pMesh mesh2){
  pEdge pa,pat;
  pPoint p0,p1,pmil,p2,p3;
  Point mil;
  double rho1,rho2,haus,d, d0,d1,dmil;
  int k,j,proj,nac1,nac2;  
  
  rho1 = 0.0;
  rho2 = 0.0;
  pmil = &mil;
  
  nac1 = 0;
  nac2 = 0;
  
  for(k=1; k<=mesh1->na;k++){
    pa = &mesh1->edge[k];
    pa->ref = 0;
  }
  
  for(k=1; k<=mesh2->na;k++){
    pa = &mesh2->edge[k];
    pa->ref = 0;
  }
  
  /* Set active edges for mesh 1 (discard bounding box) */
  for(k = 1; k<=mesh1->na;k++){
    pa = &mesh1->edge[k];
    p0 = &mesh1->point[pa->v[0]];
    p1 = &mesh1->point[pa->v[1]];
    
    if((p0->c[0]<0.01)||(p0->c[0]>0.99)|| (p0->c[1]<0.01)||(p0->c[1]>0.99))
      continue;
    pa->ref  =1;
    nac1++;
  }
  
  /* Set active edges for mesh 2 */
  for(k = 1; k<=mesh2->na;k++){
    pa = &mesh2->edge[k];
    p0 = &mesh2->point[pa->v[0]];
    p1 = &mesh2->point[pa->v[1]];
   
    if((p0->c[0]<0.01)||(p0->c[0]>0.99)|| (p0->c[1]<0.01)||(p0->c[1]>0.99))
      continue;
    
    pa->ref  =1;
    nac2++;
  }
  
  printf("Number of active edges %d %d \n", nac1, nac2);
  
  /* Compute rho(\Gamma_1,\Gamma_2)*/
  for(k=1;k<=mesh1->na;k++){
    pa = &mesh1->edge[k];
    
    if(!pa->ref) continue;
    
    p0 = &mesh1->point[pa->v[0]];
    p1 = &mesh1->point[pa->v[1]];
    pmil->c[0] = 0.5*(p0->c[0]+p1->c[0]); 
    pmil->c[1] = 0.5*(p0->c[1]+p1->c[1]);
    
    d0 = 10.0;
    d1 = 10.0;
    dmil = 10.0;
    
    for(j=1;j<=mesh2->na;j++){
      pat = &mesh2->edge[j];
      p2 = &mesh2->point[pat->v[0]];
      p3 = &mesh2->point[pat->v[1]];
      
      d = distpt_2d(p2,p3,p0,&proj);
      d0 = D_MIN(d0,d);
      
      d = distpt_2d(p2,p3,p1,&proj);
      d1 = D_MIN(d1,d);
      
      d = distpt_2d(p2,p3,pmil,&proj);
      dmil = D_MIN(dmil,d);	  
      
    }
    rho1 = D_MAX(rho1,d0);
    rho1 = D_MAX(rho1,d1);
    rho1 = D_MAX(rho1,dmil);
  }
  
  /* Compute rho(\Gamma_2,\Gamma_1) */
  for(k=1;k<=mesh2->na;k++){
    pa = &mesh2->edge[k];
    
    if(!pa->ref) continue;
    
    p0 = &mesh2->point[pa->v[0]];
    p1 = &mesh2->point[pa->v[1]];
    pmil->c[0] = 0.5*(p0->c[0]+p1->c[0]); 
    pmil->c[1] = 0.5*(p0->c[1]+p1->c[1]);
    
    d0 = 10.0;
    d1 = 10.0;
    dmil = 10.0;
    
    for(j=1;j<=mesh1->na;j++){
      pat = &mesh1->edge[j];
      
      if(!pat->ref) continue;
      
      p2 = &mesh1->point[pat->v[0]];
      p3 = &mesh1->point[pat->v[1]];
      
      d = distpt_2d(p2,p3,p0,&proj);
      d0 = D_MIN(d0,d);
      
      d = distpt_2d(p2,p3,p1,&proj);
      d1 = D_MIN(d1,d);
      
      d = distpt_2d(p2,p3,pmil,&proj);
      dmil = D_MIN(dmil,d);	  
    }
    rho2 = D_MAX(rho1,d0);
    rho2 = D_MAX(rho1,d1);
    rho2 = D_MAX(rho1,dmil);
  }
  
  haus = D_MAX(rho1,rho2);
  haus = sqrt(haus);
  
  return(haus);
}
