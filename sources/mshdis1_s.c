#include "mshdist.h"

extern char  ddb;
extern unsigned char inxt2[5];

/* Find triangles intersecting the boundary and initialize distance at their vertices */
int iniredist_s(Info info, pMesh mesh, pSol sol){
  pTria    pt;
  pPoint   p0,p1,p2;
  double  *solTmp,d;
  int     *bndy,i,j,nb,nc,i0,i1,i2,proj;
  
  nb   = 0;
  bndy = (int*)calloc(mesh->nt+1,sizeof(int));
  assert(bndy);

  /* Store intersecting background triangles in list bndy */
  for (i=1; i<=mesh->nt; i++) {
    pt = &mesh->tria[i];
    i0 = pt->v[0];
    i1 = pt->v[1];
    i2 = pt->v[2];

    if ( (sol->val[i0] * sol->val[i1] <= 0.) || (sol->val[i0] * sol->val[i2] <= 0.) || (sol->val[i1] * sol->val[i2] <= 0.)) {
      nb++;
      bndy[nb] = i;
    }
  }
  
  bndy = (int*)realloc(bndy,(nb+1)*sizeof(int));
  printf("nb= %d\n",nb);
  
  /* Temporary values stored in solTmp */
  solTmp = (double*)calloc(mesh->np+1,sizeof(double));
  assert(solTmp);

  /* Travel list bndy and compute distance at vertices of its elements */
  for (i=1; i<=nb; i++) {
    pt = &mesh->tria[bndy[i]];
    i0 = pt->v[0];
    i1 = pt->v[1];
    i2 = pt->v[2];

    p0 = &mesh->point[i0];
    p1 = &mesh->point[i1];
    p2 = &mesh->point[i2];
    
    /* Distance from i0 to the 0 level set from triangle bndy[i] */
    d = distnv0_s(mesh,sol,bndy[i],p0,&proj);
    
    if ( p0->tag == 0 ) {
      solTmp[i0] = d;
      p0->tag    = proj;
    }
    else if ( d < fabs(solTmp[i0]) ) {
      solTmp[i0] = d;
      p0->tag    = proj;
    }
          
    /* Distance from i1 to the 0 level set from triangle bndy[i] */
    d = distnv0_s(mesh,sol,bndy[i],p1,&proj);
    
    if ( p1->tag == 0 ) {
      solTmp[i1] = d;
      p1->tag    = proj;
    }
    else if ( d < fabs(solTmp[i1]) ) {
      solTmp[i1] = d;
      p1->tag    = proj;
    }
        
    /* Distance from i2 to the 0 level set from triangle bndy[i] */
    d = distnv0_s(mesh,sol,bndy[i],p2,&proj);
    
    if ( p2->tag == 0 ) {
      solTmp[i2] = d;
      p2->tag    = proj;
    }
    else if ( d < fabs(solTmp[i2]) ) {
      solTmp[i2] = d;
      p2->tag    = proj;
    }
  }
  
  /* Correction procedure, for points with tag 2 */
  nc = 0;
  for (i=1; i<=mesh->np; i++) {
    p0 = &mesh->point[i];
    if ( p0->tag < 2 )  continue;
    for (j=1; j<=nb; j++) {
      pt = &mesh->tria[bndy[j]];
      d  = distnv0_s(mesh,sol,bndy[j],p0,&proj);
      if ( proj == 1 && d < fabs(solTmp[i]) ) {
        solTmp[i] = d;
        break;
      }
    }
    p0->tag = 1;
    nc++;
  }
    
  if ( nc )   fprintf(stdout,"     %d correction(s)\n",nc);
 
  /* Tke square root */
  for (i=1; i<=mesh->np; i++){
    p0 = &mesh->point[i];
    if ( !p0->tag )
      sol->val[i] = sol->val[i] < 0.0 ? -sqrt(INIVAL_3d) : sqrt(INIVAL_3d);
    else if ( p0->tag )
      sol->val[i] = sol->val[i] < 0.0 ? -sqrt(solTmp[i]) : sqrt(solTmp[i]);
  }

  free(solTmp);
  free(bndy);
  return(1);
}

/* Initialize the signed distance function from entities existing in mesh */
int iniencdomain_s(Info info,pMesh mesh, pSol sol){
  pTria       pt,pt1;
  pPoint      p0,p1,p2,pa,pb;
  double      dd;
  int         k,l,ip0,ip,ip1,ip2,iel,jel,nb,nc,ilist,list[LONMAX],*adja,*actied;
  char        i0,i,j,j0,j1,j2,ia,ib,proj;
  
  /* reset point tags and flags */
  for(k=1; k<=mesh->np; k++){
    p0 = &mesh->point[k];
    p0->flag = 0;
    p0->tag = 0;
  }
  
  /* Large initial value */
  for (k=1; k<=sol->np; k++)
    sol->val[k] = INIVAL_3d;
  
  /* Initial distance to user-specified vertices */
  for (k=0; k<info.nsp; k++) {
    ip0 = info.sp[k];
    p0 = &mesh->point[ip0];
    pt = &mesh->tria[p0->s];
    for (i0=0; i0<3; i0++)
      if ( pt->v[i0] == ip0 ) break;
    
    ilist = boulet_2d(mesh,p0->s,i0,list);
    for (l=0; l<ilist; l++) {
      iel = list[l] / 3;
      pt = &mesh->tria[iel];
      
      for(i=0; i<3; i++) {
        ip = pt->v[i];
        p1 = &mesh->point[ip];
        dd = (p1->c[0]-p0->c[0])* (p1->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])* (p1->c[1]-p0->c[1]) + (p1->c[2]-p0->c[2])* (p1->c[2]-p0->c[2]);
        sol->val[ip] = D_MIN(sol->val[ip],dd);
        p1->tag = 1;
      }
    }
  }
  
  /* Initial distance to user-specified edges */
  for (k=0; k<info.nsa; k++) {
    
  }
  
  /* Initial distance to inner sub-phase */
  nb = 0;
  nc = 0;
  
  /* Compute nb = number of boundary edges to be considered */
  for(k=1; k<= mesh->nt; k++){
    pt = &mesh->tria[k];
    if ( !isIntDom(info,pt->ref) ) continue;
    adja = &mesh->adja[3*(k-1)+1];
    
    for(i=0; i<3; i++){
      iel = adja[i] / 3;
      /* Exclude outer boundary triangles from starting boundary */
      if( !iel ) continue;
      if ( !isIntDom(info,mesh->tria[iel].ref) ) nb++;
    }
  }
  
  if ( nb ) {
    /* Declare table of active edges */
    actied = (int*)calloc(nb+1,sizeof(int));
    assert(actied);
    
    /* Store boundary edges under the form 3*iel+ia */
    for(k=1; k<= mesh->nt; k++){
      pt = &mesh->tria[k];
      if ( !isIntDom(info,pt->ref) ) continue;
      adja = &mesh->adja[3*(k-1)+1];
      
      for(i=0; i<3; i++){
        iel = adja[i] / 3;
        if( !iel ) continue;
        if ( !isIntDom(info,mesh->tria[iel].ref) ) {
          actied[nc] = 3*k+i;
          nc++;
        }
      }
    }
    
    /* Travel boundary edges : p->flag = last edge with respect
     to which distance has been evaluated */
    for(k=0; k<nb; k++){
      iel = actied[k] / 3;
      i0 = actied[k] % 3;
      
      pt = &mesh->tria[iel];
      
      ia = inxt2[i0];
      ib = inxt2[ia];
      
      pa = &mesh->point[pt->v[ia]];
      pb = &mesh->point[pt->v[ib]];
      
      /* Travel the ball of the two vertices of the edge */
      for(j=0; j<2; j++){
        i0 =inxt2[i0];
        
        ilist = boulet_2d(mesh,iel,i0,list);
        assert(ilist);
        
        for(l=0; l<ilist; l++){
          jel = list[l] / 3;
          j0  = list[l] % 3;
          pt1 = &mesh->tria[jel];
          
          j1 = inxt2[j0];
          j2 = inxt2[j1];
          
          ip1 = pt1->v[j1];
          ip2 = pt1->v[j2];
          
          p1 = &mesh->point[ip1];
          p2 = &mesh->point[ip2];
          
          if ( p1->flag != k ) {
            dd = distpt_s(pa,pb,p1,&proj);
            if ( dd < sol->val[ip1] ) {
              sol->val[ip1] = dd;
              p1->tag = proj;
            }
            p1 ->flag = k;
          }
          
          if ( p2->flag != k ) {
            dd = distpt_s(pa,pb,p2,&proj);
            if ( dd < sol->val[ip2] ) {
              sol->val[ip2] = dd;
              p2->tag = proj;
            }
            p2->flag = k;
          }
        }
      }
    }
    
    /* Set distance to all boundary points to 0 */
    for (k=0; k<nb; k++) {
      iel = actied[k] / 3;
      i0 = actied[k] % 3;
      pt = &mesh->tria[iel];
      
      for (j=0; j<2; j++) {
        i0 = inxt2[i0];
        ip = pt->v[i0];
        sol->val[ip] = 0.0;
      }
    }
    
    /* Correction procedure */
    nc = 0;
    
    if ( !info.fini ) {
      for (k=1; k<=mesh->np; k++) {
        p0 = &mesh->point[k];
        if ( p0->tag < 2 ) continue;
        
        for (l=0; l<nb; l++) {
          iel = actied[l] / 3;
          i0 = actied[l] % 3;
          pt = &mesh->tria[iel];
          
          ia = inxt2[i0];
          ib = inxt2[ia];
          
          pa = &mesh->point[pt->v[ia]];
          pb = &mesh->point[pt->v[ib]];
          
          dd = distpt_s(pa,pb,p0,&proj);
          if ( dd < sol->val[k] ) {
            nc++;
            sol->val[k] = dd;
            if ( proj == 1 ) {
              p0->tag = 1;
              break;
            }
          }
        }
      }
    }
    if ( nc )   fprintf(stdout,"     %d correction(s)\n",nc);

  }
  
  /* Reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;
  
  /* Put sign in the function, and take sqrt */
  for (k=1; k<=mesh->nt;k++) {
    pt = &mesh->tria[k];
    for (i=0; i<3; i++) {
      ip0 = pt->v[i];
      p0 = &mesh->point[ip0];
      if ( p0->flag == 1 ) continue;
      p0->flag = 1;
      if ( isIntDom(info,pt->ref) )
        sol->val[ip0] = -sqrt(sol->val[ip0]);
      else
        sol->val[ip0] = sqrt(sol->val[ip0]);
    }
  }
  
  /* Reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;
  
  return(1);
}

/* Calculate an active value at point a from 2 accepted (positive) values v_1, v_2 at points o1, o2 */
double actival2pt_s(double o1[3], double o2[3], double a[3], double d1, double d2) {
  double       dist,nga,ng1,ng2,g1g2,alpha,ialpha,det,idet,ll,ps1,ps2,defval1,defval2;
  double       rmin,rmax,m[2][2],im[2][2],gre[2][3],Gr[3][3],g[3],aa[3],r[2];
  int          nr;
  
  dist = INIVAL_3d;
    
  /* Calculation of the gradient of basis functions */
  m[0][0]            = (o1[0]-a[0])*(o1[0]-a[0]) + (o1[1]-a[1])*(o1[1]-a[1]) + (o1[2]-a[2])*(o1[2]-a[2]);
  m[0][1] = m[1][0]  = (o1[0]-a[0])*(o2[0]-a[0]) + (o1[1]-a[1])*(o2[1]-a[1]) + (o1[2]-a[2])*(o2[2]-a[2]);
  m[1][1]            = (o2[0]-a[0])*(o2[0]-a[0]) + (o2[1]-a[1])*(o2[1]-a[1]) + (o2[2]-a[2])*(o2[2]-a[2]);
  
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
  gre[1][0] = ialpha*(-im[1][0] -im[1][1]) ;  gre[1][1] = ialpha*im[1][0] ; gre[1][2] = ialpha*im[1][1] ;
  
  /* Gr[0,1,2][i] = coordinates 0,1,2 of ith basis function */
  Gr[0][0] = gre[0][0]*(o1[0]-a[0]) + gre[1][0]*(o2[0]-a[0]);
  Gr[1][0] = gre[0][0]*(o1[1]-a[1]) + gre[1][0]*(o2[1]-a[1]);
  Gr[2][0] = gre[0][0]*(o1[2]-a[2]) + gre[1][0]*(o2[2]-a[2]);
  
  Gr[0][1] = gre[0][1]*(o1[0]-a[0]) + gre[1][1]*(o2[0]-a[0]);
  Gr[1][1] = gre[0][1]*(o1[1]-a[1]) + gre[1][1]*(o2[1]-a[1]);
  Gr[2][1] = gre[0][1]*(o1[2]-a[2]) + gre[1][1]*(o2[2]-a[2]);
  
  Gr[0][2] = gre[0][2]*(o1[0]-a[0]) + gre[1][2]*(o2[0]-a[0]);
  Gr[1][2] = gre[0][2]*(o1[1]-a[1]) + gre[1][2]*(o2[1]-a[1]);
  Gr[2][2] = gre[0][2]*(o1[2]-a[2]) + gre[1][2]*(o2[2]-a[2]);
  
  /* Local solver for Eikonal equation */
  nga  = Gr[0][0]*Gr[0][0] + Gr[1][0]*Gr[1][0] + Gr[2][0]*Gr[2][0];
  ng1  = Gr[0][1]*Gr[0][1] + Gr[1][1]*Gr[1][1] + Gr[2][1]*Gr[2][1];
  ng2  = Gr[0][2]*Gr[0][2] + Gr[1][2]*Gr[1][2] + Gr[2][2]*Gr[2][2];
  g1g2 = Gr[0][1]*Gr[0][2] + Gr[1][1]*Gr[1][2] + Gr[2][1]*Gr[2][2];
  
  aa[2] = nga;
  aa[1] = -2.0*(d1*ng1 + d2*ng2 + (d1+d2)*g1g2);
  aa[0] = d1*d1*ng1 + d2*d2*ng2 + 2.0*d1*d2*g1g2 - 1.0;
  
  nr = eqquad(aa,r);
  
  /* Only one real root */
  if ( nr == 1 ) {
    if ( r[0] >= d1 && r[0] >= d2 ) {
      /* Check the direction of the gradient in this case */
      g[0] = (d1-r[0])*Gr[0][1] + (d2-r[0])*Gr[0][2];
      g[1] = (d1-r[0])*Gr[1][1] + (d2-r[0])*Gr[1][2];
      g[2] = (d1-r[0])*Gr[2][1] + (d2-r[0])*Gr[2][2];
      ps1  = g[0]*Gr[0][1] + g[1]*Gr[1][1] + g[2]*Gr[2][1];
      ps2  = g[0]*Gr[0][2] + g[1]*Gr[1][2] + g[2]*Gr[2][2];
      
      if ( ps1 < EPS1 && ps2 < EPS1 ) dist = fabs(r[0]);
    }
  }
  
  /* Two real roots */
  else if ( nr == 2 ) {
    /* Sort the roots */
    if ( r[0] < r[1] ) {
      rmin = r[0];
      rmax = r[1];
    }
    else {
      rmin = r[1];
      rmax = r[0];
    }
        
    /* Try first rmin */
    if ( rmin >= d1 && rmin >= d2 ) {
      g[0] = (d1-rmin)*Gr[0][1] + (d2-rmin)*Gr[0][2];
      g[1] = (d1-rmin)*Gr[1][1] + (d2-rmin)*Gr[1][2];
      g[2] = (d1-rmin)*Gr[2][1] + (d2-rmin)*Gr[2][2];
      ps1  = g[0]*Gr[0][1] + g[1]*Gr[1][1] + g[2]*Gr[2][1];
      ps2  = g[0]*Gr[0][2] + g[1]*Gr[1][2] + g[2]*Gr[2][2];
      
      if ( ps1 < EPS1 && ps2 < EPS1 ) dist = fabs(rmin);
      
    }
    /* Then try rmax if no update from rmin */
    if ( ( fabs( dist - INIVAL_3d ) < EPS2 ) &&  (rmax >= d1 && rmax >= d2) ) {
      g[0] = (d1-rmax)*Gr[0][1] + (d2-rmax)*Gr[0][2];
      g[1] = (d1-rmax)*Gr[1][1] + (d2-rmax)*Gr[1][2];
      g[2] = (d1-rmax)*Gr[2][1] + (d2-rmax)*Gr[2][2];
      ps1  = g[0]*Gr[0][1] + g[1]*Gr[1][1] + g[2]*Gr[2][1];
      ps2  = g[0]*Gr[0][2] + g[1]*Gr[1][2] + g[2]*Gr[2][2];
      
      if ( ps1 < EPS1 && ps2 < EPS1 ) dist = fabs(rmax);
    }
  }
  
  /* If no other value has been assigned to dist, calculate a trial value based on both triangle edges */
  ll = (o1[0]-a[0])*(o1[0]-a[0]) + (o1[1]-a[1])*(o1[1]-a[1]) + (o1[2]-a[2])*(o1[2]-a[2]);
  ll = sqrt(ll);
  defval1 = d1 + ll;
  dist = D_MIN(fabs(defval1),dist);
  
  ll = (o2[0]-a[0])*(o2[0]-a[0]) + (o2[1]-a[1])*(o2[1]-a[1]) + (o2[2]-a[2])*(o2[2]-a[2]);
  ll = sqrt(ll);
  defval2 = d2 + ll;
  dist = D_MIN(fabs(defval2),dist);
    
  return (dist);
  
}

/* Calculate a (positive) active value at vertex i in triangle start based on the values in the other two vertices */
double actival_s(pMesh mesh,pSol sol,int start,int i) {
  pTria         pt,pt1;
  pPoint        p0,p1,p2,pa;
  double        dist,d0,d1,d2,da,ll,ps;
  int           k,ip,ipa,ip1,ip2,iel,ilist,list[LONMAX];
  char          ia,i1,i2,j;
  
  i1 = inxt2[i];
  i2 = inxt2[i1];
  
  pt = &mesh->tria[start];
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
    ps = (p1->c[0]-p0->c[0])*(p2->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])*(p2->c[1]-p0->c[1]) + (p1->c[2]-p0->c[2])*(p2->c[2]-p0->c[2]);
    
    /* Acute angle */
    if ( ps > 0.0 ) {
      ll = (p2->c[0]-p0->c[0])*(p2->c[0]-p0->c[0]) + (p2->c[1]-p0->c[1])*(p2->c[1]-p0->c[1]) + (p2->c[2]-p0->c[2])*(p2->c[2]-p0->c[2]);
      ll = sqrt(ll);
      dist = (d0 < 0.0) ? d2 -ll : d2 + ll;
      return( fabs(dist) );
    }
    
    /* Obtuse angle: look for other points in the ball of p2 */
    ilist = boulet_2d(mesh,start,i2,list);
    for (k=0; k<ilist; k++) {
      iel = list[k] / 3;
      ia = list[k] % 3;
      pt1 = &mesh->tria[iel];
      for (j=0; j<2; j++) {
        ia = inxt2[ia];
        ipa = pt1->v[ia];
        pa = &mesh->point[ipa];
        if ( pa->tag != 1 ) continue;
        da = sol->val[ipa];
        dist = actival2pt_s(&p2->c[0],&pa->c[0],&p0->c[0],fabs(d2),fabs(da));
        return(dist);
      }
    }
    

  }
  /* Else if ip2 is not accepted, calculate a trial value based on ip1 if the angle at ip0 is acute, or search for another accepted value otherwise */
  else if ( p2->tag != 1 ) {
    ps = (p1->c[0]-p0->c[0])*(p2->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])*(p2->c[1]-p0->c[1]) + (p1->c[2]-p0->c[2])*(p2->c[2]-p0->c[2]);
    
    /* Acute angle */
    if ( ps > 0.0 ) {
      ll = (p1->c[0]-p0->c[0])*(p1->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])*(p1->c[1]-p0->c[1]) + (p1->c[2]-p0->c[2])*(p1->c[2]-p0->c[2]);
      ll = sqrt(ll);
      dist = (d0 < 0.0) ? d1 -ll : d1 + ll;
      return( fabs(dist) );
    }
    
    /* Obtuse angle: look for other points in the ball of p1 */
    ilist = boulet_2d(mesh,start,i1,list);
    for (k=0; k<ilist; k++) {
      iel = list[k] / 3;
      ia = list[k] % 3;
      pt1 = &mesh->tria[iel];
      for (j=0; j<2; j++) {
        ia = inxt2[ia];
        ipa = pt1->v[ia];
        pa = &mesh->point[ipa];
        if ( pa->tag != 1 ) continue;
        da = sol->val[ipa];
        dist = actival2pt_s(&p1->c[0],&pa->c[0],&p0->c[0],fabs(d1),fabs(da));
        return(dist);
      }
    }
  }
  
  /* At this point, both values d1 and d2 are accepted */
  dist = actival2pt_s(&p1->c[0],&p2->c[0],&p0->c[0],fabs(d1),fabs(d2));
  
  return (dist);
}

/* Propagation of the signed distance function by the Fast Marching Method */
int ppgdistfmm_s(pMesh mesh,pSol sol) {
  Queue       q;
  pQueue      pq;
  pTria       pt,pt1;
  pPoint      p0,p1,p2;
  double      dist;
  int         nacc,k,l,iel,ip,ip1,ip2,ilist,list[LONMAX];
  char        i,j,j1,j2,jj;
    
  pq = &q;
  nacc = 0;
  
  /* Memory allocation for the priority queue */
  if ( !setQueue(mesh,pq) ) {
    printf("Impossible to allocate memory for priority queue. Abort program.\n");
    exit(0);
  }
  
  /* Definition of the initial set of active nodes: travel accepted nodes */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    for(i=0; i<3; i++) {
      ip = pt->v[i];
      p0 = &mesh->point[ip];
      
      if ( p0->tag != 1 ) continue;
      
      ilist = boulet_2d(mesh,k,i,list);
      
      for (l=0; l<ilist; l++) {
        iel = list[l] / 3;
        pt1 = &mesh->tria[iel];
        
        j   = list[l] % 3;
        j1  = inxt2[j];
        j2  = inxt2[j1];
        
        ip1 = pt1->v[j1];
        ip2 = pt1->v[j2];
        p1 = &mesh->point[ip1];
        p2 = &mesh->point[ip2];
        
        /* Calculate value at ip1 and put it in the queue */
        if ( p1->tag != 1 ) {
          dist = actival_s(mesh,sol,iel,j1);
          
          if ( p1->tag == 0 ) {
            insertAnod(pq,ip1,dist);
            p1->tag = 2;
          }
          else
            upAnod(pq,ip1,dist);
        }
        
        /* Calculate value at ip2 and put it in the queue */
        if ( p2->tag != 1 ) {
          dist = actival_s(mesh,sol,iel,j2);
          if ( p2->tag == 0 ) {
            insertAnod(pq,ip2,dist);
            
            p2->tag = 2;
          }
          else
            upAnod(pq,ip2,dist);
        }
      }
    }
  }
  
  /* Main loop: pop the smallest active node; it becomes accepted and the neighboring values become active */
  while ( pq->siz ) {
    ip = popAnod(pq,&dist);
    if ( !ip ) {
      printf("Problem in popping in Fast Marching Method. Abort\n");
      exit(0);
    }
    
    p0 = &mesh->point[ip];
    p0->tag = 1;
    sol->val[ip] = sol->val[ip] > 0.0 ? dist : -dist;
    
    /* Travel the ball of p0 to update the set of active nodes */
    k = p0->s;
    pt = &mesh->tria[k];
    for (i=0; i<3; i++)
      if ( pt->v[i] == ip ) break;
    assert ( i < 3 );
    
    ilist = boulet_2d(mesh,k,i,list);
    for (l=0; l<ilist; l++) {
      iel = list[l] / 3;
      j   = list[l] % 3;
      pt  = &mesh->tria[iel];
      
      for(jj=0; jj<2; jj++) {
        j   = inxt2[j];
        ip1 = pt->v[j];
        p1  = &mesh->point[ip1];
        
        /* Either insert or update active value if the point is not already accepted */
        if ( p1->tag == 1 ) continue;
        
        dist = actival_s(mesh,sol,iel,j);
        if ( p1->tag == 2 )
          upAnod(pq,ip1,dist);
        else {
          insertAnod(pq,ip1,dist);
          p1->tag = 2;
        }
      }
    }
  }
  
  /* Release memory of the queue */
  if ( !freeQueue(pq) ) {
    printf("Impossible to free priority queue. Abort program.\n");
    exit(0);
  }
  
  return(1);
}
