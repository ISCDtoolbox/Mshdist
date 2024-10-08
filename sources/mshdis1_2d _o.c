#include "mshdist.h"

extern char  ddb;
unsigned char inxt1[2] = {1,0};
extern unsigned char inxt2[5];

/* Return (squared) signed distance from point c to half space delimited by o, with (non normalized) outer normal vector n */
double distptHS_2d(double c[2],double o[2],double n[2]) {
  double dd,nn,ps;
  
  nn = n[0]*n[0] + n[1]*n[1];
  if ( nn < EPS1 )
    return(0.0);
  
  ps = (c[0]-o[0])*n[0] + (c[1]-o[1])*n[1];
  dd = ( ps > 0.0 ) ? ps*ps*nn : -ps*ps/nn;
  
  return(dd);
}

/* Return normalized unit normal vector to LS function sol in triangle iel */
int norLS_2d(pMesh mesh,pSol sol,int iel,double n[2]) {
  pTria     pt;
  pPoint    p0,p1,p2;
  double    dd,d0,d1,d2,det,idet,m[2][2],im[2][2];
  int       ip0,ip1,ip2;
  
  pt = &mesh->tria[iel];
  ip0 = pt->v[0];
  ip1 = pt->v[1];
  ip2 = pt->v[2];

  p0 = &mesh->point[ip0];
  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];
  
  d0 = sol->val[ip0];
  d1 = sol->val[ip1];
  d2 = sol->val[ip2];

  m[0][0] = p1->c[0] - p0->c[0];     m[0][1] = p1->c[1] - p0->c[1];
  m[1][0] = p2->c[0] - p0->c[0];     m[1][1] = p2->c[1] - p0->c[1];
  
  det = m[0][0]*m[1][1] - m[1][0]*m[0][1];
  if ( det < EPS1 ) return(0);

  idet = 1.0 / det;
  
  im[0][0] = idet*m[1][1];     im[0][1] = -idet*m[0][1];
  im[1][0] = -idet*m[1][0];    im[1][1] = idet*m[0][0];
  
  n[0] = im[0][0] * (d1-d0) + im[0][1] * (d2-d0);
  n[1] = im[1][0] * (d1-d0) + im[1][1] * (d2-d0);

  dd = n[0]*n[0] + n[1]*n[1];
  if ( dd < EPS1 ) return(0);
  dd = 1.0 / sqrt(dd);

  n[0] *= dd;
  n[1] *= dd;
  
  return(1);
}

/* Calculate exit point of half-line starting from c, with (normalized) direction u;
   barycentric coordinates are stored in cb;
   return: - 0: failure,
           - 1: exit through edge
           - 2: exit through vertex */
int exitPt_2d(pMesh mesh,int k,double c[2],double u[2],double cb[3]) {
  pTria   pt;
  pPoint  p0,p1,p2;
  double  det,idet,alpha,dd,ps,cb0[3],m[2][2],im[2][2],Gr[2][3];
  int     imin;
  char    i,i1,i2;
  
  pt = &mesh->tria[k];
  p0 = &mesh->point[pt->v[0]];
  p1 = &mesh->point[pt->v[1]];
  p2 = &mesh->point[pt->v[2]];

  /* Gr[*][i] = (column vector) gradient of \lambda_i */
  m[0][0] = p1->c[0] - p0->c[0];     m[0][1] = p1->c[1] - p0->c[1];
  m[1][0] = p2->c[0] - p0->c[0];     m[1][1] = p2->c[1] - p0->c[1];
  
  det = m[0][0]*m[1][1] - m[1][0]*m[0][1];
  if ( det < EPS1 ) return(0);

  idet = 1.0 / det;
  
  im[0][0] = idet*m[1][1];     im[0][1] = -idet*m[0][1];
  im[1][0] = -idet*m[1][0];    im[1][1] = idet*m[0][0];
  
  Gr[0][0] = -im[0][0] - im[0][1];    Gr[0][1] = im[0][0];    Gr[0][2] = im[0][1];
  Gr[1][0] = -im[1][0] - im[1][1];    Gr[1][1] = im[1][0];    Gr[1][2] = im[1][1];
  
  /* Barycentric coordinates of c */
  cb0[0] = 1.0 + Gr[0][0]*(c[0]-p0->c[0]) + Gr[1][0]*(c[1]-p0->c[1]);
  cb0[1] =       Gr[0][1]*(c[0]-p0->c[0]) + Gr[1][1]*(c[1]-p0->c[1]);
  cb0[2] =       Gr[0][2]*(c[0]-p0->c[0]) + Gr[1][2]*(c[1]-p0->c[1]);

  /* Search for length of exit ray */
  imin = -1;
  for (i=0; i<3; i++) {
    ps = u[0]*Gr[0][i] + u[1]*Gr[1][i];
    if ( fabs(ps) < EPS1 ) continue;
    
    dd = - cb0[i] / ps;
    if ( dd > 0.0 ) {
      if ( imin < 0 ) {
        alpha = dd;
        imin = i;
      }
      else {
        if ( dd < alpha ) {
          alpha = dd;
          imin = i;
        }
      }
    }
  }
  
  if ( imin < 0 ) return(0);
  i1 = inxt2[imin];
  i2 = inxt2[i1];
  
  /* Update barycentric coordinates */
  cb[imin]  = 0.0;
  cb[i1] = cb0[i1] + alpha*(u[0]*Gr[0][i1] + u[1]*Gr[1][i1]);
  cb[i2] = cb0[i2] + alpha*(u[0]*Gr[0][i2] + u[1]*Gr[1][i2]);

  if ( cb[i1] < EPS || cb[i2] < EPS ) return(2);
  return(1);
}

/* Return 1 if triangle k is crossed by line passing through c, with normal n,
          0 otherwise
*/
int isCrossed_2d(pMesh mesh,int k,double c[2],double n[2]) {
  pTria    pt;
  pPoint   p0,p1,p2;
  double   alpha,ps,det,idet,u[2],m[2][2],im[2][2],Gr[2][3],cb[3],cb1[3];
  char     i,i1,i2;
  
  pt = &mesh->tria[k];
  p0 = &mesh->point[pt->v[0]];
  p1 = &mesh->point[pt->v[1]];
  p2 = &mesh->point[pt->v[2]];

  /* Gr[*][i] = (column vector) gradient of \lambda_i */
  m[0][0] = p1->c[0] - p0->c[0];     m[0][1] = p1->c[1] - p0->c[1];
  m[1][0] = p2->c[0] - p0->c[0];     m[1][1] = p2->c[1] - p0->c[1];
  
  det = m[0][0]*m[1][1] - m[1][0]*m[0][1];
  if ( det < EPS1 ) return(0);

  idet = 1.0 / det;
  
  im[0][0] = idet*m[1][1];     im[0][1] = -idet*m[0][1];
  im[1][0] = -idet*m[1][0];    im[1][1] = idet*m[0][0];
  
  Gr[0][0] = -im[0][0] - im[0][1];    Gr[0][1] = im[0][0];    Gr[0][2] = im[0][1];
  Gr[1][0] = -im[1][0] - im[1][1];    Gr[1][1] = im[1][0];    Gr[1][2] = im[1][1];
  
  /* Barycentric coordinates of c */
  cb[0] = 1.0 + Gr[0][0]*(c[0]-p0->c[0]) + Gr[1][0]*(c[1]-p0->c[1]);
  cb[1] =       Gr[0][1]*(c[0]-p0->c[0]) + Gr[1][1]*(c[1]-p0->c[1]);
  cb[2] =       Gr[0][2]*(c[0]-p0->c[0]) + Gr[1][2]*(c[1]-p0->c[1]);
  
  /* Tangent vector */
  u[0] = n[1];
  u[1] = -n[0];
  
  for (i=0; i<3; i++) {
    /* Intersection point with edge i */
    ps     = u[0]*Gr[0][i] + u[1]*Gr[1][i];
    if ( fabs(ps) < EPS ) continue;
    alpha  = - cb[i] / ps;
    
    i1 = inxt2[i];
    i2 = inxt2[i1];
    cb1[i1] = cb[i1] + alpha*(u[0]*Gr[0][i1] + u[1]*Gr[1][i1]);
    cb1[i2] = cb[i2] + alpha*(u[0]*Gr[0][i2] + u[1]*Gr[1][i2]);
        
    if ( (cb1[i1] > -EPS) && (cb1[i2] > -EPS) ) break;
  }
  
  if ( i < 3 )
    return(1);
  else
    return(0);
}

/* Return 1 if triangle k is crossed by 0 level set of phi,
          0 otherwise
   Unit normal vector to 0 level set stored in n
*/
int isCrossed_LS_2d(pMesh mesh,pSol phi,int k,double n[2]) {
  pTria    pt;
  pPoint   p0,p1,p2;
  double   v0,v1,v2,dd,det,idet,m[2][2],im[2][2];
  int      ip0,ip1,ip2;
  
  pt  = &mesh->tria[k];
  ip0 = pt->v[0];
  ip1 = pt->v[1];
  ip2 = pt->v[2];
  
  p0  = &mesh->point[ip0];
  p1  = &mesh->point[ip1];
  p2  = &mesh->point[ip2];
  
  v0  = phi->val[ip0];
  v1  = phi->val[ip1];
  v2  = phi->val[ip2];
  
  if ( ( v0*v1 > EPS1 ) && ( v0*v2 > EPS1 ) && ( v1*v2 >= EPS1 ) ) return(0);

  /* Calculate tangent vector */
  m[0][0] = p1->c[0] - p0->c[0];     m[0][1] = p1->c[1] - p0->c[1];
  m[1][0] = p2->c[0] - p0->c[0];     m[1][1] = p2->c[1] - p0->c[1];
  
  det = m[0][0]*m[1][1] - m[1][0]*m[0][1];
  if ( det < EPS1 ) return(0);

  idet = 1.0 / det;
  
  im[0][0] = idet*m[1][1];     im[0][1] = -idet*m[0][1];
  im[1][0] = -idet*m[1][0];    im[1][1] = idet*m[0][0];
  
  n[0] = im[0][0]*(v1-v0) + im[0][1]*(v2-v0);
  n[1] = im[1][0]*(v1-v0) + im[1][1]*(v2-v0);
  
  dd   = n[0]*n[0] + n[1]*n[1];
  if ( dd < EPS1 ) return(0);
  
  dd = 1.0 / sqrt(dd);
  n[0] *= dd;
  n[1] *= dd;
  
  return(1);
}

/* Calculate an active value at vertex i in triangle k based on the values in the other two vertices,
       when distance is calculated tangentially to u */
double actival_tan_2d(pMesh mesh,pSol psi,int k,char i,double u[2]) {
  pTria         pt;
  pPoint        p0,p1,p2;
  double        dist,dist1,dist2,d0,d1,d2,ll;
  int           ip,ip1,ip2;
  char          i1,i2;
  
  i1 = inxt2[i];
  i2 = inxt2[i1];
  
  pt  = &mesh->tria[k];
  ip  = pt->v[i];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];
  p0  = &mesh->point[ip];
  p1  = &mesh->point[ip1];
  p2  = &mesh->point[ip2];
  
  d0  = fabs(psi->val[ip]);
  d1  = fabs(psi->val[ip1]);
  d2  = fabs(psi->val[ip2]);
  
  dist1 = INIVAL_2d;
  dist2 = INIVAL_2d;

  /* If ip1 is accepted, calculate a trial value based on ip1 */
  if ( p1->tag == 1 ) {
    ll = (p0->c[0]-p1->c[0])*u[0] + (p0->c[1]-p1->c[1])*u[1];
    ll = fabs(ll);
    dist1 = d1 + ll;
  }
  
  /* If ip2 is accepted, calculate a trial value based on ip2 */
  if ( p2->tag == 1 ) {
    ll = (p0->c[0]-p2->c[0])*u[0] + (p0->c[1]-p2->c[1])*u[1];
    ll = fabs(ll);
    dist2 = d2 + ll;
  }

  dist = D_MIN(dist1,dist2);

  return(dist);
}

/* Create adjacency relations within a line mesh: for each elt k in mesh, mesh->adja[2*(k-1)+i], i=0,1 is of the form 2*l+j, where
   - l is the neighbor to k opposite to endpoint i
   - j is the index of the vertex in l opposite to k */
int hashelt_1d(pMesh mesh) {
  pEdge   pe;
  int     k,kk,l,ip,*adja,*hcode;
  char    i,ii;
  
  /* Allocate memory: Todo: should go to loadmesh */
  hcode = (int*)calloc(mesh->np+1,sizeof(int));
  assert(hcode);
  
  mesh->adja = (int*)calloc(2*mesh->na+1,sizeof(int));
  assert(mesh->adja);
  adja = &mesh->adja[0];
  
  /* Fill adja table */
  for (k=1; k<=mesh->na; k++) {
    pe  = &mesh->edge[k];
    
    for (i=0; i<2; i++) {
      ip = pe->v[inxt1[i]];
      
      if ( !hcode[ip] )
        hcode[ip] = 2*k+i;
      else {
        l  = hcode[ip];
        kk = l / 2;
        ii = l % 2;
        
        adja[2*(k-1)+i+1]    = 2*kk + ii;
        adja[2*(kk-1)+ii+1]  = 2*k + i;
      }
    }
  }
  
  free(hcode);
  return(1);
}

/* Gives a consistent orientation to a 1d mesh and check that mesh is open */
int orimesh_1d(pMesh mesh) {
  pEdge pe,pe1;
  int   k,kk,ncc,nor,ipil,ip,iadr,*pile,*adja,*adjb;
  char  i,ii,ii1,isopen;
  
  pile = (int*)calloc(mesh->na+1,sizeof(int));
  pile[1] = 1;
  ipil    = 1;
  
  isopen = 0;
  ncc = 0;
  nor = 0;
  
  while ( ipil > 0 ) {
    ncc++;
    
    do {
      k = pile[ipil--];
      pe = &mesh->edge[k]; 
      pe->flag = ncc;
      
      adja = &mesh->adja[2*(k-1)+1];
      
      for (i=0; i<2; i++) {
        kk = adja[i] / 2;
        ii = adja[i] % 2;
        if ( !kk ) {
          isopen = 1;
          continue;
        }
        
        /* Store adjacent */
        pe1 = &mesh->edge[kk];
        if ( pe1->flag == ncc ) continue;
        
        pe1->flag = ncc;
        pile[++ipil] = kk;
        
        /* Change orientation of kk */
        if ( ii == i ) {
          ii1 = inxt1[ii];
          nor++;
          ip = pe1->v[0];
          pe1->v[0] = pe1->v[1];
          pe1->v[1] = ip;
          
          /* Change voyeur in adjacency relation of k */
          adja[i] = 2*kk + ii1;
          
          /* Change adjacencies of kk */
          adjb = &mesh->adja[2*(kk-1)+1];
          iadr = adjb[0];
          adjb[0] = adjb[1];
          adjb[1] = iadr;
        }
      }
      
    }
    while ( ipil > 0 );
      
    /* Find next unmarked edge */
    ipil = 0;
    for (kk=1; kk<=mesh->na; kk++) {
      pe = &mesh->edge[kk];
      if ( pe->flag == 0 ) {
        pile[++ipil] = kk;
        break;
      }
    }
  }
  
  /* Record */
  if ( !isopen ) {
    printf("    *** Error in mode selection: mesh is not open.\n");
    return(0);
  }
  
  printf("%d connected components, %d reoriented edges.\n",ncc,nor);
  
  free(pile);
  return(1);
}

/* Calculate value of psi at pt i in tria k as normal extension of accepted values */
double norval_2d(pMesh mesh,pSol phi,pSol psi,int k,char i) {
  pTria    pt;
  pPoint   p0,p1,p2;
  double   d0,d1,d2,v1,v2,vnor,ps,det,idet,g[2],m[2][2],im[2][2],Gr[2][3];
  int      ip,ip1,ip2;
  char     i1,i2;
  
  i1 = inxt2[i];
  i2 = inxt2[i1];
  
  pt  = &mesh->tria[k];
  ip  = pt->v[i];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];
  p0  = &mesh->point[ip];
  p1  = &mesh->point[ip1];
  p2  = &mesh->point[ip2];
  
  d0  = phi->val[ip];
  d1  = phi->val[ip1];
  d2  = phi->val[ip2];

  v1  = psi->val[ip1];
  v2  = psi->val[ip2];
  
  if ( p1->tag != 1 && p2->tag != 1 ) printf("IMPOSSIBLE %d\n",ip);
  
  if ( p1->tag != 1 ) {
    // printf("Update from one accepted value only %d\n",ip);
    return(v2);
  }
  else if ( p2->tag != 1 ) {
    // printf("Update from one accepted value only %d\n",ip);
    return(v1);
  }
  /* Gr[*][i] = (column vector) gradient of \lambda_i */
  m[0][0] = p1->c[0] - p0->c[0];     m[0][1] = p1->c[1] - p0->c[1];
  m[1][0] = p2->c[0] - p0->c[0];     m[1][1] = p2->c[1] - p0->c[1];
  
  det = m[0][0]*m[1][1] - m[1][0]*m[0][1];
  if ( det < EPS1 ) return(v1);

  idet = 1.0 / det;
  
  im[0][0] = idet*m[1][1];     im[0][1] = -idet*m[0][1];
  im[1][0] = -idet*m[1][0];    im[1][1] = idet*m[0][0];
  
  Gr[0][0] = -im[0][0] - im[0][1];    Gr[0][1] = im[0][0];    Gr[0][2] = im[0][1];
  Gr[1][0] = -im[1][0] - im[1][1];    Gr[1][1] = im[1][0];    Gr[1][2] = im[1][1];
  
  /* Gradient of phi */
  g[0] = d0*Gr[0][0] + d1*Gr[0][1] + d2*Gr[0][2];
  g[1] = d0*Gr[1][0] + d1*Gr[1][1] + d2*Gr[1][2];
  
  ps = Gr[0][0]*g[0] + Gr[1][0]*g[1];
  
  if ( fabs(ps) < EPS1 ) {
    printf("Orthogonal update %d\n",ip);
    return(v1);
  }
  vnor = -(v1*Gr[0][1]+v2*Gr[0][2])*g[0] - (v1*Gr[1][1]+v2*Gr[1][2])*g[1];
  vnor /= ps;
  
  return(vnor);
}

/* Step 1 in the definition of two LS functions for open mesh2:
   -
   -
   -
*/
int iniLS_open_2d(Info info,pMesh mesh1,pMesh mesh2,pSol sol,pSol phi,pSol psi,double *nor,pBucket bucket) {
  pTria      pt,pt1;
  pEdge      pe;
  pPoint     p0,p1,p2,pa,pb,ppt;
  double     d,dd,det,c[2],cb[3],u[2],v[2],t[2],*n;
  int        base,iadr,ilist,iball,k,l,kk,ip,ip1,ip2,iel,jel,ier,cur,nc,tag,*adja,*adja2,*list,*ball;
  char       i,j,j0,j1,j2,ia,ib,i1,voy;

  for (k=1; k<=sol->np; k++) {
    sol->val[k] = INIVAL_2d;
    phi->val[k] = INIVAL_2d;
    psi->val[k] = INIVAL_2d;
  }
  
  /* Memory allocation */
  list = (int*)calloc(mesh1->nt+1,sizeof(int));
  ball = (int*)calloc(LONMAX,sizeof(int));
  assert(list);
  assert(ball);
  
  /* Set tag = 2 to elts intersected, 1 if connected to elt intersected, 0 else */
  nc  = 0;
  
  for (k=1; k<=mesh2->na; k++){
    pe  = &mesh2->edge[k];
    p0  = &mesh2->point[pe->v[0]];
    p1  = &mesh2->point[pe->v[1]];
    iel = buckin_2d(mesh1,bucket,p0->c);
    iel = locelt_2d(mesh1,iel,p0->c,cb);
    
    if ( !iel ) {
      iel = buckin_2d(mesh1,bucket,p1->c);
      iel = locelt_2d(mesh1,iel,p1->c,cb);
    }
    
    ilist       = 1;
    list[ilist] = iel;
    base = ++mesh1->flag;
    pt   = &mesh1->tria[iel];
    pt->flag = base;
    cur      = 1;
        
    do {
      iel  = list[cur];
      pt   = &mesh1->tria[iel];

      iadr = 3*(iel-1) + 1;
      adja = &mesh1->adja[iadr];

      for (i=0; i<3; i++) {
        ia = inxt2[i];
        ib = inxt2[ia];
        pa = &mesh1->point[pt->v[ia]];
        pb = &mesh1->point[pt->v[ib]];
        ier = intersec_2d(p0,p1,pa,pb);
        jel = adja[i] / 3;
        if ( !jel ) continue;
    
        pt1 = &mesh1->tria[jel];
        if ( ier && pt1->flag < base ) {
            ilist++;
            list[ilist] = jel;
            pt1->flag   = base;
        }
      }
      cur++;
    }
    while ( cur <= ilist );

    /* List analysis */
    for (kk=1; kk<=ilist; kk++) {
      iel = list[kk];
      pt  = &mesh1->tria[iel];
      
      for (i=0; i<3; i++) {
        ip = pt->v[i];
        pa = &mesh1->point[ip];
        
        /* Orientation of pa w.r.t pe */
        u[0]  = p1->c[0] - p0->c[0];
        u[1]  = p1->c[1] - p0->c[1];
        v[0]  = pa->c[0] - p0->c[0];
        v[1]  = pa->c[1] - p0->c[1];
        det   = u[0]*v[1] - u[1]*v[0];
        
        d  = distpt_2d(p0,p1,pa,&tag);
        
        if ( d < sol->val[ip] ) {
          sol->val[ip] = d;
          phi->val[ip] = ( det > 0.0 ) ? d : -d;
          pa->tag = tag;
        }
      }
    }
  }
  
  /* Correction */
  nc = 0;
  for (k=1; k<=mesh1->np; k++){
    pa = &mesh1->point[k];
    if ( pa->tag < 2 )  continue;

    for (kk=1; kk<=mesh2->na; kk++) {
      pe = &mesh2->edge[kk];
      p0 = &mesh2->point[pe->v[0]];
      p1 = &mesh2->point[pe->v[1]];
      d  = distpt_2d(p0,p1,pa,&tag);
      if ( tag == 1 && d < sol->val[k] ) {
        sol->val[k] = d;
        break;
      }
    }
    pa->tag = 1;
    nc++;
  }
  if ( nc )   fprintf(stdout,"     %d correction(s)\n",nc);
  /* At this point, all points initialized have tag 1 */
  
  /* Travel boundary points of mesh2 */
  for (k=1; k<=mesh2->na; k++){
    pe    = &mesh2->edge[k];
    adja2 = &mesh2->adja[2*(k-1)+1];
    p0 = &mesh2->point[pe->v[0]];
    p1 = &mesh2->point[pe->v[1]];
    
    /* Tangent and normal with correct orientation */
    u[0] = p1->c[0] - p0->c[0];
    u[1] = p1->c[1] - p0->c[1];
    dd   = u[0]*u[0] + u[1]*u[1];
    dd   = sqrt(dd);
    if ( dd < EPS ) continue;
    
    u[0] /= dd;
    u[1] /= dd;
    v[0] = -u[1];
    v[1] = u[0];

    for (i=0; i<2; i++) {
      if ( adja2[i] ) continue;
      i1 = inxt1[i];
      ppt = &mesh2->point[pe->v[i1]];
      
      /* Direction for extension */
      if ( i == 0 ) {
        t[0] = u[0];
        t[1] = u[1];
      }
      else {
        t[0] = -u[0];
        t[1] = -u[1];
      }
      
      iel = buckin_2d(mesh1,bucket,ppt->c);
      iel = locelt_2d(mesh1,iel,ppt->c,cb);
      if ( !iel ) continue;
      
      pt = &mesh1->tria[iel];
      
      for (j=0; j<3; j++) {
        ip = pt->v[j];
        pa = &mesh1->point[ip];
                
        /* Update phi */
        phi->val[ip] = distptHS_2d(pa->c,ppt->c,v);
        
        /* Update normal vector */
        n    = &nor[2*(ip-1)+1];
        n[0] = v[0];
        n[1] = v[1];
                
        /* Update psi */
        psi->val[ip] = distptHS_2d(pa->c,ppt->c,t);
      }
      
      /* Find next points on extended surface for normal vectors + put a tag = 4: special trial value */
      for (j=0; j<3; j++) {
        ip = pt->v[j];
        pa = &mesh1->point[ip];
        
        n  = &nor[2*(ip-1)+1];
        dd = phi->val[ip];
        dd = dd > 0.0 ? sqrt(dd) : -sqrt(fabs(dd));
        c[0] = pa->c[0] - dd*n[0];
        c[1] = pa->c[1] - dd*n[1];
                        
        iball = boulet_2d(mesh1,iel,j,ball);
        
        for (l=0; l<iball; l++) {
          jel = ball[l] / 3;
          j0  = ball[l] % 3;
          pt1 = &mesh1->tria[jel];
          
          j1  = inxt2[j0];
          j2  = inxt2[j1];
          ip1 = pt1->v[j1];
          ip2 = pt1->v[j2];
          p1  = &mesh1->point[ip1];
          p2  = &mesh1->point[ip2];

          if ( !isCrossed_2d(mesh1,jel,c,n) ) continue;
          if ( !p1->tag ) {
            phi->val[ip1] = distptHS_2d(p1->c,c,n);
            nor[2*(ip1-1)+1] = n[0];
            nor[2*(ip1-1)+2] = n[1];
            p1->tag = 4;
          }
          if ( !p2->tag ) {
            phi->val[ip2] = distptHS_2d(p2->c,c,n);
            nor[2*(ip2-1)+1] = n[0];
            nor[2*(ip2-1)+2] = n[1];
            p2->tag = 4;
          }
        }
      }
    }
  }
  
  /* Take square roots */
  for (k=1; k<=mesh1->np; k++) {
    sol->val[k] = sqrt(sol->val[k]);
    phi->val[k] = ( phi->val[k] > 0.0 ) ? sqrt(phi->val[k]) : -sqrt(fabs(phi->val[k]));
    psi->val[k] = ( psi->val[k] > 0.0 ) ? sqrt(psi->val[k]) : -sqrt(fabs(psi->val[k]));
  }
  
  free(list);
  free(ball);
  return(1);
}

/* Step 2 in the definition of two LS functions for open mesh2 */
int ppgSolPhi_open_2d(Info info,pMesh mesh,pSol sol,pSol phi,pSol psi,double *nor) {
  Queue     q;
  pQueue    pq;
  pTria     pt,pt1;
  pPoint    p0,p1,p2;
  double    dist,dd,ps,c[2],n[2];
  int       nacc,k,l,tag,iel,ip,ip1,ip2,ilist,*list;
  char      i,j,j1,j2,jj;
  
  pq    = &q;
  nacc  = 0;
  
  /* Memory allocation */
  list = (int*)calloc(LONMAX,sizeof(int));
  assert(list);
  
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
        
        /* Normal vector and projection on the surface */
        n[0] = nor[2*(ip-1)+1];
        n[1] = nor[2*(ip-1)+2];
        c[0] = p0->c[0] - phi->val[ip]*n[0];
        c[1] = p0->c[1] - phi->val[ip]*n[1];

        j   = list[l] % 3;
        j1  = inxt2[j];
        j2  = inxt2[j1];
        
        ip1 = pt1->v[j1];
        ip2 = pt1->v[j2];
        p1 = &mesh->point[ip1];
        p2 = &mesh->point[ip2];
        
        /* Calculate value at ip1 and put it in the queue */
        if ( p1->tag != 1 ) {
          dist = actival_2d(mesh,sol,iel,j1);

          if ( p1->tag == 0 ) {
            insertAnod(pq,ip1,dist);
            p1->tag = 2;
          }
          
          else if ( p1->tag == 4 ) {
            insertAnod(pq,ip1,dist);
            p1->tag = 3;
          }
          
          else
            upAnod(pq,ip1,dist);
        }
        
        /* Calculate value at ip2 and put it in the queue */
        if ( p2->tag != 1 ) {
          dist = actival_2d(mesh,sol,iel,j2);
          if ( p2->tag == 0 ) {
            insertAnod(pq,ip2,dist);
            p2->tag = 2;
          }
          else if ( p2->tag == 4 ) {
            insertAnod(pq,ip2,dist);
            p2->tag = 3;
          }
          else
            upAnod(pq,ip2,dist);
        }
      }
    }
  }
  
  /* At this point, p0->tag = 1: accepted node;
                    p0->tag = 2: active node, not on the front;
                    p0->tag = 3; active node, on the front. */
  
  /* Main loop: pop the smallest active node; it becomes accepted and the neighboring values become active */
  while ( pq->siz ) {
    ip = popAnod(pq,&dist);
    if ( !ip ) {
      printf("Problem in popping in Fast Marching Method. Abort\n");
      exit(0);
    }
    
    p0 = &mesh->point[ip];
    tag = p0->tag;
    p0->tag = 1;
            
    /* Update sol */
    sol->val[ip] = dist;
    
    /* Update phi: - case p0->tag = 2: give sign with respect to accepted adjacent.
                   - case p0->tag = 3: nothing to do, distance is computed during tag update */
    if ( tag == 2 ) {
      ilist = boulep_2d(mesh,ip,list);
      for (l=0; l<ilist; l++) {
        ip1 = list[l];
        p1  = &mesh->point[ip1];
        if ( p1->tag == 1 ) break;
      }
      phi->val[ip] = phi->val[ip1] > 0.0 ? dist : -dist;
    }
    
    /* Update list of active nodes */
    if ( tag == 3 ) {
      /* Normal vector and projection on the surface */
      n[0] = nor[2*(ip-1)+1];
      n[1] = nor[2*(ip-1)+2];
      c[0] = p0->c[0] - phi->val[ip]*n[0];
      c[1] = p0->c[1] - phi->val[ip]*n[1];
                  
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
                
        /* If triangle is crossed by the front: - set tag = 3 to vertices,
                                                - calculate new values for phi
                                                - transmit normal vector */
        for (jj=0; jj<2; jj++) {
          j   = inxt2[j];
          ip1 = pt->v[j];
          p1  = &mesh->point[ip1];
          if ( p1->tag == 1 ) continue;
          
          dist = actival_2d(mesh,sol,iel,j);
          
          if ( p1->tag == 2 || p1->tag == 3 )
            upAnod(pq,ip1,dist);
          else {
            insertAnod(pq,ip1,dist);
            p1->tag = 2;
          }
          
          if ( isCrossed_2d(mesh,iel,c,n) ) {
            ps = (p1->c[0] - p0->c[0])*n[0] + (p1->c[1] - p0->c[1])*n[1];
            phi->val[ip1] = phi->val[ip] + ps;
            
            nor[2*(ip1-1)+1] = n[0];
            nor[2*(ip1-1)+2] = n[1];
            
            p1->tag = 3;
          }
        }
      }
    }
    
    /* Travel the ball of p0 to update the set of active nodes */
    else if ( tag == 2 ) {
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
        
        for (jj=0; jj<2; jj++) {
          j   = inxt2[j];
          ip1 = pt->v[j];
          p1  = &mesh->point[ip1];
          
          /* Either insert or update active value if the point is not already accepted */
          if ( p1->tag == 1 ) continue;
                    
          dist = actival_2d(mesh,sol,iel,j);
          if ( p1->tag == 2 || p1->tag == 3 )
            upAnod(pq,ip1,dist);
          else {
            insertAnod(pq,ip1,dist);
            p1->tag = 2;
          }
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

/* Initialize signed distance function to extended surface, and secondary LS function psi */
int inidis_open_2d(Info info,pMesh mesh,pSol phi,pSol psi) {
  Queue     q;
  pQueue    pq;
  pTria     pt,pt1;
  pPoint    p0,p1,p2;
  double    d,dist,n[2],u[2],*solTmp;
  int       k,iel,l,ip,ip0,ip1,ip2,nc,nb,ilist,*list,*bndy,*proj;
  char      i,j,jj,j1,j2;
  
  nb   = 0;
  bndy = (int*)calloc(mesh->nt+1,sizeof(int));
  assert(bndy);
  
  /* Reset point tags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].tag = 0;
  
  /* Signed distance function to the primary surface */
  /* Store intersecting triangles in list bndy */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    ip0 = pt->v[0];
    ip1 = pt->v[1];
    ip2 = pt->v[2];

    if ( (phi->val[ip0] * phi->val[ip1] <=0.0) || (phi->val[ip0] * phi->val[ip2] <=0.0) || (phi->val[ip1] * phi->val[ip2] <=0.0)){
      nb++;
      bndy[nb] = k;
    }
  }
  
  bndy = (int*)realloc(bndy,(nb+1)*sizeof(int));
  
  /* Temporary values stored in solTmp */
  solTmp = (double*)calloc(mesh->np+1,sizeof(double));
  assert(solTmp);

  /* Travel list bndy and compute distance at vertices of its elements */
  for (l=1; l<=nb; l++) {
    iel = bndy[l];
    pt = &mesh->tria[iel];
    ip0 = pt->v[0];
    ip1 = pt->v[1];
    ip2 = pt->v[2];

    p0 = &mesh->point[ip0];
    p1 = &mesh->point[ip1];
    p2 = &mesh->point[ip2];
    
    /* Distance from i0 to the 0 level set from triangle bndy[i] */
    d = distnv0_2d(mesh,phi,iel,p0,&proj);
    
    if ( p0->tag == 0 ) {
      solTmp[ip0] = d;
      p0->tag     = proj;
    }
    else if ( d < fabs(solTmp[ip0]) ) {
      solTmp[ip0] = d;
      p0->tag     = proj;
    }
          
    /* Distance from i1 to the 0 level set from triangle bndy[i] */
    d = distnv0_2d(mesh,phi,iel,p1,&proj);
    
    if ( p1->tag == 0 ) {
      solTmp[ip1] = d;
      p1->tag    = proj;
    }
    else if ( d < fabs(solTmp[ip1]) ) {
      solTmp[ip1] = d;
      p1->tag    = proj;
    }
        
    /* Distance from i2 to the 0 level set from triangle bndy[i] */
    d = distnv0_2d(mesh,phi,iel,p2,&proj);
    
    if ( p2->tag == 0 ) {
      solTmp[ip2] = d;
      p2->tag     = proj;
    }
    else if ( d < fabs(solTmp[ip2]) ) {
      solTmp[ip2] = d;
      p2->tag    = proj;
    }
  }
  
  /* Correction procedure, for points with tag = 2 */
  nc = 0;
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    if ( p0->tag < 2 )  continue;
    for (l=1; l<=nb; l++) {
      iel = bndy[l];
      pt = &mesh->tria[iel];
      d  = distnv0_2d(mesh,phi,iel,p0,&proj);
      if ( proj == 1 && d < fabs(solTmp[k]) ) {
        solTmp[k] = d;
        break;
      }
    }
    p0->tag = 1;
    nc++;
  }
    
  /* Take square roots */
  for (k=1; k<=mesh->np; k++){
    p0 = &mesh->point[k];
    if ( !p0->tag )
      phi->val[k] = phi->val[k] < 0.0 ? -sqrt(INIVAL_2d) : sqrt(INIVAL_2d);
    else if ( p0->tag )
      phi->val[k] = phi->val[k] < 0.0 ? -sqrt(solTmp[k]) : sqrt(solTmp[k]);
  }

  free(solTmp);
  free(bndy);
  
  /* Change point tags: - tag = 1: accepted node for psi
                        - tag = 2: active node for psi
                        - tag = 3: "far away" node for psi, but concerned by the calculation of psi
                        - tag = 0: "far away" node from the primary surface */
  for (k=1; k<=mesh->np; k++){
    p0 = &mesh->point[k];
    if ( fabs(psi->val[k]) < INIVAL_2d - EPS )
      p0->tag = 1;
    else if ( p0->tag == 1 )
      p0->tag = 3;
  }
  
  /* Set queue */
  pq = &q;
  if ( !setQueue(mesh,pq) ) {
    printf("Impossible to allocate memory for priority queue. Abort program.\n");
    exit(0);
  }
  
  /* Memory allocation */
  list = (int*)calloc(LONMAX,sizeof(int));
  assert(list);
    
  /* Define the initial set of active nodes for psi: travel accepted nodes */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    
    /* Catch accepted points by triangle with only accepted vertices: optimization possible */
    for(i=0; i<3; i++) {
      ip = pt->v[i];
      p0 = &mesh->point[ip];
      if ( p0->tag != 1 ) break;
    }
    
    if ( i < 3 ) continue;
    
    /* Calculate normal vector to extended surface */
    if ( !norLS_2d(mesh,psi,k,n) ) continue;
    
    for(i=0; i<3; i++) {
      ip = pt->v[i];
      p0 = &mesh->point[ip];
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
        if ( p1->tag == 2 || p1->tag == 3 ) {
          dist = psi->val[ip] + (p1->c[0]-p0->c[0])*n[0] + (p1->c[1]-p0->c[1])*n[1];
          // psi->val[ip1] = ( dist > 0.0 ) ? INIVAL_2d : -INIVAL_2d;
          psi->val[ip1] = dist;

          dist = fabs(dist);
          
          if ( p1->tag == 3 ) {
            insertAnod(pq,ip1,dist);
            p1->tag = 2;
          }
          else
            upAnod(pq,ip1,dist);
        }
        
        /* Calculate value at ip2 and put it in the queue */
        if ( p2->tag == 2 || p2->tag == 3 ) {
          dist = psi->val[ip] + (p2->c[0]-p0->c[0])*n[0] + (p2->c[1]-p0->c[1])*n[1];
          // psi->val[ip2] = ( dist > 0.0 ) ? INIVAL_2d : -INIVAL_2d;
          psi->val[ip2] = dist;
          dist = fabs(dist);
          
          if ( p2->tag == 3 ) {
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
    psi->val[ip] = ( psi->val[ip] > 0.0 ) ? dist : -dist;
        
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
      if ( !isCrossed_LS_2d(mesh,phi,iel,n) ) continue;
      
      u[0] = -n[1];
      u[1] = n[0];
      
      for(jj=0; jj<2; jj++) {
        j   = inxt2[j];
        ip1 = pt->v[j];
        p1  = &mesh->point[ip1];
        
        /* Either insert or update active value if the point is not already accepted */
        if ( !p1->tag || p1->tag == 1 ) continue;
        
        dist = actival_tan_2d(mesh,psi,iel,j,u);
        if ( p1->tag == 2 )
          upAnod(pq,ip1,dist);
        else {
           insertAnod(pq,ip1,dist);
           psi->val[ip1] = ( psi->val[ip] > 0.0 ) ? INIVAL_2d : -INIVAL_2d;
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
  
  free(list);
  return(1);
}

/* Propagation of the signed distance phi and normal extension of psi
   info:  information file;
   mesh:  mesh of computational domain;
   phi:   signed distance function, initialized at the points near the 0 LS;
   psi:   field to be extended in the normal direction, initialized at the points near the 0 LS of phi. */
int norppg_2d(Info info,pMesh mesh,pSol phi,pSol psi) {
  Queue     q;
  pQueue    pq;
  pTria     pt,pt1;
  pPoint    p0,p1,p2;
  double    dist,dold;
  int       k,l,ll,iel,ip,ip1,ip2,tag,ilist,*list;
  char      i,j,j1,j2,jj;
  
  /* Set queue */
  pq = &q;
  if ( !setQueue(mesh,pq) ) {
    printf("Impossible to allocate memory for priority queue. Abort program.\n");
    exit(0);
  }
  
  /* Memory allocation */
  list = (int*)calloc(LONMAX,sizeof(int));
  assert(list);
  
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
          dist = actival_2d(mesh,phi,iel,j1);

          if ( !p1->tag ) {
            insertAnod(pq,ip1,dist);
            p1->tag = 2;
            p1->s   = iel;
          }
          /* To do: optimization here */
          else {
            ll = pq->perm[ip1];
            dold = pq->hp[ll].d;
            if ( dist < dold )
              p1->s = iel;
            upAnod(pq,ip1,dist);
          }
        }
        
        /* Calculate value at ip2 and put it in the queue */
        if ( p2->tag != 1 ) {
          dist = actival_2d(mesh,phi,iel,j2);
          
          if ( !p2->tag ) {
            insertAnod(pq,ip2,dist);
            p2->tag = 2;
            p2->s   = iel;
          }
          /* To do: optimization here */
          else {
            ll = pq->perm[ip2];
            dold = pq->hp[ll].d;
            if ( dist < dold )
              p2->s = iel;
            upAnod(pq,ip2,dist);
          }
        }
      }
    }
  }
  
  /* At this point, p0->tag = 1: accepted node;
                    p0->tag = 2: active node.
                    p0->s contains the index of the triangle inducing the active value */
  
  /* Main loop: pop the smallest active node; it becomes accepted and the neighboring values become active */
  while ( pq->siz ) {
    ip = popAnod(pq,&dist);
    if ( !ip ) {
      printf("Problem in popping in Fast Marching Method. Abort\n");
      exit(0);
    }
    
    p0 = &mesh->point[ip];
    tag = p0->tag;
    p0->tag = 1;
            
    /* Update phi: give sign with respect to accepted adjacent. */
    ilist = boulep_2d(mesh,ip,list);
    for (l=0; l<ilist; l++) {
      ip1 = list[l];
      p1  = &mesh->point[ip1];
      if ( p1->tag == 1 ) break;
    }
    phi->val[ip] = phi->val[ip1] > 0.0 ? dist : -dist;
    
    /* Update psi by normal extension procedure */
    k = p0->s;
    pt = &mesh->tria[k];
    for (i=0; i<3; i++)
      if ( pt->v[i] == ip ) break;
    assert ( i < 3 );
    
    psi->val[ip] = norval_2d(mesh,phi,psi,k,i);

    /* Travel the ball of p0 to update the set of active nodes */
    ilist = boulet_2d(mesh,k,i,list);
    for (l=0; l<ilist; l++) {
      iel = list[l] / 3;
      j   = list[l] % 3;
      pt  = &mesh->tria[iel];
        
      for (jj=0; jj<2; jj++) {
        j   = inxt2[j];
        ip1 = pt->v[j];
        p1  = &mesh->point[ip1];
          
        /* Either insert or update active value if the point is not already accepted */
        if ( p1->tag == 1 ) continue;
                    
        dist = actival_2d(mesh,phi,iel,j);
        if ( p1->tag == 2 ) {
          ll = pq->perm[ip1];
          dold = pq->hp[ll].d;
          if ( dist < dold )
            p1->s = iel;
          upAnod(pq,ip1,dist);
        }
        else {
          insertAnod(pq,ip1,dist);
          p1->tag = 2;
          p1->s = iel;
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

/* Todo: bucksiz is #define in mshdist.c -> harmonize */
/* Generation of two level set functions for a mesh of an open (collection of) curve(s):
   info:  information file;
   mesh:  mesh of computational domain;
   mesh2: 2d line mesh of a curve; may contain multiple components;
   sol:   output unsigned distance function to mesh2;
   phi:   level set function for the primary surface, extending mesh2
   psi:   level set function for the secondary surface. */
int mshdis1_2d_o(Info info,pMesh mesh,pMesh mesh2,pSol sol,pSol phi,pSol psi) {
  pBucket  bucket;
  double   *nor;
  int      ier,bucksiz;
  
  /* Check that mesh2 is open + give consistent orientation */
  ier = hashelt_1d(mesh2);
  if ( !ier )  return(0);
  
  ier = orimesh_1d(mesh2);
  if ( !ier )  return(0);
  
  /* Initialize bucket */
  bucksiz = 16;
  bucket  = newBucket_2d(mesh,bucksiz);
  if ( !bucket )  return(0);
  
  /* (Normalized) Normal vector field to the extended surface */
  nor = (double*)calloc(2*mesh->np+1,sizeof(double));

  /* Step 1: - Initialize sol in triangles intersecting mesh2
             - Initialize phi + its sign (using orientation) in triangle intersecting mesh 2
             - Initialize psi + its sign in triangles intersecting boundary of mesh 2 */
  if ( !iniLS_open_2d(info,mesh,mesh2,sol,phi,psi,nor,bucket) ) return(0);
  
  /* Step 2: Calculate unsigned distance to mesh2 and unravel phi near \tilde S */
  if ( !ppgSolPhi_open_2d(info,mesh,sol,phi,psi,nor) ) return(0);

  free(nor);
  
  /* Step 3: Calculate psi at triangles intersecting \tilde S */
  if ( !inidis_open_2d(info,mesh,phi,psi) ) return(0);
  
  /* Step 4: propagation of phi + normal extension of psi */
  if ( !norppg_2d(info,mesh,phi,psi) ) return(0);
  
  return(1);
}
