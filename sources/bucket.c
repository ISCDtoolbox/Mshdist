#include "mshdist.h"

/* create bucket structure and store initial vertices */
pBucket newBucket_3d(pMesh mesh,int nmax) {
  pPoint   ppt;
  pBucket  bucket;
  double   dd;
  int      k,ic,ii,jj,kk;

  /* memory alloc */
  bucket = (Bucket*)malloc(sizeof(Bucket));
  assert(bucket);
  bucket->size = nmax;
  bucket->head = (int*)calloc(nmax*nmax*nmax+1,sizeof(int));
  assert(bucket->head);
  bucket->link = (int*)calloc(mesh->np+1,sizeof(int));
  assert(bucket->link);

  /* insert vertices */
  dd = nmax / (double)PRECI;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    ii = D_MAX(0,(int)(dd * ppt->c[0])-1);
    jj = D_MAX(0,(int)(dd * ppt->c[1])-1);
    kk = D_MAX(0,(int)(dd * ppt->c[2])-1);
    ic = (kk*nmax + jj)*nmax + ii;

    if ( !bucket->head[ic] )
      bucket->head[ic] = k;
    else {
      bucket->link[k]  = bucket->head[ic];
      bucket->head[ic] = k;
    }
  }

  return(bucket);
}

pBucket newBucket_2d(pMesh mesh,int nmax) {
  pPoint   ppt;
  pBucket  bucket;
  double   dd;
  int      k,ic,ii,jj;

  /* memory alloc */
  bucket = (Bucket*)malloc(sizeof(Bucket));
  assert(bucket);
  bucket->size = nmax;
  bucket->head = (int*)calloc(nmax*nmax+1,sizeof(int));
  assert(bucket->head);
  bucket->link = (int*)calloc(mesh->np+1,sizeof(int));
  assert(bucket->link);

  /* insert vertices */
  dd = nmax / (double)PRECI;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    ii = D_MAX(0,(int)(dd * ppt->c[0])-1);
    jj = D_MAX(0,(int)(dd * ppt->c[1])-1);
    ic = jj*nmax + ii;
    if ( !bucket->head[ic] )
      bucket->head[ic] = k;
    else {
      bucket->link[k]  = bucket->head[ic];
      bucket->head[ic] = k;
    }
  }

  return(bucket);
}


void freeBucket(pBucket bucket) {
  free(bucket->head);
  free(bucket->link);
  free(bucket);
}


int buckin_3d(pMesh mesh,pBucket bucket,double *c) {
  pPoint    pp1;
  double    dm,dd,d2,ux,uy,uz;
  int       i,j,k,ii,jj,kk,d,ic,ip,ip1,siz;
  int       imin,imax,jmin,jmax,kmin,kmax;

  siz = bucket->size;
  dd  = siz / (double)PRECI;
	
  ii = D_MAX(0,(int)(dd * c[0])-1);
  jj = D_MAX(0,(int)(dd * c[1])-1);
  kk = D_MAX(0,(int)(dd * c[2])-1);
  ic = (kk*siz + jj)*siz + ii;

  /* check current cell */
  if ( bucket->head[ic] ) {
    ip1 = bucket->head[ic];
    pp1 = &mesh->point[ip1];
    ux  = pp1->c[0] - c[0];
    uy  = pp1->c[1] - c[1];
    uz  = pp1->c[2] - c[2];
    dm  = ux*ux + uy*uy + uz*uz;
    ip  = ip1;

    while ( bucket->link[ip1] ) {
      ip1 = bucket->link[ip1];
      pp1 = &mesh->point[ip1];
      ux  = pp1->c[0] - c[0];
      uy  = pp1->c[1] - c[1];
      uz  = pp1->c[2] - c[2];
      d2  = ux*ux + uy*uy + uz*uz;
      if ( d2 < dm ) {
        ip = ip1;
        dm = d2;
      }
    }
    pp1 = &mesh->point[ip];
    return(pp1->s);
  }

  /* neighbours */
  d = 1;
  do {
    imin = D_MAX(0,ii-d);
    imax = D_MIN(ii+d,siz-1);
    jmin = D_MAX(0,jj-d);
    jmax = D_MIN(jj+d,siz-1);
    kmin = D_MAX(0,kk-d);
    kmax = D_MIN(kk+d,siz-1);

    for (k=kmin; k<=kmax; k++)
      for (j=jmin; j<=jmax; j++)
        for (i=imin; i<=imax; i++) {
          if ( ii == i && jj == j && kk == k )  continue;
          ic  = (k*siz + j)*siz + i;
          ip1 =  bucket->head[ic];
          if ( !ip1 )  continue;
          
          pp1 = &mesh->point[ip1];
          return(pp1->s);
        }
  }
  while ( ++d < siz / 2 );

  return(0);
}


int buckin_2d(pMesh mesh,pBucket bucket,double *c) {
  pPoint    pp1;
  double    dm,dd,d2,ux,uy;
  int       i,j,ii,jj,d,ic,ip,ip1,siz;
  int       imin,imax,jmin,jmax;
  
  siz = bucket->size;
  dd  = siz / (double)PRECI;

  ii = D_MAX(0,(int)(dd * c[0])-1);
  jj = D_MAX(0,(int)(dd * c[1])-1);
  ic = jj*siz + ii;

  /* check current cell */
  if ( bucket->head[ic] ) {
    ip1 = bucket->head[ic];
    pp1 = &mesh->point[ip1];
    ux  = pp1->c[0] - c[0];
    uy  = pp1->c[1] - c[1];
    dm  = ux*ux + uy*uy;
    ip  = ip1;

    while ( bucket->link[ip1] ) {
      ip1 = bucket->link[ip1];
      pp1 = &mesh->point[ip1];
      ux  = pp1->c[0] - c[0];
      uy  = pp1->c[1] - c[1];
      d2  = ux*ux + uy*uy;
      if ( d2 < dm ) {
        ip = ip1;
        dm = d2;
      }
    }

    pp1 = &mesh->point[ip];
    return(pp1->s);
  }

  else {
    d = 1;
    do {
      imin = ii-d;
      imax = ii+d;
      jmin = jj-d;
      jmax = jj+d;
      imin = D_MAX(0,imin);
      imax = D_MIN(imax,siz-1);
      jmin = D_MAX(0,jmin);
      jmax = D_MIN(siz-1,jmax);
      for (j=jmin; j<=jmax; j++)
        for (i=imin; i<=imax; i++) {
          if ( ii == i && jj == j )  continue;
          ic = j*siz + i;
          if ( bucket->head[ic] ) {
            ip1 = bucket->head[ic];
            pp1 = &mesh->point[ip1];
            return(pp1->s);
          }
        }
    }
    while ( ++d < siz / 2 );
  }

  return(0);
}
