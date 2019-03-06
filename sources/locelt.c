#include "mshdist.h"

#define EPST    -1.e-21 //12
#define EPSR     1.e+21 //12


/* check if p in nsdep */
int inTetra(pMesh mesh,int nsdep,double *p,double *cb) {
	pTetra   pt;
	pPoint   p0,p1,p2,p3;
	double   bx,by,bz,cx,cy,cz,dx,dy,dz,vx,vy,vz,apx,apy,apz;
	double   vol1,vol2,vol3,vol4,dd; 
	
	pt = &mesh->tetra[nsdep];
	if ( !pt->v[0] )  return(0);
	
	p0 = &mesh->point[pt->v[0]];
	p1 = &mesh->point[pt->v[1]];
	p2 = &mesh->point[pt->v[2]];
	p3 = &mesh->point[pt->v[3]];
	
	/* barycentric */
	bx  = p1->c[0] - p0->c[0];
	by  = p1->c[1] - p0->c[1];
	bz  = p1->c[2] - p0->c[2];
	cx  = p2->c[0] - p0->c[0];
	cy  = p2->c[1] - p0->c[1];
	cz  = p2->c[2] - p0->c[2];
	dx  = p3->c[0] - p0->c[0];
	dy  = p3->c[1] - p0->c[1];
	dz  = p3->c[2] - p0->c[2];
	
	/* test volume */
	vx  = cy*dz - cz*dy;
	vy  = cz*dx - cx*dz;
	vz  = cx*dy - cy*dx;
	
	dd  = bx*vx + by*vy + bz*vz;
	apx = p[0] - p0->c[0];
	apy = p[1] - p0->c[1];
	apz = p[2] - p0->c[2];
	
	/* p in 2 */
	vol2  = apx*vx + apy*vy + apz*vz;
	if ( 0.0 > vol2 )  return(0);
	
	/* p in 3 */
	vx  = by*apz - bz*apy;
	vy  = bz*apx - bx*apz;
	vz  = bx*apy - by*apx;
	vol3 = dx*vx + dy*vy + dz*vz;
	if ( 0.0 > vol3 )  return(0);
    
	/* p in 4 */
	vol4 = -cx*vx - cy*vy - cz*vz;
	if ( 0.0 > vol4 )  return(0);
	
	/* p in 1 */
	vol1 = dd - vol2 - vol3 - vol4;
	if ( 0.0 > vol1 )  return(0);
	
	dd = vol1+vol2+vol3+vol4;
	if ( dd != 0.0 )  dd = 1.0 / dd;
	cb[0] = vol1 * dd;
	cb[1] = vol2 * dd;
	cb[2] = vol3 * dd;
	cb[3] = vol4 * dd; 
	
	return(1);
}


/* find tetra containg p, starting nsdep */
int locelt_3d(pMesh mesh,int nsdep,double *p,double *cb) {
	pTetra   pt;
	pPoint   p0,p1,p2,p3;
	double   bx,by,bz,cx,cy,cz,dx,dy,dz,vx,vy,vz,apx,apy,apz;
	double   vol1,vol2,vol3,vol4,dd; 
	int     *adj,base,iadr,it,nsfin;
	
	it    = 0;
	nsfin = nsdep;
	base  = ++mesh->mark;
	do {
		if ( !nsfin )  break;
		pt = &mesh->tetra[nsfin];
		if ( !pt->v[0] )  return(0);
		if ( pt->mark == base )  break;
		pt->mark = base;
		
		iadr = 4*(nsfin-1)+1;
		adj  = &mesh->adja[iadr];
		p0 = &mesh->point[pt->v[0]];
		p1 = &mesh->point[pt->v[1]];
		p2 = &mesh->point[pt->v[2]];
		p3 = &mesh->point[pt->v[3]];
		
		/* barycentric */
		bx  = p1->c[0] - p0->c[0];
		by  = p1->c[1] - p0->c[1];
		bz  = p1->c[2] - p0->c[2];
		cx  = p2->c[0] - p0->c[0];
		cy  = p2->c[1] - p0->c[1];
		cz  = p2->c[2] - p0->c[2];
		dx  = p3->c[0] - p0->c[0];
		dy  = p3->c[1] - p0->c[1];
		dz  = p3->c[2] - p0->c[2];
		
		/* test volume */
		vx  = cy*dz - cz*dy;
		vy  = cz*dx - cx*dz;
		vz  = cx*dy - cy*dx;
		
		dd =  bx*vx + by*vy + bz*vz;
		apx = p[0] - p0->c[0];
		apy = p[1] - p0->c[1];
		apz = p[2] - p0->c[2];
		
		/* p in 2 */
		vol2  = apx*vx + apy*vy + apz*vz;
		if ( 0.0 > vol2 ) {
			nsfin = adj[1] / 4;
			continue;
		}
		
		/* p in 3 */
		vx  = by*apz - bz*apy;
		vy  = bz*apx - bx*apz;
		vz  = bx*apy - by*apx;
		vol3 = dx*vx + dy*vy + dz*vz;
		if ( 0.0 > vol3 ) {
			nsfin = adj[2] / 4;
			continue;
		}
		
		/* p in 4 */
		vol4 = -cx*vx - cy*vy - cz*vz;
		if ( 0.0 > vol4 ) {
			nsfin = adj[3] / 4;
			continue;
		}
		
		/* p in 1 */
		vol1 = dd - vol2 - vol3 - vol4;
		if ( 0.0 > vol1 ) {
			nsfin = adj[0] / 4;
			continue;
		}
		
		dd = vol1+vol2+vol3+vol4;
		if ( dd != 0.0 )  dd = 1.0 / dd;
		cb[0] = vol1 * dd;
		cb[1] = vol2 * dd;
		cb[2] = vol3 * dd;
		cb[3] = vol4 * dd; 
		
		return(nsfin);
	}
	while ( ++it <= mesh->ne );
	
	/* exhaustive search */
	base = ++mesh->mark;
	for (nsfin=1; nsfin<=mesh->ne; nsfin++) {
		pt = &mesh->tetra[nsfin];
		if ( pt->mark != base && inTetra(mesh,nsfin,p,cb) ) {
			return(nsfin);
		}
	}

	return(0);
}


/* find tetra containg p, starting nsdep ; return a positive value when tet has been found, 
   and a negative value when it is not found. In this case, return the last encountered tet, 
   with the bc of the exit point */
int nxtelt_3d(pMesh mesh,int nsdep,double *p,double *cb) {
	pTetra   pt;
	pPoint   p0,p1,p2,p3;
	double   bx,by,bz,cx,cy,cz,dx,dy,dz,vx,vy,vz,apx,apy,apz;
	double   vol1,vol2,vol3,vol4,dd; 
	int     *adj,iadr,it,nsfin,nstmp;
	
	it    = 0;
	nsfin = nsdep;
	do {
		if ( !nsfin )  break;
		pt = &mesh->tetra[nsfin];
		if ( !pt->v[0] )  return(0);
		
    nstmp = nsfin;
    
		iadr = 4*(nsfin-1)+1;
		adj  = &mesh->adja[iadr];
		p0 = &mesh->point[pt->v[0]];
		p1 = &mesh->point[pt->v[1]];
		p2 = &mesh->point[pt->v[2]];
		p3 = &mesh->point[pt->v[3]];
		
		/* barycentric */
		bx  = p1->c[0] - p0->c[0];
		by  = p1->c[1] - p0->c[1];
		bz  = p1->c[2] - p0->c[2];
		cx  = p2->c[0] - p0->c[0];
		cy  = p2->c[1] - p0->c[1];
		cz  = p2->c[2] - p0->c[2];
		dx  = p3->c[0] - p0->c[0];
		dy  = p3->c[1] - p0->c[1];
		dz  = p3->c[2] - p0->c[2];
		
		/* test volume */
		vx  = cy*dz - cz*dy;
		vy  = cz*dx - cx*dz;
		vz  = cx*dy - cy*dx;
		
		dd  = bx*vx + by*vy + bz*vz;
		apx = p[0] - p0->c[0];
		apy = p[1] - p0->c[1];
		apz = p[2] - p0->c[2];
		
		/* p in 2 */
		vol2  = apx*vx + apy*vy + apz*vz;
		if ( 0.0 > vol2 ) {
			nsfin = adj[1] / 4;
			continue;
		}
		
		/* p in 3 */
		vx  = by*apz - bz*apy;
		vy  = bz*apx - bx*apz;
		vz  = bx*apy - by*apx;
		vol3 = dx*vx + dy*vy + dz*vz;
		if ( 0.0 > vol3 ) {
			nsfin = adj[2] / 4;
			continue;
		}
		
		/* p in 4 */
		vol4 = -cx*vx - cy*vy - cz*vz;
		if ( 0.0 > vol4 ) {
			nsfin = adj[3] / 4;
			continue;
		}
		
		/* p in 1 */
		vol1 = dd - vol2 - vol3 - vol4;
		if ( 0.0 > vol1 ) {
			nsfin = adj[0] / 4;
			continue;
		}
		
		dd = vol1+vol2+vol3+vol4;
		if ( dd != 0.0 )  dd = 1.0 / dd;
		cb[0] = vol1 * dd;
		cb[1] = vol2 * dd;
		cb[2] = vol3 * dd;
		cb[3] = vol4 * dd; 
		
		return(nsfin);
	}
	while ( ++it <= mesh->ne );
  
  /* Analyze why no tetra has been found ; in this case, nstmp is the last non 0 tet */
  if( !nsfin ) {
    vol2  = apx*vx + apy*vy + apz*vz;
		vol3 = dx*vx + dy*vy + dz*vz;
    vol4 = -cx*vx - cy*vy - cz*vz;
		vol1 = dd - vol2 - vol3 - vol4;
    
    dd = 0.0;
    if( vol2 > 0.0 ) dd += vol2;
    if( vol3 > 0.0 ) dd += vol3;
    if( vol4 > 0.0 ) dd += vol4;
    if( vol1 > 0.0 ) dd += vol1;
    
    if( dd != 0.0 ) dd = 1.0 / dd;
    cb[0] = vol1 > 0.0 ? vol1*dd : 0.0;
    cb[1] = vol2 > 0.0 ? vol2*dd : 0.0;
    cb[2] = vol3 > 0.0 ? vol3*dd : 0.0;
    cb[3] = vol4 > 0.0 ? vol4*dd : 0.0;
    
    return( -nstmp );            
  }

	return(0);
}


int inTria(pMesh mesh,int nsdep,double *p,double *cb) {
  pTria      pt;
  pPoint     p0,p1,p2;
  double     ax,ay,bx,by,cx,cy;
  double     dd,aire1,aire2,aire3;
  int        isign;

  pt = &mesh->tria[nsdep];
  if ( !pt->v[0] )  return(0);

  p0 = &mesh->point[pt->v[0]];
  p1 = &mesh->point[pt->v[1]];
  p2 = &mesh->point[pt->v[2]];

  ax = p1->c[0] - p0->c[0];
  ay = p1->c[1] - p0->c[1];
  bx = p2->c[0] - p0->c[0];
  by = p2->c[1] - p0->c[1];
  dd = ax*by - ay*bx;
  isign= dd > 0.0 ? 1 : -1;

  /* barycentric */
  bx = p[0] - p1->c[0];
  by = p[1] - p1->c[1];
  cx = p[0] - p2->c[0];
  cy = p[1] - p2->c[1];
  aire1 = isign*(bx*cy - by*cx);
  if ( aire1 < 0.0 )  return(0);

  ax = p[0] - p0->c[0];
  ay = p[1] - p0->c[1];
  aire2 = isign*(cx*ay - cy*ax);
  if ( aire2 < 0.0 )  return(0);
  
  aire3 = isign*dd - aire1 - aire2;
  if ( aire3 < 0.0 )  return(0);

  aire1 = D_MAX(aire1,0.0);
  aire2 = D_MAX(aire2,0.0);
  aire3 = D_MAX(aire3,0.0);

  dd = aire1+aire2+aire3;
	if ( dd != 0.0 )  {
      dd = 1.0 / dd;
      cb[0] = aire1 * dd;
      cb[1] = aire2 * dd;
      cb[2] = aire3 * dd;
	}
  return(1);
}

/* Find triangle containing Point p, starting nsdep; return the barycentric coordinates of the point in this triangle in cb */
int locelt_2d(pMesh mesh,int nsdep,double *p,double *cb) {
  pTria     pt;
  pPoint    p0,p1,p2;
  double    ax,ay,bx,by,cx,cy;
  double    aire1,aire2,aire3,dd; 
  int      *adja,base,iadr,it,isign,nsfin;

  it    = 0;
  nsfin = nsdep;   //current triangle
  base  = ++mesh->mark;
  do {
    if ( !nsfin )  break;
    pt = &mesh->tria[nsfin];
    if ( !pt->v[0] )  return(0);
    if ( pt->mark == base )  break;  // current simplex has already been travelled

    pt->mark = base;
    iadr = 3*(nsfin-1)+1;
    adja = &mesh->adja[iadr];

    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];

    /* area of triangle */
    ax = p1->c[0] - p0->c[0];
    ay = p1->c[1] - p0->c[1];
    bx = p2->c[0] - p0->c[0];
    by = p2->c[1] - p0->c[1];
    dd = ax*by - ay*bx;
    isign = dd > 0.0 ? 1 : -1;

    /* barycentric */
    bx = p1->c[0] - p[0];
    by = p1->c[1] - p[1];
    cx = p2->c[0] - p[0];
    cy = p2->c[1] - p[1];

    /* p in half-plane lambda_0 > 0 */
    aire1 = isign*(bx*cy - by*cx);
    if ( aire1 < 0.0 ) {
      nsfin = adja[0] / 3;
      continue;
    }

    ax = p0->c[0] - p[0];
    ay = p0->c[1] - p[1];
    aire2 = isign*(cx*ay - cy*ax);
    if ( aire2 < 0.0 ) {
      nsfin = adja[1] / 3;
      continue;
    }

    aire3 = isign*dd - aire1 - aire2;
    if ( aire3 < 0.0 ) {
      nsfin = adja[2] / 3;
      continue;
    }

    aire1 = D_MAX(aire1,0.0);
    aire2 = D_MAX(aire2,0.0);
    aire3 = D_MAX(aire3,0.0);
    dd    = aire1 + aire2 + aire3;
    if ( dd != 0.0 )  dd = 1.0 / dd;
    cb[0] = aire1 * dd;
    cb[1] = aire2 * dd;
    cb[2] = aire3 * dd;

    return(nsfin);
  }
  while ( ++it <= mesh->nt );

  /* exhaustive search */
  base = ++mesh->mark;
  for (nsfin=1; nsfin<=mesh->nt; nsfin++) {
    pt = &mesh->tria[nsfin];
    if ( pt->mark != base && inTria(mesh,nsfin,p,cb) ) {
      return(nsfin);
    }
  }

  return(0);
}

/* Weak version of locelt : finds the next triangle to inspect in a locelt iteration */
int nxtelt_2d(pMesh mesh,int nsdep,double *p,double *cb) {
	pTria     pt;
	pPoint    p0,p1,p2;
	double    ax,ay,bx,by,cx,cy;
	double    epsra,aire1,aire2,aire3,dd; 
	int      *adja,iadr,it,isign,nsfin;
	
	it    = 0;
	nsfin = nsdep;
	do {
		if ( !nsfin )  break;
		pt = &mesh->tria[nsfin];
		if ( !pt->v[0] )  return(0);
		
		iadr = 3*(nsfin-1)+1;
		adja = &mesh->adja[iadr];
		
		p0 = &mesh->point[pt->v[0]];
		p1 = &mesh->point[pt->v[1]];
		p2 = &mesh->point[pt->v[2]];
		
		/* area of triangle */
		ax = p1->c[0] - p0->c[0];
		ay = p1->c[1] - p0->c[1];
		bx = p2->c[0] - p0->c[0];
		by = p2->c[1] - p0->c[1];
		dd = ax*by - ay*bx;
		isign = dd > 0.0 ? 1 : -1;
		
		/* barycentric */
		bx = p1->c[0] - p[0];
		by = p1->c[1] - p[1];
		cx = p2->c[0] - p[0];
		cy = p2->c[1] - p[1];
		
		/* p in 1 */
		aire1 = isign*(bx*cy - by*cx);
		if ( aire1 < 0.0 ) {
			nsfin = adja[0] / 3;
			continue;
		}
		
		ax = p0->c[0] - p[0];
		ay = p0->c[1] - p[1];
		aire2 = isign*(cx*ay - cy*ax);
		if ( aire2 < 0.0 ) {
			nsfin = adja[1] / 3;
			continue;
		}
		
		aire3 = isign*dd - aire1 - aire2;
		if ( aire3 < 0.0 ) {
			nsfin = adja[2] / 3;
			continue;
		}
		
		aire1 = D_MAX(aire1,0.0);
		aire2 = D_MAX(aire2,0.0);
		aire3 = D_MAX(aire3,0.0);
		dd    = aire1 + aire2 + aire3;
		if ( dd != 0.0 )  dd = 1.0 / dd;
		cb[0] = aire1 * dd;
		cb[1] = aire2 * dd;
		cb[2] = aire3 * dd;
		
		return(nsfin);
	}
	while ( ++it <= mesh->nt );
	
	return(0);
}
