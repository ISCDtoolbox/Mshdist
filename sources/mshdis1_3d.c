#include "mshdist.h"
#include "lplib3.h"

unsigned char inxt3[7] = {1,2,3,0,1,2,3};
extern unsigned char idirt[4][3];

extern hash  hTab;
char ddb;

/* Find background tetras intersecting the boundary,
 and initialize distance at their vertices */
int iniredist_3d(Info info, pMesh mesh, pSol sol) {
  pTetra  pt;
  pPoint  p0,p1,p2,p3,pa;
  double  *solTmp,d, *circum, xc[3],r,norm;
  int     *bndy,i,j,nb,nc,i0,i1,i2,i3;
  char    proj;

  nb   = 0;
  bndy = (int*)calloc(mesh->ne+1,sizeof(int));
  assert(bndy);
  
  /* Store intersecting background triangles in list bndy */
  for (i=1; i<=mesh->ne; i++) {
    pt = &mesh->tetra[i];
    i0 = pt->v[0];
    i1 = pt->v[1];
    i2 = pt->v[2];
    i3 = pt->v[3];

    if ( (sol->val[i0] * sol->val[i1] <= 0.) || (sol->val[i0] * sol->val[i2] <= 0.) ||
         (sol->val[i0] * sol->val[i3] <= 0.) || (sol->val[i1] * sol->val[i2] <= 0.) ||
         (sol->val[i1] * sol->val[i3] <= 0.) || (sol->val[i2] * sol->val[i3] <= 0.)){
      nb++;
      bndy[nb] = i;
    }
  }

  bndy = (int*)realloc(bndy,(nb+1)*sizeof(int));
  printf("nb= %d\n",nb);

  /* Store temporary values */
  solTmp = (double*)calloc(mesh->np+1,sizeof(double));
  assert(solTmp);

  /* Check list bndy and compute distance at the vertices of its elements */
  for (i=1; i<=nb; i++) {
    fflush(stdout);
    pt = &mesh->tetra[bndy[i]];
    i0 = pt->v[0];
    i1 = pt->v[1];
    i2 = pt->v[2];
    i3 = pt->v[3];
    
    p0 = &mesh->point[i0];
    p1 = &mesh->point[i1];
    p2 = &mesh->point[i2];
    p3 = &mesh->point[i3];

    /* Distance from p0 to level 0 */
    d = distnv0_3d(mesh, sol, bndy[i], p0, &proj);

  	if ( p0->tag == 0 ) {
      solTmp[i0] = d;
      p0->tag    = proj;
    }
    else if ( d < fabs(solTmp[i0]) ) {
      solTmp[i0] = d;
      p0->tag    = proj;
    }

    /* Distance from p1 to level 0 */
    d = distnv0_3d(mesh, sol, bndy[i], p1, &proj);

    if ( p1->tag == 0 ) {
      solTmp[i1] = d;
      p1->tag    = proj;
  	}
    else if ( d < fabs(solTmp[i1]) ) {
      solTmp[i1] = d;
      p1->tag    = proj;
    }

    /* Distance from p2 to level 0 */
    d = distnv0_3d(mesh, sol, bndy[i], p2, &proj);

  	if ( p2->tag == 0 ) {
      solTmp[i2] = d;
      p2->tag    = proj;
    }
    else if ( d < fabs(solTmp[i2]) ) {
      solTmp[i2] = d;
      p2->tag    = proj;
    }
    
    /* Distance from p3 to level 0 */
  	d = distnv0_3d(mesh, sol, bndy[i], p3, &proj);

    if ( p3->tag == 0 ) {
      solTmp[i3] = d;
      p3->tag    = proj;
    }

    else if ( d < fabs(solTmp[i3]) ) {
      solTmp[i3] = d;
      p3->tag    = proj;
    }
  }
  
  /* Correction procedure for points whose tag is 2, i.e. min is not achieved through an orthogonal projection */
  if ( !info.fini ) {
    circum = (double*)calloc(4*nb+1,sizeof(double));
    fprintf(stdout,"     Building circumcircle table...");
    buildcircumredis_3d(mesh,sol,bndy,nb,circum);
    printf("Done\n");
    nc = 0;
    
    for (i=1; i<=mesh->np; i++) {
      pa = &mesh->point[i];
      if ( pa->tag < 1 )  continue;
      for (j=1; j<=nb; j++) {
        xc[0] = circum[4*(j-1)+1];
        xc[1] = circum[4*(j-1)+2];
        xc[2] = circum[4*(j-1)+3];
        r = circum[4*(j-1)+4];
        norm = (pa->c[0]-xc[0])*(pa->c[0]-xc[0]) + (pa->c[1]-xc[1])*(pa->c[1]-xc[1])\
            + (pa->c[2]-xc[2])*(pa->c[2]-xc[2]);

        if ( r == 0.0 ) continue;
        if ( (r<norm) && (2* sol->val[i] *(r + norm) < (norm -r)*(norm-r)) ) continue;

	    pt = &mesh->tetra[bndy[j]];
	    d  = distnv0_3d(mesh,sol,bndy[j],pa,&proj);
	    if (/* proj == 1 &&*/ d < fabs(solTmp[i]) ) {
	      solTmp[i] = d;
	      //break;
	    }
      }
      pa->tag = 1;
      nc++;
    }
    free(circum);
  }

  if ( nc )   fprintf(stdout,"     %d correction(s)\n",nc);

   for (i=1; i<=mesh->np; i++) {
    pa = &mesh->point[i];
    if ( !pa->tag )
	  sol->val[i] = sol->val[i] < 0.0 ? -sqrt(INIVAL_3d) : sqrt(INIVAL_3d);
    else if ( pa->tag )
	  sol->val[i] = sol->val[i] < 0.0 ? -sqrt(solTmp[i]) : sqrt(solTmp[i]);
  }

  free(solTmp);
  free(bndy);
  return(1);
}

/* Initialize a (unsigned) distance function to a domain already existing in mesh, defined
   by boundary triangles of reference contained in info.sref */
int inireftrias_3d(Info info,pMesh mesh, pSol sol){
  pTetra    pt,pt1;
  pPoint    p0,p1,p2,pn0,pn1,pn2;
  int       k,l,nb,n0,n1,n2,np0,np1,np2,mins,maxs,sum,ilist,iel,jel,lball,nc;
  int       *list,ball[LONMAX+2];
  double    dd;
  char      i,j,i0,i1,i2,j0,j1,j2,ip,jp,iface,proj;

  /* Reset point tags and flags */
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    p0->flag = 0;
    p0->tag = 0;
  }

  /* Initialization of values */
  for (k=1; k<=sol->np; k++)
    sol->val[k] = INIVAL_3d;

  sol->ne = mesh->ne;

  /* Create hashing with all triangles of the mesh sharing a reference in info.sref */
  nb = hashTriaRef(info,mesh);

  /* Travel all tets of the mesh to store initialization faces under the form 4*iel + iface */
  list = (int*)calloc(nb+1,sizeof(int));
  assert(list);

  ilist = 0;
  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];

    for(i=0; i<4; i++) {
      i0 = idirt[i][0];
      i1 = idirt[i][1];
      i2 = idirt[i][2];

      n0 = pt->v[i0];
      n1 = pt->v[i1];
      n2 = pt->v[i2];

      mins = D_MIN(n0,D_MIN(n1,n2));
      maxs = D_MAX(n0,D_MAX(n1,n2));
      sum  = n0+n1+n2;

      iel = getTria(mesh,mins,maxs,sum);
      if( !iel ) continue;

      list[ilist] = 4*k+i;
      ilist++;
      assert( ilist<nb+1 );
    }
  }

  /* Make sure that all triangles with the desired ref have been reached */
  assert( ilist == nb );

  /* Free hashing structure */
  delHash(mesh);

  /* Travel boundary faces : p->flag = last triangle with respect to which distance
     has been evaluated */
  for(k=0; k<nb; k++) {
    iel = list[k] / 4;
    iface = list[k] % 4;

    pt = &mesh->tetra[iel];
    i0 = idirt[iface][0];
    i1 = idirt[iface][1];
    i2 = idirt[iface][2];

    p0 = &mesh->point[pt->v[i0]];
    p1 = &mesh->point[pt->v[i1]];
    p2 = &mesh->point[pt->v[i2]];

    for(j=0; j<3; j++) {
      ip = idirt[iface][j];

      lball = boulet_3d(mesh,iel,ip,ball);
      assert(lball);

      for(l=0; l<lball; l++) {
        jel = ball[l] / 4;
        jp = ball[l] % 4;
        pt1 = &mesh->tetra[jel];

        j0 = inxt3[jp];
        j1 = inxt3[j0];
        j2 = inxt3[j1];

        np0 = pt1->v[j0];
        np1 = pt1->v[j1];
        np2 = pt1->v[j2];

        pn0 = &mesh->point[np0];
        pn1 = &mesh->point[np1];
        pn2 = &mesh->point[np2];

        if( pn0->flag != k ) {
          dd = distpt_3d(p0,p1,p2,pn0,&proj);
          if( dd < sol->val[np0] ) {
            sol->val[np0] = dd;
            pn0->tag = proj;
          }
          pn0->flag = k;
        }
        if( pn1->flag != k ) {
          dd = distpt_3d(p0,p1,p2,pn1,&proj);
          if( dd < sol->val[np1] ) {
            sol->val[np1] = dd;
            pn1->tag = proj;
          }
          pn1->flag = k;
        }
        if( pn2->flag != k ) {
          dd = distpt_3d(p0,p1,p2,pn2,&proj);
          if( dd < sol->val[np2] ) {
            sol->val[np2] = dd;
            pn0->tag = proj;
          }
          pn2->flag = k;
        }
      }
    }
  }

  /* set distance to all boundary points to 0 */
  for(k=0; k<nb; k++){
    iel = list[k] / 4;
    iface = list[k] % 4;
    pt = &mesh->tetra[iel];

    for(j=0; j<3; j++){
      np0 = pt->v[idirt[iface][j]];
      sol->val[np0] = 0.0;
    }
  }

  /* Correction procedure */
  nc = 0;

  for(k=1; k<=mesh->np; k++){
    pn0 = &mesh->point[k];
    if(pn0->tag < 2) continue;

    for(l=0; l<nb; l++){
      iel = list[l] / 4;
      iface = list[l] % 4;

      i0 = idirt[iface][0];
      i1 = idirt[iface][1];
      i2 = idirt[iface][2];

      p0 = &mesh->point[pt->v[i0]];
      p1 = &mesh->point[pt->v[i1]];
      p2 = &mesh->point[pt->v[i2]];

      dd = distpt_3d(p0,p1,p2,pn0,&proj);
      if(dd < sol->val[k]){
        nc++;
        sol->val[k] = dd;
        if(proj == 1){
          pn0->tag = 1;
          break;
        }
      }
    }
  }
  if ( nc )   fprintf(stdout,"     %d correction(s)\n",nc);

  /* reset point flags */
  for(k=1; k<=mesh->np; k++){
    p0 = &mesh->point[k];
    p0->flag = 0;
  }

  /* Put sign in the function, and take sqrt */
  for(k=1; k<=mesh->ne;k++){
    pt = &mesh->tetra[k];
	  for(j=0;j<4;j++){
        np0 = pt->v[j];
	    p0 = &mesh->point[np0];

        if( p0->flag == 1 ) continue;
	    p0->flag = 1;
        sol->val[np0] = sqrt(sol->val[np0]);
  	}
  }

  /* reset point flags and tags */
  for(k=1; k<=mesh->np; k++){
    p0 = &mesh->point[k];
    p0->flag = 0;
    p0->tag = 0;
  }

  return(1);
}

/* Initialize a l.s. function to a domain already existing in mesh, defined by isIntDom(info,pt->ref) = 1 */
int iniencdomain_3d(Info info,pMesh mesh, pSol sol){
  pTetra  pt,pt1;
  pPoint  p0,p1,p2,pn0,pn1,pn2;
  double  dd;
  int     *list,*adja,ball[LONMAX+2];
  int     nb,k,np,nc,iel,jel,lball,np0,np1,np2,l;
  char    i,j,ip,jp,i0,i1,i2,j0,j1,j2,iface,proj;

  nb = 0;

  /* reset point tags and flags */
  for(k=1; k<=mesh->np; k++){
    p0 = &mesh->point[k];
    p0->flag = 0;
    p0->tag = 0;
  }

  /* Compute nb = number of boundary faces to be considered */
  for(k=1; k<= mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !isIntDom(info,pt->ref) ) continue;

    for(i=0; i<4; i++){
      adja = &mesh->adja[4*(k-1)+1];
      iel = adja[i] / 4;

      /* Exclude outer boundary triangles from starting boundary */
      if( !iel ) continue;
      if ( !isIntDom(info,mesh->tetra[iel].ref) ) nb++;
    }
  }

  if ( !nb ) {
    printf("****** Error - No boundary face found \n");
    return(0);
  }

  for (k=1; k<=sol->np; k++)
    sol->val[k] = INIVAL_3d;

  sol->ne = mesh->ne;

  list = (int*)calloc(nb+1,sizeof(int));
  assert(list);

  nc = 0;
  /* Store boundary faces under the form 4*iel+iface */
  for(k=1; k<= mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !isIntDom(info,pt->ref) ) continue;

    for(i=0; i<4; i++){
      adja = &mesh->adja[4*(k-1)+1];
      iel = adja[i] / 4;

      /* Exclude outer boundary triangles from starting boundary */
      if( !iel ) continue;
      if ( !isIntDom(info,mesh->tetra[iel].ref) ) {
        list[nc] = 4*k+i;
        nc++;
      }

      /* Include outer boundary triangles as starting boundary */
      /*if(!iel || mesh->tetra[iel].ref != REFINT){
        list[nc] = 4*k+i;
        nc++;
      }*/
    }
  }

  /* Travel boundary faces : p->flag = last triangle with respect
      to which distance has been evaluated */
  for(k=0; k<nb; k++){
    iel = list[k] / 4;
    iface = list[k] % 4;

    pt = &mesh->tetra[iel];

    i0 = idirt[iface][0];
    i1 = idirt[iface][1];
    i2 = idirt[iface][2];

    p0 = &mesh->point[pt->v[i0]];
    p1 = &mesh->point[pt->v[i1]];
    p2 = &mesh->point[pt->v[i2]];

    for(j=0; j<3; j++){
      ip =idirt[iface][j];
      
      /* To Do: try to take several layers of balls, using an adapted function */
      lball = boulet_3d(mesh,iel,ip,ball);
      assert(lball);

      for(l=0; l<lball; l++){
        jel = ball[l] / 4;
        jp  = ball[l] % 4;
        pt1 = &mesh->tetra[jel];

        j0 = inxt3[jp];
        j1 = inxt3[j0];
        j2 = inxt3[j1];

        np0 = pt1->v[j0];
        np1 = pt1->v[j1];
        np2 = pt1->v[j2];

        pn0 = &mesh->point[pt1->v[j0]];
        pn1 = &mesh->point[pt1->v[j1]];
        pn2 = &mesh->point[pt1->v[j2]];

        if ( pn0->flag != k ) {
          dd = distpt_3d(p0,p1,p2,pn0,&proj);
          if ( dd < sol->val[np0] ) {
            sol->val[np0] = dd;
            pn0->tag = proj;
          }
          pn0->flag = k;
        }
        
        if ( pn1->flag != k ) {
          dd = distpt_3d(p0,p1,p2,pn1,&proj);
          if ( dd < sol->val[np1] ) {
            sol->val[np1] = dd;
            pn1->tag = proj;
          }
          pn1 ->flag = k;
        }

        if ( pn2->flag != k ) {
          dd = distpt_3d(p0,p1,p2,pn2,&proj);
          if ( dd < sol->val[np2] ) {
            sol->val[np2] = dd;
            pn2->tag = proj;
          }
          pn2->flag = k;
        }
      }
    }
  }

  /* Set distance to all boundary points to 0 */
  for (k=0; k<nb; k++) {
    iel = list[k] / 4;
    iface = list[k] % 4;
    pt = &mesh->tetra[iel];

    for (j=0; j<3; j++) {
      np0 = pt->v[idirt[iface][j]];
      sol->val[np0] = 0.0;
    }
  }

  /* Correction procedure */
  nc = 0;

  if ( !info.fini ) {
    for (k=1; k<=mesh->np; k++) {
      pn0 = &mesh->point[k];
      if ( pn0->tag < 2 ) continue;

      for (l=0; l<nb; l++) {
        iel = list[l] / 4;
        iface = list[l] % 4;
        pt = &mesh->tetra[iel];

        i0 = idirt[iface][0];
        i1 = idirt[iface][1];
        i2 = idirt[iface][2];

        p0 = &mesh->point[pt->v[i0]];
        p1 = &mesh->point[pt->v[i1]];
        p2 = &mesh->point[pt->v[i2]];

        dd = distpt_3d(p0,p1,p2,pn0,&proj);
        if ( dd < sol->val[k] ) {
          nc++;
          sol->val[k] = dd;
          if ( proj == 1 ) {
            pn0->tag = 1;
            break;
          }
        }
      }
    }
  }
  if ( nc )   fprintf(stdout,"     %d correction(s)\n",nc);

  /* reset point flags */
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    p0->flag = 0;
  }

  /* Put sign in the function, and take sqrt */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    for (j=0; j<4; j++) {
      np = pt->v[j];
      p0 = &mesh->point[np];

      if ( p0->flag ==1 ) continue;
      p0->flag = 1;

      if ( isIntDom(info,pt->ref) )
	      sol->val[np] = -sqrt(sol->val[np]);
      else
        sol->val[np] = sqrt(sol->val[np]);
  	}
  }

  /* reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;


  return(1);
}

/* Initialize exact unsigned distance function at vertices of tetras intersecting mesh2 */
int inidist_3d(Info info,pMesh mesh1,pMesh mesh2,pSol sol1,pBucket bucket) {
  pTetra   pt,pt1;
  pTria    pf;
  pPoint   p1,p2, p3, pa,pb, pc;
  double   cb[4],d, *circum;
  int      *adja,*list,base,iadr,ilist,i,j,k,ia,ib, ic,nc, iel,jel,ier,cur;
  char     tag;
  double   r, xc[3], norm;

  for (k=1; k<=sol1->np; k++)  sol1->val[k] = INIVAL_3d;

  /* memory alloc */
  list = (int*)calloc(mesh1->ne+1,sizeof(int));
  assert(list);

  nc = 0;
  for (k=1; k<=mesh2->nt; k++) {
    pf  = &mesh2->tria[k];
    p1  = &mesh2->point[pf->v[0]];
    p2  = &mesh2->point[pf->v[1]];
    p3  = &mesh2->point[pf->v[2]];

    iel = buckin_3d(mesh1,bucket,p1->c);
    iel = locelt_3d(mesh1,iel,p1->c,cb);
	  if ( !iel ) {
	    iel = buckin_3d(mesh1,bucket,p2->c);
	    iel = locelt_3d(mesh1,iel,p2->c,cb);
	  }

	  if(!iel){
	    iel = buckin_3d(mesh1,bucket,p3->c);
	    iel = locelt_3d(mesh1,iel,p3->c,cb);
	  }

   	assert(iel);

    ilist       = 1;
    list[ilist] = iel;
    base = ++mesh1->flag;
    pt   = &mesh1->tetra[iel];
    pt->flag = base;
    pt->tag  = 2;
    cur      = 1;
    do {
      iel  = list[cur];
      pt   = &mesh1->tetra[iel];
      iadr = 4*(iel-1) + 1;
      adja = &mesh1->adja[iadr];
      for (i=0; i<4; i++) {
        ia = inxt3[i];
        ib = inxt3[ia];
        ic = inxt3[ib];
        pa = &mesh1->point[pt->v[ia]];
        pb = &mesh1->point[pt->v[ib]];
        pc = &mesh1->point[pt->v[ic]];
        ier = intersec_3d(p1,p2,p3,pa,pb,pc);
        jel = adja[i] / 4;
        pt1 = &mesh1->tetra[jel];
        if ( ier ) {
          pt1->tag = 2;
          if ( jel && (pt1->flag < base) ) {
            ilist++;
            list[ilist] = jel;
            pt1->flag   = base;
          }
        }
        else if ( !pt1->tag )  pt1->tag = 1;
      }
      cur++;
    }
    while ( cur <= ilist );

    /* list analysis */
    for (i=1; i<=ilist; i++) {
      iel = list[i];
      pt  = &mesh1->tetra[iel];
      for (j=0; j<4; j++) {
        ia = pt->v[j];
        pa = &mesh1->point[ia];
        pa->tag = D_MAX(pa->tag,1);
        d  = distpt_3d(p1,p2,p3,pa,&tag);
        if ( d < sol1->val[ia] ) {
          sol1->val[ia] = d;
          pa->tag = tag;
        }
      }
    }
  }
  fprintf(stdout,"     distance\n");

  /* correction */
  if ( !info.fini ) {
    circum = (double*)calloc(4*(mesh2->nt)+1,sizeof(double));
    fprintf(stdout,"     Building circumcircle table...");
    buildcircum_3d(mesh2,circum);
    fprintf(stdout,"done. \n");

    nc = 0;
    for (k=1; k<=mesh1->np; k++) {
      pa = &mesh1->point[k];
      if ( pa->tag < 1 )  continue;
      for (i=1; i<=mesh2->nt; i++) {
	    xc[0] = circum[4*(i-1)+1];
	    xc[1] = circum[4*(i-1)+2];
	    xc[2] = circum[4*(i-1)+3];
        r = circum[4*(i-1)+4];
	    norm = (pa->c[0]-xc[0])*(pa->c[0]-xc[0]) + (pa->c[1]-xc[1])*(pa->c[1]-xc[1])\
		     + (pa->c[2]-xc[2])*(pa->c[2]-xc[2]);

	    if((r<norm)&&(2* sol1->val[k] *(r + norm) < (norm -r)*(norm-r))) continue;
        pf = &mesh2->tria[i];
        p1 = &mesh2->point[pf->v[0]];
        p2 = &mesh2->point[pf->v[1]];
        p3 = &mesh2->point[pf->v[2]];
        d  = distpt_3d(p1,p2,p3,pa,&tag);
        if ( d < sol1->val[k] ) {
          sol1->val[k] = d;
        }
      }
      pa->tag = 1;
      nc++;
    }
    free(circum);
  }
  if ( nc )  fprintf(stdout,"     %d correction(s)\n",nc);

  free(list);
  return(1);

}

/* Initialize unsigned distance function to the point cloud contained in mesh2 */
int inidistpcloud_3d(pMesh mesh1,pMesh mesh2,pSol sol1,pBucket bucket) {
  return(1);
}


/* sign the implicit function in each connected component */
int sgndist_3d(Info info,pMesh mesh,pMesh mesh2,pSol sol,pBucket bucket) {
  pTetra   pt,pt1;
  pTria    pf;
  pPoint   ppt,pi,p1,p2,p3;
  double   p[3],cb[4];
  int     *pile,*adja,*pilcc,i,j,k,cc,kk,iel,base,bfin,iadr,ipile,start,ip,nc;
  char     j1,j2,j3,sgn;

  /* Memory alloc */
  pile = (int*)calloc(mesh->ne+1,sizeof(int));
  assert(pile);
  
  /* Reset flags */
  for (k=1; k<=mesh->ne; k++)  mesh->tetra[k].flag = 0;
  for (k=1; k<=mesh->np; k++)  mesh->point[k].flag = 0;
  mesh->flag = 0;
  base = ++mesh->flag;
  
  ipile = 0;

  /* Travel the exterior points (default or supplied) */
  for (k=0; k<info.nexp; k++) {
    p[0] = info.exp[3*k];
    p[1] = info.exp[3*k+1];
    p[2] = info.exp[3*k+2];
    
    iel = buckin_3d(mesh,bucket,p);
    iel = locelt_3d(mesh,iel,p,cb);
    assert(iel);
    
    ipile++;
    pile[ipile] = iel;
    pt = &mesh->tetra[iel];
    pt->flag = base;
  }
  
  start    = 1;
  while ( ipile > 0 ) {
    /* search for all elements in the current connected component */
    do {
      k = pile[ipile];
      ipile--;

      pt   = &mesh->tetra[k];
      iadr = 4*(k-1) + 1;
      adja = &mesh->adja[iadr];
      for (i=0; i<4; i++) {
        if ( !adja[i] )  continue;
        kk  = adja[i] / 4;
        pt1 = &mesh->tetra[kk];
        if ( pt1->tag < 2 && pt1->flag < base ) {
          ipile++;
          pile[ipile] = kk;
          pt1->flag   = base;
        }
      }
    }
    while ( ipile > 0 );

    /* find next component */
    ipile = 0;
    for (k=start; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( pt->flag >= 1 )  continue;
      else if ( !pt->tag ) {
        base = ++mesh->flag;
        ipile = 1;
        pile[ipile] = k;
        pt->flag    = base;
        break;
      }
    }
  }

  /* analyze components */
  fprintf(stdout,"     %d connected component(s)\n",base);
  if ( base < 2 )  return(-1);

  /* store tetrahedra intersected */
  bfin  = mesh->flag;
  sgn   = 1;
  pilcc = (int*)calloc(bfin+2,sizeof(int));
  assert(pilcc);
  base     = ++mesh->flag;
  pilcc[1] = base;
  ipile    = 0;
  do {
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( pt->flag < 1 )  continue;
      for (cc=1; cc<=bfin; cc++) {
        if ( pilcc[cc] == base && pt->flag == cc ) {
          iadr = 4*(k-1) + 1;
          adja = &mesh->adja[iadr];
          for (i=0; i<4; i++) {
            if ( !adja[i] )  continue;
            kk  = adja[i] / 4;
            pt1 = &mesh->tetra[kk];
            if ( pt1->flag == 0 ) {
              ipile++;
              pile[ipile] = kk;
              pt1->flag   = -1;
            }
          }
        }
      }
    }
    if ( !ipile )  break;

    while ( ipile > 0 ) {
      k = pile[ipile];
      ipile--;
      pt = &mesh->tetra[k];
      iadr = 4*(k-1) + 1;
      adja = &mesh->adja[iadr];
      for (i=0; i<4; i++) {
        if ( !adja[i] )  continue;
        kk  = adja[i] / 4;
        pt1 = &mesh->tetra[kk];
        if ( pt1->flag == 0 ) {
          ipile++;
          pile[ipile] = kk;
          pt1->flag   = -1;
        }
        else if ( pt1->flag > 0 && pilcc[pt1->flag] != base )
          pilcc[pt1->flag] = base + 1;
      }
    }

    /* update sign */
    base = ++mesh->flag;
    sgn  = -sgn;
    for (i=1; i<=bfin; i++) {
      if ( pilcc[i] != base )  continue;
      for (k=1; k<=mesh->ne; k++) {
        pt = &mesh->tetra[k];
        if ( pt->flag == i  && pt->tag < 2 ) {
          for (j=0; j<4; j++) {
            ppt = &mesh->point[pt->v[j]];
            ppt->flag = 1;
            sol->val[pt->v[j]] = sgn * fabs(sol->val[pt->v[j]]);
          }
        }
      }
    }
  }
  while ( 1 );
  free(pilcc);

  /* set the first component */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( pt->flag == 1 ) {
      for (j=0; j<4; j++) {
        ppt = &mesh->point[pt->v[j]];
        ppt->flag = 1;
      }
    }
  }

  /* check for elements with uninitialized vertices */
  ipile = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !pt->tag )  continue;
    for (i=0; i<4; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( !ppt->flag ) {
        ipile++;
        pile[ipile] = k;
        break;
      }
    }
  }
  fprintf(stdout,"  %%%% %d elements with a vertex uninitialized\n",ipile);

  nc = 0;
  for (i=1; i<=ipile; i++) {
    k  = pile[i];
    pt = &mesh->tetra[k];

    /* on commence par rechercher un sommet correctement initialise */
    for (j=0; j<4; j++) {
      ip = pt->v[j];
      pi = &mesh->point[ip];
      if ( pi->flag )  break;
    }
    /* unable to correct sign */
    if ( j == 4 )  continue;

    /* Search for an already initialized vertex */
    j1  = inxt3[j+0];
    j2  = inxt3[j+1];
    j3  = inxt3[j+2];
    p1  = &mesh->point[pt->v[j1]];
    p2  = &mesh->point[pt->v[j2]];
    p3  = &mesh->point[pt->v[j3]];
    if ( !p1->flag || !p2->flag || !p3->flag ) {
      for (k=1; k<=mesh2->nt; k++) {
        pf  = &mesh2->tria[k];
        if ( !p1->flag && interSegTria(mesh2,pi,p1,pf) ) {
          sol->val[pt->v[j1]] = sol->val[ip] > 0. ? -fabs(sol->val[pt->v[j1]]) : fabs(sol->val[pt->v[j1]]);
          p1->flag = 1;
          nc++;
        }
        if ( !p2->flag && interSegTria(mesh2,pi,p2,pf) ) {
          sol->val[pt->v[j2]] = sol->val[ip] > 0. ? -fabs(sol->val[pt->v[j2]]) : fabs(sol->val[pt->v[j2]]);
          p2->flag = 1;
          nc++;
        }
        if ( !p3->flag && interSegTria(mesh2,pi,p3,pf) ) {
          sol->val[pt->v[j3]] = sol->val[ip] > 0. ? -fabs(sol->val[pt->v[j3]]) : fabs(sol->val[pt->v[j3]]);
          p3->flag = 1;
          nc++;
        }
        if ( p1->flag && p2->flag && p3->flag )  break;
      }
      if ( k > mesh2->nt ) {
        if ( !p1->flag ) {
          sol->val[pt->v[j1]] = sol->val[ip] > 0. ? fabs(sol->val[pt->v[j1]]) : -fabs(sol->val[pt->v[j1]]);
          p1->flag = 1;
          nc++;
        }
        if ( !p2->flag ) {
          sol->val[pt->v[j2]] = sol->val[ip] > 0. ? fabs(sol->val[pt->v[j2]]) : -fabs(sol->val[pt->v[j2]]);
          p2->flag = 1;
          nc++;
        }
        if ( !p3->flag ) {
          sol->val[pt->v[j3]] = sol->val[ip] > 0. ? fabs(sol->val[pt->v[j3]]) : -fabs(sol->val[pt->v[j3]]);
          p3->flag = 1;
          nc++;
        }
      }
    }
  }
  fprintf(stdout,"  %%%% %d corrected vertices\n",nc);
  free(pile);

  /* take sqrt of distance */
  for (k=1; k<=mesh->np; k++) {
    mesh->point[k].flag = 0;
    sol->val[k] = sol->val[k] < 0.0 ? -sqrt(fabs(sol->val[k])) : sqrt(sol->val[k]);
  }

  return(1);
}


/* compute gradients of shape functions P1: 16 values (dx,dy,dz,cte) */
/*static void gradelt_3d(int istart,int istop,int ipth,Param *par) {
  pMesh    mesh;
  pTetra   pt;
  pPoint   p0,p1,p2,p3;
  double  *grad,m[9],mi[9];
  int      k,iadr,ier;

  mesh = par->mesh;
  grad = par->grad;

  for(k=istart; k<=istop; k++) {
    pt = &mesh->tetra[k];
    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];
    p3 = &mesh->point[pt->v[3]];

    // lambda_0
    m[0] = p1->c[0] - p0->c[0];
    m[1] = p1->c[1] - p0->c[1];
    m[2] = p1->c[2] - p0->c[2];
    m[3] = p2->c[0] - p0->c[0];
    m[4] = p2->c[1] - p0->c[1];
    m[5] = p2->c[2] - p0->c[2];
    m[6] = p3->c[0] - p0->c[0];
    m[7] = p3->c[1] - p0->c[1];
    m[8] = p3->c[2] - p0->c[2];

    ier = invmatg(m,mi);
    iadr = 16 * (k-1) + 1;
    grad[iadr + 0] = -mi[0] - mi[1] - mi[2];
    grad[iadr + 1] = -mi[3] - mi[4] - mi[5];
    grad[iadr + 2] = -mi[6] - mi[7] - mi[8];
    grad[iadr + 3] = 1.0 - grad[iadr+0]*p0->c[0] - grad[iadr+1]*p0->c[1] - grad[iadr+2]*p0->c[2];

    // lambda_1
    m[0] = p0->c[0] - p1->c[0];
    m[1] = p0->c[1] - p1->c[1];
    m[2] = p0->c[2] - p1->c[2];
    m[3] = p2->c[0] - p1->c[0];
    m[4] = p2->c[1] - p1->c[1];
    m[5] = p2->c[2] - p1->c[2];
    m[6] = p3->c[0] - p1->c[0];
    m[7] = p3->c[1] - p1->c[1];
    m[8] = p3->c[2] - p1->c[2];

    ier = invmatg(m,mi);
    grad[iadr + 4] = -mi[0] - mi[1] - mi[2];
    grad[iadr + 5] = -mi[3] - mi[4] - mi[5];
    grad[iadr + 6] = -mi[6] - mi[7] - mi[8];
    grad[iadr + 7] = 1.0 - grad[iadr+4]*p1->c[0] - grad[iadr+5]*p1->c[1] - grad[iadr+6]*p1->c[2];

    // lambda_2
    m[0] = p0->c[0] - p2->c[0];
    m[1] = p0->c[1] - p2->c[1];
    m[2] = p0->c[2] - p2->c[2];
    m[3] = p1->c[0] - p2->c[0];
    m[4] = p1->c[1] - p2->c[1];
    m[5] = p1->c[2] - p2->c[2];
    m[6] = p3->c[0] - p2->c[0];
    m[7] = p3->c[1] - p2->c[1];
    m[8] = p3->c[2] - p2->c[2];

    ier = invmatg(m,mi);
    grad[iadr + 8] = -mi[0] - mi[1] - mi[2];
    grad[iadr + 9] = -mi[3] - mi[4] - mi[5];
    grad[iadr +10] = -mi[6] - mi[7] - mi[8];
    grad[iadr +11] = 1.0 - grad[iadr+8]*p2->c[0] - grad[iadr+9]*p2->c[1] - grad[iadr+10]*p2->c[2];

    // lambda_3
    m[0] = p0->c[0] - p3->c[0];
    m[1] = p0->c[1] - p3->c[1];
    m[2] = p0->c[2] - p3->c[2];
    m[3] = p1->c[0] - p3->c[0];
    m[4] = p1->c[1] - p3->c[1];
    m[5] = p1->c[2] - p3->c[2];
    m[6] = p2->c[0] - p3->c[0];
    m[7] = p2->c[1] - p3->c[1];
    m[8] = p2->c[2] - p3->c[2];

    ier = invmatg(m,mi);
    grad[iadr +12] = -mi[0] - mi[1] - mi[2];
    grad[iadr +13] = -mi[3] - mi[4] - mi[5];
    grad[iadr +14] = -mi[6] - mi[7] - mi[8];
    grad[iadr +15] = 1.0 - grad[iadr+12]*p3->c[0] - grad[iadr+13]*p3->c[1] - grad[iadr+14]*p3->c[2];
  }
}
*/

/* Gradient of basis functions */
static void gradelt_3d(int istart,int istop,int ipth,Param *par) {
 pMesh    mesh;
 pTetra   pt;
 pPoint   p0,p1,p2,p3;
 double  *grad,m[9],mi[9],lambda[16];
 int      k,iadr,ier,i0,i1,i2,i3;

 mesh = par->mesh;
 grad = par->grad;

 for(k=istart; k<=istop; k++) {
   pt = &mesh->tetra[k];
   i0 = pt->v[0];
   i1 = pt->v[1];
   i2 = pt->v[2];
   i3 = pt->v[3];

   p0 = &mesh->point[i0];
   p1 = &mesh->point[i1];
   p2 = &mesh->point[i2];
   p3 = &mesh->point[i3];

   /* lambda_0 */
   m[0] = p1->c[0] - p0->c[0];
   m[1] = p1->c[1] - p0->c[1];
   m[2] = p1->c[2] - p0->c[2];
   m[3] = p2->c[0] - p0->c[0];
   m[4] = p2->c[1] - p0->c[1];
   m[5] = p2->c[2] - p0->c[2];
   m[6] = p3->c[0] - p0->c[0];
   m[7] = p3->c[1] - p0->c[1];
   m[8] = p3->c[2] - p0->c[2];

   ier = invmatg(m,mi);
   lambda[0] = -mi[0] - mi[1] - mi[2];
   lambda[1] = -mi[3] - mi[4] - mi[5];
   lambda[2] = -mi[6] - mi[7] - mi[8];
   lambda[3] = 1.0 - lambda[0]*p0->c[0] - lambda[1]*p0->c[1] - lambda[2]*p0->c[2];

   /* lambda_1 */
   m[0] = p0->c[0] - p1->c[0];
   m[1] = p0->c[1] - p1->c[1];
   m[2] = p0->c[2] - p1->c[2];
   m[3] = p2->c[0] - p1->c[0];
   m[4] = p2->c[1] - p1->c[1];
   m[5] = p2->c[2] - p1->c[2];
   m[6] = p3->c[0] - p1->c[0];
   m[7] = p3->c[1] - p1->c[1];
   m[8] = p3->c[2] - p1->c[2];

   ier = invmatg(m,mi);
   lambda[4] = -mi[0] - mi[1] - mi[2];
   lambda[5] = -mi[3] - mi[4] - mi[5];
   lambda[6] = -mi[6] - mi[7] - mi[8];
   lambda[7] = 1.0 - lambda[4]*p1->c[0] - lambda[5]*p1->c[1] - lambda[6]*p1->c[2];

   /* lambda_2 */
   m[0] = p0->c[0] - p2->c[0];
   m[1] = p0->c[1] - p2->c[1];
   m[2] = p0->c[2] - p2->c[2];
   m[3] = p1->c[0] - p2->c[0];
   m[4] = p1->c[1] - p2->c[1];
   m[5] = p1->c[2] - p2->c[2];
   m[6] = p3->c[0] - p2->c[0];
   m[7] = p3->c[1] - p2->c[1];
   m[8] = p3->c[2] - p2->c[2];

   ier = invmatg(m,mi);
   lambda[8] = -mi[0] - mi[1] - mi[2];
   lambda[9] = -mi[3] - mi[4] - mi[5];
   lambda[10] = -mi[6] - mi[7] - mi[8];
   lambda[11] = 1.0 - lambda[8]*p2->c[0] - lambda[9]*p2->c[1] - lambda[10]*p2->c[2];

   /* lambda_3 */
   m[0] = p0->c[0] - p3->c[0];
   m[1] = p0->c[1] - p3->c[1];
   m[2] = p0->c[2] - p3->c[2];
   m[3] = p1->c[0] - p3->c[0];
   m[4] = p1->c[1] - p3->c[1];
   m[5] = p1->c[2] - p3->c[2];
   m[6] = p2->c[0] - p3->c[0];
   m[7] = p2->c[1] - p3->c[1];
   m[8] = p2->c[2] - p3->c[2];

   ier = invmatg(m,mi);
   lambda[12] = -mi[0] - mi[1] - mi[2];
   lambda[13] = -mi[3] - mi[4] - mi[5];
   lambda[14] = -mi[6] - mi[7] - mi[8];
   lambda[15] = 1.0 - lambda[12]*p3->c[0] - lambda[13]*p3->c[1] - lambda[14]*p3->c[2];


   grad[3*(k-1)+1] = par->dtmp[i0]*lambda[0] + par->dtmp[i1]*lambda[4] + par->dtmp[i2]*lambda[8] + par->dtmp[i3]*lambda[12];
   grad[3*(k-1)+2] = par->dtmp[i0]*lambda[1] + par->dtmp[i1]*lambda[5] + par->dtmp[i2]*lambda[9] + par->dtmp[i3]*lambda[13];
   grad[3*(k-1)+3] = par->dtmp[i0]*lambda[2] + par->dtmp[i1]*lambda[6] + par->dtmp[i2]*lambda[10] + par->dtmp[i3]*lambda[14];
 }

 //return(1);
}

/* Calculate an active value at point ip from the sole value at point ip0 */
double actival1pt_3d(pMesh mesh,pSol sol,int ip0,int ip) {
  pPoint    p0,p1;
  double    d0,dist,ll;
  
  p0   = &mesh->point[ip0];
  p1  = &mesh->point[ip];
  d0 = sol->val[ip0];
  
  ll = (p1->c[0]-p0->c[0])*(p1->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])*(p1->c[1]-p0->c[1])
       + (p1->c[2]-p0->c[2])*(p1->c[2]-p0->c[2]);
  ll = sqrt(ll);
  dist = ( d0 < 0.0 ) ? d0 - ll : d0 + ll;

  return(fabs(dist));
}

/* Calculate an active value at point ip from the value v0, v1 at ip0 and ip1 */
double actival2pt_3d(pMesh mesh,pSol sol,int ip0,int ip1,int ip) {
  pPoint     p0,p1,p2;
  double     d0,d1,d2,ps0,ps1,dist,rmin,rmax,n[3],m[9],mi[9],Gr[3][3],a[3],r[2],u[3],g[3];
  int        ier,nr;
  
  p0 = &mesh->point[ip0];
  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip];
  
  d0 = sol->val[ip0];
  d1 = sol->val[ip1];
  d2 = sol->val[ip];
  
  dist = INIVAL_3d;
  
  /* Gr[i]j] is the ith coordinate of the barycentric coordinates at point j inside triangle p0p1p2 */
  m[0] = p1->c[0] - p0->c[0];   m[1] = p1->c[1] - p0->c[1];  m[2] = p1->c[2] - p0->c[2];
  m[3] = p2->c[0] - p0->c[0];   m[4] = p2->c[1] - p0->c[1];  m[5] = p2->c[2] - p0->c[2];
  m[6] = m[1]*m[5]-m[2]*m[4];
  m[7] = m[2]*m[3]-m[0]*m[5];
  m[8] = m[0]*m[4]-m[1]*m[3];
  
  ier = invmatg(m,mi);
  
  if ( !ier ) {
    d0 = actival1pt_3d(mesh,sol,ip0,ip);
    d1 = actival1pt_3d(mesh,sol,ip1,ip);
    return ( D_MIN(d0,d1) );
  }
  
  Gr[0][0] = -mi[0]-mi[1]  ;  Gr[0][1] = mi[0]  ;  Gr[0][2] = mi[1];
  Gr[1][0] = -mi[3]-mi[4]  ;  Gr[1][1] = mi[3]  ;  Gr[1][2] = mi[4];
  Gr[2][0] = -mi[6]-mi[7]  ;  Gr[2][1] = mi[6]  ;  Gr[2][2] = mi[7];
  
  /* Coefficients of the quadratic equation */
  u[0] = d0*Gr[0][0] + d1*Gr[0][1];
  u[1] = d0*Gr[1][0] + d1*Gr[1][1];
  u[2] = d0*Gr[2][0] + d1*Gr[2][1];
  
  a[2] = Gr[0][2]*Gr[0][2] + Gr[1][2]*Gr[1][2] + Gr[2][2]*Gr[2][2];
  a[1] = Gr[0][2]*u[0] + Gr[1][2]*u[1] + Gr[2][2]*u[2];
  a[1] = 2.0*a[1];
  a[0] = u[0]*u[0] + u[1]*u[1] + u[2]*u[2] - 1.0;
  
  nr = eqquad(a,r);

  if ( nr == 1 ) {
    if ( ( d2 >= 0.0 && r[0] >= d0 && r[0] >= d1 ) || ( d2 <= 0.0 && r[0] <= d0 && r[0] <= d1 ) ) {
      /* Check the direction of the gradient */
      g[0] = (d0-r[0])*Gr[0][0] + (d1-r[0])*Gr[0][1];
      g[1] = (d0-r[0])*Gr[1][0] + (d1-r[0])*Gr[1][1];
      g[2] = (d0-r[0])*Gr[2][0] + (d1-r[0])*Gr[2][1];
      ps0  = g[0]*Gr[0][0] + g[1]*Gr[1][0] + g[2]*Gr[2][0];
      ps1  = g[0]*Gr[0][1] + g[1]*Gr[1][1] + g[2]*Gr[2][1];
      
      if ( (d2 >= 0.0 && ps0 < EPS1 && ps1 < EPS1) || (d2 <= 0.0 && ps0 > -EPS1 && ps1 > -EPS1) ) dist = fabs(r[0]);
    }
  }
  else if ( nr == 2 ) {
    /* Sort the roots */
    if ( d2 >= 0.0 ) {
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
    if ( (d2 >= 0.0 && rmin >= d0 && rmin >= d1) || (d2 <= 0.0 && rmin <= d0 && rmin <= d1) ) {
      g[0] = (d0-rmin)*Gr[0][0] + (d1-rmin)*Gr[0][1];
      g[1] = (d0-rmin)*Gr[1][0] + (d1-rmin)*Gr[1][1];
      g[2] = (d0-rmin)*Gr[2][0] + (d1-rmin)*Gr[2][1];
      ps0  = g[0]*Gr[0][0] + g[1]*Gr[1][0] + g[2]*Gr[2][0];
      ps1  = g[0]*Gr[0][1] + g[1]*Gr[1][1] + g[2]*Gr[2][1];
      
      if ( ( d2 >= 0.0 && ps0 < EPS1 && ps1 < EPS1 ) || ( d2 <= 0.0 && ps0 >- EPS1 && ps1 >- EPS1 )) dist = fabs(rmin);
      
    }
    else if ( (d2 >= 0.0 && rmax >= d0 && rmax >= d1) || (d2 <= 0.0 && rmax <= d0 && rmax <= d1) ) {
      g[0] = (d0-rmax)*Gr[0][0] + (d1-rmax)*Gr[0][1];
      g[1] = (d0-rmax)*Gr[1][0] + (d1-rmax)*Gr[1][1];
      g[2] = (d0-rmax)*Gr[2][0] + (d1-rmax)*Gr[2][1];
      ps0  = g[0]*Gr[0][0] + g[1]*Gr[1][0] + g[2]*Gr[2][0];
      ps1  = g[0]*Gr[0][1] + g[1]*Gr[1][1] + g[2]*Gr[2][1];
      
      if ( ( d2 >= 0.0 && ps0 < EPS1 && ps1 < EPS1 ) || ( d2 <= 0.0 && ps0 >- EPS1 && ps1 >- EPS1 )) dist = fabs(rmax);
    }
  }
  
  /* If no other value has been assigned to dist, calculate a trial value based on both triangle edges */
  if ( fabs( dist - INIVAL_3d ) < EPS2 ) {
    d0 = actival1pt_3d(mesh,sol,ip0,ip);
    d1 = actival1pt_3d(mesh,sol,ip1,ip);
    
    dist = D_MIN(d0,d1);
  }
  
  return(dist);
}

/* Calculate a (positive) active value at vertex i in tetra k based on the values in the other three vertices */
double actival_3d(pMesh mesh,pSol sol,int k,int i) {
  pTetra            pt;
  pPoint            p0,p1,p2,p3;
  double            dist,d0,d1,d2,d3,m[9],mi[9],Gr[3][4],u[3],a[3],r[2],g[3],ps0,ps1,ps2,rmin,rmax;
  int               ip0,ip1,ip2,ip3,nr;
  char              i0,i1,i2,i3,ier;
  
  i3 = i;
  i0 = inxt3[i3];
  i1 = inxt3[i0];
  i2 = inxt3[i1];
  
  pt = &mesh->tetra[k];
  ip0 = pt->v[i0];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];
  ip3 = pt->v[i3];
  
  p0 = &mesh->point[ip0];
  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];
  p3 = &mesh->point[ip3];
  
  d0 = sol->val[ip0];
  d1 = sol->val[ip1];
  d2 = sol->val[ip2];
  d3 = sol->val[ip3];
  
  dist = INIVAL_3d;
  
  /* Case where value at ip0 is not accepted */
  if ( p0->tag != 1 ) {
    if ( p1->tag != 1 )
      return( actival1pt_3d(mesh,sol,ip2,ip3) );
    else if ( p2->tag != 1 )
      return( actival1pt_3d(mesh,sol,ip1,ip3) );
    else
      return( actival2pt_3d(mesh,sol,ip1,ip2,ip3) );
  }
  /* Case where ip0 is accepted, but not ip1 */
  else if ( p1->tag != 1 ) {
    if ( p2->tag != 1 )
      return( actival1pt_3d(mesh,sol,ip0,ip3) );
    else
      return( actival2pt_3d(mesh,sol,ip0,ip2,ip3) );
  }
  else if ( p2->tag != 1 )
    return( actival2pt_3d(mesh,sol,ip0,ip1,ip3) );

  /* At this point, ip0, ip1 and ip2 are accepted */
  /* Calculation of the barycentric coordinates */
  m[0] = p1->c[0] - p0->c[0]    ;    m[1] = p1->c[1] - p0->c[1]     ;     m[2] = p1->c[2] - p0->c[2];
  m[3] = p2->c[0] - p0->c[0]    ;    m[4] = p2->c[1] - p0->c[1]     ;     m[5] = p2->c[2] - p0->c[2];
  m[6] = p3->c[0] - p0->c[0]    ;    m[7] = p3->c[1] - p0->c[1]     ;     m[8] = p3->c[2] - p0->c[2];

  ier = invmatg(m,mi);
  if (!ier ) {
    d0 = actival2pt_3d(mesh,sol,ip1,ip2,ip3);
    d1 = actival2pt_3d(mesh,sol,ip0,ip2,ip3);
    d2 = actival2pt_3d(mesh,sol,ip0,ip1,ip3);
    
    dist = D_MIN(D_MIN(d0,d1),d2);
    return(dist);
  }

  Gr[0][0] = -mi[0]-mi[1]-mi[2] ; Gr[0][1] = mi[0] ; Gr[0][2] = mi[1] ; Gr[0][3] = mi[2];
  Gr[1][0] = -mi[3]-mi[4]-mi[5] ; Gr[1][1] = mi[3] ; Gr[1][2] = mi[4] ; Gr[1][3] = mi[5];
  Gr[2][0] = -mi[6]-mi[7]-mi[8] ; Gr[2][1] = mi[6] ; Gr[2][2] = mi[7] ; Gr[2][3] = mi[8];
  
  /* Calculate possible active values by solving quadratic equation */
  u[0] = d0*Gr[0][0] + d1*Gr[0][1] + d2*Gr[0][2];
  u[1] = d0*Gr[1][0] + d1*Gr[1][1] + d2*Gr[1][2];
  u[2] = d0*Gr[2][0] + d1*Gr[2][1] + d2*Gr[2][2];
  
  a[2] = Gr[0][3]*Gr[0][3] + Gr[1][3]*Gr[1][3] + Gr[2][3]*Gr[2][3];
  a[1] = Gr[0][3]*u[0] + Gr[1][3]*u[1] + Gr[2][3]*u[2];
  a[1] = 2.0*a[1];
  a[0] = u[0]*u[0] + u[1]*u[1] + u[2]*u[2] - 1.0;
  
  nr = eqquad(a,r);
  
  if ( nr == 1 ) {
    if ( ( d3 >= 0.0 && r[0] >= d0 && r[0] >= d1 && r[0] >= d2 ) || ( d3 <= 0.0 && r[0] <= d0 && r[0] <= d1 && r[0] <= d2 ) ) {
      /* Check the direction of the gradient */
      g[0] = (d0-r[0])*Gr[0][0] + (d1-r[0])*Gr[0][1] + (d2-r[0])*Gr[0][2];
      g[1] = (d0-r[0])*Gr[1][0] + (d1-r[0])*Gr[1][1] + (d2-r[0])*Gr[1][2];
      g[2] = (d0-r[0])*Gr[2][0] + (d1-r[0])*Gr[2][1] + (d2-r[0])*Gr[2][2];
      ps0  = g[0]*Gr[0][0] + g[1]*Gr[1][0] + g[2]*Gr[2][0];
      ps1  = g[0]*Gr[0][1] + g[1]*Gr[1][1] + g[2]*Gr[2][1];
      ps2  = g[0]*Gr[0][2] + g[1]*Gr[1][2] + g[2]*Gr[2][2];
      
      if ( (d3 >= 0.0 && ps0 < EPS1 && ps1 < EPS1 && ps2 < EPS1) || (d3 <= 0.0 && ps0 > -EPS1 && ps1 > -EPS1 && ps2 > -EPS1) )
        dist = fabs(r[0]);
    }
  }
  else if ( nr == 2 ) {
    /* Sort the roots */
    if ( d2 >= 0.0 ) {
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
    if ( ( d3 >= 0.0 && rmin >= d0 && rmin >= d1 && rmin >= d2 ) || ( d3 <= 0.0 && rmin <= d0 && rmin <= d1 && rmin <= d2 ) ) {
      /* Check the direction of the gradient */
      g[0] = (d0-rmin)*Gr[0][0] + (d1-rmin)*Gr[0][1] + (d2-rmin)*Gr[0][2];
      g[1] = (d0-rmin)*Gr[1][0] + (d1-rmin)*Gr[1][1] + (d2-rmin)*Gr[1][2];
      g[2] = (d0-rmin)*Gr[2][0] + (d1-rmin)*Gr[2][1] + (d2-rmin)*Gr[2][2];
      ps0  = g[0]*Gr[0][0] + g[1]*Gr[1][0] + g[2]*Gr[2][0];
      ps1  = g[0]*Gr[0][1] + g[1]*Gr[1][1] + g[2]*Gr[2][1];
      ps2  = g[0]*Gr[0][2] + g[1]*Gr[1][2] + g[2]*Gr[2][2];
      
      if ( (d3 >= 0.0 && ps0 < EPS1 && ps1 < EPS1 && ps2 < EPS1) || (d3 <= 0.0 && ps0 > -EPS1 && ps1 > -EPS1 && ps2 > -EPS1) )
        dist = fabs(rmin);
    }
    else if ( ( d3 >= 0.0 && rmax >= d0 && rmax >= d1 && rmax >= d2 ) || ( d3 <= 0.0 && rmax <= d0 && rmax <= d1 && rmax <= d2 ) ) {
      /* Check the direction of the gradient */
      g[0] = (d0-rmax)*Gr[0][0] + (d1-rmax)*Gr[0][1] + (d2-rmax)*Gr[0][2];
      g[1] = (d0-rmax)*Gr[1][0] + (d1-rmax)*Gr[1][1] + (d2-rmax)*Gr[1][2];
      g[2] = (d0-rmax)*Gr[2][0] + (d1-rmax)*Gr[2][1] + (d2-rmax)*Gr[2][2];
      ps0  = g[0]*Gr[0][0] + g[1]*Gr[1][0] + g[2]*Gr[2][0];
      ps1  = g[0]*Gr[0][1] + g[1]*Gr[1][1] + g[2]*Gr[2][1];
      ps2  = g[0]*Gr[0][2] + g[1]*Gr[1][2] + g[2]*Gr[2][2];
      
      if ( (d3 >= 0.0 && ps0 < EPS1 && ps1 < EPS1 && ps2 < EPS1) || (d3 <= 0.0 && ps0 > -EPS1 && ps1 > -EPS1 && ps2 > -EPS1) )
        dist = fabs(rmax);
    }
  }
  
  /* If no other value has been assigned, calculate a default value */
  if ( fabs( dist - INIVAL_3d ) < EPS2 ) {
    d0 = actival2pt_3d(mesh,sol,ip1,ip2,ip3);
    d1 = actival2pt_3d(mesh,sol,ip0,ip2,ip3);
    d2 = actival2pt_3d(mesh,sol,ip0,ip1,ip3);
    
    dist = D_MIN(D_MIN(d0,d1),d2);
  }
  
  return(dist);
}

/* Propagation of the signed distance function by the Fast Marching Method */
int ppgdistfmm_3d(pMesh mesh,pSol sol) {
  Queue            q;
  pQueue           pq;
  pTetra           pt,pt1;
  pPoint           p0,p1,p2,p3;
  double           dist;
  int              nacc,k,l,iel,ip,ip1,ip2,ip3,ilist,*list;
  char             i,j,jj,j1,j2,j3;
  
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
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    for (i=0; i<4; i++) {
      ip = pt->v[i];
      p0 = &mesh->point[ip];
      
      if ( p0->tag != 1 ) continue;
      
      ilist = boulet_3d(mesh,k,i,list);
      
      for (l=0; l<ilist; l++) {
        iel = list[l] / 4;
        pt1 = &mesh->tetra[iel];
        
        j  = list[l] % 4;
        
        /* Calculate value at the three other points and put it in the queue */
        for (jj=0; jj<3; jj++) {
          j1   = inxt3[j];
          ip1  = pt1->v[j1];
          p1   = &mesh->point[ip1];
          
          if ( p1->tag == 1 ) continue;
          dist = actival_3d(mesh,sol,iel,j1);
            
          if ( p1->tag == 0 ) {
            insertAnod(pq,ip1,dist);
            p1->tag = 2;
          }
          else
            upAnod(pq,ip1,dist);
        }
      }
    }
  }
  
  /* Main loop: pop the smallest active node; it becomes accepted and the neighboring values become active */
  while ( pq->siz ) {
    ip = popAnod(pq,&dist);
    if ( !ip ) {
      printf("Problem in popping in Fast Marching Method. Abort.\n");
      exit(0);
    }
    
    p0 = &mesh->point[ip];
    p0->tag = 1;
    sol->val[ip] = sol->val[ip] > 0.0 ? dist : -dist;
    
    /* Travel the ball of p0 to update the set of active nodes */
    k   = p0->s;
    pt  = &mesh->tetra[k];
    for (i=0; i<4; i++)
      if ( pt->v[i] == ip ) break;
    assert( i < 4 );
    
    ilist = boulet_3d(mesh,k,i,list);
    for (l=0; l<ilist; l++) {
      iel = list[l] / 4;
      j   = list[l] % 4;
      pt  = &mesh->tetra[iel];
      
      for(jj=0; jj<3; jj++) {
        j    = inxt3[j];
        ip1  = pt->v[j];
        p1   = &mesh->point[ip1];
        
        /* Either insert or update active value if the point is not already accepted */
        if ( p1->tag == 1 ) continue;
        
        dist = actival_3d(mesh,sol,iel,j);
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

/* compute new distance at points in // */
static void tmpdist_3d(int istart,int istop,int ipth,Param *par) {
  pMesh     mesh;
  pSol      sol;
  pTetra    pt,pt1;
  pPoint    p0,p1,ppt[4];
  double   *grad,*dtmp,dx,dy,dz,ddx,ddy,ddz,dd,p[3],cb[4],drec,ldt,px,py,pz,dp;
  int       j,k,base,i0,iel,ip;

  mesh = par->mesh;
  sol  = par->sol;
  grad = par->grad;
  dtmp = par->dtmp;

  base = mesh->flag;

  for (k=istart; k<=istop; k++) {
    pt = &mesh->tetra[k];

    /* pt->flag = -1 when gradient is close to 1 in pt */
    if( pt->flag == -1 ) continue;

    dx = grad[3*(k-1)+1];
    dy = grad[3*(k-1)+2];
    dz = grad[3*(k-1)+3];
    dd = dx*dx + dy*dy + dz*dz;

    /* Don't move when gradient is close to 1 in pt */
    if( fabs(dd - 1.) <0.001 ) {
        pt->flag =-1;
        continue;
    }

    /* Vector to track characteristic line */
    if ( dd < EPS1 )  continue;
    dd = par->dt / sqrt(dd);
    dx *= dd;
    dy *= dd;
    dz *= dd;

    for (j=0; j<4; j++) {
      ip = pt->v[j];
      p0 = &mesh->point[ip];

      if ( p0->tag > 0 )  continue;
      else if (( fabs(sol->val[ip]) < 0.7 * par->dtfin )&&( p0->flag && p0->flag < mesh->mark-2 ))  continue;
      else if ( fabs(sol->val[ip]) < 0.4 * par->dtfin )  continue;

      if ( fabs(sol->val[ip]) < EPS1 )
        continue;
      else if ( sol->val[ip] > 0.0 ) {
        p[0] = p0->c[0] - dx;
        p[1] = p0->c[1] - dy;
        p[2] = p0->c[2] - dz;
      }
      else {
        p[0] = p0->c[0] + dx;
        p[1] = p0->c[1] + dy;
        p[2] = p0->c[2] + dz;
      }
      if ( p[0] <= EPS1 || p[0] >= 1.0-EPS1 )  continue;
      if ( p[1] <= EPS1 || p[1] >= 1.0-EPS1 )  continue;
      if ( p[2] <= EPS1 || p[2] >= 1.0-EPS1 )  continue;

      /* find enclosing tetra, k is guessed */
      iel = nxtelt_3d(mesh,k,p,cb);
      if ( iel < 1 ) {
        if( abs(iel) > 0 ) {
          iel = abs(iel);

          /* P1 Lagrange interpolation */
          pt1  = &mesh->tetra[iel];
          i0   = pt1->v[0];
          p1   = &mesh->point[i0];

          ddx = grad[3*(iel-1)+1];
          ddy = grad[3*(iel-1)+2];
          ddz = grad[3*(iel-1)+3];

          /* Coordinates of the exit point */
          ppt[0] = &mesh->point[pt1->v[0]];
          ppt[1] = &mesh->point[pt1->v[1]];
          ppt[2] = &mesh->point[pt1->v[2]];
          ppt[3] = &mesh->point[pt1->v[3]];

          px = ppt[0]->c[0]*cb[0] + ppt[1]->c[0]*cb[1] + ppt[2]->c[0]*cb[2] + ppt[3]->c[0]*cb[3];
          py = ppt[0]->c[1]*cb[0] + ppt[1]->c[1]*cb[1] + ppt[2]->c[1]*cb[2] + ppt[3]->c[1]*cb[3];
          pz = ppt[0]->c[2]*cb[0] + ppt[1]->c[2]*cb[1] + ppt[2]->c[2]*cb[2] + ppt[3]->c[2]*cb[3];

          drec = sol->val[i0] + ddx*(px-p1->c[0]) + ddy*(py-p1->c[1]) + ddz*(pz-p1->c[2]);

          /* Local time step, based on the physical distance between the exit and initial points */
          ldt = (px-p0->c[0])*(px-p0->c[0]) + (py-p0->c[1])*(py-p0->c[1]) + (pz-p0->c[2])*(pz-p0->c[2]);
          ldt = sqrt(ldt);

          //dp = drec + ddx*(p[0]-px) + ddy*(p[1]-py) + ddz*(p[2]-pz);

          // Si ldt est vraiment tout petit, faire quelque chose, du style tagger le point,
          // point, à la fin de l'itération, lui mettre la moyenne des points autour qui
          // n'ont pas été tagge. ou, extrapoler a partir de x + dt, ou d'une diretion rentrante.
          if( ldt < EPS1 ) {
            /*if ( dtmp[ip] > 0.0 )
              dtmp[ip] = D_MIN(dtmp[ip],drec + ddx*(p[0]-px) + ddx*(p[1]-py) + ddx*(p[2]-pz) + par->dt);
            else
              dtmp[ip] = D_MAX(dtmp[ip],drec + ddx*(p[0]-px) + ddx*(p[1]-py) + ddx*(p[2]-pz) - par->dt);
            */
            continue;
          }
          else {
            if ( (dtmp[ip] > 0.0) && (drec + ldt > 0.0) ) {
              dtmp[ip] = D_MIN(dtmp[ip],drec + ldt);
            }
            else if ( (dtmp[ip] < 0.0) && (drec - ldt < 0.0) ) {
              dtmp[ip] = D_MAX(dtmp[ip],drec - ldt);
            }
            continue;
          }
        }
        else {
          for (iel=1; iel<=mesh->ne; iel++) {
            pt1 = &mesh->tetra[iel];
            if ( inTetra(mesh,iel,p,cb) )  break;
          }
          if ( iel > mesh->ne ) {
            continue;
          }
        }
      }

      /* P1 Lagrange interpolation */
      pt1  = &mesh->tetra[iel];
      i0   = pt1->v[0];
      p0   = &mesh->point[i0];

      ddx = grad[3*(iel-1)+1];
      ddy = grad[3*(iel-1)+2];
      ddz = grad[3*(iel-1)+3];

      drec = sol->val[i0] + ddx*(p[0]-p0->c[0]) + ddy*(p[1]-p0->c[1]) + ddz*(p[2]-p0->c[2]);

      if ( (dtmp[ip] > 0.0) && (drec + par->dt > 0.0) ) {
      dtmp[ip] = D_MIN(dtmp[ip],drec + par->dt);
    }
    else if ( (dtmp[ip] < 0.0) && (drec - par->dt < 0.0) ) {
      dtmp[ip] = D_MAX(dtmp[ip],drec - par->dt);
    }
    }
  }
}


/* update the distance values in // */
static void upddist_3d(int istart,int istop,int ipth,Param *par) {
  pMesh    mesh;
  pSol     sol;
  pPoint   p0;
  double   res;
  int      k;

  mesh = par->mesh;
  sol  = par->sol;

  /* update sol */
  res = 0.0;
  for (k=istart; k<=istop; k++) {
    p0 = &mesh->point[k];
    if ( p0->tag )  continue;
    else if (( fabs(sol->val[k]) < 0.7 * par->dtfin )&&( p0->flag && p0->flag < mesh->mark-2 ))  continue;
    else if ( fabs(sol->val[k]) < 0.4 * par->dtfin )  continue;

    if ( fabs(par->dtmp[k]) < EPS1 )  continue;
    else if ( (sol->val[k] < 0.0) && (par->dtmp[k] < 0.0) ) {
      if ( par->dtmp[k] > sol->val[k] ) {
        res        += par->dtmp[k] - sol->val[k];
        sol->val[k] = par->dtmp[k];
        p0->flag    = mesh->mark;
      }
    }
    else if ( (sol->val[k] > 0.0) && (par->dtmp[k] > 0.0) ) {
      if ( par->dtmp[k] < sol->val[k] ) {
        res        += sol->val[k] - par->dtmp[k];
        sol->val[k] = par->dtmp[k];
        p0->flag    = mesh->mark;
      }
    }
  }

  par->res[ipth] += res;
}

/* expand distance function to the entire domain by solving the Eikonal eq. using characteristics */
int ppgdist_3d(Info info,pMesh mesh,pSol sol) {
  Param     par;
  double    res0;
  int       it,i;
  FILE     *out;
  int k;
  pTetra pt;

  /* memory alloc */
  par.grad = (double*)calloc(3*(mesh->ne)+1,sizeof(double));
  assert(par.grad);
  par.dtmp = (double*)malloc((mesh->np+1)*sizeof(double));
  assert(par.dtmp);
  memcpy(&par.dtmp[0],&sol->val[0],(mesh->np+1)*sizeof(double));

  par.mesh  = mesh;
  par.sol   = sol;
  par.dt    = info.dt;
	par.dtfin = 0.0;
	par.res   = (double*)calloc(info.ncpu,sizeof(double));
	assert(par.res);

  if ( info.ncpu > 1 ) {
    info.libpid = InitParallel(info.ncpu);
    assert(info.libpid);
    info.typ[0] = NewType(info.libpid,mesh->ne);
    info.typ[1] = NewType(info.libpid,mesh->np);
    //LaunchParallel(info.libpid,info.typ[0],0,(void *)gradelt_3d,(void *)&par);
  }
  /*else {
    gradelt_3d(1,mesh->ne,0,&par);
  }*/

  /* Tets intersecting boundary don't move*/
  for(k=1; k<=mesh->ne; k++) {
	  pt = &mesh->tetra[k];
	  if((fabs(sol->val[pt->v[0]])<sqrt(INIVAL_3d))&&(fabs(sol->val[pt->v[1]])<sqrt(INIVAL_3d))\
	     &&(fabs(sol->val[pt->v[2]])<sqrt(INIVAL_3d))&&(fabs(sol->val[pt->v[3]])<sqrt(INIVAL_3d)))
		  pt->flag = -1;
  }

  /* main cvg loop */
  it   = 1;
  res0 = 0.0;
	if ( info.ddebug )  out = fopen("residual","w");
  do {
  	//if(it<=60) par.dt = 0.1 * info.dt;
	  if ((it>60)&&(it<100)) par.dt = info.dt;
	  if ( it == 100 )  par.dt = 4.0 * info.dt;
	  //if(it >40) par.dt += 5e-5;

	  if ( info.ncpu > 1 ) {
	    LaunchParallel(info.libpid,info.typ[0],0,(void *)gradelt_3d,(void *)&par);
	  }
	  else {
	    gradelt_3d(1,mesh->ne,0,&par);
  	}

    mesh->mark++;
		memset(par.res,0,info.ncpu*sizeof(double));

    if ( info.ncpu > 1 ) {
      LaunchParallel(info.libpid,info.typ[0],0,(void *)tmpdist_3d,(void *)&par);
      LaunchParallel(info.libpid,info.typ[1],0,(void *)upddist_3d,(void *)&par);
			for (i=1; i<info.ncpu; i++)  par.res[0] += par.res[i];
    }
    else {
      tmpdist_3d(1,mesh->ne,0,&par);
      upddist_3d(1,mesh->np,0,&par);
    }

    if ( it == 1 )  res0 = par.res[0];
    else if ( par.res[0] < info.res * res0 )  break;

    fprintf(stdout,"     %9.7f  %8d\r",par.res[0]/res0,it);  fflush(stdout);
		if ( info.ddebug ) fprintf(out,"%E\n",par.res[0]);
		par.dtfin += par.dt;
  }
	while ( ++ it < info.maxit );

  if ( info.ddebug )  fclose(out);
  fprintf(stdout,"     Residual %E after %d iterations\n",par.res[0] / res0,it);
  free(par.grad);
  free(par.dtmp);

  //corrGrad_3d(mesh,sol);

  if ( info.ncpu > 1 ) {
    FreeType(info.libpid,info.typ[0]);
    FreeType(info.libpid,info.typ[1]);
    StopParallel(info.libpid);
  }

  return(1);
}
