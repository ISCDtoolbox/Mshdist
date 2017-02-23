#include "mshdist.h"

unsigned char inxt2[5] = {1,2,0,1,2};
extern Info  info;

/* Only used when generating a signed distance function
      from a domain explicitely discretized in current mesh */
typedef struct {
  int   ia,ib,num,nxt;
} Hedge;


Hedge *tab;
int    hnext;

/* Hash mesh edges for creating P2 nodes */
int hashEdge_2d(pMesh mesh) {
  pEdge    pa;
  int      k,ia,ib,kb,key;
  
  tab   = (Hedge*)calloc(2*mesh->na+1,sizeof(Hedge));
  hnext = mesh->na;
  
  for (k=1; k<=mesh->na; k++) {
    pa = &mesh->edge[k];
	ia = pa->v[0];
	ib = pa->v[1];
	key = D_MIN(ia,ib);
	key = 1 + key % mesh->na;
	kb  = D_MAX(ia,ib);
	
	if ( tab[key].ia == 0 ) {
	  tab[key].ia = D_MIN(ia,ib);
	  tab[key].ib = D_MAX(ia,ib);
	  tab[key].num = k;
	}
	else {
	  while ( tab[key].nxt ) {
		key = tab[key].nxt;
	  }
	  hnext++; 
	  assert(hnext < 2*mesh->na);
	  tab[key].nxt  = hnext;
	  tab[hnext].ia = D_MIN(ia,ib);
	  tab[hnext].ib = D_MAX(ia,ib);
	  tab[hnext].num = k;
    }
  }
  return(1);
}


int getEdge(pMesh mesh,int ia,int ib) {
  int   key,kb;
  
  key = D_MIN(ia,ib);
  key = 1 + key % mesh->na;
  kb  = D_MAX(ia,ib);
  
  if ( tab[key].ia == 0 ) {
	return(0);
  }
  else {
	while ( (tab[key].ia != D_MIN(ia,ib) || tab[key].ib != D_MAX(ia,ib)) && tab[key].nxt ) {  //PROBLEME ICI !!!!!
	  key = tab[key].nxt;
	}
	if ( tab[key].ia == D_MIN(ia,ib) && tab[key].ib ==D_MAX(ia,ib) )
	  return(tab[key].num);
	else
	  return(0);
  }
  return(0);
}

/* find background triangles intersecting the boundary, 
 and initialize distance at their vertices */
int iniredist_2d(pMesh mesh, pSol sol){
  pTria    pt;
  pPoint   p0,p1,p2,pa;
  double  *solTmp,d;
  int     *bndy,i,j,nb,nc,i0,i1,i2,proj;
	
  nb   = 0;	
  bndy = (int*)calloc(mesh->nt+1,sizeof(int));
  assert(bndy);	  

  /* store intersecting background triangles in list bndy */
  for (i=1; i<=mesh->nt; i++) {
	pt = &mesh->tria[i];
	i0 = pt->v[0]; 
	i1 = pt->v[1]; 
	i2 = pt->v[2]; 

	if ( (sol->val[i0] * sol->val[i1] <=0.)||(sol->val[i0] * sol->val[i2] <=0.)|| \
		 (sol->val[i1] * sol->val[i2] <=0.)){
	  nb++;
	  bndy[nb] = i;
	}	
  }
  
  bndy = (int*)realloc(bndy,(nb+1)*sizeof(int));	
	printf("nb= %d\n",nb);
  
  /* Temporary values stored in solTmp */	
  solTmp = (double*)calloc(mesh->np+1,sizeof(double));
  assert(solTmp);

  /* Check list bndy and compute distance at vertices of its elements */
  for (i=1; i<=nb; i++) {
	  pt = &mesh->tria[bndy[i]];
	  i0 = pt->v[0]; 
	  i1 = pt->v[1]; 
	  i2 = pt->v[2];

	  p0 = &mesh->point[i0]; 
	  p1 = &mesh->point[i1]; 
	  p2 = &mesh->point[i2]; 
	  
	  d = distnv0_2d(mesh, sol, bndy[i], p0, &proj);
	  
	  if (p0->tag == 0) {
	    solTmp[i0] = d;
	    p0->tag    = proj;
	  }
	  else if (d < fabs(solTmp[i0])) {
	    solTmp[i0] = d;
	    p0->tag    = proj;			
	  }
	  d = distnv0_2d(mesh, sol, bndy[i], p1, &proj);
	  
	  if (p1->tag == 0) {
	    solTmp[i1] = d;
	    p1->tag    = proj;
	  }
	  else if (d < fabs(solTmp[i1])) {
      solTmp[i1] = d;
	    p1->tag    = proj;			
	  }
    
    d = distnv0_2d(mesh, sol, bndy[i], p2, &proj);
	  
  	if (p2->tag == 0) {
	    solTmp[i2] = d;
	    p2->tag    = proj;
	  }
	  else if (d < fabs(solTmp[i2])) {
	    solTmp[i2] = d;
	    p2->tag    = proj;			
	  }
  }
  
  /* correction procedure, for points whose tag is 2 */
  nc = 0;
  for (i=1; i<=mesh->np; i++) {
	  pa = &mesh->point[i];
	  if ( pa->tag < 2 )  continue;
	  for (j=1; j<=nb; j++) {
	    pt = &mesh->tria[bndy[j]];
	    d  = distnv0_2d(mesh,sol,bndy[j],pa,&proj);
	    if ( proj == 1 && d < fabs(solTmp[i]) ) {
		    solTmp[i] = d;
		    break;
	    }
	  }
	  pa->tag = 1;
	  nc++;
  }
  if ( nc )   fprintf(stdout,"     %d correction(s)\n",nc);
 
  for (i=1; i<=mesh->np; i++){
	  pa = &mesh->point[i];
	  if ( !pa->tag ) 
	    sol->val[i] = sol->val[i] < 0.0 ? -sqrt(INIVAL_2d) : sqrt(INIVAL_2d);
	  else if ( pa->tag ) 
	    sol->val[i] = sol->val[i] < 0.0 ? -sqrt(solTmp[i]) : sqrt(solTmp[i]);
	}

  free(solTmp);
  free(bndy);
  return(1);
}

/* Initialize the exact unsigned distance function 
   to the boundary of mesh2 at vertices of the triangles of mesh1 intersecting this contour */
int inidist_2d(pMesh mesh1,pMesh mesh2,pSol sol1,pBucket bucket) {
  pTria      pt,pt1;
  pEdge      pe;
  pPoint     p1,p2,pa,pb;
  double     cb[3],d;
  int        *adja,*list,base,iadr,ilist,i,j,k,ia,ib,iel,jel,ier,cur,nc,tag,npp,i0,i1;

  for (k=1; k<=sol1->np; k++)
    sol1->val[k] = INIVAL_2d;
  sol1->nt = mesh1->nt; 	
  
  /* Memory allocation */
  list = (int*)calloc(mesh1->nt+1,sizeof(int));
  assert(list);
  sol1->ref = (int*)calloc(mesh1->nt+1,sizeof(int));
  assert(sol1->ref);	

  /* set all references to -1 by default */	
  for (k=1; k<=sol1->nt; k++)
    sol1->ref[k] = -1;

  /* set tag=2 if elt intersected, 1 if connected to elt intersected, 0 else */
  nc  = 0;
  npp = 0;
  for (k=1; k<=mesh2->na; k++){
    pe  = &mesh2->edge[k];
	  p1  = &mesh2->point[pe->v[0]];
    p2  = &mesh2->point[pe->v[1]];
    iel = buckin(mesh1,bucket,p1->c);
    iel = locelt(mesh1,iel,p1->c,cb);
	  
	  if(!iel){
	    iel = buckin(mesh1,bucket,p2->c);
		  iel = locelt(mesh1,iel,p2->c,cb);
	  }
	  
    ilist       = 1;
    list[ilist] = iel;
    base = ++mesh1->flag;
    pt   = &mesh1->tria[iel];
    pt->flag = base;
    pt->tag  = 2;
    cur      = 1;
    
    do {
      iel  = list[cur];
      pt   = &mesh1->tria[iel];
	  
    /* Intersecting triangle inherits reference from the underlying edge of mesh2 */
	    if(sol1->ref[iel] == -1)  
	      sol1->ref[iel] = pe->ref;
	    else if((sol1->ref[iel] >= 0)&&(sol1->ref[iel] != pe->ref)){
		    if ( info.ddebug )  fprintf(stdout,"         WARNING : conflict in setting reference in triangle %d \n", iel);
		    sol1->ref[iel] = pe->ref;
		    npp++;
	    }

      iadr = 3*(iel-1) + 1;
      adja = &mesh1->adja[iadr];

      for (i=0; i<3; i++) {
        ia = inxt2[i];
        ib = inxt2[ia];
        pa = &mesh1->point[pt->v[ia]];      
        pb = &mesh1->point[pt->v[ib]];       
        ier = intersec_2d(p1,p2,pa,pb);
        jel = adja[i] / 3;
		    if(!jel ) continue;
		
        pt1 = &mesh1->tria[jel];
        if ( ier ) {
          pt1->tag = 2;
          if ( jel && pt1->flag < base ) {
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
      pt  = &mesh1->tria[iel];
      for (j=0; j<3; j++) {
        ia = pt->v[j];
        pa = &mesh1->point[ia];
        d  = distpt_2d(p1,p2,pa,&tag);
        if ( d < sol1->val[ia] ) {
          sol1->val[ia] = d;
          pa->tag = tag;
        }
      }
    }
  }
  fprintf(stdout,"     distance\n");
	
  /* correction */
  nc = 0;
  for (k=1; k<=mesh1->np; k++){
    pa = &mesh1->point[k];
    if ( pa->tag < 2 )  continue;

    /* possible optimization here */
    for (i=1; i<=mesh2->na; i++) {
      pe = &mesh2->edge[i];
      p1 = &mesh2->point[pe->v[0]];      
      p2 = &mesh2->point[pe->v[1]];       
      d  = distpt_2d(p1,p2,pa,&tag);
      if ( tag == 1 && d < sol1->val[k] ) {
        sol1->val[k] = d;
        break;
      }
    }
    pa->tag = 1;
    nc++;
  }
  if ( nc )   fprintf(stdout,"     %d correction(s)\n",nc);

  if ( npp )	
	fprintf(stdout,"  ##  WARNING : conflict in setting reference in %d triangles\n",npp);
	
  printf("NB D ARETES %d\n", mesh1->na);

  /* if bbbc mode has been actived, take references from the boundary box */
  if(info.bbbc){
    for(k=1;k<=mesh1->nt;k++){
	    sol1->ref[k] =0;
	  }
	
	  ier = hashEdge_2d(mesh1);
	  assert(ier);

	  for(k=1;k<=mesh1->nt;k++){
	    pt = &mesh1->tria[k];
	  
	    for(j=0;j<3;j++){
	      i0 = pt->v[inxt2[j]];
		    i1 = pt->v[inxt2[inxt2[j]]];

        ia = getEdge(mesh1,i0,i1);
        if(!ia) continue;
		    pe = &mesh1->edge[ia];
		    sol1->ref[k] = pe->ref; 
	    }
	  }
  }	
  
  free(list);
  return(1);
}

/* Initialize the sign of the implicit function in each connected component */
int sgndist_2d(pMesh mesh,pMesh mesh2,pSol sol,pBucket bucket) {
  pEdge    pe;
  pTria    pt,pt1;
  pPoint   ppt,pi,p1,p2,pa,pb;
  double   p[2],cb[3];
  int     *pile,*adja,*pilcc,i,j,k,cc,kk,iel,base,bfin,iadr,ipile,start,ip,nc;
  char     j1,j2,sgn;
  
  /* memory alloc */
  pile = (int*)calloc(mesh->nt+1,sizeof(int));
  assert(pile);

  printf("coords %E, %E ", mesh->point[1].c[0], mesh->point[1].c[1]); 

  /* identify triangle close to lower corner (boundary) */
  p[0] = 0.05;
  p[1] = 0.05;
  
/*  //for gripping mechanism
  p[0] = 0.31;
  p[1] = 0.62;
  //p[1] = 0.02;
*/    
  iel = buckin(mesh,bucket,p);
  iel = locelt(mesh,iel,p,cb);
  assert(iel);
  printf("ELEMENT DE DEPART POUR LE SIGNE %d \n", iel);

  /* reset colors */
  for (k=1; k<=mesh->nt; k++)  mesh->tria[k].flag = 0;
  for (k=1; k<=mesh->np; k++)  mesh->point[k].flag = 0;
  mesh->flag = 0;

  ipile       = 1;
  pile[ipile] = iel;
  base = ++mesh->flag;
  pt   = &mesh->tria[iel];
  pt->flag = base;
  start    = 1;
  while ( ipile > 0 ) {
    /* search for all elements in the current connected component */
    do {
      k = pile[ipile];
      ipile--;

      pt = &mesh->tria[k];
      iadr = 3*(k-1) + 1;
      adja = &mesh->adja[iadr];
      for (i=0; i<3; i++) {
        if ( !adja[i] )  continue;
        kk  = adja[i] / 3;
        pt1 = &mesh->tria[kk];
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
    for (k=start; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( pt->flag >= 1 )  continue;      
      else if (pt->tag <2) { 
        base = ++mesh->flag;
        ipile = 1;
        pile[ipile] = k;
        pt->flag    = base;
        break;                  
      }
    }
  }

  /* At this point, each non-boundary element is associated to a flag */
  
  /* analyze components */
  if ( fabs(info.imprim) > 3 )  fprintf(stdout,"     %d connected component(s)\n",base);
  if ( base < 2 )  return(-1);

  if ( info.ddebug ) {
		FILE  *out;
		
		out = fopen("cc.sol","w");
		fprintf(out,"MeshVersionFormatted 1\n\nDimension 2\n\n");
		fprintf(out,"SolAtTriangles\n%d\n 1 1 \n",mesh->nt);
		for (k=1; k<=mesh->nt; k++)
			fprintf(out,"%d\n",mesh->tria[k].flag);
		fclose(out);
	}
	
	
  /* store triangles intersected */
  bfin  = mesh->flag;   
  sgn   = 1;
  pilcc = (int*)calloc(bfin+1,sizeof(int));
  assert(pilcc);
  base     = ++mesh->flag;   
  pilcc[1] = base;
  ipile    = 0;
  do {
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( pt->flag < 1 )  continue; 
      for (cc=1; cc<=bfin; cc++) {   
        if ( pilcc[cc] == base && pt->flag == cc ) {   
          iadr = 3*(k-1) + 1;
          adja = &mesh->adja[iadr];
          for (i=0; i<3; i++) {
            if ( !adja[i] )  continue;
            kk  = adja[i] / 3;
            pt1 = &mesh->tria[kk];
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
    /* At this point, pile stores all boundary trias which neighbour the current component */
    
    /* Complete boundary pile */
    while ( ipile > 0 ) {
      k = pile[ipile];
      ipile--;
      pt = &mesh->tria[k];
      iadr = 3*(k-1) + 1;
      adja = &mesh->adja[iadr];
      for (i=0; i<3; i++) {
        if ( !adja[i] )  continue;
        kk  = adja[i] / 3;
        pt1 = &mesh->tria[kk];
        if ( pt1->flag == 0 ) {
          ipile++;
          pile[ipile] = kk;
          pt1->flag   = -1;
        }
        else if ( pt1->flag > 0 && pilcc[pt1->flag] != base )
          pilcc[pt1->flag] = base + 1; 
      }
    }

    if ( info.ddebug ) {
      fprintf(stdout,"     component list: ");
      for (cc=1; cc<=bfin; cc++)
        if ( pilcc[cc] == base )  printf(" %d  ",cc);
      printf("\n");
    }

    /* update sign */
    base = ++mesh->flag;
    sgn = -sgn;
    for (i=1; i<=bfin; i++) {
      if ( pilcc[i] != base )  continue; 
      for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( pt->flag == i && pt->tag < 2 ) {
          for (j=0; j<3; j++) {
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
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( pt->flag == 1 ) {
      for (j=0; j<3; j++) {
        ppt = &mesh->point[pt->v[j]];
        ppt->flag = 1;
      }
    }
  }

  /* Treatment of uninitialized triangles */
  ipile = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->tag )  continue;
    for (j=0; j<3; j++) {
      ppt = &mesh->point[pt->v[j]];   
      if ( !ppt->flag ) {
        ipile++;
        pile[ipile] = k;
        break;
      }
    }
  }
  fprintf(stdout,"  %%%% %d elements with a vertex uninitialized\n",ipile);
  
  /* Assignation of a sign to the vertices of uninitialized elements by comparing them with an
     already defined vertex */ 
  nc = 0;
	
  for (i=1; i<=ipile; i++) {
    k  = pile[i];
    pt = &mesh->tria[k]; 
    
    /* Search for an already initialized vertex */
    for (j=0; j<3; j++) {
      ip = pt->v[j];
      pi = &mesh->point[ip];
      if (( pi->flag )&&(fabs(sol->val[ip]>EPS)))  break;
    }
   
    /* unable to correct sign */
    if ( j == 3 )  continue;

    /* Search for uninitialized vertices in tria */
    j1  = inxt2[j+0];
    j2  = inxt2[j+1];
    p1  = &mesh->point[pt->v[j1]];
    p2  = &mesh->point[pt->v[j2]];
	
    if ( !p1->flag || !p2->flag ) {
      for (k=1; k<=mesh2->na; k++) {
        pe  = &mesh2->edge[k];
        pa  = &mesh2->point[pe->v[0]];
        pb  = &mesh2->point[pe->v[1]];
        if ( !p1->flag && intersec_2d(pi,p1,pa,pb) ) {
		      sol->val[pt->v[j1]] = sol->val[ip] > 0 ? -fabs(sol->val[pt->v[j1]]) : fabs(sol->val[pt->v[j1]]);
          p1->flag = 1; 
		      nc++;
        }
        if ( !p2->flag && intersec_2d(pi,p2,pa,pb) ) {
		      sol->val[pt->v[j2]] = sol->val[ip] > 0 ? -fabs(sol->val[pt->v[j2]]) : fabs(sol->val[pt->v[j2]]);
		      p2->flag = 1;
		      nc++;
        }
        if ( p1->flag && p2->flag )  break;
      }
      
      if ( k > mesh2->na ) {
        if ( !p1->flag ) {
          sol->val[pt->v[j1]] = sol->val[ip] > 0 ? fabs(sol->val[pt->v[j1]]) : -fabs(sol->val[pt->v[j1]]);
          p1->flag = 1;
					nc++;
					
        } 
        if ( !p2->flag ) {
          sol->val[pt->v[j2]] = sol->val[ip] > 0 ? fabs(sol->val[pt->v[j2]]) : -fabs(sol->val[pt->v[j2]]);
          p2->flag = 1;
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
  
  // case when initializing a l-set from full domain
  /*for (k=1; k<=mesh->np; k++) {
    mesh->point[k].flag = 0;
    sol->val[k] = sol->val[k] < 0.0 ? sqrt(fabs(sol->val[k])) : -sqrt(sol->val[k]);
  }*/

  return(1);
}

/* Initialize a (unsigned) distance function to a domain already existing in mesh, defined 
 by boundary triangles of reference contained in info.sref */
int inireftrias_2d(pMesh mesh, pSol sol){
  
  return(1);
}

/* Initialize a signed distance function to a domain already existing in mesh,
    defined by pt->ref = REFINT */
int iniencdomain_2d(pMesh mesh, pSol sol){
  int ier,*actiedg,nb,k,i0,i1,ia, *adja,iel,jel;
  int ilist, *list, base,cur,nc,j,iadr,i,ib,tag; 
  pTria pt,pt1;
  pEdge pe;
  pPoint pa,pb,p1,p2;
  double d;
  
  /* Hash edges of the mesh */
  ier = hashEdge_2d(mesh);
  assert(ier);
  
  nb = 0;
  actiedg = (int*)calloc(mesh->na+1,sizeof(int));

  for(k=1;k<=mesh->nt;k++){
    pt = &mesh->tria[k];

    if ( pt->ref != REFINT ) continue;
    for(j=0; j<3; j++){
      i0 = inxt2[j];
      i1 = inxt2[i0];
      iadr = 3*(k-1) + 1;
      adja = &mesh->adja[iadr];
      iel = adja[j]/3; //note iel can be 0, if pt is near an external boundary... Such edges must be accepted
      
      if( /*!iel || ((mesh->tria[iel]).ref != REFINT )*/ iel && ((mesh->tria[iel]).ref != REFINT)){
        ia = getEdge(mesh, pt->v[i0], pt->v[i1]);
        assert(ia); 
    
        if(((mesh->edge[ia]).ref == REFDIR)|| ((mesh->edge[ia]).ref == REFDIRSWP) || ((mesh->edge[ia]).ref == REFSYM)) continue; //MODIF for moving boundary conditions
        
        nb++;
        actiedg[nb] = ia;
        (mesh->edge[ia]).flag = k; // starting triangle
      }
    }
  }
  
  printf("Number of active edges : %d \n", nb);
  
  for (k=1; k<=sol->np; k++) 
    sol->val[k] = INIVAL_2d;
  
  sol->nt = mesh->nt; 	
  /* memory alloc */
  list = (int*)calloc(mesh->nt+1,sizeof(int));
  assert(list);
    
  /* travel the list of boundary edges */	
  for (k=1; k<=nb; k++) {
    pe  = &mesh->edge[actiedg[k]];
	  p1 = &mesh->point[pe->v[0]];
	  p2 = &mesh->point[pe->v[1]];
    iel = pe->flag; 
    assert(iel);
	
    ilist       = 1;
    list[ilist] = iel;
    base = ++mesh->flag;
    pt   = &mesh->tria[iel];
    pt->flag = base;
    pt->tag  = 2;
    cur      = 1;
    do {
      iel  = list[cur];
      pt   = &mesh->tria[iel];
	  
      iadr = 3*(iel-1) + 1;
      adja = &mesh->adja[iadr];
      
      for (i=0; i<3; i++) {
        ia = inxt2[i];
        ib = inxt2[ia];
        pa = &mesh->point[pt->v[ia]];      
        pb = &mesh->point[pt->v[ib]];       
        ier = intersec_2d(p1,p2,pa,pb);
        jel = adja[i] / 3;
        pt1 = &mesh->tria[jel];
        if ( ier ) {
          pt1->tag = 2;
          if ( jel && pt1->flag < base ) {
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
      pt  = &mesh->tria[iel];
      for (j=0; j<3; j++) {
        ia = pt->v[j];
        pa = &mesh->point[ia];
        d  = distpt_2d(p1,p2,pa,&tag);
        if ( d < sol->val[ia] ) {
          sol->val[ia] = d;
          pa->tag = tag;
        }
      }
    }
  }
  fprintf(stdout,"     distance\n");
  
  /* Set all points on the boundary to exactly 0 */
  for (k=1; k<=nb; k++) {
    pe  = &mesh->edge[actiedg[k]];
	  i0 = pe->v[0];
	  i1 = pe->v[1];
	  sol->val[i0] = 0.0;
	  sol->val[i1] = 0.0;
  }	  
  
  /* correction */
  nc = 0;
  for (k=1; k<=mesh->np; k++) {
    pa = &mesh->point[k];
    if ( pa->tag < 2 )  continue;
	
    /* possible optimization here */
    for (i=1; i<=nb; i++) {
      pe = &mesh->edge[actiedg[i]];
      p1 = &mesh->point[pe->v[0]];      
      p2 = &mesh->point[pe->v[1]];       
      d  = distpt_2d(p1,p2,pa,&tag);
      if ( tag == 1 && d < sol->val[k] ) {
        sol->val[k] = d;
        break;
      }
    }
    pa->tag = 1;
    nc++;
  }
  if ( nc )   fprintf(stdout,"     %d correction(s)\n",nc);
  
  /* Put sign in the function, and take sqrt */ 
  for(k=1; k<=mesh->nt;k++){
    pt = &mesh->tria[k];
	  for(j=0;j<3;j++){
	    ia = pt->v[j];
	    pa = &mesh->point[ia];
	    if(pa->flag ==1) continue;
	    pa->flag =1;
	    if(pt->ref == REFINT){
	      sol->val[ia] = -sqrt(sol->val[ia]);
      }
	    else{
        sol->val[ia] = sqrt(sol->val[ia]);
	    } 
  	}
  }
  
  /* reset flags */
  for(k =1;k<=mesh->np;k++){
    pa = &mesh->point[k];
	  pa->flag = 0;
  }
  
  free(actiedg);
  return(1);
}

/* compute gradients of shape functions P1: 9 values (dx,dy,cte) */
static void gradelt_2d(int istart,int istop,int ipth,Param *par) {
  pMesh     mesh;
  pTria     pt;
  pPoint    p0,p1,p2;
	double   *grad,v0x,v0y,v1x,v1y,v2x,v2y,det,dx,dy;
  int      k,iadr;

  mesh = par->mesh;
  grad = par->grad;

  for(k=istart; k<=istop; k++) {
    pt = &mesh->tria[k]; 
    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];
    
    /* lambda_0 */
    v1x = p1->c[0] - p0->c[0]; 
    v1y = p1->c[1] - p0->c[1];
    v2x = p2->c[0] - p0->c[0];
    v2y = p2->c[1] - p0->c[1];
    det = v1x*v2y - v2x*v1y;
    dx  = dy = 0.0;
    if ( fabs(det) > EPS1 ) {
      det = 1.0 / det; 
      dx = (-v2y + v1y) * det;
      dy = ( v2x - v1x) * det;
    }
    iadr = 9 * (k-1) + 1;
    grad[iadr + 0] = dx;    
    grad[iadr + 1] = dy;
    grad[iadr + 2] = 1.0 - dx*(p0->c[0]) - dy*(p0->c[1]);

    /* lambda_1 */
    v0x = -v1x;
    v0y = -v1y;
    v2x = p2->c[0] - p1->c[0];
    v2y = p2->c[1] - p1->c[1];
    det = v0x*v2y - v2x*v0y;
    dx  = dy = 0.0; 
    if ( fabs(det) > EPS1 ) {
      det = 1.0 / det;
      dx  = (-v2y + v0y) * det;
      dy  = ( v2x - v0x) * det;
    }
    grad[iadr + 3] = dx;    
    grad[iadr + 4] = dy;
    grad[iadr + 5] = 1.0 - dx*(p1->c[0]) - dy*(p1->c[1]);
    
    /* lambda_2 */
    v0x = p0->c[0] - p2->c[0]; 
    v0y = p0->c[1] - p2->c[1];
    v1x = p1->c[0] - p2->c[0];
    v1y = p1->c[1] - p2->c[1];
    det = v0x*v1y - v1x*v0y;
    dx = dy = 0.0;
    if ( fabs(det) > EPS1 ) {
      det = 1.0 / det;
      dx  = (-v1y + v0y) * det;
      dy  = ( v1x - v0x) * det;
    }
    grad[iadr + 6] = dx;
    grad[iadr + 7] = dy;
    grad[iadr + 8] = 1.0 - dx*(p2->c[0]) - dy*(p2->c[1]);
  }
}

/* compute new distance at points */
static void tmpdist_2d(int istart,int istop,int ipth,Param *par) {
  pMesh     mesh;
  pSol      sol;
  pTria     pt,pt1;
  pPoint    p0;
  double   *grad,*dtmp,dx,dy,dd,p[2],cb[3],drec;
  int       j,k,base,iadr,iel,ip;

  mesh = par->mesh;
  sol  = par->sol;
  grad = par->grad;
  dtmp = par->dtmp;

  base = mesh->flag;

  for (k=istart; k<=istop; k++) {
    pt   = &mesh->tria[k]; 
	
	/*test : say pt->flag = -1 when gradient is close to 1 in pt*/
	if(pt->flag ==-1) continue;  
	  
    iadr = 9 * (k-1) + 1;
    dx = sol->val[pt->v[0]] * grad[iadr+0] + sol->val[pt->v[1]] * grad[iadr+3] + sol->val[pt->v[2]] * grad[iadr+6];
    dy = sol->val[pt->v[0]] * grad[iadr+1] + sol->val[pt->v[1]] * grad[iadr+4] + sol->val[pt->v[2]] * grad[iadr+7];
    dd = dx*dx + dy*dy;
  
	/*test : don't move when gradient is close to 1 in pt*/
	  if(fabs(dd -1.) <0.005) {
		  pt->flag =-1;
		  continue;
	  }
	  
	  
    /* vector to track characteristic line */
    if ( dd < EPS1 )  continue;
    dd = par->dt / sqrt(dd);
    dx *= dd;
    dy *= dd;
  
    for (j=0; j<3; j++) {
      ip = pt->v[j];
      p0 = &mesh->point[ip];
      if ( p0->tag > 0 )  continue;
      else if (( fabs(sol->val[ip]) < 0.8 * par->dtfin )&&( p0->flag && p0->flag < mesh->mark-2 )) continue;
      else if ( fabs(sol->val[ip]) < 0.1 * par->dtfin )  continue;
  
      /* follow characteristic line */
      if ( fabs(sol->val[ip]) < EPS1 )
        continue;
		
      else if ( sol->val[ip] > 0.0 ) {
        p[0] = p0->c[0] - dx; 
        p[1] = p0->c[1] - dy;
      }
      else {
        p[0] = p0->c[0] + dx; 
        p[1] = p0->c[1] + dy;
      }
      if ( p[0] <= EPS1 || p[0] >= 1.0-EPS1 )  continue;
      if ( p[1] <= EPS1 || p[1] >= 1.0-EPS1 )  continue;
  
      /* find enclosing triangle, k is guessed */ 
      iel = nxtelt(mesh,k,p,cb);
      if ( iel < 1 ) {
        fprintf(stdout,"--> exhaustive search: %d",k); fflush(stdout);
        for (iel=1; iel<=mesh->nt; iel++) {
          pt1 = &mesh->tria[iel];
          if ( inTria(mesh,iel,p,cb) )  break;
        }
        if ( iel > mesh->nt ) {
          fprintf(stdout,"\n  ## Oops...no simplex found (%d).\n",k);
          continue;
        }
        else
          fprintf(stdout," %d\n",iel);
      }
  
      /* P1 Lagrange interpolation */
      pt1  = &mesh->tria[iel]; 
      iadr = 9 * (iel-1) + 1;
      drec = sol->val[pt1->v[0]] * (p[0]*grad[iadr+0] + p[1]*grad[iadr+1] + grad[iadr+2])
           + sol->val[pt1->v[1]] * (p[0]*grad[iadr+3] + p[1]*grad[iadr+4] + grad[iadr+5])
           + sol->val[pt1->v[2]] * (p[0]*grad[iadr+6] + p[1]*grad[iadr+7] + grad[iadr+8]);   
  
      if ( (dtmp[ip] > 0.0) && (drec + par->dt > 0.0) ) {
        dtmp[ip] = D_MIN(dtmp[ip],drec + par->dt);
      }
      else if ( (dtmp[ip] < 0.0) && (drec - par->dt < 0.0) ) {
        dtmp[ip] = D_MAX(dtmp[ip],drec - par->dt);
      }
    }
  }
		
}

/* Proceed to update of distance table from an iteration to the other */
static void upddist_2d(int istart,int istop,int ipth,Param *par) {
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
    else if (( fabs(sol->val[k]) < 0.8 * par->dtfin )&&( p0->flag && p0->flag < mesh->mark-2 ))  continue;
    else if ( fabs(sol->val[k]) < 0.1 * par->dtfin )  continue;

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


/* Expand distance function to the domain by solving
   Eikonal equation using method of characteristics */
int ppgdist_2d(pMesh mesh,pSol sol) {
	Param     par;
	double    res0;
	int       it,i,j;
	FILE     *out;
	int k;
	pTria pt;

  /*char numit[10];
  char *ptr;*/
  
  /* Memory allocation */
  par.grad = (double*)calloc(9*(mesh->nt)+1,sizeof(double));
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
    info.typ[0] = NewType(info.libpid,mesh->nt);
    info.typ[1] = NewType(info.libpid,mesh->np);
    LaunchParallel(info.libpid,info.typ[0],0,(void *)gradelt_2d,(void *)&par);
  }
  else {
    gradelt_2d(1,mesh->nt,0,&par);
  }
   /*test : say triangles intersecting boundary don't move*/
	for(k = 1; k<= mesh->nt; k++){
		pt = &mesh->tria[k]; 
		if((fabs(sol->val[pt->v[0]])<INIVAL_2d)&&(fabs(sol->val[pt->v[1]])<INIVAL_2d) &&(fabs(sol->val[pt->v[2]])<INIVAL_2d)){
			pt->flag = -1;
		}
	}
		
  /* main cvg loop */
  it   = 1;
  res0 = 0.0;

	if ( info.ddebug )  out = fopen("residual","w");
  do {
	if(it<=60) par.dt = 0.1 * info.dt;
	if ((it>60)&&(it<100)) par.dt = info.dt;
	if ( it == 100 )  par.dt = 4.0 * info.dt;
	//if(it>40) par.dt += info.dt/4.; 
	mesh->mark++;
	  
	  /*for(j=1; j<= mesh->np; j++){
		  if(sol->val[j]>0.0) par.dtmp[j] = INIVAL_2d;
		  if(sol->val[j]==0.0) par.dtmp[j] = 0.0;
		  if(sol->val[j]<0.0) par.dtmp[j] = -INIVAL_2d;
		  
	  }*/
		memset(par.res,0,info.ncpu*sizeof(double));
    if ( info.ncpu > 1 ) {
      LaunchParallel(info.libpid,info.typ[0],0,(void *)tmpdist_2d,(void *)&par);
      LaunchParallel(info.libpid,info.typ[1],0,(void *)upddist_2d,(void *)&par);
			for (i=1; i<info.ncpu; i++)  par.res[0] += par.res[i];
    }
    else {
      tmpdist_2d(1,mesh->nt,0,&par);
      upddist_2d(1,mesh->np,0,&par);
    }

    if ( it == 1 )  res0 = par.res[0];
    else if ( par.res[0] < info.res * res0 )  break;
    fprintf(stdout,"     %9.7f  %8d\r",par.res[0]/res0,it);  fflush(stdout);
		if ( info.ddebug ) fprintf(out,"%E\n",par.res[0]);
		par.dtfin += par.dt;
    
    /*strcat(sol->name,".");
    sprintf(numit,"%d",it);
    strcat(sol->name,numit);
    if ( !saveSol(sol) )     return(1);
    numit[0] ='\0';
    ptr = strstr(sol->name,".");
    ptr[0] = '\0';*/
    
    
  }  
  while ( ++ it < info.maxit );
	if ( info.ddebug )  fclose(out);

  fprintf(stdout,"     Residual %E after %d iterations\n",par.res[0] / res0,it);
  free(par.grad);
  free(par.dtmp);

  if ( info.ncpu > 1 ) {
    FreeType(info.libpid,info.typ[0]);
    FreeType(info.libpid,info.typ[1]);
    StopParallel(info.libpid);
  }
  
  return(1);
}

/* Compute L^1,L^2 and L^\infty errors when compared to interpolation of real sdf */
int errdist(pMesh mesh, pMesh mesh2, pSol sol){
  int k,l,proj,i0,i1,i2;
  pTria pt;
  pEdge pe;
  pPoint p0,p1,p2;
  double errLInfty,errL1,errL2,*dist,area;
  
  errLInfty = 0.0;
  errL1     = 0.0;
  errL2     = 0.0;
  
  dist = (double*)calloc(mesh->np+1,sizeof(double));
  assert(dist);
  
  /* Computation of exact unsigned distance function to mesh2 */
  for(k=1;k<=mesh->np;k++){
    dist[k] = INIVAL_2d;
  }
  
  for(k=1;k<=mesh->np;k++){
    p0 = &mesh->point[k];
    for(l=1;l<=mesh2->na;l++){
	    pe  = &mesh2->edge[l];
	    p1  = &mesh2->point[pe->v[0]];
	    p2  = &mesh2->point[pe->v[1]];
	  
	    dist[k] = D_MIN(dist[k],distpt_2d(p1,p2,p0,&proj));	
	  } 
  }

  for(k=1;k<=mesh->np;k++){
    dist[k] = sqrt(dist[k]);
  }
  
  /* Computation of L^\infty error */ 
  for(k=1;k<=mesh->np;k++){
    errLInfty = D_MAX(errLInfty,fabs(fabs(sol->val[k]) - dist[k]));
  } 
  
  /* Computation of L^1 and L^2 errors */
  for(k=1;k<=mesh->nt;k++){
    pt = &mesh->tria[k];
	  i0 = pt->v[0];
	  i1 = pt->v[1];
	  i2 = pt->v[2];
	
  	p0 = &mesh->point[i0];
  	p1 = &mesh->point[i1];
  	p2 = &mesh->point[i2];
	
	  area = 0.5*fabs((p1->c[0]-p0->c[0])*(p2->c[1]-p0->c[1]) - \
	           (p1->c[1]-p0->c[1])*(p2->c[0]-p0->c[0]));
	
	  /* Exact quadrature formula for \mathbb{P}^1 functions */
	  if((sol->val[i0]*sol->val[i1] <= 0.0) || (sol->val[i0]*sol->val[i2] <= 0.0))
	    continue;
	
	  errL1 += area/3.0*fabs(dist[i0]-fabs(sol->val[i0]) + dist[i1]-fabs(sol->val[i1]) \
	           + dist[i2]-fabs(sol->val[i2]));  
			   
  	/* Exact quadrature formula for \mathbb{P}^2 functions */	
	  errL2 += (area/3.0)*( (0.5*dist[i0] + 0.5*dist[i1] - 0.5*fabs(sol->val[i0]) - 0.5*fabs(sol->val[i1]))* \  
	    (0.5*dist[i0] + 0.5*dist[i1] - 0.5*fabs(sol->val[i0]) - 0.5*fabs(sol->val[i1])) + \
			  (0.5*dist[i0] + 0.5*dist[i2] - 0.5*fabs(sol->val[i0]) - 0.5*fabs(sol->val[i2]))* \  
				  	   (0.5*dist[i0] + 0.5*dist[i2] - 0.5*fabs(sol->val[i0]) - 0.5*fabs(sol->val[i2])) +\
			  (0.5*dist[i1] + 0.5*dist[i2] - 0.5*fabs(sol->val[i1]) - 0.5*fabs(sol->val[i2]))* \  
				  	(0.5*dist[i1] + 0.5*dist[i2] - 0.5*fabs(sol->val[i1]) - 0.5*fabs(sol->val[i2])));
  }
 
  errL2 = sqrt(errL2);
  printf("Linfty error with respect to the exact signed distance function : \t %f\n",errLInfty);
  printf("L1 error with respect to the exact signed distance function : \t %f\n",errL1);  
  printf("L2 error with respect to the exact signed distance function : \t %f\n",errL2);  

  free(dist);
  return(1);
}
