#include "mshdist.h"

extern unsigned char inxt2[5];
extern unsigned char inxt3[7];


/* Store in list the points in the ball of ip; return length of the ball */
int boulep_2d(pMesh mesh,int ip,int *list) {
  pTria     pt;
  pPoint    p0;
  int       start,k,ipl,ilist,*adja;
  char      i,i1,i2;
  
  ilist = 0;
  
  p0 = &mesh->point[ip];
  start = p0->s;
  k = start;
  pt = &mesh->tria[k];
  
  for (i=0; i<3; i++)
    if ( pt->v[i] == ip ) break;
  
  assert ( i < 3 );
  i1 = inxt2[i];
  ipl = pt->v[i1];
  
  do {
    pt = &mesh->tria[k];
    i1 = inxt2[i];
    i2 = inxt2[i1];
    
    if ( ilist >= LONMAX ) {
      printf(" **** Problem function boulep_2d; point %d, more than %d points in the ball; abort.\n",ip,LONMAX);
      exit(0);
    }
    
    list[ilist] = pt->v[i2];
    adja = &mesh->adja[3*(k-1)+1];
    k = adja[i1] / 3;
    i = adja[i1] % 3;
    ilist++;
  }
  while ( k && k != start );
  
  /* If travel of the ball ends because an external boundary has been met, add the first point in the 
   starting triangle to the ball */
  if ( !k ) {
    if ( ilist >= LONMAX ) {
      printf(" **** Problem function boulep_2d; point %d, more than %d points in the ball; abort.\n",ip,LONMAX);
      exit(0);
    }
    list[ilist] = ipl;
    ilist++;
  }
  
  return(ilist);
}

/* Store in list the triangles in the ball of point ip in triangle start; return length of the ball */
int boulet_2d(pMesh mesh,int start,int ip,int *list) {
  pTria     pt;
  pPoint    p0;
  int       k,ilist,*adja;
  char      i,i1,i2;
  
  ilist = 0;
  k = start;
  i = ip;
  
  /* First loop, in the trigonometric sense */
  do {
    list[ilist] = 3*k+i;
    ilist++;
    if ( ilist >= LONMAX ) {
      fprintf(stdout," **** Problem function boulet_2d; more than %d points in a ball; abort.\n",LONMAX);
      exit(0);
    }
    
    i1 = inxt2[i];
    adja = &mesh->adja[3*(k-1)+1];
    k = adja[i1] / 3;
    i2 = adja[i1] % 3;
    i = inxt2[i2];
  }
  while( k && k!= start );
  
  /* If the loop has ended because start has been reached, return */
  if ( k == start ) return(ilist);
  
  /* Travel in the converse sense */
  k   = start;
  i   = ip;
  
  while( k ) {
    i2 = inxt2[inxt2[i]];
    adja = &mesh->adja[3*(k-1)+1];
    k = adja[i2] / 3;
    if ( k == 0 ) break;
    i1 = adja[i2] % 3;
    i = inxt2[inxt2[i1]];
    
    list[ilist] = 3*k+i;
    ilist++;
    if ( ilist >= LONMAX ) {
      fprintf(stdout," **** Problem function boulet_2d; more than %d points in a ball; abort.\n",LONMAX);
      exit(0);
    }
  }
  
  return(ilist);
}

/* Return volumic ball (i.e. filled with tetrahedra) of point ip in tetra start. 
Results are stored under the form 4*kel + jel , kel = number of the tetra, jel = local 
index of p within kel */
int boulet_3d (pMesh mesh, int start, int ip, int * list){
  pTetra pt,pt1;
  int    nump,ilist,base,cur,k,k1, *adja;
  char   j,l,i;  
  
  base = ++mesh->base;
  ilist = 0;
   
  pt = &mesh->tetra[start];
  nump = pt->v[ip];
  
  /* Store initial tetrahedron */
  pt->flag = base;
  list[ilist] = 4*start + ip;
  ilist++;
  
  /* Explore list and travel by adjacency through elements sharing p */
  cur = 0;
  while(cur < ilist){
    k = list[cur]/4;
    i = list[cur]%4; // index of point p in tetra k
	  adja = &mesh->adja[4*(k-1)+1];
    	
	  for(l = 0 ; l < 3 ; l++){
	    i = inxt3[i];
	    k1 = adja[i] / 4;
	    if(!k1) continue;

	    pt1 = &mesh->tetra[k1];
	  
	    if(pt1->flag == base) continue;
	    pt1->flag = base; 
	  
	    for(j=0 ; j<4 ; j++){
	      if(pt1->v[j] == nump)
		    break;
	    }
	    assert(j<4);
	  
	    /* overflow */
	    assert ( ilist <= LONMAX-3 );
	    list[ilist] = 4*k1+j; 
	    ilist++;	  	  
	  }
    cur++;
  }
  return(ilist); 
}
