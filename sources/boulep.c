#include "mshdist.h"

extern unsigned char inxt3[7];

/* Return volumic ball (i.e. filled with tetrahedra) of point ip in tetra start. 
Results are stored under the form 4*kel + jel , kel = number of the tetra, jel = local 
index of p within kel */
int boulep (pMesh mesh, int start, int ip, int * list){
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
