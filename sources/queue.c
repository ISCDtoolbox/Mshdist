#include "mshdist.h"

extern unsigned char inxt2[5];
extern unsigned char inxt3[7];
extern char  ddb;


/* Create the queue structure for the Fast Marching method */
int setQueue(pMesh mesh,pQueue pq) {
  
  pq->sizmax   = mesh->np+1;
  pq->siz      = 0;
  
  /* memory allocation */
  pq->hp     = (Anod*)calloc(mesh->np+1,sizeof(Anod));
  pq->perm   = (int*)calloc(mesh->np+1,sizeof(int));
  
  return(1);
}

/* Free the memory for the queue */
int freeQueue(pQueue pq) {
  
  free(pq->hp);
  free(pq->perm);
  pq->sizmax = pq->siz = 0;
  
  return(1);
}

/* Left leaf of a node with index $k$ */
static inline int leftLeaf(int k) {
  return(2*k);
}

/* Right leaf of a node with index $k$ */
static inline int rightLeaf(int k) {
  return(2*k+1);
}

/* Root of a node with index $k$ */
static inline int root(int k) {
  return(k/2);
}

/* Update the value at node ip */
int upAnod(pQueue pq,int ip,double val) {
  pAnod    pa;
  int      k;
  
  k  = pq->perm[ip];
  pa = &pq->hp[k];
  
  if ( val < pa->d ) {
    pa->d = val;
    if ( !upPrio(pq,ip,val) ) {
      printf("Func. upAnod: impossible to change priority of %d. Abort\n",ip);
      exit(0);
    }
  }
  
  return(1);
}

/* Insert node ip */
int insertAnod(pQueue pq,int ip,double val) {
  pAnod     pa;
  int       k;
  
  k = ++pq->siz;
  if ( pq->siz >= pq->sizmax -1 ) {
    fprintf(stdout,"Func. insertAnod: impossible to insert new element. Abort\n");
    exit(0);
  }
  
  pq->perm[ip]  = k;
  pa            = &pq->hp[k];
  pa->indp      = ip;
  pa->d         = val;
  
  if ( !upPrio(pq,ip,val) ) {
    fprintf(stdout,"Func. insertAnod: impossible to change priority of %d. Abort\n",ip);
    exit(0);
  }
  
  return(1);
}

/* Increase the priority associated to node ip,
   which is possibly wrong, assuming the rest of the heap is well sorted */
int upPrio(pQueue pq,int ip,double val) {
  pAnod    par,pa,pa0;
  double   dr;
  int      k,kr,ipr,*key;
  
  key  = &pq->perm[0];
  pa0  = &pq->hp[0];
  
  k = key[ip];
  
  while ( k >= 2 ) {
    kr  = root(k);
    par = &pq->hp[kr];
    ipr = par->indp;
    dr  = par->d;
    
    if ( val >= dr ) break;
    
    /* Swap the active nodes between k and kr */
    pa  = &pq->hp[k];
    memcpy(pa0,par,sizeof(Anod));
    memcpy(par,pa,sizeof(Anod));
    memcpy(pa,pa0,sizeof(Anod));
    
    /* Update keys in the perm table */
    key[ipr] = k;
    key[ip] = kr;
    
    k = kr;
  }
  
  return(1);
}

/* Decrease the priority associated to node ip, 
   which is possibly wrong, assuming the rest of the heap is well sorted */
int downPrio(pQueue pq,int ip,double val) {
  pAnod   pa,pal,par,pa0;
  double  lval,rval;
  int     k,kl,kr,ipl,ipr,*key;
  
  key  = &pq->perm[0];
  pa0  = &pq->hp[0];
  k  = key[ip];
  kl = leftLeaf(k);
  kr = rightLeaf(k);
  
  while ( kl <= pq->siz ) {
        
    pa   = &pq->hp[k];
    pal  = &pq->hp[kl];
    lval = pal->d;
    ipl  = pal->indp;
    
    if ( kr > pq->siz ) {
      /* swap pa and pal */
      if ( val > lval ) {
        memcpy(pa0,pal,sizeof(Anod));
        memcpy(pal,pa,sizeof(Anod));
        memcpy(pa,pa0,sizeof(Anod));
      
        key[ip]  = kl;
        key[ipl] = k;
        k = kl;
      }
      else
        break;
    }
    else {
      par  = &pq->hp[kr];
      rval = par->d;
      ipr  = par->indp;
      
      /* swap pa and par */
      if ( rval < lval && val > rval ) {
        memcpy(pa0,par,sizeof(Anod));
        memcpy(par,pa,sizeof(Anod));
        memcpy(pa,pa0,sizeof(Anod));
        
        key[ip]  = kr;
        key[ipr] = k;
        k = kr;
      }
      /* swap pa and pal */
      else if ( lval <= rval && val > lval ) {
        memcpy(pa0,pal,sizeof(Anod));
        memcpy(pal,pa,sizeof(Anod));
        memcpy(pa,pa0,sizeof(Anod));
        
        key[ip]  = kl;
        key[ipl] = k;
        k = kl;
      }
      else
        break;
    }
    
    kl = leftLeaf(k);
    kr = rightLeaf(k);
  }
  
  return(1);
}

/* Get the value of the smallest active node */
int popAnod(pQueue pq,double *d) {
  pAnod          pa,paf;
  double         nval;
  int            ip,ipn,k;
  
  k = pq->siz;
  if ( !k ) return(0);
  
  /* Get the data from the root */
  pa   = &pq->hp[1];
  *d   = pa->d;
  ip   = pa->indp;
  
  /* Erase the root form the perm table */
  pq->perm[ip] = 0;
  
  /* If there is only the root, delete it */
  if ( k == 1 ) {
    memset(pa,0,sizeof(Anod));
    pq->siz = 0;
  }
  else {
    /* Copy last leaf in the root */
    paf  = &pq->hp[k];
    memcpy(pa,paf,sizeof(Anod));
    ipn  = pa->indp;
    nval = pa->d;
    pq->perm[ipn] = 1;
    
    memset(paf,0,sizeof(Anod));
    pq->siz--;
    
    if ( !downPrio(pq,ipn,nval) ) {
      printf("Func. popAnod: impossible to change priority of %d. Abort\n",1);
      exit(0);
    }
  }
  
  return(ip);
}

/* Check the order in the priority queue */
int checkHeap(pQueue pq) {
  pAnod    pa,pal;
  int      k,kr,kl;
  
  for (k=1; k<=pq->siz; k++) {
    pa = &pq->hp[k];
    kl = leftLeaf(k);
    kr = rightLeaf(k);
    if ( kl <= pq->siz ) {
      pal = &pq->hp[kl];
      if ( pal->d < pa->d ) {
        printf("Problem between numbers %d %d : %f vs %f \n",k,kl,pal->d,pa->d);
        return(0);
      }
    }
    if ( kr <= pq->siz ) {
      pal = &pq->hp[kr];
      if ( pal->d < pa->d ) {
        printf("Problem between numbers %d %d: %f vs %f \n",k,kr,pal->d,pa->d);
        return(0);
      }
    }
  }
  
  return(1);
}


