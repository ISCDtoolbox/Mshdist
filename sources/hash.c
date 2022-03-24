#include "mshdist.h"

#define KA     31
#define KB     57
#define KC     79
#define KTA     7
#define KTB    11

#define HA     7
#define HB    11
#define HC    13


unsigned char idir[5]     = {0,1,2,0,1};
unsigned char idirt[4][3] = { {1,2,3}, {0,3,2}, {0,1,3}, {0,2,1} };
hash hTab;

/* Identify whether ref corresponds to a reference of an interior subdomain */
inline int isIntDom(Info info,int ref) {
  int k;
  
  for (k=0; k<info.nintel; k++)
    if ( info.intel[k] == ref ) return(1);
  
  return(0);
}

/* Identify whether ref corresponds to a starting triangle */
inline int isStartTri(Info info,int ref) {
  int k;
  
  if ( info.nst ) {
    for (k=0; k<info.nst; k++)
      if ( info.st[k] == ref ) return(1);
  }
  
  return(0);
}

/* Identify whether ref corresponds to a starting edges */
inline int isStartEdg(Info info,int ref) {
  int k;
  
  if ( info.nsa ) {
    for (k=0; k<info.nsa; k++)
      if ( info.sa[k] == ref ) return(1);
  }
  
  return(0);
}

/* Identify whether ref corresponds to a starting vertex */
inline int isStartVer(Info info,int ref) {
  int k;
  
  if ( info.nsp ) {
    for (k=0; k<info.nsp; k++)
      if ( info.sp[k] == ref ) return(1);
  }
  
  return(0);
}

/* Create adjacency table in a 3d mesh */
int hashelt_3d(pMesh mesh) {
  pTetra    pt,pt1;
  int       k,kk,pp,l,ll,mins,mins1,maxs,maxs1,sum,sum1,iadr;
  int      *hcode,*link,inival,hsize;
  unsigned char   i,ii,i1,i2,i3;
  unsigned int    key;

  /* memory alloc */
  hcode = (int*)calloc(mesh->ne+1,sizeof(int));
  assert(hcode);
  link  = mesh->adja;
  hsize = mesh->ne;

  /* init */
  inival = 2147483647;
  for (k=0; k<=mesh->ne; k++)
    hcode[k] = -inival;

  /* build hash table */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !pt->v[0] )  continue;
    for (i=0; i<4; i++) {
      i1 = idirt[i][0];
      i2 = idirt[i][1];
      i3 = idirt[i][2];
      mins = D_MIN(pt->v[i1],pt->v[i2]);
      mins = D_MIN(mins,pt->v[i3]);
      maxs = D_MAX(pt->v[i1],pt->v[i2]);
      maxs = D_MAX(maxs,pt->v[i3]);

      /* compute key */
      sum = pt->v[i1] + pt->v[i2] + pt->v[i3];
      key = KA*mins + KB*maxs + KC*sum;
      key = key % hsize + 1;

      /* insert */
      iadr = 4*(k-1) + i+1;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }

  /* set adjacency */
  for (l=4*mesh->ne; l>0; l--) {
    if ( link[l] >= 0 )  continue;
    k = ((l-1) >> 2) + 1;
    i = (l-1) % 4;
    i1 = idirt[i][0];
    i2 = idirt[i][1];
    i3 = idirt[i][2];
    pt = &mesh->tetra[k];

    sum  = pt->v[i1] + pt->v[i2] + pt->v[i3];
    mins = D_MIN(pt->v[i1],pt->v[i2]);
    mins = D_MIN(mins,pt->v[i3]);
    maxs = D_MAX(pt->v[i1],pt->v[i2]);
    maxs = D_MAX(maxs,pt->v[i3]);

    /* accross link */
    ll = -link[l];
    pp = 0;
    link[l] = 0;
    while ( ll != inival ) {
      kk = ((ll-1) >> 2) + 1;
      ii = (ll-1) % 4;
      i1 = idirt[ii][0];
      i2 = idirt[ii][1];
      i3 = idirt[ii][2];
      pt1  = &mesh->tetra[kk];
      sum1 = pt1->v[i1] + pt1->v[i2] + pt1->v[i3];
      if ( sum1 == sum ) {
        mins1 = D_MIN(pt1->v[i1],pt1->v[i2]);
        mins1 = D_MIN(mins1,pt1->v[i3]);
        if ( mins1 == mins ) {
          maxs1 = D_MAX(pt1->v[i1],pt1->v[i2]);
          maxs1 = D_MAX(maxs1,pt1->v[i3]);
          if ( maxs1 == maxs ) {
            /* adjacent found */
            if ( pp != 0 )  link[pp] = link[ll];
            link[l] = 4*kk + ii;
            link[ll]= 4*k + i;
            break;
          }
        }
      }
      pp = ll;
      ll = -link[ll];
    }
  }

  free(hcode);
  return(1);
}

/* Create adjacency table in a 2d mesh: mesh->adja[3*(k-1)+1+i], i=0,1,2 = neighbors of k through vertex i */
int hashelt_2d(pMesh mesh) {
  pTria     pt,pt1;
  int       k,kk,pp,l,ll,mins,mins1,maxs,maxs1,iadr;
  int      *hcode,*link,inival,hsize;
  unsigned char   i,ii,i1,i2;
  unsigned int    key;
  
  /* memory alloc */
  hcode = (int*)calloc(mesh->nt+1,sizeof(int));
  assert(hcode);
  link  = mesh->adja;
  hsize = mesh->nt;

  /* init */
  inival = 2147483647;
  for (k=0; k<=mesh->nt; k++)
    hcode[k] = -inival;

  /* build hash table */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] )  continue;
    
    for (i=0; i<3; i++) {
      i1 = idir[i+1];
      i2 = idir[i+2];
      if ( pt->v[i1] < pt->v[i2] ) {
        mins = pt->v[i1];
        maxs = pt->v[i2];
      }
      else {
        mins = pt->v[i2];
        maxs = pt->v[i1];
      }

      /* compute key */
      key = KTA*mins + KTB*maxs;
      key = key % hsize + 1;

      /* insert */
      iadr = 3*(k-1) + i+1;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }

  /* set adjacency */
  for (l=3*mesh->nt; l>0; l--) {
    if ( link[l] >= 0 )  continue;
    k = (l-1) / 3 + 1;
    i = (l-1) % 3;
    i1 = idir[i+1];
    i2 = idir[i+2];
    pt = &mesh->tria[k];

    mins = D_MIN(pt->v[i1],pt->v[i2]);
    maxs = D_MAX(pt->v[i1],pt->v[i2]);

    /* accross link */
    ll = -link[l];
    pp = 0;
    link[l] = 0;
    while ( ll != inival ) {
      kk = (ll-1) / 3 + 1;
      ii = (ll-1) % 3;
      i1 = idir[ii+1];
      i2 = idir[ii+2];
      pt1  = &mesh->tria[kk];
      if ( pt1->v[i1] < pt1->v[i2] ) {
        mins1 = pt1->v[i1];
        maxs1 = pt1->v[i2];
      }
      else {
        mins1 = pt1->v[i2];
        maxs1 = pt1->v[i1];
      }
      
      if ( mins1 == mins  && maxs1 == maxs ) {
        /* adjacent found */
        if ( pp != 0 )  link[pp] = link[ll];
        link[l] = 3*kk + ii;
        link[ll]= 3*k + i;
        break;
      }
      pp = ll;
      ll = -link[ll];
    }
  }

  free(hcode);
  return(1);
}

/* Store all triangles of mesh with refs contained in info.sref in hash table hTab ; return number of
   successfully hashed triangles */
int hashTriaRef(Info info,pMesh mesh){ 
  pTria    ptt;
  hTria    *tab,*ph;
  int      k,n0,n1,n2,mins,maxs,sum,key,hnxt,l,nb;
  char     ier; 
  
  hTab.thsiz = (int)(0.51*mesh->nt);
  hTab.thmax = (int)(1.51*mesh->nt);
  hnxt = hTab.thsiz;
  
  nb = 0;
  
  hTab.ttab = (hTria*)calloc((int)(1.51*mesh->nt),sizeof(hTria));
  assert(hTab.ttab);
  tab = hTab.ttab;
  
  for(k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    
    for(l=1; l<=info.nsref; l++) {
      ier = (ptt->ref == info.sref[l]);
      if(ier) break;
    }
    if(!ier) continue;
    
    nb++;
    n0 = ptt->v[0];
    n1 = ptt->v[1];
    n2 = ptt->v[2];
    
    mins = D_MIN(n0,D_MIN(n1,n2));
    maxs = D_MAX(n0,D_MAX(n1,n2));
    sum = n0+n1+n2;
    
    key = ( KA*mins + KB*maxs )%hTab.thsiz;
    
    /* Uninitialized entry */
    if( tab[key].mins == 0 ){
      ph = &tab[key];
      ph->mins = mins;
      ph->maxs = maxs;
      ph->s  = sum;
      ph->k    = k; 
      ph->nxt  = 0;
    }
    else {
      while( tab[key].nxt ) {
        key = tab[key].nxt;
      }
      assert( hnxt < hTab.thmax ); 
      tab[key].nxt = hnxt;
      key = tab[key].nxt;
      ph = &tab[key];
      ph->mins = mins;
      ph->maxs = maxs;
      ph->s  = sum;
      ph->k    = k; 
      ph->nxt  = 0;
      hnxt++;
    }
  }
  
  return(nb);
}

/* Find triangle with keys mins, maxs, sum, and delete entry ; return 0 if not found */
int getTria(pMesh mesh,int mins,int maxs,int sum){
  hTria    *tab;
  int      key;
  
  tab = hTab.ttab;
  key = ( KA*mins + KB*maxs )%hTab.thsiz;
  
  if( tab[key].mins == 0 ){
    return(0);
  }
  else{
    while( (tab[key].mins != mins || tab[key].maxs != maxs || tab[key].s != sum)\
          && tab[key].nxt ) {
      key = tab[key].nxt;
    }
    
    /* If key is found, and positive, return it, else, if key is negative (i.e. it has already been 
       travelled ) return 0 */
    if( tab[key].mins == mins && tab[key].maxs == maxs && tab[key].s == sum ) {
      if( tab[key].k > 0 ) {
        tab[key].k *= -1;
        return(abs(tab[key].k));
      }
      else {
        return(0);
      }
    }
    else
      return(0);  
  } 
}

/* Free triangle hashing */
void delHash(pMesh mesh) {
  
  free(hTab.ttab);
  hTab.thsiz = 0;
  hTab.ttab = 0;
  
}

/* Extract points of a 3d mesh which are part of the surface triangulation, renumber them so that they are contiguous, and pack the values of sol accordingly */
int pack_s(pMesh mesh,pSol sol,int *perm) {
  pTria        ptt;
  pPoint       p0,p1;
  int          k,npc;
  char         i;
  
  /* Reset flag field */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;
  
  /* Identify points of the surface triangulation */
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !ptt->v[0] ) continue;
    for (i=0; i<3; i++){
      p0 = &mesh->point[ptt->v[i]];
      p0->flag = 1;
    }
  }
  
  /* Compress points: p1->flag contains the original index */
  npc = 0;
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    if ( p0->flag ) {
      npc++;
      p1 = &mesh->point[npc];
      memcpy(p1,p0,sizeof(Point));
      p1->flag = k;
    }
  }
  
  mesh->np = npc;
  
  /* Compress solution */
  sol->npi = sol->np;
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    sol->val[k] = sol->val[p0->flag];
  }
  
  sol->np = npc;
  
  /* Store permutation; perm[k] = original index of (new) point k */
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    perm[k] = p0->flag;
  }
  
  return(1);
}

/* Restore values of sol at the position pointed by perm */
int unpack_s(pMesh mesh,pSol sol,int *perm) {
  int k,ki;

  for (k=sol->np; k>=1; k--) {
    ki = perm[k];
    sol->val[ki] = sol->val[k];
  }
  
  sol->np = sol->npi;
  return(1);
}
