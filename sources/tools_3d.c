#include "mshdist.h"

extern unsigned char inxt3[7];
extern char ddb;


/* Compute volume of tetra */
inline double volume(double *a,double *b,double *c,double *d) {
  double  ax,ay,az,bx,by,bz,vol;
  
  ax = b[0] - a[0];
  ay = b[1] - a[1];
  az = b[2] - a[2];
  
  bx = c[0] - a[0];
  by = c[1] - a[1];
  bz = c[2] - a[2];
  
  vol = (d[0]-a[0]) * (ay*bz - az*by) + (d[1]-a[1]) * (az*bx - ax*bz) \
  + (d[2]-a[2]) * (ax*by - ay*bx);
  
  return(fabs(vol) / 6.0);
}

/* invert 3x3 non-symmetric matrix */
int invmatg(double m[9],double mi[9]) {
  double  aa,bb,cc,det,vmin,vmax,maxx;
  int     k;
  
  /* check ill-conditionned matrix */
  vmin = vmax = fabs(m[0]);
  for (k=1; k<9; k++) {
    maxx = fabs(m[k]);
    if ( maxx < vmin )  vmin = maxx;
    else if ( maxx > vmax )  vmax = maxx;
  }
  if ( vmax == 0.0 )  return(0);
  
  /* compute sub-dets */
  aa = m[4]*m[8] - m[5]*m[7];
  bb = m[5]*m[6] - m[3]*m[8];
  cc = m[3]*m[7] - m[4]*m[6];
  det = m[0]*aa + m[1]*bb + m[2]*cc;
  if ( fabs(det) < EPS1 )  return(0);
  det = 1.0f / det;
  
  mi[0] = aa*det;
  mi[3] = bb*det;
  mi[6] = cc*det;
  mi[1] = (m[2]*m[7] - m[1]*m[8])*det;
  mi[4] = (m[0]*m[8] - m[2]*m[6])*det;
  mi[7] = (m[1]*m[6] - m[0]*m[7])*det;
  mi[2] = (m[1]*m[5] - m[2]*m[4])*det;
  mi[5] = (m[2]*m[3] - m[0]*m[5])*det;
  mi[8] = (m[0]*m[4] - m[1]*m[3])*det;
  
  return(1);
}

/* Store in circum the 3 coordinates of circumcenter of each surface triangle of mesh
 and the associated (squared) circumradius. circum is a table of length 4* mesh->ntria +1*/
int buildcircum_3d(pMesh mesh, double *circum){
  int i;
  pTria pt;
  pPoint p0,p1,p2;
  double *circ;	
	
  for(i=1;i<=mesh->nt;i++){
    pt = &mesh->tria[i];
    p0 = &mesh->point[pt->v[0]];  
    p1 = &mesh->point[pt->v[1]];  
    p2 = &mesh->point[pt->v[2]];  
    circ = &circum[4*(i-1)+1];
    
    circumcoords(p0,p1,p2,circ);  
	  
  }
	
  return(1);
}

/* Compute [p,q,r,s] := det(q-p, r-p, s-p) */
double detOrient_3d(pPoint p,pPoint q,pPoint r,pPoint s) {
  double m11, m12, m13, m21, m22, m23, m31, m32, m33;
  
  m11 = q->c[0] - p->c[0]; 
  m21 = q->c[1] - p->c[1]; 
  m31 = q->c[2] - p->c[2]; 
  
  m12 = r->c[0] - p->c[0]; 
  m22 = r->c[1] - p->c[1]; 
  m32 = r->c[2] - p->c[2]; 
  
  m13 = s->c[0] - p->c[0]; 
  m23 = s->c[1] - p->c[1]; 
  m33 = s->c[2] - p->c[2];
  
  return(m11*m22*m33 + m12*m23*m31 + m13*m21*m32 - m31*m22*m13 - m32*m23*m11 - m33*m21*m12); 
}


/* Compute det(q-p, r-p, v)*/
double determinant3Pts1Vct_3d(pPoint p, pPoint q, pPoint r, pPoint v) {
  double m11, m21, m31, m12, m22, m32;
  
  m11 = q->c[0] - p->c[0]; 
  m21 = q->c[1] - p->c[1]; 
  m31 = q->c[2] - p->c[2]; 
  
  m12 = r->c[0] - p->c[0]; 
  m22 = r->c[1] - p->c[1]; 
  m32 = r->c[2] - p->c[2];
  
  return(m11*m22*v->c[2] + m12*v->c[1]*m31 + v->c[0]*m21*m32 - m31*m22*v->c[0] - m32*v->c[1]*m11 - v->c[2]*m21*m12);
}


/* Store in (circ[0],circ[1],circ[2]) the coordinates of the circumcenter of tria p0p1p2
 and in circ[3] its squared circumradius */
int circumcoords(pPoint p0, pPoint p1, pPoint p2, double *circ){
  double r;  
  
  circ[0] = 0.33333*(p0->c[0]+p1->c[0]+p2->c[0]);
  circ[1] = 0.33333*(p0->c[1]+p1->c[1]+p2->c[1]); 
  circ[2] = 0.33333*(p0->c[2]+p1->c[2]+p2->c[2]);
	
  circ[3] = (p0->c[0]-circ[0])*(p0->c[0]-circ[0])+(p0->c[1]-circ[1])*(p0->c[1]-circ[1])+\
  (p0->c[2]-circ[2])*(p0->c[2]-circ[2]);  
	
  r =  (p1->c[0]-circ[0])*(p1->c[0]-circ[0])+(p1->c[1]-circ[1])*(p1->c[1]-circ[1])+\
	(p1->c[2]-circ[2])*(p1->c[2]-circ[2]); 	
  
  circ[3] = D_MAX(circ[3], r);	
	
  r =  (p2->c[0]-circ[0])*(p2->c[0]-circ[0])+(p2->c[1]-circ[1])*(p2->c[1]-circ[1])+\
	(p2->c[2]-circ[2])*(p2->c[2]-circ[2]);	
	
  circ[3] = D_MAX(circ[3], r);	
  return(1);
	
  /*int     ier;
   double  a,b,c,norm, det1,det2,det3;
   double  m[9],mi[9];	
   
   //computation of the circumradius
   
   a = (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0]) + (p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1])\
   + (p2->c[2]-p1->c[2])*(p2->c[2]-p1->c[2]);
   b = (p1->c[0]-p0->c[0])*(p1->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])*(p1->c[1]-p0->c[1])\
   + (p1->c[2]-p0->c[2])*(p1->c[2]-p0->c[2]);
   c = (p2->c[0]-p0->c[0])*(p2->c[0]-p0->c[0]) + (p2->c[1]-p0->c[1])*(p2->c[1]-p0->c[1])\
   + (p2->c[2]-p0->c[2])*(p2->c[2]-p0->c[2]);  
   
   det1 = (p1->c[1]-p0->c[1])*(p2->c[2]-p0->c[2]) - (p1->c[2]-p0->c[2])*(p2->c[1]-p0->c[1]) ;   
   det2 = (p1->c[2]-p0->c[2])*(p2->c[0]-p0->c[0]) - (p1->c[0]-p0->c[0])*(p2->c[2]-p0->c[2]) ;   
   det3 = (p1->c[0]-p0->c[0])*(p2->c[1]-p0->c[1]) - (p1->c[1]-p0->c[1])*(p2->c[0]-p0->c[0]) ;   
   norm = det1*det1 + det2*det2 + det3*det3;
   
   assert(norm >0.0);
   
   circ[3] = 0.25*a*b*c/norm;
   
   //computation of the coordinates of the circumcenter
   m[0] = p1->c[0] - p0->c[0];  
   m[1] = p1->c[1] - p0->c[1];  
   m[2] = p1->c[2] - p0->c[2];  
   m[3] = p2->c[0] - p0->c[0];  
   m[4] = p2->c[1] - p0->c[1];  
   m[5] = p2->c[2] - p0->c[2];
   m[6] = (p1->c[1]-p0->c[1])*(p2->c[2]-p0->c[2]) - (p1->c[2]-p0->c[2])*(p2->c[1]-p0->c[1]);  
   m[7] = -(p1->c[0]-p0->c[0])*(p2->c[2]-p0->c[2]) + (p1->c[2]-p0->c[2])*(p2->c[0]-p0->c[0]);  
   m[8] = (p1->c[0]-p0->c[0])*(p2->c[1]-p0->c[1]) - (p1->c[1]-p0->c[1])*(p2->c[0]-p0->c[0]);
   ier = invmatg(m, mi);
   //assert(ier);
   if(!ier){
   circ[0] = 0.33*(p0->c[0] + p1->c[0] + p2->c[0]); 
   circ[1] = 0.33*(p0->c[1] + p1->c[1] + p2->c[1]);  
   circ[2] = 0.33*(p0->c[2] + p1->c[2] + p2->c[2]);
   printf("coucou \n"); 
   return(1);   
   }	
   
   a = 0.5*((p1->c[0])*(p1->c[0]) + (p1->c[1])*(p1->c[1]) + (p1->c[2])*(p1->c[2]) \
   - (p0->c[0])*(p0->c[0]) - (p0->c[1])*(p0->c[1]) - (p0->c[2])*(p0->c[2]));
   
   b = 0.5*((p2->c[0])*(p2->c[0]) + (p2->c[1])*(p2->c[1]) + (p2->c[2])*(p2->c[2]) \
   - (p0->c[0])*(p0->c[0]) - (p0->c[1])*(p0->c[1]) - (p0->c[2])*(p0->c[2]));
   
   c = m[6]*p0->c[0] + m[7]*p0->c[1] + m[8]*p0->c[2];  
   
   circ[0] = mi[0]*a + mi[1]*b + mi[2]*c; 
   circ[1] = mi[3]*a + mi[4]*b + mi[5]*c;  
   circ[2] = mi[6]*a + mi[7]*b + mi[8]*c;  
   
   return(1);*/
}


/* check intersection triangles p1q1r1 with p2q2e2 
 return 0: no intersection, 1: intersection */
int intersec_3d(pPoint p1,pPoint q1,pPoint r1,pPoint p2,pPoint q2,pPoint r2 ) {
  Point    prodVect; 
  pPoint   p1Prime, q1Prime, r1Prime,  p2Prime, q2Prime, r2Prime;
  double   orientp2, orientq2, orientr2, orientp1, orientq1, orientr1;
  int      renum1, renum2;
  
  // Un premier test permet de positionner un triangle par rapport au plan decrit par l autre
  // On commence par regarder les points du triangle p2q2r2 par rapport a ceux de p1q1r1
  
  orientp2 = detOrient_3d(p1, q1 , r1, p2); 
  orientq2 = detOrient_3d(p1, q1 , r1, q2);
  orientr2 = detOrient_3d(p1, q1 , r1, r2);
  
  // de meme, on oriente les points de T1 par rapport au plan de T2
  
  orientp1 = detOrient_3d(p2, q2 , r2, p1); 
  orientq1 = detOrient_3d(p2, q2 , r2, q1);
  orientr1 = detOrient_3d(p2, q2 , r2, r1);
  
  // si p2, q2 et r2 sont strictement du meme cote du plan decrit par T1, aucune intersection n est possible
  
  if(((orientp2>0.0)&&(orientq2>0.0)&&(orientr2>0.0))||((orientp2<0.0)&&(orientq2<0.0)&&(orientr2<0.0)))  //prendre une marge a epsilon pres ici
  { return(0);}
  
  // et symetriquement
  
  if(((orientp1>0.0)&&(orientq1>0.0)&&(orientr1>0.0))||((orientp1<0.0)&&(orientq1<0.0)&&(orientr1<0.0)))
  { return(0);}  
  
  //Deux autres cas particuliers sont a traiter : tout d abord celui ou les deux triangles sont coplanaires...
  
  if((fabs(orientp2)<0.0)&&(fabs(orientq2)<0.0)&&(fabs(orientr2)<0.0))   //prendre une marge aussi.
  {
    return(1);       //attention, bien sur, il n est pas toujours vrai de dire que si les deux triangles sont coplanaires, ils s intersectent necessairement !
    // mais dnas le contexte dans lequel on va utiliser cette fonction, ce sera toujours le cas
  }
  
  // et le cas ou seul l'un des trois points p1,q1, r1 appartient au plan de T2 et ou les deux autres sont du meme cote de ce plan (et son symetrique)
  // Dans ce cas, il y aura intersection ssi le point en question appartient au triangle, ceci se testant au moyen de trois tests d orientation
  
  if((fabs(orientp1)< EPS1)&&(orientq1*orientr1>=0.0))  //prendre une marge a epsilon pres pour le test ==
  {
    //On calcule le produit vectoriel p2q2 /wedge p2r2 : il y a surement des optimisations a faire avec le calcul de determinants d avant
    
    prodVect.c[0] = (q2->c[1] - p2->c[1])*(r2->c[2] - p2->c[2]) - (q2->c[2] - p2->c[2])*(r2->c[1] - p2->c[1]);
    prodVect.c[1] = (q2->c[2] - p2->c[2])*(r2->c[0] - p2->c[0]) - (q2->c[0] - p2->c[0])*(r2->c[2] - p2->c[2]);
    prodVect.c[2] = (q2->c[0] - p2->c[0])*(r2->c[1] - p2->c[1]) - (q2->c[1] - p2->c[1])*(r2->c[0] - p2->c[0]);
    
    if((determinant3Pts1Vct_3d(p1, p2, q2, & prodVect) >=0.0)&&(determinant3Pts1Vct_3d(p1, q2, r2, & prodVect) >=0.0)&&(determinant3Pts1Vct_3d(p1, r2, p2, & prodVect) >=0.0))
    {return(1);}
    
    else {return(0);}
    
  }
  
  if((fabs(orientq1) < EPS1)&&(orientp1*orientr1>=0.0))  
  {
    //exactement la meme chose ici, sauf que l on regarde le point q1
    
    prodVect.c[0] = (q2->c[1] - p2->c[1])*(r2->c[2] - p2->c[2]) - (q2->c[2] - p2->c[2])*(r2->c[1] - p2->c[1]);
    prodVect.c[1] = (q2->c[2] - p2->c[2])*(r2->c[0] - p2->c[0]) - (q2->c[0] - p2->c[0])*(r2->c[2] - p2->c[2]);
    prodVect.c[2] = (q2->c[0] - p2->c[0])*(r2->c[1] - p2->c[1]) - (q2->c[1] - p2->c[1])*(r2->c[0] - p2->c[0]);
    
    if((determinant3Pts1Vct_3d(q1, p2, q2, & prodVect) >=0.0)&&(determinant3Pts1Vct_3d(q1, q2, r2, & prodVect) >=0.0)&&(determinant3Pts1Vct_3d(q1, r2, p2, & prodVect) >=0.0))
    {return 1;}
    
    else {return 0;}
  } 
  
  if((fabs(orientr1) < EPS1)&&(orientq1*orientp1>=0.0))  
  {   
    // et le point r1
    prodVect.c[0] = (q2->c[1] - p2->c[1])*(r2->c[2] - p2->c[2]) - (q2->c[2] - p2->c[2])*(r2->c[1] - p2->c[1]);
    prodVect.c[1] = (q2->c[2] - p2->c[2])*(r2->c[0] - p2->c[0]) - (q2->c[0] - p2->c[0])*(r2->c[2] - p2->c[2]);
    prodVect.c[2] = (q2->c[0] - p2->c[0])*(r2->c[1] - p2->c[1]) - (q2->c[1] - p2->c[1])*(r2->c[0] - p2->c[0]);
    
    if((determinant3Pts1Vct_3d(r1, p2, q2, & prodVect) >=0.0)&&(determinant3Pts1Vct_3d(r1, q2, r2, & prodVect) >=0.0)&&(determinant3Pts1Vct_3d(r1, r2, p2, & prodVect) >=0.0))
    {return 1;}
    
    else {return 0;}      
  }
  
  // on fait la meme chose pour les points du triangle T2
  
  if((fabs(orientp2) < EPS1)&&(orientq2*orientr2>=0.0))  //prendre une marge a epsilon pres pour le test ==
  {
    //On calcule le produit vectoriel p1q1 /wedge p1r1 : il y a surement des optimisations a faire avec le calcul de determinants d avant
    
    prodVect.c[0] = (q1->c[1] - p1->c[1])*(r1->c[2] - p1->c[2]) - (q1->c[2] - p1->c[2])*(r1->c[1] - p1->c[1]);
    prodVect.c[1] = (q1->c[2] - p1->c[2])*(r1->c[0] - p1->c[0]) - (q1->c[0] - p1->c[0])*(r1->c[2] - p1->c[2]);
    prodVect.c[2] = (q1->c[0] - p1->c[0])*(r1->c[1] - p1->c[1]) - (q1->c[1] - p1->c[1])*(r1->c[0] - p1->c[0]);
    
    if((determinant3Pts1Vct_3d(p2, p1, q1, & prodVect) >=0.0)&&(determinant3Pts1Vct_3d(p2, q1, r1, & prodVect) >=0.0)&&(determinant3Pts1Vct_3d(p2, r1, p1, & prodVect) >=0.0))
    {return(1);}
    
    else {return(0);}      
  }
  
  if((fabs(orientq2) < EPS1)&&(orientp2*orientr2>=0.0))  
  {
    prodVect.c[0] = (q1->c[1] - p1->c[1])*(r1->c[2] - p1->c[2]) - (q1->c[2] - p1->c[2])*(r1->c[1] - p1->c[1]);
    prodVect.c[1] = (q1->c[2] - p1->c[2])*(r1->c[0] - p1->c[0]) - (q1->c[0] - p1->c[0])*(r1->c[2] - p1->c[2]);
    prodVect.c[2] = (q1->c[0] - p1->c[0])*(r1->c[1] - p1->c[1]) - (q1->c[1] - p1->c[1])*(r1->c[0] - p1->c[0]);
    
    if((determinant3Pts1Vct_3d(q2, p1, q1, & prodVect) >=0.0)&&(determinant3Pts1Vct_3d(q2, q1, r1, & prodVect) >=0.0)&&(determinant3Pts1Vct_3d(q2, r1, p1, & prodVect) >=0.0))
    {return(1);}
    
    else {return(0);}
    
  } 
  
  if((fabs(orientr2) < EPS1)&&(orientq2*orientp2>=0.0))  
  {
    prodVect.c[0] = (q1->c[1] - p1->c[1])*(r1->c[2] - p1->c[2]) - (q1->c[2] - p1->c[2])*(r1->c[1] - p1->c[1]);
    prodVect.c[1] = (q1->c[2] - p1->c[2])*(r1->c[0] - p1->c[0]) - (q1->c[0] - p1->c[0])*(r1->c[2] - p1->c[2]);
    prodVect.c[2] = (q1->c[0] - p1->c[0])*(r1->c[1] - p1->c[1]) - (q1->c[1] - p1->c[1])*(r1->c[0] - p1->c[0]);
    
    if((determinant3Pts1Vct_3d(r2, p1, q1, & prodVect) >=0.0)&&(determinant3Pts1Vct_3d(r2, q1, r1, & prodVect) >=0.0)&&(determinant3Pts1Vct_3d(r2, r1, p1, & prodVect) >=0.0))
    {return(1);}
    
    return(0);
    
  } 
  
  // On a donc ecarte les cas ou les deux triangles sont coplanaires, ou l'un des deux est strictement d un seul cote du plan 
  // decrit par l autre, et ou l un des deux a un point dans le plan de l autre, et les deux autres points sont du meme cote.
  // Maintenant, on renumerote les points de sorte que p1 et q1 soient seuls du cote du plan de l autre triangle
  // les entiers renum serviront plus loin, pour marquer dans quel cas de figure on se trouve sans avoir besoin de refaire des tests compliques
  
  p1Prime = p2Prime = q1Prime = q2Prime = r1Prime = r2Prime = 0;
  
  if(((orientp1<0.0)&&(orientq1>=0.0)&&(orientr1>=0.0))||((orientp1>0.0)&&(orientq1<=0.0)&&(orientr1<=0.0)))  //prendre une marge a epsilon pres pour le test >
  {
    p1Prime = p1; 
    renum1 =1;
    
  }
  
  if(((orientq1<0.0)&&(orientp1>=0.0)&&(orientr1>=0.0))||((orientq1>0.0)&&(orientp1<=0.0)&&(orientr1<=0.0)))  //prendre une marge a epsilon pres pour le test >
  {
    p1Prime = q1; 
    renum1 =2 ;
    
  } 
  
  if(((orientr1<0.0)&&(orientq1>=0.0)&&(orientp1>=0.0))||((orientr1>0.0)&&(orientq1<=0.0)&&(orientp1<=0.0)))  //prendre une marge a epsilon pres pour le test >
  {
    p1Prime = r1; 
    renum1 = 3;
  }
  
  // on fait la meme chose pour les points du triangle T2
  
  if(((orientp2<0.0)&&(orientq2>=0.0)&&(orientr2>=0.0))||((orientp2>0.0)&&(orientq2<=0.0)&&(orientr2<=0.0)))  //prendre une marge a epsilon pres pour le test >
  {
    p2Prime = p2; 
    renum2 = 1;
    
  }
  
  if(((orientq2<0.0)&&(orientp2>=0.0)&&(orientr2>=0.0))||((orientq2>0.0)&&(orientp2<=0.0)&&(orientr2<=0.0)))  //prendre une marge a epsilon pres pour le test >
  {
    p2Prime = q2; 
    renum2 = 2;
  } 
  
  if(((orientr2<0.0)&&(orientq2>=0.0)&&(orientp2>=0.0))||((orientr2>0.0)&&(orientq2<=0.0)&&(orientp2<=0.0)))  //prendre une marge a epsilon pres pour le test >
  {
    p2Prime = r2; 
    renum2 = 3;
  }
  
  // Il faut alors pratiquer encore d autres tests et renumerotations pour s assurer que p2q2r2 tournent dans le sens trigo par rapport a p1 et reciproquement
  
  if((renum1 ==1)&&(renum2==1))
  {
    if(orientp2 >=0.0)
    {
      q1Prime = q1;
      r1Prime = r1;
    }
    
    else
    {
      q1Prime = r1;
      r1Prime = q1;
    }
    
    if(orientp1 >=0.0)
    {
      q2Prime = q2;
      r2Prime = r2;
    }
    
    else
    {
      q2Prime = r2;
      r2Prime = q2;
    }
  }
  
  if((renum1 ==1)&&(renum2==2))
  {
    if(orientq2 >=0.0)
    {
      q1Prime = q1;
      r1Prime = r1;
    }
    
    else
    {
      q1Prime = r1;
      r1Prime = q1;
    }
    
    if(orientp1 >=0.0)   //dans ce cas, il faut preserver le sens de rotation des points p2q2r2 dans l ecriture avec des "prime"
    {
      q2Prime = r2;
      r2Prime = p2;
    }
    
    else
    {
      q2Prime = p2;
      r2Prime = r2;
    }
  }
  
  if((renum1 ==1)&&(renum2==3))
  {
    if(orientr2 >=0.0)
    {
      q1Prime = q1;
      r1Prime = r1;
    }
    
    else
    {
      q1Prime = r1;
      r1Prime = q1;
    }
    
    if(orientp1 >=0.0)
    {
      q2Prime = p2;
      r2Prime = q2;
    }
    
    else
    {
      q2Prime = q2;
      r2Prime = p2;
    }
  }
  
  if((renum1 ==2)&&(renum2==1))
  {
    if(orientp2 >=0.0)
    {
      q1Prime = r1;
      r1Prime = p1;
    }
    
    else
    {
      q1Prime = p1;
      r1Prime = r1;
    }
    
    if(orientq1 >=0.0)
    {
      q2Prime = q2;
      r2Prime = r2;
    }
    
    else
    {
      q2Prime = r2;
      r2Prime = q2;
    }
  }
  
  if((renum1 ==2)&&(renum2==2))
  {
    if(orientq2 >=0.0)
    {
      q1Prime = r1;
      r1Prime = p1;
    }
    
    else
    {
      q1Prime = p1;
      r1Prime = r1;
    }
    
    if(orientq1 >=0.0)   //dans ce cas, il faut preserver le sens de rotation des points p2q2r2 dans l ecriture avec des "prime"
    {
      q2Prime = r2;
      r2Prime = p2;
    }
    
    else
    {
      q2Prime = p2;
      r2Prime = r2;
    }
  }
  
  if((renum1 ==2)&&(renum2==3))
  {
    if(orientr2 >=0.0)
    {
      q1Prime = r1;
      r1Prime = p1;
    }
    
    else
    {
      q1Prime = p1;
      r1Prime = r1;
    }
    
    if(orientq1 >=0.0)
    {
      q2Prime = p2;
      r2Prime = q2;
    }
    
    else
    {
      q2Prime = q2;
      r2Prime = p2;
    }
  }
  
  if((renum1 ==3)&&(renum2==1))
  {
    if(orientp2 >=0.0)
    {
      q1Prime = p1;
      r1Prime = q1;
    }
    
    else
    {
      q1Prime = q1;
      r1Prime = p1;
    }
    
    if(orientr1 >=0.0)
    {
      q2Prime = q2;
      r2Prime = r2;
    }
    
    else
    {
      q2Prime = r2;
      r2Prime = q2;
    }
  }
  
  if((renum1 ==3)&&(renum2==2))
  {
    if(orientq2 >=0.0)
    {
      q1Prime = p1;
      r1Prime = q1;
    }
    
    else
    {
      q1Prime = q1;
      r1Prime = p1;
    }
    
    if(orientr1 >=0.0)   //dans ce cas, il faut preserver le sens de rotation des points p2q2r2 dans l ecriture avec des "prime"
    {
      q2Prime = r2;
      r2Prime = p2;
    }
    
    else
    {
      q2Prime = p2;
      r2Prime = r2;
    }
  }
  
  if((renum1 ==3)&&(renum2==3))
  {
    if(orientr2 >=0.0)
    {
      q1Prime = p1;
      r1Prime = q1;
    }
    
    else
    {
      q1Prime = q1;
      r1Prime = p1;
    }
    
    if(orientr1 >=0.0)
    {
      q2Prime = p2;
      r2Prime = q2;
    }
    
    else
    {
      q2Prime = q2;
      r2Prime = p2;
    }
  }
  
  /* Conclusion */
  assert(p1Prime);
  assert(p2Prime);
  assert(q1Prime);
  assert(q2Prime);
  assert(r1Prime);
  assert(r2Prime);
  
  return ((detOrient_3d(p1Prime, q1Prime, p2Prime, q2Prime)<=0.0)&&(detOrient_3d(p1Prime, r1Prime, r2Prime, p2Prime)<=0.0));   //idem, prendre une marge avec un eps>0
}

/* return distance from point pa to segment (p1p2)
 vertex tag set to 2 if distance is realized by p1 or p2 */
double distpt_23d(pPoint p1,pPoint p2,pPoint pa) {
  double   a,b,c,d,dd,ux,uy,vx,vy,wx,wy,xp,yp;
  
  a = p1->c[1] - p2->c[1];
  b = p2->c[0] - p1->c[0];
  c = -b*p1->c[1] - a*p1->c[0];
  d = INIVAL_3d;
  
  dd = a*a + b*b;
  if ( dd < EPS1 ) {
    d = (pa->c[0]-p1->c[0])*(pa->c[0]-p1->c[0]) + (pa->c[1]-p1->c[1])*(pa->c[1]-p1->c[1]);
    return(d);
  }
  xp =  b*b * pa->c[0] - a*b * pa->c[1] - a*c;
  yp = -a*b * pa->c[0] + a*a * pa->c[1] - b*c;
  dd = 1.0 / dd;
  xp *= dd;
  yp *= dd;
  
  ux = xp - p1->c[0];
  uy = yp - p1->c[1];
  vx = xp - p2->c[0];
  vy = yp - p2->c[1];
  wx = p2->c[0] - p1->c[0];
  wy = p2->c[1] - p1->c[1];
  
  if ( fabs(b) < EPS1 ) {
    if ( uy*wy <= 0.0 ) {
      d = (pa->c[0]-p1->c[0])*(pa->c[0]-p1->c[0]) + (pa->c[1]-p1->c[1])*(pa->c[1]-p1->c[1]);
    }
    else if ( vy*wy <= 0.0 )
      d = (pa->c[0]-xp)*(pa->c[0]-xp) + (pa->c[1]-yp)*(pa->c[1]-yp);
    else {
      d = (pa->c[0]-p2->c[0])*(pa->c[0]-p2->c[0]) + (pa->c[1]-p2->c[1])*(pa->c[1]-p2->c[1]);
    }
  }
  else {
    if ( ux*wx <= 0.0 ) {
      d = (pa->c[0]-p1->c[0])*(pa->c[0]-p1->c[0]) + (pa->c[1]-p1->c[1])*(pa->c[1]-p1->c[1]);
    }
    else if ( vx*wx <= 0.0 )
      d = (pa->c[0]-xp)*(pa->c[0]-xp) + (pa->c[1]-yp)*(pa->c[1]-yp);
    else {
      d = (pa->c[0]-p2->c[0])*(pa->c[0]-p2->c[0]) + (pa->c[1]-p2->c[1])*(pa->c[1]-p2->c[1]);
    }
  }
  
  return(d);
}

/* Compute (squared) distance from point pa to segment (p0p1) in 3d */
double distPtSeg_3d(pPoint p0, pPoint p1, pPoint pa){
  double lx1,ly1,lz1,lxa,lya,lza,longp0p1,cosAlpha,sinAlpha;
  double m11,m12,m13,m21,m22,m23,m31,m32,m33;	
  double p0XTemp,p0YTemp,p0ZTemp,p1XTemp,p1YTemp,p1ZTemp,paXTemp,paYTemp,paZTemp;
  double aa,bb,ab,ll,l;	
  
  lx1 = p1->c[0] - p0->c[0];
  ly1 = p1->c[1] - p0->c[1];  
  lz1 = p1->c[2] - p0->c[2];
  
  lxa = pa->c[0] - p0->c[0];
  lya = pa->c[1] - p0->c[1];
  lza = pa->c[2] - p0->c[2];
  longp0p1 = sqrt(lx1*lx1 + ly1*ly1 + lz1*lz1); 
  
  /* Case when points are too close from each other */
  if(longp0p1 <EPS1) 
    return((pa->c[0]-p0->c[0])*(pa->c[0]-p0->c[0])\
        + (pa->c[1]-p0->c[1])*(pa->c[1]-p0->c[1]) + (pa->c[2]-p0->c[2])*(pa->c[2]-p0->c[2])); 	
  
  cosAlpha = lz1/longp0p1; 
  sinAlpha = sqrt(1.0-cosAlpha*cosAlpha); 
	
  /* Apply translation of vector p0p1 */
  p0XTemp = 0.0;
  p0YTemp = 0.0;
  p0ZTemp = 0.0; 
  p1XTemp = lx1; 
  p1YTemp = ly1; 
  p1ZTemp = lz1;
  paXTemp = lxa;
  paYTemp = lya;	
  paZTemp = lza;
	
  aa = lx1*lx1;
  bb = ly1*ly1; 
  ab = lx1*ly1;
  ll = aa + bb;
  l = sqrt(ll);	
	
  /* Special case : vector p0p1 is already vertical */	
  if(ll<EPS1){       
    if(p1ZTemp <= 0.0){     
      p1ZTemp = -p1ZTemp;
      paZTemp = -paZTemp;  
    }
  }	
	
  else{      	  
    m11 = (aa*cosAlpha + bb)/ll;
    m12 = (ab*cosAlpha - ab )/ll;
    m13 = -lx1*sinAlpha/l;
	  
    m21 = (ab * cosAlpha - ab)/ll;
    m22 = (bb*cosAlpha + aa)/ll;
    m23 = - ly1 * sinAlpha/l; 
	  
    m31 = (lx1* sinAlpha)/l;
    m32 = (ly1 * sinAlpha)/l;
    m33 = cosAlpha;   
	  
	  
    p1XTemp = m11*lx1 + m12*ly1 + m13 * lz1;  
    p1YTemp = m21*lx1 + m22*ly1 + m23 * lz1;
    p1ZTemp = longp0p1;   
	  
    paXTemp = m11*lxa + m12*lya + m13 * lza;
    paYTemp = m21*lxa + m22*lya + m23 * lza;  
    paZTemp = m31*lxa + m32*lya + m33 * lza;
  }	
  
  /* At this point, p0Temp = 0.° , p1Temp is in axis 0z, making projection of paTemp easier */
  if(paZTemp<0.0){
    return((paXTemp)*(paXTemp)+(paYTemp)*(paYTemp)+(paZTemp)*(paZTemp));  
  }
	
  if(paZTemp>longp0p1){
    return((paXTemp)*(paXTemp)+(paYTemp)*(paYTemp)+(paZTemp-longp0p1)*(paZTemp-longp0p1));  
  }	
	
  else{
    return((paXTemp)*(paXTemp)+(paYTemp)*(paYTemp));  
  }
}

/* Compute squared distance from pq to tria p0p1p2, 
 proj = 1: projection onto face, =2: distance to vertex or edge */
double distpt_3d(pPoint p0,pPoint p1,pPoint p2,pPoint pq,char *proj) {
  Point pointPlan0, pointPlan1, pointPlan2, pointPlana; 
  double lx1, ly1, lz1, lx2, ly2, lz2, lxq, lyq, lzq, longp0p1, cosAlpha, sinAlpha, aa, bb, ab, ll, l; 
  double p0XTemp, p0YTemp, p0ZTemp, p1XTemp, p1YTemp, p1ZTemp, p2XTemp, p2YTemp, p2ZTemp, pqXTemp, pqYTemp, pqZTemp; // pour stocker temporairement
  double m11, m12, m13, m21, m22, m23, m31, m32, m33, d01, d12, d02, dTmp;
  double zone01, zone02, zone12;
  
  *proj = 1;
  lx1 = p1->c[0] - p0->c[0]; 
  ly1 = p1->c[1] - p0->c[1];  
  lz1 = p1->c[2] - p0->c[2];  
  
  lx2 = p2->c[0] - p0->c[0]; 
  ly2 = p2->c[1] - p0->c[1];  
  lz2 = p2->c[2] - p0->c[2];    
  
  lxq = pq->c[0] - p0->c[0]; 
  lyq = pq->c[1] - p0->c[1]; 
  lzq = pq->c[2] - p0->c[2];
    
  longp0p1 = sqrt(lx1*lx1 + ly1*ly1 + lz1*lz1);
  
  cosAlpha = lz1/longp0p1; 
  sinAlpha = sqrt(1.0-cosAlpha*cosAlpha); 
    
  p0XTemp = 0.0;
  p0YTemp = 0.0;
  p0ZTemp = 0.0; 
  p1XTemp = lx1; 
  p1YTemp = ly1; 
  p1ZTemp = lz1;
  p2XTemp = lx2; 
  p2YTemp = ly2; 
  p2ZTemp = lz2;  
  pqXTemp = lxq; 
  pqYTemp = lyq; 
  pqZTemp = lzq;
  
  /* Apply rotation with angle alpha */
  
  aa = lx1*lx1;
  bb = ly1*ly1; 
  ab = lx1*ly1;
  ll = aa + bb;
  l = sqrt(ll);
  
  if ( ll<EPS1 ){
    if ( p1ZTemp <= 0.0 ){
      p1ZTemp = -p1ZTemp;
      p2ZTemp = -p2ZTemp; 
      pqZTemp = -pqZTemp;   
    }
  }
  
  else{         
    m11 = (aa*cosAlpha + bb)/ll;
    m12 = (ab*cosAlpha - ab )/ll;
    m13 = -lx1*sinAlpha/l;
    
    m21 = (ab * cosAlpha - ab)/ll;
    m22 = (bb*cosAlpha + aa)/ll;
    m23 = - ly1 * sinAlpha/l; 
    
    m31 = (lx1* sinAlpha)/l;
    m32 = (ly1 * sinAlpha)/l;
    m33 = cosAlpha;   
        
    p1XTemp = m11*lx1 + m12*ly1 + m13 * lz1;   
    p1YTemp = m21*lx1 + m22*ly1 + m23 * lz1;
    p1ZTemp = longp0p1; 
    
    p2XTemp = m11*lx2 + m12*ly2 + m13 * lz2;
    p2YTemp = m21*lx2 + m22*ly2 + m23 * lz2;  
    p2ZTemp = m31*lx2 + m32*ly2 + m33 * lz2;  
    
    pqXTemp = m11*lxq + m12*lyq + m13 * lzq;
    pqYTemp = m21*lxq + m22*lyq + m23 * lzq;  
    pqZTemp = m31*lxq + m32*lyq + m33 * lzq;
  }
  
  /* Apply second rotation to put point p2 in plane (yz), in half plane y > 0*/
  /* p2 is almost on edge p0p1*/
  if ( (p2XTemp* p2XTemp + p2YTemp* p2YTemp) < EPS1 )
    return( distPtSeg_3d(p0,p1,pq) );
  
  cosAlpha = p2YTemp/(sqrt(p2XTemp* p2XTemp + p2YTemp* p2YTemp));
  sinAlpha = sqrt(1.0-cosAlpha * cosAlpha);
  if(p2XTemp <=0.0) {sinAlpha = -sinAlpha; }  
      
  lx2 = 0.0;      
  ly2 = sinAlpha * p2XTemp + cosAlpha * p2YTemp; 
  lz2 = p2ZTemp;       
      
  lxq = cosAlpha * pqXTemp - sinAlpha * pqYTemp;    
  lyq = sinAlpha * pqXTemp + cosAlpha * pqYTemp; 
  lzq = pqZTemp; 
  
  zone01 = lyq; 
  zone02 = ly2*lzq - lz2*lyq; 
  zone12 = lyq*(lz2-longp0p1) - ly2*(lzq - longp0p1); 
    
  if((zone01>=0.0)&&(zone02 >=0.0)&&(zone12>=0.0))
  {
    return (lxq*lxq);    
  }
  
  else
  {
    pointPlan0.c[0] = 0.0;
    pointPlan0.c[1] = 0.0; 
    
    pointPlan1.c[0] = 0.0; 
    pointPlan1.c[1] = longp0p1;
    
    pointPlan2.c[0] = ly2; 
    pointPlan2.c[1] = lz2;     
    
    pointPlana.c[0] = lyq; 
    pointPlana.c[1] = lzq;
    
    d01 = distpt_23d(&pointPlan0, &pointPlan1, &pointPlana);
    d02 = distpt_23d(&pointPlan0, &pointPlan2, &pointPlana);
    d12 = distpt_23d(&pointPlan1, &pointPlan2, &pointPlana);
    
    dTmp = D_MIN(d01, D_MIN(d02, d12)); 
    *proj = 2;
    
    return(dTmp + lxq*lxq);
  }
  
}

/* Computes (squared) distance from point pa to plane p0p1p2 */
double distptplan(pPoint p0, pPoint p1, pPoint p2, pPoint pa){
  double norm,ps,pv[3];
  
  pv[0] = (p1->c[1]-p0->c[1])*(p2->c[2]-p0->c[2]) -(p1->c[2]-p0->c[2])*(p2->c[1]-p0->c[1]);
  pv[1] = (p1->c[2]-p0->c[2])*(p2->c[0]-p0->c[0]) -(p1->c[0]-p0->c[0])*(p2->c[2]-p0->c[2]);	
  pv[2] = (p1->c[0]-p0->c[0])*(p2->c[1]-p0->c[1]) -(p1->c[1]-p0->c[1])*(p2->c[0]-p0->c[0]);
	
  norm = pv[0]*pv[0] + pv[1]*pv[1] + pv[2]*pv[2];
  assert(norm >0.0);	
  ps = pv[0]*(pa->c[0]-p0->c[0]) + pv[1]*(pa->c[1]-p0->c[1]) + pv[2]*(pa->c[2]-p0->c[2]);
  ps = ps*ps;
	
  return(ps/norm);
  
}

int interSegTria(pMesh mesh,pPoint pa,pPoint pb,pTria pt) {
  pPoint   p1,p2,p3;
  Point    v;
  double   d1,d2,d3,dd;
	
  p1 = &mesh->point[pt->v[0]];
  p2 = &mesh->point[pt->v[1]];
  p3 = &mesh->point[pt->v[2]];
  
  d1 = detOrient_3d(p1,p2,p3,pa);
  d2 = detOrient_3d(p1,p2,p3,pb);
  if ( d1*d2 > 0.0 )  return(0);
  
  v.c[0] = pb->c[0] - pa->c[0];
  v.c[1] = pb->c[1] - pa->c[1];
  v.c[2] = pb->c[2] - pa->c[2];
  d1 = determinant3Pts1Vct_3d(p1,p2,p3,&v);
  d2 = determinant3Pts1Vct_3d(p1,p2,pb,&v);
  d3 = determinant3Pts1Vct_3d(p1,pb,p3,&v);
  dd = d2 * d3;
  if ( dd <= 0.0 || d2*d1 <= 0.0 )  return(0);
  
  d1 = determinant3Pts1Vct_3d(p2,p3,p1,&v);
  d2 = determinant3Pts1Vct_3d(p2,p3,pb,&v);
  d3 = determinant3Pts1Vct_3d(p2,pb,p1,&v);
  dd = d2 * d3;
  if ( dd <= 0.0 || d2*d1 <= 0.0 )  return(0);
  
  return(1);
}

/* Computes distance from pa to level 0 of sol in tetra ntetra */
double distnv0_3d(pMesh mesh, pSol sol, int ntetra, pPoint pa,char *proj){
  pTetra      pt;
  pPoint      p0,p1,p2,p3,p0p,p1p,p2p,p3p;
  double      v0,v1,v2,v3,v0p,v1p,v2p,v3p,lambda,d1,d2,d3,d4,dd;
  int         i0,i1,i2,i3;
  Point       p,q,r,t;
  int         nzeros,nplus,nmoins;
  char        projTmp;	
  
  pt = &mesh->tetra[ntetra];
  i0 = pt->v[0];	
  i1 = pt->v[1];	
  i2 = pt->v[2];	
  i3 = pt->v[3];	
  p0 = &mesh->point[i0];
  p1 = &mesh->point[i1];
  p2 = &mesh->point[i2];
  p3 = &mesh->point[i3];
  v0 = sol->val[i0];
  v1 = sol->val[i1];
  v2 = sol->val[i2];
  v3 = sol->val[i3];
  nplus =0;
  nmoins = 0;
  nzeros =0;
  if(v0>0.0)  ++nplus;	
  if(v0==0.0) ++nzeros;	
  if(v0<0.0)  ++nmoins;
  if(v1>0.0)  ++nplus;	
  if(v1==0.0) ++nzeros;	
  if(v1<0.0)  ++nmoins;
  if(v2>0.0)  ++nplus;	
  if(v2==0.0) ++nzeros;	
  if(v2<0.0)  ++nmoins; 
  if(v3>0.0)  ++nplus;	
  if(v3==0.0) ++nzeros;	
  if(v3<0.0)  ++nmoins;
  
  /* wrong configuration : flat level set function */ 
  if ( !((v0 != 0.)||(v1 != 0.)||(v2 != 0.)||(v3!=0.)) )
		printf("ERROR distnv0_3d : flat ls function ; ntetra = %d: %E %E %E %E\n",ntetra,v0,v1,v2,v3);
	assert((v0 != 0.)||(v1 != 0.)||(v2 != 0.)||(v3!=0.));	
	
  /* case lv 0 in tetra ntetra is given by 3 0 vertices */
	
  if((v0==0.)&&(v1==0.)&&(v2==0.)){
    return(distpt_3d(p0, p1, p2, pa, proj));	  
  }
  
  if((v0==0.)&&(v1==0.)&&(v3==0.)){
    return(distpt_3d(p0, p1, p3, pa, proj));	  
  }	
	
  if((v0==0.)&&(v2==0.)&&(v3==0.)){
    return(distpt_3d(p0, p2, p3, pa, proj));	  
  }	
	
  if((v1==0.)&&(v2==0.)&&(v3==0.)){
    return(distpt_3d(p1, p2, p3, pa, proj));	  
  }	
	
  /* case 2 vertices of tetra are on level set 0 */
	
  if((v0==0.)&&(v1==0.)){
    if(v2*v3>0.){
      *proj =2;
      return(distPtSeg_3d(p0,p1,pa));	
    }
    else{
      lambda = v2/(v2-v3);
      p.c[0] = p2->c[0] + lambda*(p3->c[0] - p2->c[0]);	
      p.c[1] = p2->c[1] + lambda*(p3->c[1] - p2->c[1]);	
      p.c[2] = p2->c[2] + lambda*(p3->c[2] - p2->c[2]);	
      return(distpt_3d(p0,p1,&p,pa,proj));	
    }
  }
	
  if((v0==0.)&&(v2==0.)){
    if(v1*v3>0.){
      *proj =2;
      return(distPtSeg_3d(p0,p2,pa));	
    }
    else{
      lambda = v1/(v1-v3);
      p.c[0] = p1->c[0] + lambda*(p3->c[0] - p1->c[0]);	
      p.c[1] = p1->c[1] + lambda*(p3->c[1] - p1->c[1]);	
      p.c[2] = p1->c[2] + lambda*(p3->c[2] - p1->c[2]);	
      return(distpt_3d(p0,p2,&p,pa,proj));	
    }
  }	
  
  if((v0==0.)&&(v3==0.)){
    if(v1*v2>0.){
      *proj =2;
      return(distPtSeg_3d(p0,p3,pa));	
    }
    else{
      lambda = v1/(v1-v2);
      p.c[0] = p1->c[0] + lambda*(p2->c[0] - p1->c[0]);	
      p.c[1] = p1->c[1] + lambda*(p2->c[1] - p1->c[1]);	
      p.c[2] = p1->c[2] + lambda*(p2->c[2] - p1->c[2]);	
      return(distpt_3d(p0,p3,&p,pa,proj));	
    }
  }	
	
  if((v1==0.)&&(v2==0.)){
    if(v0*v3>0.){
      *proj =2;
      return(distPtSeg_3d(p1,p2,pa));	
    }
    else{
      lambda = v0/(v0-v3);
      p.c[0] = p0->c[0] + lambda*(p3->c[0] - p0->c[0]);	
      p.c[1] = p0->c[1] + lambda*(p3->c[1] - p0->c[1]);	
      p.c[2] = p0->c[2] + lambda*(p3->c[2] - p0->c[2]);	
      return(distpt_3d(p1,p2,&p,pa,proj));	
    }
  }	
	
  if((v1==0.)&&(v3==0.)){
    if(v0*v2>0.){
      *proj =2;
      return(distPtSeg_3d(p1,p3,pa));	
    }
    else{
      lambda = v0/(v0-v2);
      p.c[0] = p0->c[0] + lambda*(p2->c[0] - p0->c[0]);	
      p.c[1] = p0->c[1] + lambda*(p2->c[1] - p0->c[1]);	
      p.c[2] = p0->c[2] + lambda*(p2->c[2] - p0->c[2]);	
      return(distpt_3d(p1,p3,&p,pa,proj));	
    }
  }
  
  if((v2==0.)&&(v3==0.)){
    if(v0*v1>0.){
      *proj =2;
      return(distPtSeg_3d(p2,p3,pa));	
    }
    else{
      lambda = v0/(v0-v1);
      p.c[0] = p0->c[0] + lambda*(p1->c[0] - p0->c[0]);	
      p.c[1] = p0->c[1] + lambda*(p1->c[1] - p0->c[1]);	
      p.c[2] = p0->c[2] + lambda*(p1->c[2] - p0->c[2]);	
      return(distpt_3d(p2,p3,&p,pa,proj));	
    }
  }			
  
  if(nzeros ==1){
    if(v0 ==0.0) p0p = p0;
    if(v1 ==0.0) p0p = p1;
    if(v2 ==0.0) p0p = p2;
    if(v3 ==0.0) p0p = p3;
    v0p = 0.0;
	  
    if(nplus ==1){
      if(v0 >0.0) { p1p = p0; v1p = v0; }
      if(v1 >0.0) { p1p = p1; v1p = v1; }
      if(v2 >0.0) { p1p = p2; v1p = v2; }
      if(v3 >0.0) { p1p = p3; v1p = v3; }
      
      if((v0 <0.0)&&(v1 <0.0)) { p2p = p0; v2p = v0; p3p = p1; v3p = v1;}
      if((v0 <0.0)&&(v2 <0.0)) { p2p = p0; v2p = v0; p3p = p2; v3p = v2;}
      if((v0 <0.0)&&(v3 <0.0)) { p2p = p0; v2p = v0; p3p = p3; v3p = v3;}
      if((v1 <0.0)&&(v2 <0.0)) { p2p = p1; v2p = v1; p3p = p2; v3p = v2;}
      if((v1 <0.0)&&(v3 <0.0)) { p2p = p1; v2p = v1; p3p = p3; v3p = v3;}
      if((v2 <0.0)&&(v3 <0.0)) { p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
    }
	  
    else{
      /* Return squared distance to the point where the ls function is 0.0 */
      if (nmoins != 1) {
        if(v0 ==0.0) p0p = p0;
        if(v1 ==0.0) p0p = p1;
        if(v2 ==0.0) p0p = p2;
        if(v3 ==0.0) p0p = p3;
        dd = (pa->c[0]-p0p->c[0])*(pa->c[0]-p0p->c[0]);
        dd += (pa->c[1]-p0p->c[1])*(pa->c[1]-p0p->c[1]);
        dd += (pa->c[2]-p0p->c[2])*(pa->c[2]-p0p->c[2]);
        return(dd);
      }
      
      if(v0 <0.0) { p1p = p0; v1p = v0; }
      if(v1 <0.0) { p1p = p1; v1p = v1; }
      if(v2 <0.0) { p1p = p2; v1p = v2; }
      if(v3 <0.0) { p1p = p3; v1p = v3; }
		  
      if((v0 >0.0)&&(v1 >0.0)) { p2p = p0; v2p = v0; p3p = p1; v3p = v1;}
      if((v0 >0.0)&&(v2 >0.0)) { p2p = p0; v2p = v0; p3p = p2; v3p = v2;}
      if((v0 >0.0)&&(v3 >0.0)) { p2p = p0; v2p = v0; p3p = p3; v3p = v3;}
      if((v1 >0.0)&&(v2 >0.0)) { p2p = p1; v2p = v1; p3p = p2; v3p = v2;}
      if((v1 >0.0)&&(v3 >0.0)) { p2p = p1; v2p = v1; p3p = p3; v3p = v3;}
      if((v2 >0.0)&&(v3 >0.0)) { p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
    }
    
    /* now, p0p has value 0, p1p has different sign from p2p and p3p, which have same sign */
    lambda = v1p/(v1p-v2p);
    p.c[0] = p1p->c[0] + lambda*(p2p->c[0] - p1p->c[0]);	
    p.c[1] = p1p->c[1] + lambda*(p2p->c[1] - p1p->c[1]);	
    p.c[2] = p1p->c[2] + lambda*(p2p->c[2] - p1p->c[2]);	
	  
    lambda = v1p/(v1p-v3p);
    q.c[0] = p1p->c[0] + lambda*(p3p->c[0] - p1p->c[0]);	
    q.c[1] = p1p->c[1] + lambda*(p3p->c[1] - p1p->c[1]);	
    q.c[2] = p1p->c[2] + lambda*(p3p->c[2] - p1p->c[2]);
	  
    return(distpt_3d(p0p, &p, &q, pa, proj));  
    
  }
  
  /* general case : 2 possible configurations */	
  assert(nzeros ==0);
  
  if(nplus ==1){
    assert(nmoins == 3);
    if(v0 >0.0) { p0p = p0; v0p = v0; p1p = p1; v1p = v1; p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
    if(v1 >0.0) { p0p = p1; v0p = v1; p1p = p0; v1p = v0; p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
    if(v2 >0.0) { p0p = p2; v0p = v2; p1p = p0; v1p = v0; p2p = p1; v2p = v1; p3p = p3; v3p = v3;}
    if(v3 >0.0) { p0p = p3; v0p = v3; p1p = p0; v1p = v0; p2p = p1; v2p = v1; p3p = p2; v3p = v2;}  
	  
    lambda = v0p/(v0p-v1p);
    p.c[0] = p0p->c[0] + lambda*(p1p->c[0] - p0p->c[0]);	
    p.c[1] = p0p->c[1] + lambda*(p1p->c[1] - p0p->c[1]);	
    p.c[2] = p0p->c[2] + lambda*(p1p->c[2] - p0p->c[2]);	
	  
    lambda = v0p/(v0p-v2p);
    q.c[0] = p0p->c[0] + lambda*(p2p->c[0] - p0p->c[0]);	
    q.c[1] = p0p->c[1] + lambda*(p2p->c[1] - p0p->c[1]);	
    q.c[2] = p0p->c[2] + lambda*(p2p->c[2] - p0p->c[2]); 
	  
    lambda = v0p/(v0p-v3p);
    r.c[0] = p0p->c[0] + lambda*(p3p->c[0] - p0p->c[0]);	
    r.c[1] = p0p->c[1] + lambda*(p3p->c[1] - p0p->c[1]);	
    r.c[2] = p0p->c[2] + lambda*(p3p->c[2] - p0p->c[2]);
    
    return(distpt_3d(&p, &q, &r, pa, proj));  
    
  }
	
  if(nmoins ==1){
    assert(nplus == 3);
    if(v0 <0.0) { p0p = p0; v0p = v0; p1p = p1; v1p = v1; p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
    if(v1 <0.0) { p0p = p1; v0p = v1; p1p = p0; v1p = v0; p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
    if(v2 <0.0) { p0p = p2; v0p = v2; p1p = p0; v1p = v0; p2p = p1; v2p = v1; p3p = p3; v3p = v3;}
    if(v3 <0.0) { p0p = p3; v0p = v3; p1p = p0; v1p = v0; p2p = p1; v2p = v1; p3p = p2; v3p = v2;}  
		
    lambda = v0p/(v0p-v1p);
    p.c[0] = p0p->c[0] + lambda*(p1p->c[0] - p0p->c[0]);	
    p.c[1] = p0p->c[1] + lambda*(p1p->c[1] - p0p->c[1]);	
    p.c[2] = p0p->c[2] + lambda*(p1p->c[2] - p0p->c[2]);	
		
    lambda = v0p/(v0p-v2p);
    q.c[0] = p0p->c[0] + lambda*(p2p->c[0] - p0p->c[0]);	
    q.c[1] = p0p->c[1] + lambda*(p2p->c[1] - p0p->c[1]);	
    q.c[2] = p0p->c[2] + lambda*(p2p->c[2] - p0p->c[2]); 
		
    lambda = v0p/(v0p-v3p);
    r.c[0] = p0p->c[0] + lambda*(p3p->c[0] - p0p->c[0]);	
    r.c[1] = p0p->c[1] + lambda*(p3p->c[1] - p0p->c[1]);	
    r.c[2] = p0p->c[2] + lambda*(p3p->c[2] - p0p->c[2]);
		
    return(distpt_3d(&p, &q, &r, pa, proj));  
  }	
  
  if(nplus ==2){
    assert(nmoins ==2);
    *proj =2; //default value, unless...
    if((v0 >0.0)&&(v1>0.0)) { p0p = p0; v0p = v0; p1p = p1; v1p = v1; p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
    if((v0 >0.0)&&(v2>0.0)) { p0p = p0; v0p = v0; p1p = p2; v1p = v2; p2p = p1; v2p = v1; p3p = p3; v3p = v3;}
    if((v0 >0.0)&&(v3>0.0)) { p0p = p0; v0p = v0; p1p = p3; v1p = v3; p2p = p1; v2p = v1; p3p = p2; v3p = v2;}
    if((v1 >0.0)&&(v2>0.0)) { p0p = p1; v0p = v1; p1p = p2; v1p = v2; p2p = p0; v2p = v0; p3p = p3; v3p = v3;}
    if((v1 >0.0)&&(v3>0.0)) { p0p = p1; v0p = v1; p1p = p3; v1p = v3; p2p = p0; v2p = v0; p3p = p2; v3p = v2;}
    if((v2 >0.0)&&(v3>0.0)) { p0p = p2; v0p = v2; p1p = p3; v1p = v3; p2p = p0; v2p = v0; p3p = p1; v3p = v1;}
    
    /* p, q and r,t go together (same face) */
    lambda = v0p/(v0p-v2p);
    p.c[0] = p0p->c[0] + lambda*(p2p->c[0] - p0p->c[0]);	
    p.c[1] = p0p->c[1] + lambda*(p2p->c[1] - p0p->c[1]);	
    p.c[2] = p0p->c[2] + lambda*(p2p->c[2] - p0p->c[2]);	
	  
    lambda = v0p/(v0p-v3p);
    q.c[0] = p0p->c[0] + lambda*(p3p->c[0] - p0p->c[0]);	
    q.c[1] = p0p->c[1] + lambda*(p3p->c[1] - p0p->c[1]);	
    q.c[2] = p0p->c[2] + lambda*(p3p->c[2] - p0p->c[2]);  
	  
    lambda = v1p/(v1p-v2p);
    r.c[0] = p1p->c[0] + lambda*(p2p->c[0] - p1p->c[0]);	
    r.c[1] = p1p->c[1] + lambda*(p2p->c[1] - p1p->c[1]);	
    r.c[2] = p1p->c[2] + lambda*(p2p->c[2] - p1p->c[2]);	
	  
    lambda = v1p/(v1p-v3p);
    t.c[0] = p1p->c[0] + lambda*(p3p->c[0] - p1p->c[0]);	
    t.c[1] = p1p->c[1] + lambda*(p3p->c[1] - p1p->c[1]);	
    t.c[2] = p1p->c[2] + lambda*(p3p->c[2] - p1p->c[2]);  
	  
    d1 = distpt_3d(&p, &q, &r, pa, &projTmp);
    if(projTmp ==1){*proj = 1; return(d1);}
    d2 = distpt_3d(&p, &q, &t, pa, &projTmp);
    if(projTmp ==1){*proj = 1; return(d2);}
    d3 = distpt_3d(&r, &t, &p, pa, &projTmp);
    if(projTmp ==1){*proj = 1; return(d3);}
    d4 = distpt_3d(&r, &t, &q, pa, &projTmp);
    if(projTmp ==1){*proj = 1; return(d4);}
    
    return(D_MIN(D_MIN(D_MIN(d1,d2),d3),d4));
	  
  }
	
  return(1);
}

/* Computes distance from pa to level 0 of sol in tetra ntetra */
double distnv0approx_3d(pMesh mesh, pSol sol, int ntetra, pPoint pa){
	pTetra      pt;
	pPoint      p0,p1,p2,p3,p0p,p1p,p2p,p3p;
	double      v0,v1,v2,v3,v0p,v1p,v2p,v3p,lambda,d1,d2,d3,d4;	
	int         i0,i1,i2,i3;
	Point       p,q,r,t;
	int         nzeros,nplus,nmoins;
  
	pt = &mesh->tetra[ntetra];
	i0 = pt->v[0];	
	i1 = pt->v[1];	
	i2 = pt->v[2];	
	i3 = pt->v[3];	
	p0 = &mesh->point[i0];
	p1 = &mesh->point[i1];
	p2 = &mesh->point[i2];
	p3 = &mesh->point[i3];
	v0 = sol->val[i0];
	v1 = sol->val[i1];
	v2 = sol->val[i2];
	v3 = sol->val[i3];
	nplus =0;
	nmoins = 0;
	nzeros =0;
	if(v0>0.0)  ++nplus;	
	if(v0==0.0) ++nzeros;	
	if(v0<0.0)  ++nmoins;
	if(v1>0.0)  ++nplus;	
	if(v1==0.0) ++nzeros;	
	if(v1<0.0)  ++nmoins;
	if(v2>0.0)  ++nplus;	
	if(v2==0.0) ++nzeros;	
	if(v2<0.0)  ++nmoins; 
	if(v3>0.0)  ++nplus;	
	if(v3==0.0) ++nzeros;	
	if(v3<0.0)  ++nmoins;
	
	// wrong configuration : flat level set function 
	if ( !((v0 != 0.)||(v1 != 0.)||(v2 != 0.)||(v3!=0.)) )
		printf("ERROR distnv0_3d : flat ls function ; ntetra = %d: %E %E %E %E\n",ntetra,v0,v1,v2,v3);
	assert((v0 != 0.)||(v1 != 0.)||(v2 != 0.)||(v3!=0.));	
	
	/* case lv 0 in tetra ntetra is given by 3 0 vertices */
	
	if((v0==0.)&&(v1==0.)&&(v2==0.)){
		return(distptplan(p0, p1, p2, pa));	  
	}
	
	if((v0==0.)&&(v1==0.)&&(v3==0.)){
		return(distptplan(p0, p1, p3, pa));	  
	}	
	
	if((v0==0.)&&(v2==0.)&&(v3==0.)){
		return(distptplan(p0, p2, p3, pa));	  
	}	
	
	if((v1==0.)&&(v2==0.)&&(v3==0.)){
		return(distptplan(p1, p2, p3, pa));	  
	}	
	
	/* case 2 vertices of tetra are on level set 0 */
	
	if((v0==0.)&&(v1==0.)){
		if(v2*v3>0.){
			return(distPtSeg_3d(p0,p1,pa));	
		}
		else{
			lambda = v2/(v2-v3);
			p.c[0] = p2->c[0] + lambda*(p3->c[0] - p2->c[0]);	
			p.c[1] = p2->c[1] + lambda*(p3->c[1] - p2->c[1]);	
			p.c[2] = p2->c[2] + lambda*(p3->c[2] - p2->c[2]);	
			return(distptplan(p0,p1,&p,pa));	
		}
	}
	
	if((v0==0.)&&(v2==0.)){
		if(v1*v3>0.){
			return(distPtSeg_3d(p0,p2,pa));	
		}
		else{
			lambda = v1/(v1-v3);
			p.c[0] = p1->c[0] + lambda*(p3->c[0] - p1->c[0]);	
			p.c[1] = p1->c[1] + lambda*(p3->c[1] - p1->c[1]);	
			p.c[2] = p1->c[2] + lambda*(p3->c[2] - p1->c[2]);	
			return(distptplan(p0,p2,&p,pa));	
		}
	}	
  
	if((v0==0.)&&(v3==0.)){
		if(v1*v2>0.){
			return(distPtSeg_3d(p0,p3,pa));	
		}
		else{
			lambda = v1/(v1-v2);
			p.c[0] = p1->c[0] + lambda*(p2->c[0] - p1->c[0]);	
			p.c[1] = p1->c[1] + lambda*(p2->c[1] - p1->c[1]);	
			p.c[2] = p1->c[2] + lambda*(p2->c[2] - p1->c[2]);	
			return(distptplan(p0,p3,&p,pa));	
		}
	}	
	
	if((v1==0.)&&(v2==0.)){
		if(v0*v3>0.){
			return(distPtSeg_3d(p1,p2,pa));	
		}
		else{
			lambda = v0/(v0-v3);
			p.c[0] = p0->c[0] + lambda*(p3->c[0] - p0->c[0]);	
			p.c[1] = p0->c[1] + lambda*(p3->c[1] - p0->c[1]);	
			p.c[2] = p0->c[2] + lambda*(p3->c[2] - p0->c[2]);	
			return(distptplan(p1,p2,&p,pa));	
		}
	}	
	
	if((v1==0.)&&(v3==0.)){
		if(v0*v2>0.){
			return(distPtSeg_3d(p1,p3,pa));	
		}
		else{
			lambda = v0/(v0-v2);
			p.c[0] = p0->c[0] + lambda*(p2->c[0] - p0->c[0]);	
			p.c[1] = p0->c[1] + lambda*(p2->c[1] - p0->c[1]);	
			p.c[2] = p0->c[2] + lambda*(p2->c[2] - p0->c[2]);	
			return(distptplan(p1,p3,&p,pa));	
		}
	}
	
	if((v2==0.)&&(v3==0.)){
		if(v0*v1>0.){
			return(distPtSeg_3d(p2,p3,pa));	
		}
		else{
			lambda = v0/(v0-v1);
			p.c[0] = p0->c[0] + lambda*(p1->c[0] - p0->c[0]);	
			p.c[1] = p0->c[1] + lambda*(p1->c[1] - p0->c[1]);	
			p.c[2] = p0->c[2] + lambda*(p1->c[2] - p0->c[2]);	
			return(distptplan(p2,p3,&p,pa));	
		}
	}			
	
	if(nzeros ==1){
		if(v0 ==0.0) p0p = p0;
		if(v1 ==0.0) p0p = p1;
		if(v2 ==0.0) p0p = p2;
		if(v3 ==0.0) p0p = p3;
		v0p = 0.0;
		
		if(nplus ==1){
			if(v0 >0.0) { p1p = p0; v1p = v0; }
			if(v1 >0.0) { p1p = p1; v1p = v1; }
			if(v2 >0.0) { p1p = p2; v1p = v2; }
			if(v3 >0.0) { p1p = p3; v1p = v3; }
			
			if((v0 <0.0)&&(v1 <0.0)) { p2p = p0; v2p = v0; p3p = p1; v3p = v1;}
			if((v0 <0.0)&&(v2 <0.0)) { p2p = p0; v2p = v0; p3p = p2; v3p = v2;}
			if((v0 <0.0)&&(v3 <0.0)) { p2p = p0; v2p = v0; p3p = p3; v3p = v3;}
			if((v1 <0.0)&&(v2 <0.0)) { p2p = p1; v2p = v1; p3p = p2; v3p = v2;}
			if((v1 <0.0)&&(v3 <0.0)) { p2p = p1; v2p = v1; p3p = p3; v3p = v3;}
			if((v2 <0.0)&&(v3 <0.0)) { p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
		}
		
		else{
			assert(nmoins ==1);	
			if(v0 <0.0) { p1p = p0; v1p = v0; }
			if(v1 <0.0) { p1p = p1; v1p = v1; }
			if(v2 <0.0) { p1p = p2; v1p = v2; }
			if(v3 <0.0) { p1p = p3; v1p = v3; }
			
			if((v0 >0.0)&&(v1 >0.0)) { p2p = p0; v2p = v0; p3p = p1; v3p = v1;}
			if((v0 >0.0)&&(v2 >0.0)) { p2p = p0; v2p = v0; p3p = p2; v3p = v2;}
			if((v0 >0.0)&&(v3 >0.0)) { p2p = p0; v2p = v0; p3p = p3; v3p = v3;}
			if((v1 >0.0)&&(v2 >0.0)) { p2p = p1; v2p = v1; p3p = p2; v3p = v2;}
			if((v1 >0.0)&&(v3 >0.0)) { p2p = p1; v2p = v1; p3p = p3; v3p = v3;}
			if((v2 >0.0)&&(v3 >0.0)) { p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
		}
		
		/* Now, p0p has value 0, p1p has different sign from p2p and p3p, which have same sign */
		lambda = v1p/(v1p-v2p);
		p.c[0] = p1p->c[0] + lambda*(p2p->c[0] - p1p->c[0]);	
		p.c[1] = p1p->c[1] + lambda*(p2p->c[1] - p1p->c[1]);	
		p.c[2] = p1p->c[2] + lambda*(p2p->c[2] - p1p->c[2]);	
		
		lambda = v1p/(v1p-v3p);
		q.c[0] = p1p->c[0] + lambda*(p3p->c[0] - p1p->c[0]);	
		q.c[1] = p1p->c[1] + lambda*(p3p->c[1] - p1p->c[1]);	
		q.c[2] = p1p->c[2] + lambda*(p3p->c[2] - p1p->c[2]);
		
		return(distptplan(p0p, &p, &q, pa));  
		
	}
	
	/* general case : 2 possible configurations	*/
	assert(nzeros ==0);
	
	if(nplus ==1){
		assert(nmoins == 3);
		if(v0 >0.0) { p0p = p0; v0p = v0; p1p = p1; v1p = v1; p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
		if(v1 >0.0) { p0p = p1; v0p = v1; p1p = p0; v1p = v0; p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
		if(v2 >0.0) { p0p = p2; v0p = v2; p1p = p0; v1p = v0; p2p = p1; v2p = v1; p3p = p3; v3p = v3;}
		if(v3 >0.0) { p0p = p3; v0p = v3; p1p = p0; v1p = v0; p2p = p1; v2p = v1; p3p = p2; v3p = v2;}  
		
		lambda = v0p/(v0p-v1p);
		p.c[0] = p0p->c[0] + lambda*(p1p->c[0] - p0p->c[0]);	
		p.c[1] = p0p->c[1] + lambda*(p1p->c[1] - p0p->c[1]);	
		p.c[2] = p0p->c[2] + lambda*(p1p->c[2] - p0p->c[2]);	
		
		lambda = v0p/(v0p-v2p);
		q.c[0] = p0p->c[0] + lambda*(p2p->c[0] - p0p->c[0]);	
		q.c[1] = p0p->c[1] + lambda*(p2p->c[1] - p0p->c[1]);	
		q.c[2] = p0p->c[2] + lambda*(p2p->c[2] - p0p->c[2]); 
		
		lambda = v0p/(v0p-v3p);
		r.c[0] = p0p->c[0] + lambda*(p3p->c[0] - p0p->c[0]);	
		r.c[1] = p0p->c[1] + lambda*(p3p->c[1] - p0p->c[1]);	
		r.c[2] = p0p->c[2] + lambda*(p3p->c[2] - p0p->c[2]);
		
		return(distptplan(&p, &q, &r, pa));  
		
	}
	
	if(nmoins ==1){
		assert(nplus == 3);
		if(v0 <0.0) { p0p = p0; v0p = v0; p1p = p1; v1p = v1; p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
		if(v1 <0.0) { p0p = p1; v0p = v1; p1p = p0; v1p = v0; p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
		if(v2 <0.0) { p0p = p2; v0p = v2; p1p = p0; v1p = v0; p2p = p1; v2p = v1; p3p = p3; v3p = v3;}
		if(v3 <0.0) { p0p = p3; v0p = v3; p1p = p0; v1p = v0; p2p = p1; v2p = v1; p3p = p2; v3p = v2;}  
		
		lambda = v0p/(v0p-v1p);
		p.c[0] = p0p->c[0] + lambda*(p1p->c[0] - p0p->c[0]);	
		p.c[1] = p0p->c[1] + lambda*(p1p->c[1] - p0p->c[1]);	
		p.c[2] = p0p->c[2] + lambda*(p1p->c[2] - p0p->c[2]);	
		
		lambda = v0p/(v0p-v2p);
		q.c[0] = p0p->c[0] + lambda*(p2p->c[0] - p0p->c[0]);	
		q.c[1] = p0p->c[1] + lambda*(p2p->c[1] - p0p->c[1]);	
		q.c[2] = p0p->c[2] + lambda*(p2p->c[2] - p0p->c[2]); 
		
		lambda = v0p/(v0p-v3p);
		r.c[0] = p0p->c[0] + lambda*(p3p->c[0] - p0p->c[0]);	
		r.c[1] = p0p->c[1] + lambda*(p3p->c[1] - p0p->c[1]);	
		r.c[2] = p0p->c[2] + lambda*(p3p->c[2] - p0p->c[2]);
		
		return(distptplan(&p, &q, &r, pa));  
	}	
  
	if(nplus ==2){
		assert(nmoins ==2);
		if((v0 >0.0)&&(v1>0.0)) { p0p = p0; v0p = v0; p1p = p1; v1p = v1; p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
		if((v0 >0.0)&&(v2>0.0)) { p0p = p0; v0p = v0; p1p = p2; v1p = v2; p2p = p1; v2p = v1; p3p = p3; v3p = v3;}
		if((v0 >0.0)&&(v3>0.0)) { p0p = p0; v0p = v0; p1p = p3; v1p = v3; p2p = p1; v2p = v1; p3p = p2; v3p = v2;}
		if((v1 >0.0)&&(v2>0.0)) { p0p = p1; v0p = v1; p1p = p2; v1p = v2; p2p = p0; v2p = v0; p3p = p3; v3p = v3;}
		if((v1 >0.0)&&(v3>0.0)) { p0p = p1; v0p = v1; p1p = p3; v1p = v3; p2p = p0; v2p = v0; p3p = p2; v3p = v2;}
		if((v2 >0.0)&&(v3>0.0)) { p0p = p2; v0p = v2; p1p = p3; v1p = v3; p2p = p0; v2p = v0; p3p = p1; v3p = v1;}
		
		/* p, q and r,t go together (same face) */
		lambda = v0p/(v0p-v2p);
		p.c[0] = p0p->c[0] + lambda*(p2p->c[0] - p0p->c[0]);	
		p.c[1] = p0p->c[1] + lambda*(p2p->c[1] - p0p->c[1]);	
		p.c[2] = p0p->c[2] + lambda*(p2p->c[2] - p0p->c[2]);	
		
		lambda = v0p/(v0p-v3p);
		q.c[0] = p0p->c[0] + lambda*(p3p->c[0] - p0p->c[0]);	
		q.c[1] = p0p->c[1] + lambda*(p3p->c[1] - p0p->c[1]);	
		q.c[2] = p0p->c[2] + lambda*(p3p->c[2] - p0p->c[2]);  
		
		lambda = v1p/(v1p-v2p);
		r.c[0] = p1p->c[0] + lambda*(p2p->c[0] - p1p->c[0]);	
		r.c[1] = p1p->c[1] + lambda*(p2p->c[1] - p1p->c[1]);	
		r.c[2] = p1p->c[2] + lambda*(p2p->c[2] - p1p->c[2]);	
		
		lambda = v1p/(v1p-v3p);
		t.c[0] = p1p->c[0] + lambda*(p3p->c[0] - p1p->c[0]);	
		t.c[1] = p1p->c[1] + lambda*(p3p->c[1] - p1p->c[1]);	
		t.c[2] = p1p->c[2] + lambda*(p3p->c[2] - p1p->c[2]);  
		
		d1 = distptplan(&p, &q, &r, pa);
		d2 = distptplan(&p, &q, &t, pa);
		d3 = distptplan(&r, &t, &p, pa);
		d4 = distptplan(&r, &t, &q, pa);
		
		return(D_MIN(D_MIN(D_MIN(d1,d2),d3),d4));
	}
	return(1);
}


/* Store in circum (length 4* nb +1)
   - the 3 coordinates of circumcenter of the surface triangle defined by the
     0 level of sol in each triangle in bndy
   - the associated (squared) circumradius.
*/
int buildcircumredis_3d(pMesh mesh, pSol sol, int *bndy, int nb, double *circum){
  pTetra      pt;
  pPoint      p0,p1,p2,p3,p0p,p1p,p2p,p3p;
  Point       p,q,r,t;
  double      v0,v1,v2,v3,v0p,v1p,v2p,v3p,*circ,lambda;
  int         i, i0,i1,i2,i3,nzeros,nplus,nmoins;

  /* Travel tetras crossed by the 0 lv of sol */
  for(i=1; i<=nb; i++){
    pt = &mesh->tetra[bndy[i]];
    i0 = pt->v[0];	
    i1 = pt->v[1];	
    i2 = pt->v[2];	
    i3 = pt->v[3];
        
    p0 = &mesh->point[i0];
    p1 = &mesh->point[i1];
    p2 = &mesh->point[i2];
    p3 = &mesh->point[i3];
    
    v0 = sol->val[i0];
    v1 = sol->val[i1];
    v2 = sol->val[i2];
    v3 = sol->val[i3];
    
    nplus =0;
    nmoins = 0;
    nzeros =0;
    
    if ( v0 > 0.0 )  ++nplus;
    if ( v0 == 0.0 ) ++nzeros;
    if ( v0 < 0.0 )  ++nmoins;
    if ( v1 > 0.0 )  ++nplus;
    if ( v1 == 0.0 ) ++nzeros;
    if ( v1 < 0.0 )  ++nmoins;
    if ( v2 > 0.0 )  ++nplus;
    if ( v2 == 0.0 ) ++nzeros;
    if ( v2 < 0.0 )  ++nmoins;
    if ( v3 > 0.0 )  ++nplus;
    if ( v3 == 0.0 ) ++nzeros;
    if ( v3 < 0.0 )  ++nmoins;
    circ = &circum[4*(i-1)+1];  
	  
    /* Case where lv 0 in bndy[i] is given by 3 0 vertices */
    if ( (v0 == 0.) && ( v1 == 0. ) && ( v2 == 0. ) ) {
      circumcoords(p0, p1, p2, circ);
      continue;
    }
	  
    if ( ( v0 == 0. ) && ( v1 == 0. ) && ( v3 == 0. ) ) {
      circumcoords(p0, p1, p3, circ);
      continue;
    }	
	  
    if ( (v0 == 0. ) && ( v2 == 0. ) && ( v3 == 0. ) ) {
      circumcoords(p0, p2, p3, circ);
      continue;
    }	
	  
    if ( ( v1 == 0. ) && ( v2 == 0. ) && ( v3 == 0. ) ) {
      circumcoords(p1, p2, p3, circ);
      continue;
    }
    
    /* Case where 2 vertices of bndy[i] are on lv 0 */
    if ( ( v0 == 0. ) && ( v1 == 0. ) ) {
      if ( v2*v3 >0. ) {
        circ[0] = 0.0;	
        circ[1] = 0.0;	
        circ[2] = 0.0;	
        circ[3] = 0.0;	
      }
      else{
        lambda = v2/(v2-v3);
        p.c[0] = p2->c[0] + lambda*(p3->c[0] - p2->c[0]);	
        p.c[1] = p2->c[1] + lambda*(p3->c[1] - p2->c[1]);	
        p.c[2] = p2->c[2] + lambda*(p3->c[2] - p2->c[2]);	
        circumcoords(p0,p1,&p,circ);
        continue;	 
      }
    }
	  
    if ( ( v0 == 0. ) && ( v2 == 0. ) ) {
      if ( v1*v3 >0. ) {
        circ[0] = 0.0;	
        circ[1] = 0.0;	
        circ[2] = 0.0;	
        circ[3] = 0.0;	
      }
      else {
        lambda = v1/(v1-v3);
        p.c[0] = p1->c[0] + lambda*(p3->c[0] - p1->c[0]);	
        p.c[1] = p1->c[1] + lambda*(p3->c[1] - p1->c[1]);	
        p.c[2] = p1->c[2] + lambda*(p3->c[2] - p1->c[2]);	
        circumcoords(p0,p2,&p,circ);
        continue;   
      }
    }	
	  
    if ( ( v0 == 0. ) && ( v3 == 0. ) ) {
      if ( v1*v2 >0. ) {
        circ[0] = 0.0;	
        circ[1] = 0.0;	
        circ[2] = 0.0;	
        circ[3] = 0.0;
      }
      else{
        lambda = v1/(v1-v2);
        p.c[0] = p1->c[0] + lambda*(p2->c[0] - p1->c[0]);	
        p.c[1] = p1->c[1] + lambda*(p2->c[1] - p1->c[1]);	
        p.c[2] = p1->c[2] + lambda*(p2->c[2] - p1->c[2]);	
        circumcoords(p0,p3,&p,circ);
        continue;    
      }
    }	
	  
    if ( ( v1 == 0. ) && ( v2 == 0. ) ) {
      if ( v0*v3 > 0. ) {
        circ[0] = 0.0;	
        circ[1] = 0.0;	
        circ[2] = 0.0;	
        circ[3] = 0.0;	
      }
      else {
        lambda = v0/(v0-v3);
        p.c[0] = p0->c[0] + lambda*(p3->c[0] - p0->c[0]);	
        p.c[1] = p0->c[1] + lambda*(p3->c[1] - p0->c[1]);	
        p.c[2] = p0->c[2] + lambda*(p3->c[2] - p0->c[2]);	
        circumcoords(p1,p2,&p,circ);
        continue;  
      }
    }	
	  
    if ( ( v1 == 0. ) && ( v3 == 0. ) ) {
      if ( v0*v2 > 0. ) {
        circ[0] = 0.0;	
        circ[1] = 0.0;	
        circ[2] = 0.0;	
        circ[3] = 0.0;	
      }
      else {
        lambda = v0/(v0-v2);
        p.c[0] = p0->c[0] + lambda*(p2->c[0] - p0->c[0]);	
        p.c[1] = p0->c[1] + lambda*(p2->c[1] - p0->c[1]);	
        p.c[2] = p0->c[2] + lambda*(p2->c[2] - p0->c[2]);	
        circumcoords(p1,p3,&p,circ);
        continue;  
      }
    }
	  
    if ( ( v2 == 0. ) && ( v3 == 0. ) ) {
      if ( v0*v1 > 0. ) {
        circ[0] = 0.0;	
        circ[1] = 0.0;	
        circ[2] = 0.0;	
        circ[3] = 0.0;		
      }
      else {
        lambda = v0/(v0-v1);
        p.c[0] = p0->c[0] + lambda*(p1->c[0] - p0->c[0]);	
        p.c[1] = p0->c[1] + lambda*(p1->c[1] - p0->c[1]);	
        p.c[2] = p0->c[2] + lambda*(p1->c[2] - p0->c[2]);	
        circumcoords(p2,p3,&p,circ);
        continue;  
      }
    }
            
    /* Case where only one pt of the 0 lv set if vertex of bndy[i] */
    if ( nzeros == 1 ) {
      if ( v0 == 0.0 ) p0p = p0;
      else if ( v1 == 0.0 ) p0p = p1;
      else if ( v2 == 0.0 ) p0p = p2;
      else p0p = p3;
      v0p = 0.0;
		  
      if ( nplus == 1 ) {
        if ( v0 >0.0 ) { p1p = p0; v1p = v0; }
        else if ( v1 >0.0 ) { p1p = p1; v1p = v1; }
        else if ( v2 >0.0 ) { p1p = p2; v1p = v2; }
        else { p1p = p3; v1p = v3; }
			  
        if ( ( v0 < 0.0 ) && ( v1 < 0.0 ) ) { p2p = p0; v2p = v0; p3p = p1; v3p = v1;}
        else if ( ( v0 < 0.0 ) && ( v2 < 0.0 ) ) { p2p = p0; v2p = v0; p3p = p2; v3p = v2;}
        else if ( ( v0 < 0.0 ) && ( v3 < 0.0 ) ) { p2p = p0; v2p = v0; p3p = p3; v3p = v3;}
        else if ( ( v1 < 0.0 ) && ( v2 < 0.0 ) ) { p2p = p1; v2p = v1; p3p = p2; v3p = v2;}
        else if ( ( v1 < 0.0 ) && ( v3 < 0.0 ) ) { p2p = p1; v2p = v1; p3p = p3; v3p = v3;}
        else { p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
      }
		  
      else {
        if ( nmoins != 1 ){
          circ[0] = p0p->c[0];
          circ[1] = p0p->c[1];
          circ[2] = p0p->c[2];
          circ[3] = 0.0;
        }
        
        if ( v0 < 0.0 ) { p1p = p0; v1p = v0; }
        else if ( v1 < 0.0 ) { p1p = p1; v1p = v1; }
        else if ( v2 < 0.0 ) { p1p = p2; v1p = v2; }
        else { p1p = p3; v1p = v3; }
        
        if ( ( v0 > 0.0 ) && ( v1 > 0.0 ) ) { p2p = p0; v2p = v0; p3p = p1; v3p = v1;}
        else if ( ( v0 > 0.0 ) && ( v2 > 0.0 ) ) { p2p = p0; v2p = v0; p3p = p2; v3p = v2;}
        else if ( ( v0 > 0.0 ) && ( v3 > 0.0 ) ) { p2p = p0; v2p = v0; p3p = p3; v3p = v3;}
        else if ( ( v1 > 0.0 ) && ( v2 > 0.0 ) ) { p2p = p1; v2p = v1; p3p = p2; v3p = v2;}
        else if ( ( v1 > 0.0 ) && ( v3 > 0.0 ) ) { p2p = p1; v2p = v1; p3p = p3; v3p = v3;}
        else { p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
      }
		  
      /* p0p has value 0, p1p has different sign from p2p and p3p, which have same sign */
      if ( fabs(v1p-v2p) < EPS1 || fabs(v1p-v3p) < EPS1 ) {
        memset(circ,0.0,4*sizeof(double));
        continue;
      }
      
      lambda = v1p/(v1p-v2p);
      p.c[0] = p1p->c[0] + lambda*(p2p->c[0] - p1p->c[0]);	
      p.c[1] = p1p->c[1] + lambda*(p2p->c[1] - p1p->c[1]);	
      p.c[2] = p1p->c[2] + lambda*(p2p->c[2] - p1p->c[2]);	
		  
      lambda = v1p/(v1p-v3p);
      q.c[0] = p1p->c[0] + lambda*(p3p->c[0] - p1p->c[0]);	
      q.c[1] = p1p->c[1] + lambda*(p3p->c[1] - p1p->c[1]);	
      q.c[2] = p1p->c[2] + lambda*(p3p->c[2] - p1p->c[2]);
		  
      circumcoords(p0p, &p, &q, circ);
      continue;
		  
    }  
    
    /* General case : 2 possible configurations	*/
    assert( nzeros == 0 );
	  
    if ( nplus == 1 ) {
      assert( nmoins == 3 );
      if ( v0 > 0.0 ) { p0p = p0; v0p = v0; p1p = p1; v1p = v1; p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
      else if ( v1 > 0.0 ) { p0p = p1; v0p = v1; p1p = p0; v1p = v0; p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
      else if ( v2 > 0.0 ) { p0p = p2; v0p = v2; p1p = p0; v1p = v0; p2p = p1; v2p = v1; p3p = p3; v3p = v3;}
      else { p0p = p3; v0p = v3; p1p = p0; v1p = v0; p2p = p1; v2p = v1; p3p = p2; v3p = v2;}
      
      if ( fabs(v0p-v1p) < EPS1 || fabs(v0p-v2p) < EPS1 || fabs(v0p-v3p) < EPS1 ) {
        memset(circ,0.0,4*sizeof(double));
        continue;
      }
      
      lambda = v0p/(v0p-v1p);
      p.c[0] = p0p->c[0] + lambda*(p1p->c[0] - p0p->c[0]);	
      p.c[1] = p0p->c[1] + lambda*(p1p->c[1] - p0p->c[1]);	
      p.c[2] = p0p->c[2] + lambda*(p1p->c[2] - p0p->c[2]);	
		  
      lambda = v0p/(v0p-v2p);
      q.c[0] = p0p->c[0] + lambda*(p2p->c[0] - p0p->c[0]);	
      q.c[1] = p0p->c[1] + lambda*(p2p->c[1] - p0p->c[1]);	
      q.c[2] = p0p->c[2] + lambda*(p2p->c[2] - p0p->c[2]); 
		  
      lambda = v0p/(v0p-v3p);
      r.c[0] = p0p->c[0] + lambda*(p3p->c[0] - p0p->c[0]);	
      r.c[1] = p0p->c[1] + lambda*(p3p->c[1] - p0p->c[1]);	
      r.c[2] = p0p->c[2] + lambda*(p3p->c[2] - p0p->c[2]);
		  
      circumcoords(&p, &q, &r, circ);
      continue;	
		  
	  }
	  
    if ( nmoins == 1 ) {
      assert ( nplus == 3 );
      if (v0 < 0.0 ) { p0p = p0; v0p = v0; p1p = p1; v1p = v1; p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
      else if (v1 < 0.0 ) { p0p = p1; v0p = v1; p1p = p0; v1p = v0; p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
      else if (v2 < 0.0 ) { p0p = p2; v0p = v2; p1p = p0; v1p = v0; p2p = p1; v2p = v1; p3p = p3; v3p = v3;}
      else { p0p = p3; v0p = v3; p1p = p0; v1p = v0; p2p = p1; v2p = v1; p3p = p2; v3p = v2;}
		  
      if ( fabs(v0p-v1p) < EPS1 || fabs(v0p-v2p) < EPS1 || fabs(v0p-v3p) < EPS1 ) {
        memset(circ,0.0,4*sizeof(double));
        continue;
      }
      
      lambda = v0p/(v0p-v1p);
      p.c[0] = p0p->c[0] + lambda*(p1p->c[0] - p0p->c[0]);	
      p.c[1] = p0p->c[1] + lambda*(p1p->c[1] - p0p->c[1]);	
      p.c[2] = p0p->c[2] + lambda*(p1p->c[2] - p0p->c[2]);	
		  
      lambda = v0p/(v0p-v2p);
      q.c[0] = p0p->c[0] + lambda*(p2p->c[0] - p0p->c[0]);	
      q.c[1] = p0p->c[1] + lambda*(p2p->c[1] - p0p->c[1]);	
      q.c[2] = p0p->c[2] + lambda*(p2p->c[2] - p0p->c[2]); 
		  
      lambda = v0p/(v0p-v3p);
      r.c[0] = p0p->c[0] + lambda*(p3p->c[0] - p0p->c[0]);	
      r.c[1] = p0p->c[1] + lambda*(p3p->c[1] - p0p->c[1]);	
      r.c[2] = p0p->c[2] + lambda*(p3p->c[2] - p0p->c[2]);
		  
      circumcoords(&p, &q, &r, circ); 
      continue;
    }	
	  
    if ( nplus == 2 ) {
      assert ( nmoins == 2 );
      if ( (v0 > 0.0 ) && ( v1 > 0.0 ) ) { p0p = p0; v0p = v0; p1p = p1; v1p = v1; p2p = p2; v2p = v2; p3p = p3; v3p = v3;}
      else if ( (v0 > 0.0 ) && ( v2 > 0.0 ) ) { p0p = p0; v0p = v0; p1p = p2; v1p = v2; p2p = p1; v2p = v1; p3p = p3; v3p = v3;}
      else if ( (v0 > 0.0 ) && ( v3 > 0.0 ) ) { p0p = p0; v0p = v0; p1p = p3; v1p = v3; p2p = p1; v2p = v1; p3p = p2; v3p = v2;}
      else if ( (v1 > 0.0 ) && ( v2 > 0.0 ) ) { p0p = p1; v0p = v1; p1p = p2; v1p = v2; p2p = p0; v2p = v0; p3p = p3; v3p = v3;}
      else if ( (v1 > 0.0 ) && ( v3 > 0.0 ) ) { p0p = p1; v0p = v1; p1p = p3; v1p = v3; p2p = p0; v2p = v0; p3p = p2; v3p = v2;}
      else { p0p = p2; v0p = v2; p1p = p3; v1p = v3; p2p = p0; v2p = v0; p3p = p1; v3p = v1;}
		  
      /* p, q and r,t go together (same face) */
      if ( fabs(v0p-v2p) < EPS1 || fabs(v0p-v3p) < EPS1 || fabs(v1p-v2p) < EPS1 || fabs(v1p-v3p) < EPS1 ) {
        memset(circ,0.0,4*sizeof(double));
        continue;
      }
      
      lambda = v0p/(v0p-v2p);
      p.c[0] = p0p->c[0] + lambda*(p2p->c[0] - p0p->c[0]);	
      p.c[1] = p0p->c[1] + lambda*(p2p->c[1] - p0p->c[1]);	
      p.c[2] = p0p->c[2] + lambda*(p2p->c[2] - p0p->c[2]);	
		  
      lambda = v0p/(v0p-v3p);
      q.c[0] = p0p->c[0] + lambda*(p3p->c[0] - p0p->c[0]);	
      q.c[1] = p0p->c[1] + lambda*(p3p->c[1] - p0p->c[1]);	
      q.c[2] = p0p->c[2] + lambda*(p3p->c[2] - p0p->c[2]);  
      
      lambda = v1p/(v1p-v2p);
      r.c[0] = p1p->c[0] + lambda*(p2p->c[0] - p1p->c[0]);	
      r.c[1] = p1p->c[1] + lambda*(p2p->c[1] - p1p->c[1]);	
      r.c[2] = p1p->c[2] + lambda*(p2p->c[2] - p1p->c[2]);	
      
      lambda = v1p/(v1p-v3p);
      t.c[0] = p1p->c[0] + lambda*(p3p->c[0] - p1p->c[0]);	
      t.c[1] = p1p->c[1] + lambda*(p3p->c[1] - p1p->c[1]);	
      t.c[2] = p1p->c[2] + lambda*(p3p->c[2] - p1p->c[2]);  
		  
      circ[0] = 0.25*(p.c[0] + q.c[0] + r.c[0] + t.c[0]);
      circ[1] = 0.25*(p.c[1] + q.c[1] + r.c[1] + t.c[1]);
      circ[2] = 0.25*(p.c[2] + q.c[2] + r.c[2] + t.c[2]);
      
      circ[3] = (p.c[0] - circ[0])*(p.c[0] - circ[0]) + (p.c[1] - circ[1])*(p.c[1] - circ[1])\
		  + (p.c[2] - circ[2])*(p.c[2] - circ[2]);
      
      circ[3] = D_MAX(circ[3], (q.c[0] - circ[0])*(q.c[0] - circ[0]) + \
                      (q.c[1] - circ[1])*(q.c[1] - circ[1]) + (q.c[2] - circ[2])*(q.c[2] - circ[2]));
      
      circ[3] = D_MAX(circ[3], (r.c[0] - circ[0])*(r.c[0] - circ[0]) + \
                      (r.c[1] - circ[1])*(r.c[1] - circ[1]) + (r.c[2] - circ[2])*(r.c[2] - circ[2]));
      
      circ[3] = D_MAX(circ[3], (t.c[0] - circ[0])*(t.c[0] - circ[0]) + \
                      (t.c[1] - circ[1])*(t.c[1] - circ[1]) + (t.c[2] - circ[2])*(t.c[2] - circ[2]));
		  
      continue;
    }
  }
  
  return(1);
}

/* Analyze gradient of sol, and truncate inconsistent values */
int corrGrad_3d(pMesh mesh, pSol sol) {
  pTetra    pt,pt1; 
  pPoint    p0,p1,p2,p3;
  double    *grd,*area,m[9],im[9],vol,u0u1,u0u2,u0u3,gr[3],dd,dl,mgrd2;
  int       ier,k,np[4],l,ilist,list[LONMAX+2];
  char      i,j,jp;
  
  mgrd2 = MAXGRD*MAXGRD;
  
  /* Reset flag field to points and tetrahedra */
  for(k=1; k<= mesh->ne; k++) {
    pt = &mesh->tetra[k]; 
    pt->flag = 0;
    
    for(i=0; i<4; i++) {
      mesh->point[pt->v[i]].flag = 0;
    }
  }
  
  grd = (double*)calloc(3*mesh->np+1,sizeof(double));
  area = (double*)calloc(mesh->np+1,sizeof(double));
    
  /* Compute an approximate gradient at each point */
  for(k=1; k<= mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !pt->v[0] ) continue; 
    
    np[0] = pt->v[0];
    np[1] = pt->v[1];
    np[2] = pt->v[2];
    np[3] = pt->v[3]; 
    
    p0 = &mesh->point[np[0]];
    p1 = &mesh->point[np[1]];
    p2 = &mesh->point[np[2]];
    p3 = &mesh->point[np[3]];
    
    m[0] = p1->c[0] - p0->c[0];  m[1] = p1->c[1] - p0->c[1]; m[2] = p1->c[2] - p0->c[2]; 
    m[3] = p2->c[0] - p0->c[0];  m[4] = p2->c[1] - p0->c[1]; m[5] = p2->c[2] - p0->c[2]; 
    m[6] = p3->c[0] - p0->c[0];  m[7] = p3->c[1] - p0->c[1]; m[8] = p3->c[2] - p0->c[2]; 
    
    ier = invmatg(m,im);
    assert(ier); 
    
    u0u1 = sol->val[np[1]] - sol->val[np[0]];
    u0u2 = sol->val[np[2]] - sol->val[np[0]];
    u0u3 = sol->val[np[3]] - sol->val[np[0]];
    
    vol = volume(p0->c,p1->c,p2->c,p3->c);
    
    gr[0] = vol*(im[0]*u0u1 + im[1]*u0u2 + im[2]*u0u3);
    gr[1] = vol*(im[3]*u0u1 + im[4]*u0u2 + im[5]*u0u3);
    gr[2] = vol*(im[6]*u0u1 + im[7]*u0u2 + im[8]*u0u3);
    
    for(i=0; i<4; i++) {
      grd[3*(np[i]-1)+1] += gr[0];
      grd[3*(np[i]-1)+2] += gr[1];
      grd[3*(np[i]-1)+3] += gr[2];
      
      area[np[i]] += vol; 
    }
  }
  
  for(k=1; k<=mesh->np; k++) {
    grd[3*(k-1)+1] =  grd[3*(k-1)+1] / area[k];
    grd[3*(k-1)+2] =  grd[3*(k-1)+2] / area[k];
    grd[3*(k-1)+3] =  grd[3*(k-1)+3] / area[k];
  }
  
  /* Identify those points where the gradient of u is very steep, and correct them */
  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !pt->v[0] ) continue; 
  
    for(i=0; i<4; i++) {
      np[0] = pt->v[i]; 
      p0 = &mesh->point[np[0]];
      
      if ( p0->flag ) continue;
      dd = grd[3*(np[0]-1)+1]*grd[3*(np[0]-1)+1] + grd[3*(np[0]-1)+2]*grd[3*(np[0]-1)+2]\
           + grd[3*(np[0]-1)+3]*grd[3*(np[0]-1)+3];
          
      if ( dd < mgrd2 ) {
        p0->flag = 1; 
        continue;
      }
      
      ilist = boulet_3d(mesh,k,i,list);
      assert(ilist);
      dd = 0.0;
      dl = (double)ilist;
      
      for(l=0; l<ilist; l++) {
        pt1 = &mesh->tetra[list[l]/4];
        jp = list[l] % 4;
        for(j=0; j<3; j++) {
          jp = inxt3[jp];
          dd += ATHIRD*sol->val[pt1->v[jp]];
        }
      }
      dd = dd/dl;
      sol->val[np[0]] = dd;
      p0->flag = 1;
    }
  }
 
 
  free(area);
  free(grd);
  return(1);
}
