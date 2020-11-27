#include "mshdist.h"

extern unsigned char inxt2[5];

/* Return (squared) distance from point pa to segment (p1,p2);
   proj = 2 if distance is realized by p1 or p2,
          1 if it is realized by the orthogonal projection of pa on (p1p2) */
double distpt_s(pPoint p0,pPoint p1,pPoint pa,int *proj) {
  double   ux,uy,uz,p1p0,pap0,ps,lambda;
  
  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];
  p1p0 = ux*ux + uy*uy + uz*uz;
  pap0 = (pa->c[0]-p0->c[0])*(pa->c[0]-p0->c[0]) + (pa->c[1]-p0->c[1])*(pa->c[1]-p0->c[1]) + (pa->c[2]-p0->c[2])*(pa->c[2]-p0->c[2]);

  /* If p0p1 is too short, return squared distance to p0 */
  if ( p1p0 < EPS1 )
    return ( pap0 );
  
  ps = (pa->c[0]-p0->c[0])*ux + (pa->c[1]-p0->c[1])*uy + (pa->c[2]-p0->c[2])*uz;
  lambda = ps / p1p0;
  
  /* Closest point is p0 */
  if ( lambda < 0.0 ) {
    *proj = 2;
    return (pap0);
  }
  /* Closest point is p1 */
  else if ( lambda > 1.0 ) {
    *proj = 2;
    return ( (pa->c[0]-p1->c[0])*(pa->c[0]-p1->c[0]) + (pa->c[1]-p1->c[1])*(pa->c[1]-p1->c[1]) + (pa->c[2]-p1->c[2])*(pa->c[2]-p1->c[2]) );
  }
  /* Closest point is q = p0 + lambda*(p1-p0) \in [p0,p1] */
  else {
    *proj = 1;
    return ( fabs(pap0 - lambda*lambda*p1p0) );
  }
}