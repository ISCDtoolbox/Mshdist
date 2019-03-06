#include "mshdist.h"
#include <math.h>

extern unsigned char inxt3[7];
extern char ddb;


/* Generate holes in the pcb example */
int genHolesPCB_2d(pMesh mesh,pSol sol) {
  pPoint   p0;
  double   o[2],r,dd,dj;
  int      k,j;
  
  /* Initialization */
  for(k=1; k<=mesh->np; k++)
    sol->val[k] = 10.0;
  
  
  for (j=0; j<=4; j++) {
    dj  = (double)j;
    o[0] = 0.0 + dj*0.25;
    o[1] = 1.0;
    r = 0.05;
  
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(dd,sol->val[k]);
    }
  }
  
  for (j=1; j<=4; j++) {
    dj  = (double)j;
    o[0] = 0.0 + dj*0.2;
    o[1] = 0.75;
    r = 0.05;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(dd,sol->val[k]);
    }
  }
  
  for (j=0; j<=4; j++) {
    dj  = (double)j;
    o[0] = 0.0 + dj*0.25;
    o[1] = 0.5;
    r = 0.05;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(dd,sol->val[k]);
    }
  }
  
  for (j=1; j<=4; j++) {
    dj  = (double)j;
    o[0] = 0.0 + dj*0.2;
    o[1] = 0.25;
    r = 0.05;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(dd,sol->val[k]);
    }
  }
  
  for (j=0; j<=4; j++) {
    dj  = (double)j;
    o[0] = 0.0 + dj*0.25;
    o[1] = 0.0;
    r = 0.05;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(dd,sol->val[k]);
    }
  }

  return(1);
}

/* Generate two holes in a unit square */
int gen2Holes_2d(pMesh mesh,pSol sol) {
  pPoint   p0;
  double   o[2],r,dd;
  int      k;
  
  /* Initialization */
  for(k=1; k<=mesh->np; k++)
    sol->val[k] = 10.0;
  
  o[0] = 0.5;
  o[1] = 0.5;
  r = 0.475;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
    dd = sqrt(dd) - r;
    sol->val[k] = D_MIN(dd,sol->val[k]);
  }
  
  /*o[0] = 0.75;
  o[1] = 0.5;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
    dd = sqrt(dd)-r;
    sol->val[k] = D_MIN(dd,sol->val[k]);
  }*/
  
  return(1);
}

/* Generate an implicit function corresponding to the shape optimization cantilever test-case */
int genHolesCanti_2d(pMesh mesh,pSol sol) {
  pPoint    p0;
  double    o[2],dj,dd,r;
  int       k;
  char      j;
  
  r = 0.1;
  
  /* Initialization */
  for(k=1; k<=mesh->np; k++)
    sol->val[k] = 10.0;
    
  for(j=0; j<4; j++) {
    dj  = (double)j;
    o[0] = 0.4 + dj*0.4;
    o[1] = 0.75;  
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(sol->val[k],dd);
    }
  }  
  
  for(j=0; j<3; j++) {
    dj  = (double)j;
    o[0] = 0.5 + dj*0.5;
    o[1] = 0.5;  
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(sol->val[k],dd);
    }
  }  
  
  for(j=0; j<4; j++) {
    dj  = (double)j;
    o[0] = 0.4 + dj*0.4;
    o[1] = 0.25;  
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(sol->val[k],dd);
    }
  }  
  
  
  for(k=1; k<=mesh->np; k++) {
    sol->val[k] *= -1.0;
  }

  return(1);
}

/* Generate an implicit function corresponding to the shape optimization radiator test-case */
int genHolesRadia_2d(pMesh mesh,pSol sol) {
  pPoint    p0;
  double    o[2],dj,dd,r;
  int       k;
  char      j;
  
  r = 0.05;
  
  /* Initialization */
  for(k=1; k<=mesh->np; k++)
    sol->val[k] = 10.0;
  
  for(j=0; j<3; j++) {
    dj  = (double)j;
    o[0] = 0.25 + dj*0.25;
    o[1] = 1.0;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(sol->val[k],dd);
    }
  }
  
  for(j=0; j<4; j++) {
    dj  = (double)j;
    o[0] = 0.2 + dj*0.2;
    o[1] = 0.833333;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(sol->val[k],dd);
    }
  }
  
  for(j=0; j<3; j++) {
    dj  = (double)j;
    o[0] = 0.25 + dj*0.25;
    o[1] = 0.666666;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(sol->val[k],dd);
    }
  }

  for(j=0; j<4; j++) {
    dj  = (double)j;
    o[0] = 0.2 + dj*0.2;
    o[1] = 0.5;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(sol->val[k],dd);
    }
  }
  
  for(j=0; j<3; j++) {
    dj  = (double)j;
    o[0] = 0.25 + dj*0.25;
    o[1] = 0.333333;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(sol->val[k],dd);
    }
  }
  
  for(j=0; j<4; j++) {
    dj  = (double)j;
    o[0] = 0.2 + dj*0.2;
    o[1] = 0.166666;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(sol->val[k],dd);
    }
  }
  
  for(j=0; j<3; j++) {
    dj  = (double)j;
    o[0] = 0.25 + dj*0.25;
    o[1] = 0.0;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(sol->val[k],dd);
    }
  }
  
  for(k=1; k<=mesh->np; k++)
    sol->val[k] *= -1.0;
  
  return(1);
}

/* Generate an implicit function corresponding to the shape optimization mast test-case */
int genHolesMast_2d(pMesh mesh,pSol sol) {
  pPoint    p0;
  double    o[2],dj,dd,r;
  int       k;
  char      j;
  
  r = 5.0;
  
  /* Initialization */
  for(k=1; k<=mesh->np; k++)
    sol->val[k] = 10.0;
  
  for(j=0; j<3; j++) {
    dj  = (double)j;
    o[0] = -10.0 + dj*30.0;
    o[1] = 120.0;  
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(sol->val[k],dd);
    }
  }  
  
  for(j=0; j<3; j++) {
    dj  = (double)j;
    o[0] = -10.0 + dj*30.0;
    o[1] = 90.0;   
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(sol->val[k],dd);
    }
  }  
  
  for(j=0; j<3; j++) {
    dj  = (double)j;
    o[0] = 20.0;
    o[1] = 20.0 + dj*20.0;  
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      dd = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1]);
      dd = sqrt(dd) - r;
      sol->val[k] = D_MIN(sol->val[k],dd);
    }
  }  
  
  
  for(k=1; k<=mesh->np; k++) {
    sol->val[k] *= -1.0;
  }
  
  return(1);
}

/* Generate an implicit function associated to a circle */
int gen1Hole_2d(pMesh mesh,pSol sol) {
  pPoint    p0;
  double    o[2],dd,r;
  int       k;
  
  r = 0.3; 
  o[0] = 0.0;
  o[1] = 0.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    dd = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]);
    sol->val[k] = sqrt(dd) - r;
  }
  
  return(1);
}



