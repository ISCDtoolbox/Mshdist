#include "mshdist.h"
#include <math.h>

extern unsigned char inxt3[7];
extern char ddb;

/* Generate an analytic function (cube) */
int anafuncsq(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     dd;
  int        k;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    dd = (p0->c[0]-0.5)*(p0->c[0]-0.5) + (p0->c[1]-0.5)*(p0->c[1]-0.5) + (p0->c[2]-0.5)*(p0->c[2]-0.5);
    dd = sqrt(dd) -0.2;
    if ( fabs(dd) < 0.05 ) sol->val[k] = dd;
    else if ( dd < 0.0 ) sol->val[k] = -INIVAL_3d;
    else sol->val[k] = INIVAL_3d;
  }
  
  return(1);
}

/* Generate an analytic function (8 branches test-case ) */
int anafunc(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     dd;
  int        k;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    dd = tanh((p0->c[2]-2.0)*3)*tanh((p0->c[1])*3); //pow(p0->c[0]+5.0,4)*pow(p0->c[0]-4.0,4)*pow(p0->c[1],6)*
    sol->val[k] = dd;
  }

  return(1);
}

/* Generate an analytic function (buddha ) */
int anafuncbuddha(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     dd;
  int        k;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    dd = 100.0*tanh(fabs(p0->c[2])*0.01)*tanh(fabs(p0->c[0])*0.01); //pow(p0->c[0]+5.0,4)*pow(p0->c[0]-4.0,4)*pow(p0->c[1],6)*
    sol->val[k] = dd;
  }
  
  return(1);
}


/* Generate an analytic function (saddle) */
int anafuncsaddle(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     dd;
  int        k;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    dd = p0->c[2] - 0.5*( (2.0*p0->c[0]-1.0)*(2.0*p0->c[0]-1.0) \
         - (2.0*p0->c[1]-1.0)*(2.0*p0->c[1]-1.0) +1.0);
    sol->val[k] = dd;
  }
  
  return(1);
}

/* Generate an analytic function (sphere for Lagrangian motion) */
int anafuncspherelag(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     dd;
  int        k;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    dd = (p0->c[0]-0.5)*(p0->c[0]-0.5) + (p0->c[1]-0.5)*(p0->c[1]-0.5) + (p0->c[2]-0.5)*(p0->c[2]-0.5) -0.3*0.3;
    sol->val[k] = dd;
  }
  
  return(1);
}

/* Generate an analytic function (torus of radii R and r, of axis z) */
int anafunctorus(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     dd,R,r;
  int        k;
  
  R = 0.3;
  r = 0.15;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    dd = (p0->c[2]-0.5)*(p0->c[2]-0.5) + (R - sqrt((p0->c[0]-0.5)*(p0->c[0]-0.5) \
       + (p0->c[1]-0.5)*(p0->c[1]-0.5)))*(R - sqrt((p0->c[0]-0.5)*(p0->c[0]-0.5) + (p0->c[1]-0.5)*(p0->c[1]-0.5))) - r*r;
    sol->val[k] = dd;
  }
  
  return(1);
}

/* Generate an analytic function (helix) */
int anafunchelix(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     dd,r,rho,theta,phi,pi;
  int        k;
  
  pi = 3.14159;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];

    r = (p0->c[0]-0.5)*(p0->c[0]-0.5) + (p0->c[1]-0.5)*(p0->c[1]-0.5);
    r = sqrt(r);
    
    rho = (p0->c[0]-0.5)*(p0->c[0]-0.5) + (p0->c[1]-0.5)*(p0->c[1]-0.5) + (p0->c[2]-0.5)*(p0->c[2]-0.5);
    rho = sqrt(rho);
    
    phi = acos((p0->c[2]-0.5)/rho);
    
    theta = ( p0->c[1] > 0.5 ) ? acos( (p0->c[1]-0.5) / r) : (2*pi - acos( (p0->c[1]-0.5) / r));
        
    
    sol->val[k] = dd;
  }
  
  return(1);
}

/* Generate an implicit function corresponding to m*n*o holes */
int genholes_3d(pMesh mesh,pSol sol,int m,int n, int o){
  pPoint    p0;
  int       nx,ny,nz,k;
  double    x,y,z,r,d;
pTetra pt; 
char j;  
int nc;
  
  /* radius of balls, computed so they do not intersect each other */
  r = D_MIN(0.3/(m-1.0), D_MIN(0.3/(n-1.0),0.3/(o-1.0)));
  printf("Le rayon : %f \n",r);
  
  /* Initialization */
  for(k=1; k<=mesh->np; k++)
    sol->val[k] = 1.0;
    
  for(k=1; k<=mesh->np; k++){
    p0 = &mesh->point[k];
    
    /* Multiply by the signed distance function associated to each hole */
    for(nx = 0; nx<m; nx++){
      for(ny = 0; ny < n; ny++){
        for(nz = 0; nz< o; nz++){
          x = (double) nx/(m-1);
          y = (double) ny/(n-1);   
          z = (double) nz/(o-1);
                    
          d = (p0->c[0] - x)*(p0->c[0] - x) + (p0->c[1] - y)*(p0->c[1] - y) \
            +(p0->c[2] - z)*(p0->c[2] - z); 
            
          d = sqrt(d); 
          d -= r;
              
          sol->val[k] = D_MIN(d,sol->val[k]); 
        }
      }
    } 
  }
  
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    nc = 0; 
    for(j=0;j<4;j++){
      if(fabs(sol->val[pt->v[j]]) < 1.e-6){
        nc++;
      }
      if(nc == 4){
        printf("tetra %d : 4 vertices to 0 \n",k);
        exit(0);
      }
    }
  }

  return(1);
}

/* Generate an implicit function corresponding to the shape optimization cantilever test-case */
int genHolesCanti_3d(pMesh mesh,pSol sol) {
  pPoint    p0;
  double    a,b,c,x,y,z,d,op[3],dj;
  int       m,n,o,nx,ny,nz,k;
  char      j;
  
  /* Half-radii of the ellipsoids */
  a = 0.3;
  b = 0.8; 
  c = 0.4;
  
  /* Number of holes in each direction */
  m = 2; 
  n = 3; 
  o = 2;
  
  /* Initialization */
  for(k=1; k<=mesh->np; k++)
    sol->val[k] = 10.0;
  
  /* 4 edges */
  for(j=0; j<4; j++) {
    dj = (double)j;
    
    op[0] = 2.4;
    op[1] = 0.0 + j*1.666667;
    op[2] = 3.0;
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      d = (p0->c[0] - op[0])*(p0->c[0] - op[0]) + (p0->c[1] - op[1])*(p0->c[1] - op[1]) 
          + (p0->c[2] - op[2])*(p0->c[2] - op[2]); 
      d = sqrt(d) - 0.2;      
      sol->val[k] = D_MIN(d,sol->val[k]); 
    }
  } 
  
  for(j=0; j<4; j++) {
    dj = (double)j;
    
    op[0] = 0.0;
    op[1] = 0.0 + j*1.666667;
    op[2] = 3.0;
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      d = (p0->c[0] - op[0])*(p0->c[0] - op[0]) + (p0->c[1] - op[1])*(p0->c[1] - op[1]) 
          + (p0->c[2] - op[2])*(p0->c[2] - op[2]); 
      d = sqrt(d) - 0.2;      
      sol->val[k] = D_MIN(d,sol->val[k]); 
    }
  } 
  
  for(j=0; j<4; j++) {
    dj = (double)j;
    
    op[0] = 2.4;
    op[1] = 0.0 + j*1.666667;
    op[2] = 0.0;
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      d = (p0->c[0] - op[0])*(p0->c[0] - op[0]) + (p0->c[1] - op[1])*(p0->c[1] - op[1]) 
        + (p0->c[2] - op[2])*(p0->c[2] - op[2]); 
      d = sqrt(d) - 0.2;      
      sol->val[k] = D_MIN(d,sol->val[k]); 
    }
  } 
  
  for(j=0; j<4; j++) {
    dj = (double)j;
    
    op[0] = 0.0;
    op[1] = 0.0 + j*1.666667;
    op[2] = 0.0;
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      d = (p0->c[0] - op[0])*(p0->c[0] - op[0]) + (p0->c[1] - op[1])*(p0->c[1] - op[1]) 
        + (p0->c[2] - op[2])*(p0->c[2] - op[2]);  
      d = sqrt(d) - 0.2;      
      sol->val[k] = D_MIN(d,sol->val[k]); 
    }
  }
  
  for(j=1; j<4; j++) {
    dj = (double)j;
    
    op[0] = 1.2;
    op[1] = 0.0 + j*1.666667;
    op[2] = 1.5;
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      d = (p0->c[0] - op[0])*(p0->c[0] - op[0]) + (p0->c[1] - op[1])*(p0->c[1] - op[1]) 
        + (p0->c[2] - op[2])*(p0->c[2] - op[2]);  
      d = sqrt(d) - 0.2;      
      sol->val[k] = D_MIN(d,sol->val[k]); 
    }
  } 
   
  /*for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
        
    for(nx=0; nx<m; nx++) {
      for(ny=0; ny<n; ny++) {
        for(nz=0; nz<o; nz++) {
          x = m == 1 ? 1.2 : 2.4*(double) nx/(m-1);
          y = n == 1 ? 2.5 : 5.0*(double) ny/(n-1);   
          z = o == 1 ? 1.5 : 3.0*(double) nz/(o-1);
                     
          d = (p0->c[0] - x)*(p0->c[0] - x)/(a*a) + (p0->c[1] - y)*(p0->c[1] - y)/(b*b) \
              + (p0->c[2] - z)*(p0->c[2] - z)/(c*c) - 1.0; 
                                      
          sol->val[k] = D_MIN(d,sol->val[k]); 
        }
      }
    }  
  }*/
  
  for(k=1; k<=mesh->np; k++) {
    sol->val[k] *= -1.0;
  }

  return(1);
}

/* Generate an implicit function corresponding to the shape optimization 
   perturbed short cantilever test-case */
int genHolesSCantia_3d(pMesh mesh,pSol sol) {
  pPoint    p0;
  double    r,o[3],d,di,dj;
  int       k,n;
  char      i,j;
  
  r = 0.2;
  
  /* Initialization */
  for(k=1; k<=mesh->np; k++)
    sol->val[k] = 1000.0;
  
  /* Left side */
  for(i=0; i<3; i++) {
    di = (double)i;
    for(j=0; j<3; j++) {
      dj = (double)j;
      o[0] = 0.5*di;
      o[1] = 0.0;
      o[2] = 1.0*dj;   
    
      for(k=1; k<=mesh->np; k++) {
        p0 = &mesh->point[k];
        d = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1])\
            + (p0->c[2] - o[2])*(p0->c[2] - o[2]);  
        d = sqrt(d) - r;
        sol->val[k] = D_MIN(sol->val[k],d);
      }
    } 
  } 
  
  /* Middle holes */
  for(i=0; i<3; i++) {
    di = (double)i;
    for(j=0; j<3; j++) {
      if ( j==1 ) continue; 
      dj = (double)j;
      o[0] = 0.5*di;
      o[1] = 0.5;
      o[2] = 1.0*dj;   
      
      for(k=1; k<=mesh->np; k++) {
        p0 = &mesh->point[k];
        d = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1])\
        + (p0->c[2] - o[2])*(p0->c[2] - o[2]);  
        d = sqrt(d) - r;
        sol->val[k] = D_MIN(sol->val[k],d);
      }
    } 
  } 
     
  /* Right side */
  for(i=0; i<3; i++) {
    di = (double)i;
    for(j=0; j<3; j++) {
      dj = (double)j;
      o[0] = 0.5*di;
      o[1] = 1.0;
      o[2] = 1.0*dj;   
      
      for(k=1; k<=mesh->np; k++) {
        p0 = &mesh->point[k];
        d = (p0->c[0] - o[0])*(p0->c[0] - o[0]) + (p0->c[1] - o[1])*(p0->c[1] - o[1])\
        + (p0->c[2] - o[2])*(p0->c[2] - o[2]);  
        d = sqrt(d) - r;
        sol->val[k] = D_MIN(sol->val[k],d);
      }
    } 
  }      
  
  for(k=1; k<=mesh->np; k++) {
    sol->val[k] *= -1.0;
  }
  
  return(1);
}

/* Generate an implicit function corresponding to the shape optimization cantilever test-case */
int genHolesCantia_3d(pMesh mesh,pSol sol) {
  pPoint    p0;
  double    r,x,y,z,d;
  int       k,n;
  
  r = 0.4;
  
  /* Initialization */
  for(k=1; k<=mesh->np; k++)
    sol->val[k] = 1000.0;
  
  /* First edge */
  for(n=0; n<3; n++) {
    x = 2.4;
    y = 1.25 + 1.25*n;
    z = 3.0;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
    
      d = (p0->c[0] - x)*(p0->c[0] - x) + (p0->c[1] - y)*(p0->c[1] - y) \
        + (p0->c[2] - z)*(p0->c[2] - z) - r*r; 
          
      sol->val[k] = D_MIN(d,sol->val[k]); 
    }
  }
  
  /* Second edge */
  for(n=0; n<3; n++) {
    x = 0.0;
    y = 1.25 + 1.25*n;
    z = 3.0;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      
      d = (p0->c[0] - x)*(p0->c[0] - x) + (p0->c[1] - y)*(p0->c[1] - y) \
      + (p0->c[2] - z)*(p0->c[2] - z) - r*r; 
      
      sol->val[k] = D_MIN(d,sol->val[k]); 
    }
  }
  
  /* Third edge */
  for(n=0; n<3; n++) {
    x = 2.4;
    y = 1.25 + 1.25*n;
    z = 0.0;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      
      d = (p0->c[0] - x)*(p0->c[0] - x) + (p0->c[1] - y)*(p0->c[1] - y) \
      + (p0->c[2] - z)*(p0->c[2] - z) - r*r; 
      
      sol->val[k] = D_MIN(d,sol->val[k]); 
    }
  }
  
  /* Fourth edge */
  for(n=0; n<3; n++) {
    x = 0.0;
    y = 1.25 + 1.25*n;
    z = 0.0;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      
      d = (p0->c[0] - x)*(p0->c[0] - x) + (p0->c[1] - y)*(p0->c[1] - y) \
      + (p0->c[2] - z)*(p0->c[2] - z) - r*r; 
      
      sol->val[k] = D_MIN(d,sol->val[k]); 
    }
  }
  
  /* First face */
  for(n=0; n<4; n++) {
    x = 1.2;
    y = 1.666666*n;
    z = 3.0;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      
      d = (p0->c[0] - x)*(p0->c[0] - x) + (p0->c[1] - y)*(p0->c[1] - y) \
      + (p0->c[2] - z)*(p0->c[2] - z) - r*r; 
      
      sol->val[k] = D_MIN(d,sol->val[k]); 
    }
  }
  
  /* Second face */
  for(n=0; n<4; n++) {
    x = 1.2;
    y = 1.666666*n;
    z = 0.0;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      
      d = (p0->c[0] - x)*(p0->c[0] - x) + (p0->c[1] - y)*(p0->c[1] - y) \
      + (p0->c[2] - z)*(p0->c[2] - z) - r*r; 
      
      sol->val[k] = D_MIN(d,sol->val[k]); 
    }
  }
  
  /* Third face */
  for(n=0; n<4; n++) {
    x = 0.0;
    y = 1.666666*n;
    z = 1.5;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      
      d = (p0->c[0] - x)*(p0->c[0] - x) + (p0->c[1] - y)*(p0->c[1] - y) \
      + (p0->c[2] - z)*(p0->c[2] - z) - r*r; 
      
      sol->val[k] = D_MIN(d,sol->val[k]); 
    }
  }
  
  /* Fourth face */
  for(n=0; n<4; n++) {
    x = 2.4;
    y = 1.666666*n;
    z = 1.5;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      
      d = (p0->c[0] - x)*(p0->c[0] - x) + (p0->c[1] - y)*(p0->c[1] - y) \
      + (p0->c[2] - z)*(p0->c[2] - z) - r*r; 
      
      sol->val[k] = D_MIN(d,sol->val[k]); 
    }
  }
  
  /* Inner holes */
  for(n=1; n<3; n++) {
    x = 1.2;
    y = 1.666666*n;
    z = 1.5;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      
      d = (p0->c[0] - x)*(p0->c[0] - x) + (p0->c[1] - y)*(p0->c[1] - y) \
      + (p0->c[2] - z)*(p0->c[2] - z) - r*r; 
      
      sol->val[k] = D_MIN(d,sol->val[k]); 
    }
  }
  
  
  for(k=1; k<=mesh->np; k++) {
    sol->val[k] *= -1.0;
  }
  
  return(1);
}

/* Generate a hole for Neumann boundary conditions in the cantilever test-case */
int holeCl_3d(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     r,d,c[3]; 
  int        k;
  
  r = 0.15;
  c[0] = 1.2;
  c[1] = 0.0;
  c[2] = 1.5;
  
  for(k=1;k<=mesh->np;k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
         + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    sol->val[k] = sqrt(d) - r;
  }

  return(1);
}

/* Generate holes for initialization in the Mast test-case */
int genHolesMast_3d(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     d,o[3],a,b,c,dj; 
  int        k,j;
  
  for (k=1; k<=mesh->np; k++) 
    sol->val[k] = 1000.0;
  
  for (j=0; j<4; j++) {
    dj = (double)j;
    o[0] = 20.0;
    o[1] = 20.0;
    o[2] = 0.0 + dj*30.0;
    a = 5.0;
    
    /* Define successive holes */
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      d =  (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
           + (p0->c[2]-o[2])*(p0->c[2]-o[2]) - a*a;
      sol->val[k] = D_MIN(sol->val[k],d);
    }
  }
  
  o[0] = 20.0;
  o[1] = 0.0;
  o[2] = 40.0;
  a = 5.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
         + (p0->c[2]-o[2])*(p0->c[2]-o[2]) - a*a;
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  o[0] = 20.0;
  o[1] = 40.0;
  o[2] = 40.0;
  a = 5.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]) - a*a;
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  o[0] = 0.0;
  o[1] = 20.0;
  o[2] = 40.0;
  a = 5.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
       + (p0->c[2]-o[2])*(p0->c[2]-o[2]) - a*a;
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  o[0] = 40.0;
  o[1] = 20.0;
  o[2] = 40.0;
  a = 5.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
       + (p0->c[2]-o[2])*(p0->c[2]-o[2]) - a*a;
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  /* Other series of lateral holes */
  o[0] = 20.0;
  o[1] = -20.0;
  o[2] = 100.0;
  a = 5.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]) - a*a;
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  o[0] = 20.0;
  o[1] = 60.0;
  o[2] = 100.0;
  a = 5.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]) - a*a;
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  o[0] = -20.0;
  o[1] = 20.0;
  o[2] = 100.0;
  a = 5.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]) - a*a;
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  o[0] = 60.0;
  o[1] = 20.0;
  o[2] = 100.0;
  a = 5.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]) - a*a;
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  /*o[0] = 40.0;
  o[1] = 20.0;
  o[2] = 126.0;
  a = 5.0;
  b = 25.0;
  c = 10.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
         + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c) - 1.0;
    sol->val[k] = d;
  }
    
  o[0] = 0.0;
  o[1] = 20.0;
  o[2] = 126.0;
  a = 5.0;
  b = 25.0;
  c = 10.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c) - 1.0;
    sol->val[k] = D_MIN(sol->val[k],d);
  }  
    
  o[0] = 40.0;
  o[1] = 20.0;
  o[2] = 40.0;
  a = 5.0;
  b = 8.0;
  c = 30.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c) - 1.0;
    sol->val[k] = D_MIN(sol->val[k],d);
  }    
  
  o[0] = 0.0;
  o[1] = 20.0;
  o[2] = 40.0;
  a = 5.0;
  b = 8.0;
  c = 30.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c) - 1.0;
    sol->val[k] = D_MIN(sol->val[k],d);
  }      
  
  o[0] = 20.0;
  o[1] = 0.0;
  o[2] = 40.0;
  a = 5.0;
  b = 8.0;
  c = 30.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c) - 1.0;
    sol->val[k] = D_MIN(sol->val[k],d);
  }    
  
  o[0] = 20.0;
  o[1] = 40.0;
  o[2] = 40.0;
  a = 5.0;
  b = 8.0;
  c = 30.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c) - 1.0;
    sol->val[k] = D_MIN(sol->val[k],d);
  }*/
  
  /* Change of signs */
  for(k=1; k<=mesh->np; k++) {
    sol->val[k] *= -1.0;
  }
  
  return(1);
}

/* Generate holes for Neumann boundary conditions in the Mast test-case */
int holeClMast_3d(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     r,d,c[3]; 
  int        k;
  
  for(k=1; k<=mesh->np; k++) 
    sol->val[k] = 1000.0;  
  
  /*r = 10.0;
  c[0] = 20.0;
  c[1] = -20.0;
  c[2] = 80.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
      + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    sol->val[k] = sqrt(d) - r;
  }
  
  r = 10.0;
  c[0] = 20.0;
  c[1] = 60.0;
  c[2] = 80.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
      + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
  }*/
    
  r = 5.0;
  c[0] = 60.0;
  c[1] = -20.0;
  c[2] = 80.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
    + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
  }
  
  c[0] = -20.0;
  c[1] = -20.0;
  c[2] = 80.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
    + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
  }
  
  c[0] = 60.0;
  c[1] = 60.0;
  c[2] = 80.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
    + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
  }  
  
  c[0] = -20.0;
  c[1] = 60.0;
  c[2] = 80.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
    + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
  }
  
  c[0] = 40.0;
  c[1] = 0.0;
  c[2] = 0.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
    + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
  }
  
  c[0] = 0.0;
  c[1] = 0.0;
  c[2] = 0.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
    + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
  }
  
  c[0] = 40.0;
  c[1] = 40.0;
  c[2] = 0.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
    + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
  }  
  
  c[0] = 0.0;
  c[1] = 40.0;
  c[2] = 0.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
    + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
  }
  
  /* Change signs */
  for(k=1; k<=mesh->np; k++) 
    sol->val[k] *=-1.0;
  
  return(1);
}

/* Generate holes for Dirichlet and Neumann boundary conditions in the Bridge test-case */
int holeClBridge_3d(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     d,c; 
  int        k;
  
  d = 20.0;
  
  /* First Dirichlet boundary condition */
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    sol->val[k] = p0->c[1] - d;
  }
  
  d = 180.0;
  
  /* Second Dirichlet boundary condition */
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    sol->val[k] = D_MIN(sol->val[k],d - p0->c[1]);
  }
  
  d = 15.0;
  c = 100.0;
  
  /* Neumann boundary conditions */
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    sol->val[k] = D_MIN(sol->val[k],fabs(p0->c[1]-c) - d);
  }
  
  return(1);
}

/* Generate holes for Dirichlet and Neumann boundary conditions in the Bridge test-case */
int holeClBridge2_3d(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     dd,dd1,n[3],c0[3],c1[3]; 
  int        k;
  
  /* First Dirichlet boundary condition */
  n[0] = 0.0; 
  n[1] = 5.0;
  n[2] = -6.0;
  
  dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  dd = sqrt(dd); 
  if ( dd > EPS2 ) {
    dd = 1.0 / dd;
    n[0] *= dd; 
    n[1] *= dd; 
    n[2] *= dd;  
  } 
  
  c0[0] = 20.0;
  c0[1] = 30.0; 
  c0[2] = 50.0 - 30.0/60.0*50.0;
    
  c1[0] = 20.0;
  c1[1] = 50.0; 
  c1[2] = 50.0 - 50.0/60.0*50.0;

  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    dd = -( (p0->c[0] - c0[0])*n[0] + (p0->c[1] - c0[1])*n[1] + (p0->c[1] - c0[1])*n[1]);
    dd1 = ( (p0->c[0] - c1[0])*n[0] + (p0->c[1] - c1[1])*n[1] + (p0->c[1] - c1[1])*n[1]);
    dd = D_MAX(dd,dd1);
    sol->val[k] = dd;
  }
    
  /* Second Dirichlet boundary condition */
  n[0] = 0.0; 
  n[1] = 5.0;
  n[2] = 6.0;
  
  dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  dd = sqrt(dd); 
  if ( dd > EPS2 ) {
    dd = 1.0 / dd;
    n[0] *= dd; 
    n[1] *= dd; 
    n[2] *= dd;  
  } 
  
  c0[0] = 20.0;
  c0[1] = 150.0; 
  c0[2] = 50.0 - 50.0/60.0*50.0;
  
  c1[0] = 20.0;
  c1[1] = 170.0; 
  c1[2] = 50.0-30.0/60.0*50.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    dd = -( (p0->c[0] - c0[0])*n[0] + (p0->c[1] - c0[1])*n[1] + (p0->c[1] - c0[1])*n[1]);
    dd1 = ( (p0->c[0] - c1[0])*n[0] + (p0->c[1] - c1[1])*n[1] + (p0->c[1] - c1[1])*n[1]);
    dd = D_MAX(dd,dd1);
    sol->val[k] = D_MIN(sol->val[k],dd);
  }

  
  return(1);
}


/* Generate holes for initialization in the Bridge test-case */
int genHolesBridge_3d(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     h,r,dd,c[3]; 
  int        k;
  
  /*h = 20.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    sol->val[k] = p0->c[2] - h;  
  }*/
    
  for(k=1; k<=mesh->np; k++) 
    sol->val[k] = 1000.0;
    
  r = 10.0;
  
  c[0] = 0.0;
  c[1] = 100.0;
  c[2] = 25.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    dd = (p0->c[0]-c[0])*(p0->c[0]-c[0]) +  (p0->c[1]-c[1])*(p0->c[1]-c[1])
         + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    dd = sqrt(dd)-r;
    sol->val[k] = D_MIN(sol->val[k],dd); 
  }
  
  c[0] = 40.0;
  c[1] = 100.0;
  c[2] = 25.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    dd = (p0->c[0]-c[0])*(p0->c[0]-c[0]) +  (p0->c[1]-c[1])*(p0->c[1]-c[1])
    + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    dd = sqrt(dd)-r;
    sol->val[k] = D_MIN(sol->val[k],dd); 
  }
  
  c[0] = 20.0;
  c[1] = 100.0;
  c[2] = 0.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    dd = (p0->c[0]-c[0])*(p0->c[0]-c[0]) +  (p0->c[1]-c[1])*(p0->c[1]-c[1])
    + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    dd = sqrt(dd)-r;
    sol->val[k] = D_MIN(sol->val[k],dd); 
  }
    
  for(k=1; k<=mesh->np; k++) 
    sol->val[k] *= -1.0;  
        
  return(1);
}

/* Generate holes for boundary conditions in the Starfish test-case */
int holeClStarfish_3d(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     r,d,R,pi,dn,c[3]; 
  int        k,n;
  
  pi = 3.141592653;
  
  r = 0.3;
  c[0] = 2.0;
  c[1] = 2.0;
  c[2] = 0.0;
  
  /* First (Neumann) boundary condition */
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
    + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    sol->val[k] = sqrt(d) - r;
  }
  
  /* Dirichlet boundary conditions */
  for(n=0; n<5; n++) {
    R = 1.6;
    r = 0.1;
    
    dn = (double)n;
    c[0] = 2.0 + R*cos(2.0*pi*dn / 5.0);
    c[1] = 2.0 + R*sin(2.0*pi*dn / 5.0);
    c[2] = 0.0;  
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
      + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
      sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
    }
  
  }
  
  return(1);
}

/* Generate initialization for the Starfish test-case */
int genHolesStarfish_3d(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     h; 
  int        k;
  
  h = 0.5;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    sol->val[k] = p0->c[2] - h;  
  }
  
  return(1);
}

/* Generate holes for boundary conditions in the crane test-case */
int holeClCrane(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     r,d,w,c[3]; 
  int        k;
  
  r = 5.0;
  c[0] = 15.0;
  c[1] = 150.0;
  c[2] = 125.0;
  
  /* First Neumann boundary condition */
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
       + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    sol->val[k] = sqrt(d) - r;
  }
    
  /* Second Neumann boundary condition */
  c[0] = 15.0;
  c[1] = -30.0;
  c[2] = 90.0;
  w = 20.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    sol->val[k] = D_MIN(sol->val[k],fabs(p0->c[1]-c[1])-w);
  }
  
  /* Dirichlet boundary conditions */
  r = 10.0;
  c[0] = 0.0;
  c[1] = 0.0;
  c[2] = 0.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
       + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
  }
  
  c[0] = 30.0;
  c[1] = 0.0;
  c[2] = 0.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
       + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
  }
  
  c[0] = 0.0;
  c[1] = 30.0;
  c[2] = 0.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
       + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
  }
  
  c[0] = 30.0;
  c[1] = 30.0;
  c[2] = 0.0;
  
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
       + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
    sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
  }
  
  return(1);
}

/* Generate initial guess in the crane test-case */
int holeCraneIni(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     r,d,dj,c[3]; 
  int        k,j;
  
  for(k=1; k<=mesh->np; k++) 
    sol->val[k] = 1000.0;  

  r = 5.0;
  
  /* Lower ranges of hole */
  for(j=0; j<=2; j++) {
    dj = (double)j;
    c[0] = 0.0;
    c[1] = 40.0 + dj*40.0;
    c[2] = 110.0;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
         + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
      sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
    }
  }
  
  
  for(j=0; j<=2; j++) {
    dj = (double)j;
    c[0] = 30.0;
    c[1] = 40.0 + dj*40.0;
    c[2] = 110.0;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
      + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
      sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
    }
  }
  
  /* Upper ranges of hole */
  for(j=0; j<=3; j++) {
    dj = (double)j;
    c[0] = 0.0;
    c[1] = 0.0 + dj*40.0;
    c[2] = 140.0;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
      + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
      sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
    }
  }
  
  
  for(j=0; j<=3; j++) {
    dj = (double)j;
    c[0] = 30.0;
    c[1] = 0.0 + dj*40.0;
    c[2] = 140.0;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
      + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
      sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
    }
  }
  
  /* Holes on the vertical part*/
  for(j=0; j<=1; j++) {
    dj = (double)j;
    c[0] = 0.0;
    c[1] = 15.0;
    c[2] = 40.0 + dj*40.0;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
      + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
      sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
    }
  }
  
  
  for(j=0; j<=1; j++) {
    dj = (double)j;
    c[0] = 15.0;
    c[1] = 0.0;
    c[2] = 40.0 + dj*40.0;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
      + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
      sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
    }
  }
  
  for(j=0; j<=1; j++) {
    dj = (double)j;
    c[0] = 30.0;
    c[1] = 15.0;
    c[2] = 40.0 + dj*40.0;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
      + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
      sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
    }
  }
  
  
  for(j=0; j<=1; j++) {
    dj = (double)j;
    c[0] = 15.0;
    c[1] = 30.0;
    c[2] = 40.0 + dj*40.0;
    
    for(k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      d =  (p0->c[0]-c[0])*(p0->c[0]-c[0]) + (p0->c[1]-c[1])*(p0->c[1]-c[1]) \
      + (p0->c[2]-c[2])*(p0->c[2]-c[2]);
      sol->val[k] = D_MIN(sol->val[k],sqrt(d) - r);
    }
  }

  for(k=1; k<=mesh->np; k++) 
    sol->val[k] *= -1.0; 
  
  return(1);
}

/* Generate boundary conditions in the gripping mechanism test-case */
int holeClGrip(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     w;
  int        k;
  
  w = 0.1; 
  
  /* Neumann boundary conditions */
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    sol->val[k] = fabs(p0->c[2]-1.0) - 0.1;
  }  
    
  /* Dirichlet boundary conditions */
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    sol->val[k] = D_MIN(sol->val[k],(fabs(p0->c[2]-1.0) - 0.8)*(fabs(p0->c[0]-0.0) - 0.1));
  }
  
  return(1);
}

/* Generate boundary conditions in the gripping mechanism test-case */
int holeGripIni(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     o[3],a,b,c,d;
  int        k;
  
  for(k=1; k<=mesh->np; k++) 
    sol->val[k] = 1000.0;  
    
  /* Holes on first edge */
  /*o[0] = 0.0;
  o[1] = 0.0;
  o[2] = 1.5;
  
  a = 0.1;
  b = 0.1;
  c = 0.3;
    
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
        + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }  
  
  o[0] = 0.0;
  o[1] = 0.0;
  o[2] = 0.5;
  
  a = 0.1;
  b = 0.1;
  c = 0.3;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  } */ 
  
  /* Holes on second edge */
  /*o[0] = 0.0;
  o[1] = 1.0;
  o[2] = 1.5;
  
  a = 0.1;
  b = 0.1;
  c = 0.3;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }  
  
  o[0] = 0.0;
  o[1] = 1.0;
  o[2] = 0.5;
  
  a = 0.1;
  b = 0.1;
  c = 0.3;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  } */
  
  /* Holes on third edge */
  o[0] = 1.0;
  o[1] = 0.0;
  o[2] = 1.5;
  
  a = 0.1;
  b = 0.1;
  c = 0.3;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  } 
  
  o[0] = 1.0;
  o[1] = 0.0;
  o[2] = 0.5;
  
  a = 0.1;
  b = 0.1;
  c = 0.3;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  } 
  
  /* Holes on fourth edge */
  o[0] = 1.0;
  o[1] = 1.0;
  o[2] = 1.5;
  
  a = 0.1;
  b = 0.1;
  c = 0.3;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }   
  
  o[0] = 1.0;
  o[1] = 1.0;
  o[2] = 0.5;
  
  a = 0.1;
  b = 0.1;
  c = 0.3;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  } 

  /* Holes on first wall */
  o[0] = 0.5;
  o[1] = 0.0;
  o[2] = 2.0;
  
  a = 0.2;
  b = 0.1;
  c = 0.2;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }  
  
  o[0] = 0.5;
  o[1] = 0.0;
  o[2] = 1.0;
  
  a = 0.2;
  b = 0.1;
  c = 0.2;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  } 
  
  o[0] = 0.5;
  o[1] = 0.0;
  o[2] = 0.0;
  
  a = 0.2;
  b = 0.1;
  c = 0.2;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  } 
  
  /* Holes on second wall */
  o[0] = 0.5;
  o[1] = 1.0;
  o[2] = 2.0;
  
  a = 0.2;
  b = 0.1;
  c = 0.2;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }  
  
  o[0] = 0.5;
  o[1] = 1.0;
  o[2] = 1.0;
  
  a = 0.2;
  b = 0.1;
  c = 0.2;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  } 
  
  o[0] = 0.5;
  o[1] = 1.0;
  o[2] = 0.0;
  
  a = 0.2;
  b = 0.1;
  c = 0.2;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  } 
  
  /* Holes in the middle */
  o[0] = 0.5;
  o[1] = 0.5;
  o[2] = 0.0;
  
  a = 0.2;
  b = 0.2;
  c = 0.2;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  } 
  
  o[0] = 0.5;
  o[1] = 0.5;
  o[2] = 2.0;
  
  a = 0.2;
  b = 0.2;
  c = 0.2;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  } 
  
  o[0] = 0.5;
  o[1] = 0.5;
  o[2] = 1.0;
  
  a = 0.2;
  b = 0.2;
  c = 0.3;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  } 
  
  /* Last series of holes */
  o[0] = 1.0;
  o[1] = 0.5;
  o[2] = 0.0;
  
  a = 0.1;
  b = 0.3;
  c = 0.1;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  } 
  
  o[0] = 1.0;
  o[1] = 0.5;
  o[2] = 2.0;
  
  a = 0.1;
  b = 0.3;
  c = 0.1;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0])/(a*a) + (p0->c[1]-o[1])*(p0->c[1]-o[1])/(b*b) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2])/(c*c);
    d = sqrt(d)-1.0;    
    sol->val[k] = D_MIN(sol->val[k],d);
  } 
  
  for(k=1; k<=mesh->np; k++) 
    sol->val[k] *= -1.0; 
                        
  
  return(1);
}

/* Generate boundary conditions in the L-Beam test-case */
int holeLBeamIni(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     o[3],d,dn,r;
  int        k,n;
  
  for(k=1; k<=mesh->np; k++) 
    sol->val[k] = 1000.0;  
  
  /* Holes on the outer corner edge */
  o[0] = 0.0;
  o[1] = 0.0;
  o[2] = 0.0;
  
  r = 0.1;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
    d = sqrt(d)-r;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  o[0] = 0.0;
  o[1] = 0.5;
  o[2] = 0.0;
  
  r = 0.1;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
    d = sqrt(d)-r;    
    sol->val[k] = D_MIN(sol->val[k],d);
  } 
  
  o[0] = 0.0;
  o[1] = 1.0;
  o[2] = 0.0;
  
  r = 0.1;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
    d = sqrt(d)-r;    
    sol->val[k] = D_MIN(sol->val[k],d);
  } 
    
  /* Holes on the left face */
  for(n=0; n<3; n++) {
    dn = (double)n;
    o[0] = 0.5+dn*0.5;
    o[1] = 0.0;
    o[2] = 0.33333333;
    
    r = 0.1;
    
    for(k=1; k<=mesh->np; k++) { 
      p0 = &mesh->point[k]; 
      d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
      + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
      d = sqrt(d)-r;    
      sol->val[k] = D_MIN(sol->val[k],d);
    }
    
    o[0] = 0.5+dn*0.5;
    o[1] = 0.0;
    o[2] = 0.66666666;
    
    r = 0.1;
    
    for(k=1; k<=mesh->np; k++) { 
      p0 = &mesh->point[k]; 
      d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
      + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
      d = sqrt(d)-r;    
      sol->val[k] = D_MIN(sol->val[k],d);
    }
  } 

  /* Holes on the right face */
  for(n=0; n<3; n++) {
    dn = (double)n;
    o[0] = 0.5+dn*0.5;
    o[1] = 1.0;
    o[2] = 0.33333333;
  
    r = 0.1;
  
    for(k=1; k<=mesh->np; k++) { 
      p0 = &mesh->point[k]; 
      d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
        + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
      d = sqrt(d)-r;    
      sol->val[k] = D_MIN(sol->val[k],d);
    }
    
    o[0] = 0.5+dn*0.5;
    o[1] = 1.0;
    o[2] = 0.66666666;
    
    r = 0.1;
    
    for(k=1; k<=mesh->np; k++) { 
      p0 = &mesh->point[k]; 
      d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
      + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
      d = sqrt(d)-r;    
      sol->val[k] = D_MIN(sol->val[k],d);
    }
  }  
  
  /* Internal holes */
  for(n=0; n<3; n++) {
    dn = (double)n;
    o[0] = 0.5+dn*0.5;
    o[1] = 0.5;
    o[2] = 0.33333333;
    
    r = 0.1;
    
    for(k=1; k<=mesh->np; k++) { 
      p0 = &mesh->point[k]; 
      d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
      + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
      d = sqrt(d)-r;    
      sol->val[k] = D_MIN(sol->val[k],d);
    }
    
    o[0] = 0.5+dn*0.5;
    o[1] = 0.5;
    o[2] = 0.66666666;
    
    r = 0.1;
    
    for(k=1; k<=mesh->np; k++) { 
      p0 = &mesh->point[k]; 
      d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
      + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
      d = sqrt(d)-r;    
      sol->val[k] = D_MIN(sol->val[k],d);
    }
  }  
  

  for(k=1; k<=mesh->np; k++) 
    sol->val[k] *= -1.0; 
  
  
  return(1);
}


/* Generate boundary conditions in the LBB-Beam test-case */
int holeCl_LBBBeam(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     o[3],d,dn,r;
  int        k,n;
  
  for(k=1; k<=mesh->np; k++) 
    sol->val[k] = 1000.0;  
  
  /* Holes on the bottom corners for Dirichlet conditions */
  o[0] = 0.0;
  o[1] = 0.0;
  o[2] = 0.0;
  
  r = 5.0;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
    d = sqrt(d)-r;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  o[0] = 40.0;
  o[1] = 0.0;
  o[2] = 0.0;
  
  r = 5.0;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
    d = sqrt(d)-r;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }
    
  o[0] = 0.0;
  o[1] = 200.0;
  o[2] = 0.0;
  
  r = 5.0;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
    d = sqrt(d)-r;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  o[0] = 40.0;
  o[1] = 200.0;
  o[2] = 0.0;
  
  r = 5.0;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
    d = sqrt(d)-r;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  /* Upper part of the plane for Neumann conditions */
  o[0] = 20.0;
  o[1] = 100.0;
  o[2] = 60.0;
  
  r = 10.0;

  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = fabs(p0->c[1]-o[1]) -r;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  for(k=1; k<=mesh->np; k++) 
    sol->val[k] *= -1.0; 
  
  return(1);
}

/* Generate initial shape in the LBB-Beam test-case */
int holeLBBBeamIni(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     o[3],d,dn,r;
  int        k,n;
  
  for(k=1; k<=mesh->np; k++) 
    sol->val[k] = 1000.0;  
  
  /* Inner holes */
  for(n=1; n<4; n++) {
    dn = (double)n;
    o[0] = 20.0;
    o[1] = 0.0 + dn*50.0;
    o[2] = 30.0;
    
    r = 5.0;
    
    for(k=1; k<=mesh->np; k++) { 
      p0 = &mesh->point[k]; 
      d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
      + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
      d = sqrt(d)-r;    
      sol->val[k] = D_MIN(sol->val[k],d);
    }    
  }
  
  /* Holes on the boundary edges */
  for(n=0; n<4; n++) {
    dn = (double)n;
    o[0] = 0.0;
    o[1] = 40.0 + dn*40.0;
    o[2] = 0.0;
    
    r = 5.0;
    
    for(k=1; k<=mesh->np; k++) { 
      p0 = &mesh->point[k]; 
      d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
      + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
      d = sqrt(d)-r;    
      sol->val[k] = D_MIN(sol->val[k],d);
    }    
  }
  
  for(n=0; n<4; n++) {
    dn = (double)n;
    o[0] = 40.0;
    o[1] = 40.0 + dn*40.0;
    o[2] = 0.0;
    
    r = 5.0;
    
    for(k=1; k<=mesh->np; k++) { 
      p0 = &mesh->point[k]; 
      d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
      + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
      d = sqrt(d)-r;    
      sol->val[k] = D_MIN(sol->val[k],d);
    }    
  }
  
  for(n=0; n<4; n++) {
    dn = (double)n;
    o[0] = 0.0;
    o[1] = 40.0 + dn*40.0;
    o[2] = 60.0;
    
    r = 5.0;
    
    for(k=1; k<=mesh->np; k++) { 
      p0 = &mesh->point[k]; 
      d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
      + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
      d = sqrt(d)-r;    
      sol->val[k] = D_MIN(sol->val[k],d);
    }    
  }
  
  for(n=0; n<4; n++) {
    dn = (double)n;
    o[0] = 40.0;
    o[1] = 40.0 + dn*40.0;
    o[2] = 60.0;
    
    r = 5.0;
    
    for(k=1; k<=mesh->np; k++) { 
      p0 = &mesh->point[k]; 
      d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
      + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
      d = sqrt(d)-r;    
      sol->val[k] = D_MIN(sol->val[k],d);
    }    
  }
      
  for(k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    sol->val[k] = p0->c[2]-30.0;
  }
  
  return(1);
}

/* Generate boundary conditions in the Chair test-case */
int holeCl_Chair(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     o[3],d,dn,r;
  int        k,n;
  
  for(k=1; k<=mesh->np; k++) 
    sol->val[k] = 1000.0;  
  
  /* Holes on the bottom corners for Dirichlet conditions */
  o[0] = 0.0;
  o[1] = 0.0;
  o[2] = 0.0;
  
  r = 0.08;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
    d = sqrt(d)-r;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  o[0] = 0.7;
  o[1] = 0.0;
  o[2] = 0.0;
  
  r = 0.08;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
    d = sqrt(d)-r;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  o[0] = 0.7;
  o[1] = 0.5;
  o[2] = 0.0;
  
  r = 0.08;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
    d = sqrt(d)-r;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  o[0] = 0.0;
  o[1] = 0.5;
  o[2] = 0.0;
  
  r = 0.08;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
    d = sqrt(d)-r;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  /* First Load case */
  o[0] = 0.3;
  o[1] = 0.25;
  o[2] = 0.5;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = o[0]-p0->c[0];    
    if ( fabs(o[2]-p0->c[2]) < 0.01 ) 
      sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  /* Second Load case */
  o[0] = 0.0;
  o[1] = 0.25;
  o[2] = 0.6;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = o[2]-p0->c[2];    
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  for(k=1; k<=mesh->np; k++) 
    sol->val[k] *= -1.0; 
  
  return(1);
}

/* Generate initial shape in the chair test-case */
int holeChairIni(pMesh mesh, pSol sol) {
  pPoint     p0;
  double     o[3],d,dn,r;
  int        k,n;
  
  for(k=1; k<=mesh->np; k++) 
    sol->val[k] = 1000.0;  
  
  /* Bottom hole */
  o[0] = 0.35;
  o[1] = 0.25;
  o[2] = 0.0;
  
  r = 0.1;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
    d = sqrt(d)-r;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }
  
  /* Left-hand hole */
  o[0] = 0.35;
  o[1] = 0.5;
  o[2] = 0.25;
  
  r = 0.1;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
    d = sqrt(d)-r;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }  
  
  /* Right-hand hole */
  o[0] = 0.35;
  o[1] = 0.0;
  o[2] = 0.25;
  
  r = 0.1;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
    d = sqrt(d)-r;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }  
  
  /* Front hole */
  o[0] = 0.7;
  o[1] = 0.25;
  o[2] = 0.25;
  
  r = 0.1;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
    d = sqrt(d)-r;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }  
  
  /* Back hole */
  o[0] = 0.0;
  o[1] = 0.25;
  o[2] = 0.5;
  
  r = 0.1;
  
  for(k=1; k<=mesh->np; k++) { 
    p0 = &mesh->point[k]; 
    d = (p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
    + (p0->c[2]-o[2])*(p0->c[2]-o[2]);
    d = sqrt(d)-r;    
    sol->val[k] = D_MIN(sol->val[k],d);
  }  
  
  for(k=1; k<=mesh->np; k++) 
    sol->val[k] *= -1.0; 
  
  
  return(1);
}


