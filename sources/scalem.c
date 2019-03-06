#include "mshdist.h"


/* warning: rectangular bounding box */
int scaleMesh(Info *info,pMesh mesh1,pMesh mesh2,pSol sol1) {
  pPoint    ppt;
  double    dd,deltb,delta;
  int       i,k;

  /* find bounding box */
	for (i=0; i<mesh1->dim; i++) {
		info->min1[i] = info->min2[i] =  1.e30;
		info->max1[i] = info->max2[i] = -1.e30;
	}
  
	/* Get the size of mesh1 in every direction */ 	
	for (k=1; k<=mesh1->np; k++) {
		ppt = &mesh1->point[k];
		for (i=0; i<mesh1->dim; i++) {
			if ( ppt->c[i] > info->max1[i] )  info->max1[i] = ppt->c[i];
			if ( ppt->c[i] < info->min1[i] )  info->min1[i] = ppt->c[i];
		}
	}
  
	delta = 0.0;
	for (i=0; i<mesh1->dim; i++) {
		dd = fabs(info->max1[i]-info->min1[i]);
		if ( dd > delta )  delta = dd;
	}
	if ( delta < EPS1 ) {
		fprintf(stdout,"  ## Unable to scale mesh\n");
		return(0);
	}

	/* Normalize coordinates of mesh1 */
	dd = 1.0 / delta;
	for (k=1; k<=mesh1->np; k++) {
		ppt = &mesh1->point[k];
		for (i=0; i<mesh1->dim; i++)
			ppt->c[i] = dd * (ppt->c[i] - info->min1[i]);
	}
  
  /* Store scaling */
  for (i=0; i<mesh1->dim; i++)
    info->delta1[i] =  (info->max1[i]-info->min1[i]);  /* Bug fix 27/12/2015: old dd* */
  
  /* scale initial solution */
  if ( info->option == 2 ) {
    for (k=1; k<=sol1->np; k++) {
	    sol1->val[k] *= dd;
  	}
    return(1);
  }
  
  /* scale starting points */
  if ( info->nexp > 0 ) {
    for (k=0; k<info->nexp; k++) {
      for (i=0; i<mesh1->dim; i++) {
        info->exp[mesh1->dim*k+i] = dd*(info->exp[mesh1->dim*k+i] - info->min1[i]);
      }
    }
  }
  /* Assign a default value */
  else {
    info->nexp = 1;
    for (i=0; i<mesh1->dim; i++)
      info->exp[i] = 0.01;
  }

  if ( info->option == 3 )
    return(1);
  
  /* 2nd mesh is put to scale SIZE with respect to mesh 1 */
  /* store in info->max/min the max/min coordinate of all the points of mesh2 */
  for (k=1; k<=mesh2->np; k++) {
    ppt = &mesh2->point[k];
    for (i=0; i<mesh2->dim; i++) {
      if ( ppt->c[i] > info->max2[i] )  info->max2[i] = ppt->c[i];
      if ( ppt->c[i] < info->min2[i] )  info->min2[i] = ppt->c[i];
    }
  }
  
  /* In the noscale situation, use the same operation for mesh2 as for mesh1 */
  if ( info->noscale ) {
    dd = 1.0 / delta;
    for (k=1; k<=mesh2->np; k++) {
      ppt = &mesh2->point[k];
      for (i=0; i<mesh2->dim; i++)
        ppt->c[i] = dd * (ppt->c[i] - info->min1[i]);
    }
    return(1);
  }
  
  /* Get the size of mesh2 in every direction */ 	
  deltb = 0.0;
  for (i=0; i<mesh2->dim; i++) {
    info->delta2[i] = fabs(info->max2[i]-info->min2[i]);
    if ( info->delta2[i] > deltb )  deltb = info->delta2[i];
  }
  if ( deltb < EPS1 ) {
    fprintf(stdout,"  ## Unable to scale mesh2\n");
    return(0);
  }
  
  /* centering */	
  /*deltb = 0.0;
  for (i=0; i<mesh2->dim; i++) {
  	info->delta2[i] = (info->max2[i]-info->min2[i]);
		dd = info->delta2[i] / info->delta1[i];
		if ( dd > deltb )  deltb = dd;
  }*/
  

    /*printf("min %f %f %f  max %f %f %f\n",info->min1[0],info->min1[1],info->min1[2],info->max1[0],info->max1[1],info->max1[2]);
  	printf("min %f %f %f  max %f %f %f\n",info->min2[0],info->min2[1],info->min2[2],info->max2[0],info->max2[1],info->max2[2]);
    printf("delta1 %f %f %f \n",info->delta1[0],info->delta1[1],info->delta1[2]);
    printf("delta2 %f %f %f \n",info->delta2[0],info->delta2[1],info->delta2[2]);
    printf("deltb = %f\n",deltb);*/
  
  /* Resize mesh2 to SIZE*mesh1, and locate it at the centre of mesh1 */
  dd = info->size / deltb;
  for (k=1; k<=mesh2->np; k++) {
    ppt = &mesh2->point[k];
    for (i=0; i<mesh2->dim; i++) {
      ppt->c[i]  = dd * (ppt->c[i]-info->min2[i]);
      ppt->c[i] += 0.5*(1.0-dd*info->delta2[i]);
      //ppt->c[i] -= dd * 0.5*info->delta2[i];
      //ppt->c[i] +=      0.5*info->delta1[i];
    }
  }
  
  return(1);
}

/* Unscale solution with respect to the stored information in scaleMesh*/
int unscaleSol(Info info,pSol sol) {
  double     dd;
  int        k;
  char       i;
  
  /* de-normalize solution */
  dd = 0.0;
  for (i=0; i<sol->dim; i++){
    if ( info.delta1[i] > dd )  dd = info.delta1[i];
  }
  
  dd = dd / (double)PRECI;  
  for (k=1; k<=sol->np; k++) {
    sol->val[k] *= dd;
  }
  
  return(1);
}



