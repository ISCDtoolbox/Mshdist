#ifndef _DISTEXTERNS_H
#define _DISTEXTERNS_H

#include "mshdist.h"

extern pBucket (*newBucket)(pMesh ,int );
extern int     (*closept)(pMesh ,double *);
extern int     (*inidist)(Info info,pMesh ,pMesh ,pSol ,pBucket );
extern int     (*inidistpcloud)(pMesh ,pMesh ,pSol ,pBucket );
extern int     (*iniredist)(Info ,pMesh ,pSol );
extern int     (*iniencdomain)(Info ,pMesh ,pSol );
extern int     (*inireftrias)(Info ,pMesh, pSol);
extern int     (*ppgdist)(Info info,pMesh mesh, pSol sol);
extern int     (*ppgdistfmm)(pMesh mesh, pSol sol);

#endif
