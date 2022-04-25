#include "mshdist.h"

extern pBucket (*newBucket)(pMesh ,int ) = NULL;
extern int     (*closept)(pMesh ,double *) = NULL;
extern int     (*inidist)(Info info,pMesh ,pMesh ,pSol ,pBucket ) = NULL;
extern int     (*inidistpcloud)(pMesh ,pMesh ,pSol ,pBucket ) = NULL;
extern int     (*iniredist)(Info ,pMesh ,pSol ) = NULL;
extern int     (*iniencdomain)(Info ,pMesh ,pSol ) = NULL;
extern int     (*inireftrias)(Info ,pMesh, pSol) = NULL;
extern int     (*ppgdist)(Info info,pMesh mesh, pSol sol) = NULL;
extern int     (*ppgdistfmm)(pMesh mesh, pSol sol) = NULL;
