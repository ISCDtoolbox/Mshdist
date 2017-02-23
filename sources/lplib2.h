

/*----------------------------------------------------------*/
/*															*/
/*						LPLIB	V2.94						*/
/*															*/
/*----------------------------------------------------------*/
/*															*/
/*	Description:		Handles scheduling & threads launch	*/
/*	Author:				Loic MARECHAL						*/
/*	Creation date:		feb 25 2008							*/
/*	Last modification:	oct 13 2008							*/
/*															*/
/*----------------------------------------------------------*/


/*----------------------------------------------------------*/
/* User available procedures' prototypes					*/
/*----------------------------------------------------------*/

int InitParallel(int);
void StopParallel(int);
int NewType(int, int);
void FreeType(int, int);
int BeginDependency(int, int, int);
void AddDependency(int, int, int);
void EndDependency(int, float [2]);
float LaunchParallel(int, int, int, void *, void *);
int LaunchPipeline(int, void *, void *, int, int *);
void WaitPipeline(int);
int ParallelMemClear(int , void *, long);
