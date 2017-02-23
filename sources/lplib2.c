

/*----------------------------------------------------------*/
/*															*/
/*						LPLIB	V2.97						*/
/*															*/
/*----------------------------------------------------------*/
/*															*/
/*	Description:		Handles threads, scheduling			*/
/*						& dependencies						*/
/*	Author:				Loic MARECHAL						*/
/*	Creation date:		feb 25 2008							*/
/*	Last modification:	jul 30 2009							*/
/*															*/
/*----------------------------------------------------------*/


#ifndef SERIAL

/*----------------------------------------------------------*/
/* Includes													*/
/*----------------------------------------------------------*/

#include <unistd.h>
#include <signal.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include "lplib2.h"


/*----------------------------------------------------------*/
/* Defines													*/
/*----------------------------------------------------------*/

#define MaxLibPar 10
#define MaxPth 128
#define MaxTyp 100
#define DefNmbSmlBlk 64
#define DefNmbDepBlk 256
#define MaxTotPip 65536
#define MaxPipDep 100

enum ParCmd {RunBigWrk, RunSmlWrk, ClrMem, EndPth};


/*----------------------------------------------------------*/
/* Structures' prototypes									*/
/*----------------------------------------------------------*/

typedef struct WrkSct
{
	int BegIdx, EndIdx, NmbDep, *DepWrdTab;
	int *DepIdxTab;
	struct WrkSct *pre, *nex;
}WrkSct;

typedef struct
{
	int NmbLin, NmbSmlWrk, NmbBigWrk, SmlWrkSiz, BigWrkSiz, DepWrkSiz, NmbDepWrd, *DepWrdMat, *DepIdxMat;
	char *RunDepTab;
	WrkSct *SmlWrkTab, *BigWrkTab;
}TypSct;

typedef struct
{
	int idx;
	char *ClrAdr;
	WrkSct *wrk;
	pthread_mutex_t mtx;
	pthread_cond_t cnd;
	pthread_t pth;
	struct ParSct *par;
}PthSct;

typedef struct PipSct
{
	int idx, NmbDep, DepTab[ MaxPipDep ];
	void *prc, *arg;
	pthread_t pth;
	struct ParSct *par;
}PipSct;

typedef struct ParSct
{
	int NmbCpu, WrkCpt, NmbPip, PenPip, RunPip, NmbTyp, BufMax, BufCpt, req, cmd, ClrLinSiz, *PipWrd;
	double sta[2];
	void (*prc)(int, int, int, void *), *arg;
	pthread_cond_t ParCnd, PipCnd;
	pthread_mutex_t ParMtx, PipMtx;
	pthread_t PipPth;
	PthSct *PthTab;
	TypSct *TypTab, *CurTyp, *DepTyp, *typ1, *typ2;
	WrkSct *NexWrk, *BufWrk[ MaxPth / 4 ];
}ParSct;


/*----------------------------------------------------------*/
/* Private procedures' prototypes							*/
/*----------------------------------------------------------*/

static int SetBit(int *, int);
static int GetBit(int *, int);
static int AndWrd(WrkSct *, char *);
static void SetWrd(WrkSct *, char *);
static void ClrWrd(WrkSct *, char *);
int CmpWrk(const void *, const void *);
static void *PipHdl(void *);
static void *PthHdl(void *);
static WrkSct *NexWrk(ParSct *, int);


/*----------------------------------------------------------*/
/* Global variables											*/
/*----------------------------------------------------------*/

ParSct *ParTab[ MaxLibPar+1 ];
int IniLibPar = 0;


/*----------------------------------------------------------*/
/* Init structures, scheduler and launch threads			*/
/*----------------------------------------------------------*/

int InitParallel(int NmbCpu)
{
	int i, ParIdx;
	ParSct *par = NULL;
	PthSct *pth;

	/* Check the number of requested cpu and clear the main par table at first call */

	if(NmbCpu > MaxPth)
		return(0);

	if(!IniLibPar)
	{
		IniLibPar = 1;

		for(i=1;i<=MaxLibPar;i++)
			ParTab[i] = NULL;
	}

	/* Allocate and build main parallel structure */

	for(ParIdx=1; ParIdx<=MaxLibPar; ParIdx++)
		if(!ParTab[ ParIdx ])
		{
			par = ParTab[ ParIdx ] = calloc(1, sizeof(ParSct));
			break;
		}

	if(!par)
		return(0);

	if(!(par->PthTab = calloc(NmbCpu, sizeof(PthSct))))
		return(0);

	if(!(par->TypTab = calloc((MaxTyp + 1), sizeof(TypSct))))
		return(0);

	if(!(par->PipWrd = calloc(MaxTotPip/32, sizeof(int))))
		return(0);

	par->NmbCpu = NmbCpu;
	par->WrkCpt = par->NmbPip = par->PenPip = par->RunPip = 0;

	/* Set the size of WP buffer */

	if(NmbCpu >= 4)
		par->BufMax = NmbCpu / 4;
	else
		par->BufMax = 1;

	pthread_mutex_init(&par->ParMtx, NULL);
	pthread_mutex_init(&par->PipMtx, NULL);
	pthread_cond_init(&par->ParCnd, NULL);
	pthread_cond_init(&par->PipCnd, NULL);

	/* Launch pthreads */

	for(i=0;i<par->NmbCpu;i++)
	{
		pth = &par->PthTab[i];
		pth->idx = i;
		pth->par = par;
		pthread_mutex_init(&pth->mtx, NULL);
		pthread_cond_init(&pth->cnd, NULL);
		pthread_create(&pth->pth, NULL, PthHdl, (void *)pth);
	}

	/* Wait for all threads to be up and wainting */

	pthread_mutex_lock(&par->ParMtx);

	while(par->WrkCpt < par->NmbCpu)
		pthread_cond_wait(&par->ParCnd, &par->ParMtx);

	pthread_mutex_unlock(&par->ParMtx);

	return(ParIdx);
}


/*----------------------------------------------------------*/
/* Stop all threads and free memories						*/
/*----------------------------------------------------------*/

void StopParallel(int ParIdx)
{
	int i;
	PthSct *pth;
	ParSct *par;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return;

	/* Send stop to all threads */

	pthread_mutex_lock(&par->ParMtx);
	par->cmd = EndPth;
	pthread_mutex_unlock(&par->ParMtx);

	/* Wait for all threads to complete */

	for(i=0;i<par->NmbCpu;i++)
	{
		pth = &par->PthTab[i];
		pthread_mutex_lock(&pth->mtx);
		pthread_cond_signal(&pth->cnd);
		pthread_mutex_unlock(&pth->mtx);
		pthread_join(pth->pth, NULL);
	}

	pthread_mutex_destroy(&par->ParMtx);
	pthread_cond_destroy(&par->ParCnd);

	WaitPipeline(ParIdx);

	pthread_mutex_destroy(&par->PipMtx);
	pthread_cond_destroy(&par->PipCnd);

	/* Free memories */

	for(i=1;i<=MaxTyp;i++)
		if(par->TypTab[i].NmbLin)
			FreeType(ParIdx, i);

	free(par->PthTab);
	free(par->TypTab);
	free(par->PipWrd);
	free(par);

	ParTab[ ParIdx ] = NULL;
}


/*----------------------------------------------------------*/
/* Launch the loop prc on typ1 element depending on typ2	*/
/*----------------------------------------------------------*/

float LaunchParallel(int ParIdx, int TypIdx1, int TypIdx2, void *prc, void *PtrArg)
{
	int i;
	PthSct *pth;
	ParSct *par;
	TypSct *typ1, *typ2 = NULL;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return(-1.);

	/* Check bounds */

	if( (TypIdx1 < 1) || (TypIdx1 > MaxTyp) || (TypIdx2 < 0) || (TypIdx2 > MaxTyp) || (TypIdx1 == TypIdx2) )
		return(-1.);

	typ1 =  &par->TypTab[ TypIdx1 ];

	if(TypIdx2)
	{
		/* Lock acces to global parameters */
    
		pthread_mutex_lock(&par->ParMtx);

		par->cmd = RunSmlWrk;
		par->prc = (void (*)(int, int, int, void *))prc;
		par->arg = PtrArg;
		par->typ1 = typ1;
		par->typ2 = typ2 = &par->TypTab[ TypIdx2 ];
		par->NexWrk = typ1->SmlWrkTab;
		par->BufCpt = 0;
		par->WrkCpt = 0;
		par->sta[0] = par->sta[1] = 0.;
    	par->req = 0;

		/* Clear running wp */

		for(i=0;i<par->NmbCpu;i++)
			par->PthTab[i].wrk = NULL;

		memset(typ1->RunDepTab, 0, typ1->NmbDepWrd * 32 * sizeof(char));

		/* Build a linked list of wp */

		for(i=0;i<par->typ1->NmbSmlWrk;i++)
		{
			typ1->SmlWrkTab[i].pre = &typ1->SmlWrkTab[ i-1 ];
			typ1->SmlWrkTab[i].nex = &typ1->SmlWrkTab[ i+1 ];
		}

		typ1->SmlWrkTab[0].pre = typ1->SmlWrkTab[ typ1->NmbSmlWrk - 1 ].nex = NULL;

		/* Start main loop : wake up threads and wait for completion or blocked threads */

		do
		{
			/* Search for some idle threads */

			par->req = 0;

			for(i=0;i<par->NmbCpu;i++)
			{
				pth = &par->PthTab[i];

				if(pth->wrk)
					continue;

				if(!(pth->wrk = NexWrk(par, i)))
				{
					par->req = 1;
					break;
				}

				/* Wake up the thread and provide it with a WP list */

				pthread_mutex_lock(&pth->mtx);
				pthread_cond_signal(&pth->cnd);
				pthread_mutex_unlock(&pth->mtx);
			}

			/* If every WP are done : exit the parallel loop */

			if(par->WrkCpt == typ1->NmbSmlWrk)
				break;

			/* Otherwise, wait for a blocked thread */

			pthread_cond_wait(&par->ParCnd, &par->ParMtx);
		}while(1);

		pthread_mutex_unlock(&par->ParMtx);

		/* Return the average speedup */

		return(par->sta[1] / par->sta[0]);
	}
	else
	{
		/* Lock acces to global parameters */
	
		pthread_mutex_lock(&par->ParMtx);

		par->cmd = RunBigWrk;
		par->prc = (void (*)(int, int, int, void *))prc;
		par->arg = PtrArg;
		par->typ1 = typ1;
		par->typ2 = NULL;
		par->WrkCpt = 0;

		for(i=0;i<typ1->NmbBigWrk;i++)
		{
			pth = &par->PthTab[i];
			pth->wrk = &typ1->BigWrkTab[i];
		}

		for(i=0;i<typ1->NmbBigWrk;i++)
		{
			pth = &par->PthTab[i];
			pthread_mutex_lock(&pth->mtx);
			pthread_cond_signal(&pth->cnd);
			pthread_mutex_unlock(&pth->mtx);
		}

		pthread_cond_wait(&par->ParCnd, &par->ParMtx);

		pthread_mutex_unlock(&par->ParMtx);

		/* Return the average speedup */

		return(par->NmbCpu);
	}
}


/*----------------------------------------------------------*/
/* Pthread handler, waits for job, does it, then signal end	*/
/*----------------------------------------------------------*/

static void *PthHdl(void *ptr)
{
	PthSct *pth = (PthSct *)ptr;
	ParSct *par = pth->par;

	/* Send prap'n ready to the scheduler */ 

	pthread_mutex_lock(&par->ParMtx);
	par->WrkCpt++;
	pthread_cond_signal(&par->ParCnd);
	pthread_mutex_unlock(&par->ParMtx);

	/* Enter main loop until StopParallel is send */

	pthread_mutex_lock(&pth->mtx);

	do
	{
		/* Wait for a wake-up signal from the main loop */

		pthread_cond_wait(&pth->cnd, &pth->mtx);

		switch(par->cmd)
		{
			case RunBigWrk :
			{
				/* Launch a single big wp and signal completion to the scheduler */

				par->prc(pth->wrk->BegIdx, pth->wrk->EndIdx, pth->idx, par->arg);

				pthread_mutex_lock(&par->ParMtx);
				par->WrkCpt++;

				if(par->WrkCpt >= par->typ1->NmbBigWrk)
					pthread_cond_signal(&par->ParCnd);

				pthread_mutex_unlock(&par->ParMtx);
			}break;

			case RunSmlWrk :
			{
				do
				{
					/* Run the WP */

					par->prc(pth->wrk->BegIdx, pth->wrk->EndIdx, pth->idx, par->arg);

					/* Locked acces to global parameters : update WP count, tag WP done and signal the main loop */

					pthread_mutex_lock(&par->ParMtx);

					par->WrkCpt++;

					if(!(pth->wrk = NexWrk(par, pth->idx)))
					{
						par->req = 1;
						pthread_cond_signal(&par->ParCnd);
						pthread_mutex_unlock(&par->ParMtx);
						break;
					}

					if(par->req)
						pthread_cond_signal(&par->ParCnd);

					pthread_mutex_unlock(&par->ParMtx);
				}while(1);
			}break;

			case ClrMem :
			{
				/* Clear memory and signal completion to the scheduler */

				memset(pth->ClrAdr, 0, par->ClrLinSiz);

				pthread_mutex_lock(&par->ParMtx);
				par->WrkCpt++;
				pthread_cond_signal(&par->ParCnd);
				pthread_mutex_unlock(&par->ParMtx);
			}break;

			case EndPth :
			{
				/* Destroy the thread mutex and condition and call for join */

				pthread_mutex_unlock(&pth->mtx);
				pthread_mutex_destroy(&pth->mtx);
				pthread_cond_destroy(&pth->cnd);
				return(NULL);
			}break;
		}
	}while(1);

	return(NULL);
}


/*----------------------------------------------------------*/
/* Get the next WP to be computed							*/
/*----------------------------------------------------------*/

static WrkSct *NexWrk(ParSct *par, int PthIdx)
{
	int i;
	PthSct *pth = &par->PthTab[ PthIdx ];
	WrkSct *wrk;

	/* Update stats */

	par->sta[0]++;

	for(i=0;i<par->NmbCpu;i++)
		if(par->PthTab[i].wrk)
			par->sta[1]++;

	/* Remove previous work's tags */

	if(pth->wrk)
		ClrWrd(pth->wrk, par->typ1->RunDepTab);

	/* If the wp's buffer is empty search for some new compatible wp to fill in */

	if(!par->BufCpt)
	{
		wrk = par->NexWrk;

		while(wrk)
		{
			/* Check for dependencies */
    
			if(!AndWrd(wrk, par->typ1->RunDepTab))
			{
				par->BufWrk[ par->BufCpt++ ] = wrk;

				/* Unlink wp */

				if(wrk->pre)
					wrk->pre->nex = wrk->nex;
				else
					par->NexWrk = wrk->nex;

				if(wrk->nex)
					wrk->nex->pre = wrk->pre;

				/* Add new work's tags */

				SetWrd(wrk, par->typ1->RunDepTab);

				if(par->BufCpt == par->BufMax)
					break;
			}

			wrk = wrk->nex;
		}
	}

	/* Return the next available wp in buffer and unlink it from the todo list */

	if(par->BufCpt)
		return(par->BufWrk[ --par->BufCpt ]);
	else
		return(NULL);
}


/*----------------------------------------------------------*/
/* Allocate a new kind of elements and set work-packages	*/
/*----------------------------------------------------------*/

int NewType(int ParIdx, int NmbLin)
{
	int i, TypIdx=0, idx;
	TypSct *typ;
	ParSct *par;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return(0);

	if(NmbLin < par->NmbCpu)
		return(0);

	/* Search for a free type structure */

	for(i=1;i<=MaxTyp;i++)
		if(!par->TypTab[i].NmbLin)
		{
			TypIdx = i;
			break;
		}

	if(!TypIdx)
		return(0);

	typ = &par->TypTab[ TypIdx ];
	typ->NmbLin = NmbLin;

	/* Compute the size of small work-packages */

	if(NmbLin >= DefNmbSmlBlk * par->NmbCpu)
		typ->SmlWrkSiz = NmbLin / (DefNmbSmlBlk * par->NmbCpu);
	else
		typ->SmlWrkSiz = NmbLin / par->NmbCpu;

	typ->NmbSmlWrk = NmbLin / typ->SmlWrkSiz;

	if(NmbLin != typ->NmbSmlWrk * typ->SmlWrkSiz)
		typ->NmbSmlWrk++;

	if(!(typ->SmlWrkTab = calloc(typ->NmbSmlWrk , sizeof(WrkSct))))
		return(0);

	/* Set small work-packages */

	idx = 0;

	for(i=0;i<typ->NmbSmlWrk;i++)
	{
		typ->SmlWrkTab[i].BegIdx = idx + 1;
		typ->SmlWrkTab[i].EndIdx = idx + typ->SmlWrkSiz;
		idx += typ->SmlWrkSiz;
	}

	typ->SmlWrkTab[ typ->NmbSmlWrk - 1 ].EndIdx = NmbLin;

	/* Compute the size of big work-packages */

	typ->BigWrkSiz = NmbLin / par->NmbCpu;
	typ->NmbBigWrk = par->NmbCpu;

	if(!(typ->BigWrkTab = calloc(typ->NmbBigWrk , sizeof(WrkSct))))
		return(0);

	/* Set big work-packages */

	idx = 0;

	for(i=0;i<typ->NmbBigWrk;i++)
	{
		typ->BigWrkTab[i].BegIdx = idx + 1;
		typ->BigWrkTab[i].EndIdx = idx + typ->BigWrkSiz;
		idx += typ->BigWrkSiz;
	}

	typ->BigWrkTab[ typ->NmbBigWrk - 1 ].EndIdx = NmbLin;

	return(TypIdx);
}


/*----------------------------------------------------------*/
/* Add this kind of element to the free-list				*/
/*----------------------------------------------------------*/

void FreeType(int ParIdx, int TypIdx)
{
	TypSct *typ;
	ParSct *par;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return;

	/* Check bounds and free mem */

	if( (TypIdx < 1) || (TypIdx > MaxTyp) )
		return;

	typ = &par->TypTab[ TypIdx ];

	if(typ->SmlWrkTab)
		free(typ->SmlWrkTab);

	if(typ->BigWrkTab)
		free(typ->BigWrkTab);

	if(typ->DepIdxMat)
		free(typ->DepIdxMat);

	if(typ->RunDepTab)
		free(typ->RunDepTab);

	if(typ->DepWrdMat)
		free(typ->DepWrdMat);

	memset(typ, 0, sizeof(TypSct));
}


/*----------------------------------------------------------*/
/* Allocate a dependency matrix linking both types			*/
/*----------------------------------------------------------*/

int BeginDependency(int ParIdx, int TypIdx1, int TypIdx2)
{
	int i;
	ParSct *par;
	TypSct *typ1, *typ2;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return(0);

	/* Check bounds */

	par->CurTyp = typ1 = &par->TypTab[ TypIdx1 ];
	par->DepTyp = typ2 = &par->TypTab[ TypIdx2 ];

	if( (TypIdx1 < 1) || (TypIdx1 > MaxTyp) || (TypIdx2 < 1) || (TypIdx2 > MaxTyp) || (typ1 == typ2) \
	|| !typ1->NmbLin || !typ2->NmbLin)
	{
		return(0);
	}

	/* Compute dependency table's size */

	if(typ2->NmbLin >= DefNmbDepBlk * par->NmbCpu)
		typ1->DepWrkSiz = typ2->NmbLin / (DefNmbDepBlk * par->NmbCpu);
	else
		typ1->DepWrkSiz = typ2->NmbLin / par->NmbCpu;

	typ1->NmbDepWrd = typ2->NmbLin / (typ1->DepWrkSiz * 32);

	if(typ2->NmbLin != typ1->NmbDepWrd * typ1->DepWrkSiz * 32)
		typ1->NmbDepWrd++;

	/* Allocate a global dependency table */

	if(!(typ1->DepWrdMat = calloc(typ1->NmbSmlWrk * typ1->NmbDepWrd, sizeof(int))))
		return(0);

	/* Then spread sub-tables among WP */

	for(i=0;i<typ1->NmbSmlWrk;i++)
	{
		typ1->SmlWrkTab[i].NmbDep = 0;
		typ1->SmlWrkTab[i].DepWrdTab = &typ1->DepWrdMat[ i * typ1->NmbDepWrd ];
	}

	/* Allocate a running tags table */

	if(!(typ1->RunDepTab = calloc(typ1->NmbDepWrd * 32, sizeof(char))))
		return(0);

	return(1);
}


/*----------------------------------------------------------*/
/* Type1 element idx1 depends on type2 element idx2			*/
/*----------------------------------------------------------*/

void AddDependency(int ParIdx, int idx1, int idx2)
{
	WrkSct *wrk;
	ParSct *par;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return;

	/* Set and count dependency bit */

	wrk = &par->CurTyp->SmlWrkTab[ (idx1-1) / par->CurTyp->SmlWrkSiz ];

	if(!SetBit(wrk->DepWrdTab, (idx2-1) / par->CurTyp->DepWrkSiz ))
		wrk->NmbDep++;
}


/*----------------------------------------------------------*/
/* Sort wp depending on their number of dependencies		*/
/*----------------------------------------------------------*/

void EndDependency(int ParIdx, float DepSta[2])
{
	int i, j, idx=0, NmbDep, NmbDepBit, TotNmbDep = 0;
	ParSct *par;
	TypSct *typ1, *typ2;
	WrkSct *wrk;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return;

	/* Compute average number of collisions */

	DepSta[1] = 0.;
	typ1 = par->CurTyp;
	typ2 = par->DepTyp;

	for(i=0;i<typ1->NmbSmlWrk;i++)
	{
		TotNmbDep += typ1->SmlWrkTab[i].NmbDep;

		if(typ1->SmlWrkTab[i].NmbDep > DepSta[1])
			DepSta[1] = typ1->SmlWrkTab[i].NmbDep;
	}

	DepSta[0] = TotNmbDep;

	/* Allocate a global dependency index table */

	if(!(typ1->DepIdxMat = calloc(TotNmbDep, sizeof(int))))
		return;

	/* Then spread and fill the sub-tables among WP */

	NmbDep = typ1->NmbDepWrd * 32;

	for(i=0;i<typ1->NmbSmlWrk;i++)
	{
		wrk = &typ1->SmlWrkTab[i];
		wrk->DepIdxTab = &typ1->DepIdxMat[ idx ];
		idx += wrk->NmbDep;
		wrk->NmbDep = 0;

		for(j=0;j<NmbDep;j++)
			if(GetBit(wrk->DepWrdTab, j))
				wrk->DepIdxTab[ wrk->NmbDep++ ] = j;
	}

	/* Compute stats */

	NmbDepBit = typ2->NmbLin / typ1->DepWrkSiz;

	if(typ2->NmbLin - NmbDepBit * typ1->DepWrkSiz)
		NmbDepBit++;

	DepSta[0] = 100 * DepSta[0] / (typ1->NmbSmlWrk * NmbDepBit);
	DepSta[1] = 100 * DepSta[1] / NmbDepBit;

	/* Sort WP from highest collision number to the lowest */

	qsort(typ1->SmlWrkTab, typ1->NmbSmlWrk, sizeof(WrkSct), CmpWrk);

	par->CurTyp = NULL;
}


/*----------------------------------------------------------*/
/* Test and set a bit in a multibyte number					*/
/*----------------------------------------------------------*/

static int SetBit(int *tab, int idx)
{
	int res = ( tab[ idx >> 5 ] & (1 << (idx & 31)) );
	tab[ idx >> 5 ] |= 1 << (idx & 31);
	return(res);
}


/*----------------------------------------------------------*/
/* Test a bit in a multibyte number							*/
/*----------------------------------------------------------*/

static int GetBit(int *tab, int idx)
{
	return( tab[ idx >> 5 ] & (1 << (idx & 31)) );
}


/*----------------------------------------------------------*/
/* Check wether two WP share common resources -> locked		*/
/*----------------------------------------------------------*/

static int AndWrd(WrkSct *wrk, char *wrd)
{
	int i;

	for(i=0;i<wrk->NmbDep;i++)
		if(wrd[ wrk->DepIdxTab[i] ])
			return(1);

	return(0);
}

static void SetWrd(WrkSct *wrk, char *wrd)
{
	int i;

	for(i=0;i<wrk->NmbDep;i++)
		wrd[ wrk->DepIdxTab[i] ] = 1;
}

static void ClrWrd(WrkSct *wrk, char *wrd)
{
	int i;

	for(i=0;i<wrk->NmbDep;i++)
		wrd[ wrk->DepIdxTab[i] ] = 0;
}


/*----------------------------------------------------------*/
/* Compare two workpackages number of bits					*/
/*----------------------------------------------------------*/

int CmpWrk(const void *ptr1, const void *ptr2)
{
	WrkSct *w1, *w2;

	w1 = (WrkSct *)ptr1;
	w2 = (WrkSct *)ptr2;

	if(w1->NmbDep > w2->NmbDep)
		return(-1);
	else if(w1->NmbDep < w2->NmbDep)
		return(1);
	else
		return(0);
}


/*----------------------------------------------------------*/
/* Launch the loop prc on typ1 element depending on typ2	*/
/*----------------------------------------------------------*/

int ParallelMemClear(int ParIdx, void *PtrArg, long siz)
{
	char *tab = (char *)PtrArg;
	int i;
	PthSct *pth;
	ParSct *par;

	/* Get and check lib parallel instance, adresse and size */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) || !tab || (siz < par->NmbCpu) )
		return(0);

	/* Lock acces to global parameters */
	
	pthread_mutex_lock(&par->ParMtx);

	par->cmd = ClrMem;
	par->ClrLinSiz = siz / par->NmbCpu;
	par->WrkCpt = 0;

	/* Spread the buffer among each thread and wake then up */

	for(i=0;i<par->NmbCpu;i++)
	{
		pth = &par->PthTab[i];
		pth->ClrAdr = &tab[ i * par->ClrLinSiz ];

		pthread_mutex_lock(&pth->mtx);
		pthread_cond_signal(&pth->cnd);
		pthread_mutex_unlock(&pth->mtx);
	}

	/* Wait for each thread to complete */

	while(par->WrkCpt < par->NmbCpu)
		pthread_cond_wait(&par->ParCnd, &par->ParMtx);

	pthread_mutex_unlock(&par->ParMtx);

	return(1);
}


/*----------------------------------------------------------*/
/* Wait for a condition, launch and detach a user procedure	*/
/*----------------------------------------------------------*/

int LaunchPipeline(int ParIdx, void *prc, void *PtrArg, int NmbDep, int *DepTab)
{
	int i;
	PipSct *NewPip=NULL;
	ParSct *par;

	/* Get and check lib parallel instance and the number of pipes and dependencies */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) \
	||	(NmbDep > MaxPipDep) || (par->NmbPip >= MaxTotPip) )
	{
		return(0);
	}

	/* Allocate and setup a new pipe */

	if(!(NewPip = malloc(sizeof(PipSct))))
		return(0);

	NewPip->prc = prc;
	NewPip->arg = PtrArg;
	NewPip->par = par;
	NewPip->NmbDep = NmbDep;

	for(i=0;i<NmbDep;i++)
		NewPip->DepTab[i] = DepTab[i];

	/* Lock pipe mutex, increment pipe counter and launch the pipe regardless dependencies */

	pthread_mutex_lock(&par->PipMtx);
	NewPip->idx = ++par->NmbPip;
	par->PenPip++;
	pthread_create(&NewPip->pth, NULL, PipHdl, (void *)NewPip);
	pthread_mutex_unlock(&par->PipMtx);

	return(NewPip->idx);
}


/*----------------------------------------------------------*/
/* Thread handler launching and waitint for user's procedure*/
/*----------------------------------------------------------*/

static void *PipHdl(void *ptr)
{
	int RunFlg, i;
	PipSct *pip = (PipSct *)ptr;
	ParSct *par = pip->par;
	void (*prc)(void *);

	/* Wait for conditions to be met */

	do
	{
		pthread_mutex_lock(&par->PipMtx);

		if(par->RunPip < par->NmbCpu)
		{
			RunFlg = 1;

			for(i=0;i<pip->NmbDep;i++)
				if(!GetBit(par->PipWrd, pip->DepTab[i]))
				{
					RunFlg = 0;
					break;
				}
		}

		if(!RunFlg)
		{
			pthread_mutex_unlock(&par->PipMtx);
			usleep(1000);
		}
	}while(!RunFlg);

	/* Execute the user's procedure and set the flag to 2 (done) */

	prc = (void (*)(void *))pip->prc;
	par->RunPip++;

	pthread_mutex_unlock(&par->PipMtx);

	prc(pip->arg);

	pthread_mutex_lock(&par->PipMtx);
	SetBit(par->PipWrd, pip->idx);
	par->PenPip--;
	par->RunPip--;
	free(pip);
	pthread_mutex_unlock(&par->PipMtx);

	return(NULL);
}


/*----------------------------------------------------------*/
/* Wait for all pipelined procedures to complete			*/
/*----------------------------------------------------------*/

void WaitPipeline(int ParIdx)
{
	int PenPip;
	ParSct *par;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return;

	do
	{
		pthread_mutex_lock(&par->PipMtx);
		PenPip = par->PenPip;
		pthread_mutex_unlock(&par->PipMtx);
		usleep(1000);
	}while(PenPip);
}


#else


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libparallel2.h"

#define MaxLibPar 10
#define MaxPth 128
#define MaxTyp 100

typedef struct
{
	int NmbLin;
}TypSct;

typedef struct ParSct
{
	int NmbTyp;
	float sta[2];
	void (*prc)(int, int, int, void *), *arg;
	TypSct TypTab[ MaxTyp+1 ];
}ParSct;

ParSct *ParTab[ MaxLibPar+1 ];
int IniLibPar = 0;


int InitParallel(int NmbCpu)
{
	int i, ParIdx;
	ParSct *par = NULL;

	if(NmbCpu > MaxPth)
		return(0);

	if(!IniLibPar)
	{
		IniLibPar = 1;

		for(i=1;i<=MaxLibPar;i++)
			ParTab[i] = NULL;
	}

	/* Allocate and build main parallel structure */

	for(ParIdx=1; ParIdx<=MaxLibPar; ParIdx++)
		if(!ParTab[ ParIdx ])
		{
			par = ParTab[ ParIdx ] = calloc(1, sizeof(ParSct));
			break;
		}

	if(!par)
		return(0);

	return(ParIdx);
}


void StopParallel(int ParIdx)
{
	int i;
	ParSct *par;

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return;

	for(i=1;i<=MaxTyp;i++)
		FreeType(ParIdx, i);

	free(par);
	ParTab[ ParIdx ] = NULL;
}


int NewType(int ParIdx, int NmbLin)
{
	int i, TypIdx=0;
	TypSct *typ;
	ParSct *par;

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return(0);

	if(NmbLin <= 0)
		return(0);

	for(i=1;i<=MaxTyp;i++)
		if(!par->TypTab[i].NmbLin)
		{
			TypIdx = i;
			break;
		}

	if(!TypIdx)
		return(0);

	typ = &par->TypTab[ TypIdx ];
	typ->NmbLin = NmbLin;

	return(TypIdx);
}


void FreeType(int ParIdx, int TypIdx)
{
}


int BeginDependency(int ParIdx, int TypIdx1, int TypIdx2)
{
	return(0);
}


void AddDependency(int ParIdx, int idx1, int idx2)
{
}


void EndDependency(int ParIdx, float DepSta[2])
{
}


float LaunchParallel(int ParIdx, int typ1, int typ2, void *prc, void *PtrArg)
{
	ParSct *par;
	void (*UsrPrc)(int, int, int, void *) = (void (*)(int, int, int, void *))prc;

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return(-1.);

	if( (typ1 < 1) || (typ1 > MaxTyp) || (typ2 < 0) || (typ2 > MaxTyp) || (typ1 == typ2) )
		return(-1.);

	UsrPrc(1, par->TypTab[ typ1 ].NmbLin, 0, PtrArg);

	return(1.);
}


int ParallelMemClear(int ParIdx, void *PtrArg, long siz)
{
	memset(PtrArg, 0, siz);
	return(1);
}

#endif
