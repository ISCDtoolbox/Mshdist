#include "mshdist.h"
#include "mshdistexterns.h"

#include "compil.date"

#define BUCKSIZ       16


/* Exceptions */
static void excfun(int sigid) {
  fprintf(stdout,"\n Unexpected error:");  fflush(stdout);
  switch(sigid) {
    case SIGABRT:
      fprintf(stdout,"  Abnormal stop\n");  exit(1);
    case SIGFPE:
      fprintf(stdout,"  Floating-point exception\n"); exit(1);
    case SIGILL:
      fprintf(stdout,"  Illegal instruction\n"); exit(1);
    case SIGSEGV:
      fprintf(stdout,"  Segmentation fault\n");  exit(1);
    case SIGTERM:
    case SIGINT:
      fprintf(stdout,"  Program killed\n");  exit(1);
  }
  exit(1);
}

/* Manual */
static void usage(char *prog) {
  fprintf(stdout,"usage: %s [-v[n]] [-h] file1[.mesh] [file2[.mesh]] options\n",prog);

  fprintf(stdout,"\n** Generic options :\n");
  fprintf(stdout,"-d      Turn on debug mode\n");
  fprintf(stdout,"-h      Print this message\n");
  fprintf(stdout,"-dt     Time stepping (hmin)\n");
  fprintf(stdout,"-it n   Max number of iterations\n");
  fprintf(stdout,"-ncpu n Use n CPUs\n");
  fprintf(stdout,"-r res  Residual\n");
  fprintf(stdout,"-v [n]  Tune level of verbosity\n");

  exit(1);
}

/* Parse arguments on the command line */
static int parsar(int argc,char *argv[],Info *info,pMesh mesh1,pSol sol1,pMesh mesh2) {
  int      i,ier;
  char    *ptr;

  i = 1;

  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {

    case '?':
      usage(argv[0]);
      break;

	  case 'b':
		if ( !strcmp(argv[i],"-bbbc") ) {
		  info->bbbc = 1;
		}
		break;

      case 'd':
        if ( !strcmp(argv[i],"-dt") ) {
          ++i;
          if ( i < argc && isdigit(argv[i][0]) )
            info->dt = atof(argv[i]);
          else
            --i;
        }
		else if ( !strcmp(argv[i],"-dom") ) {
		  info->option = 3;
		  ++i;
		  if ( i < argc && isdigit(argv[i][0]) )
			info->ref = atoi(argv[i]);
		  else{
			--i;
		  }
		}
		else if ( !strcmp(argv[i],"-d") )
          info->ddebug = 1;
		break;
      
      /* Calculation of the distance field with the Fast marching method */
      case 'f':
        if ( !strcmp(argv[i],"-fmm") )
          info->fmm = 1;
        else if ( !strcmp(argv[i],"-fini") )
          info->fini = 1;
        break;

      /* Calculate Hausdorff distance */
	  case 'h':
		if ( !strcmp(argv[i],"-hausdorff") ) {
		  info->hausdorff = 1;
		}
        else
          usage(argv[0]);
		break;

      case 'i':
        if ( ++i < argc ) {
          if ( isdigit(argv[i][0]) ) {
            info->maxit = atoi(argv[i]);
            info->maxit = D_MAX(10,info->maxit);
          }
          else
            i--;
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          usage(argv[0]);
        }
        break;
          
      case 'n':
        if ( !strcmp(argv[i],"-ncpu") ) {
          ++i;
          if ( i < argc && isdigit(argv[i][0]) )
            info->ncpu = atoi(argv[i]);
          else
            --i;
        }
        else if ( !strcmp(argv[i],"-noscale") )
          info->noscale = 1;
          
        break;
          
      /* Generate unsigned distance with respect to a point cloud */
      case 'p':
        if ( !strcmp(argv[i],"-pcloud") )
          info->pcloud = 1;
        break;

      case 'r':
        if ( ++i < argc ) {
          if ( isdigit(argv[i][0]) )
            info->res = atof(argv[i]);
          else
            i--;
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          usage(argv[0]);
        }
        break;

      /* Generate a particular, analytical implicit function */
      case 's':
        if ( !strcmp(argv[i],"-specdist") ) {
          info->specdist = 1;
          info->option = 3;
        }
        else if ( !strcmp(argv[i],"-scale") ) {
          ++i;
          if ( i < argc && isdigit(argv[i][0]) )
            info->size = atof(argv[i]);
          else
            i--;
        }
        else if ( !strcmp(argv[i],"-sref") ) {
            info->startref = 1;
            info->option = 3;
        }
        else if ( !strcmp(argv[i],"-surf") )
          info->dsurf = 1;
          
        break;

      case 'v':
        if ( ++i < argc ) {
          if ( argv[i][0] == '-' || isdigit(argv[i][0]) )
            info->imprim = atoi(argv[i]);
          else
            i--;
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          usage(argv[0]);
        }
        break;

      default:
        fprintf(stderr,"  Unrecognized option %s\n",argv[i]);
        usage(argv[0]);
        break;
      }
    }

    else {
      if ( mesh1->name == NULL ) {
        mesh1->name = argv[i];
        if ( info->imprim == -99 )  info->imprim = 5;
      }
      else if ( mesh2->name == NULL )
        mesh2->name = argv[i];
      else {
        fprintf(stdout,"  Argument %s ignored\n",argv[i]);
        usage(argv[0]);
      }
    }
    i++;
  }

  /* Check parameters */
  if ( mesh1->name == NULL  || info->imprim == -99 ) {
    fprintf(stdout,"\n  -- PRINT (0 10(advised) -10) ?\n");
    fflush(stdin);
    ier = fscanf(stdin,"%d",&i);
    info->imprim = i;
  }

  if ( mesh1->name == NULL ) {
    mesh1->name = (char *)calloc(128,sizeof(char));
    assert(mesh1->name);
    fprintf(stdout,"  -- MESH1 BASENAME ?\n");
    fflush(stdin);
    ier = fscanf(stdin,"%s",mesh1->name);
  }

  sol1->name = (char *)calloc(128,sizeof(char));
  assert(sol1->name);
  strcpy(sol1->name,mesh1->name);
  ptr = strstr(sol1->name,".mesh");
  if ( ptr ) *ptr = '\0';

  /* Option 3 prevails */
  if ( ( mesh2->name == NULL ) && (info->option != 3) )
	  info->option = 2;
  
  return(1);
}

/* Read file DEFAULT.mshdist */
static int parsop(Info *info,pMesh mesh) {
  float      fbuf;
  int        ret,i,k,l;
  char       *ptr,data[256];
  FILE       *in;

  strcpy(data,mesh->name);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  strcat(data,".mshdist");
  in = fopen(data,"r");

  if ( !in ) {
    sprintf(data,"%s","DEFAULT.mshdist");
    in = fopen(data,"r");
    if ( !in )  return(1);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /* Read parameters */
  while ( !feof(in) ) {
    ret = fscanf(in,"%s",data);
    if ( !ret || feof(in) )  break;
    for (i=0; i<strlen(data); i++) data[i] = tolower(data[i]);

    /* in mode -dom: read interior triangles; if none, default reference is REFINT */
    if ( !strcmp(data,"interiordomains") ) {
      fscanf(in,"%d",&info->nintel);
      info->intel = (int*)calloc(info->nintel,sizeof(int));
      assert ( info->nintel );
      for (k=0; k<info->nintel; k++)
        fscanf(in,"%d",&info->intel[k]);
    }

    /* in mode -dom: read starting triangles (useful in 3d only) */
    if ( !strcmp(data,"starttrias") ) {
      fscanf(in,"%d",&info->nst);
      info->st = (int*)calloc(info->nst,sizeof(int));
      assert ( info->nst );
      for (k=0; k<info->nst; k++)
        fscanf(in,"%d",&info->st[k]);
    }

    /* in mode -dom: read starting edges */
    if ( !strcmp(data,"startedges") ) {
      fscanf(in,"%d",&info->nsa);
      info->sa = (int*)calloc(info->nsa,sizeof(int));
      assert ( info->nsa );
      for (k=0; k<info->nsa; k++)
        fscanf(in,"%d",&info->sa[k]);
    }

    /* in mode -dom: read starting vertices */
    if ( !strcmp(data,"startver") ) {
      fscanf(in,"%d",&info->nsp);
      info->sp = (int*)calloc(info->nsp,sizeof(int));
      assert ( info->nsp );
      for (k=0; k<info->nsp; k++)
        fscanf(in,"%d",&info->sp[k]);
    }
    /* In generating distance, specify one or several exterior points for the sign */
    if ( !strcmp(data,"exteriorpoints") ) {
      fscanf(in,"%d",&info->nexp);
      info->exp = (double*)calloc(mesh->dim*info->nexp,sizeof(double));
      for (k=0; k<info->nexp; k++) {
        for (l=0; l<mesh->dim; l++) {
          fscanf(in,"%f",&fbuf);
          info->exp[mesh->dim*k+l] = fbuf;
        }
      }
    }

    /* Read references to initialize distance function (used with option startref) */
    /* Only available in 3D -> to update to match the 2d version */
    if ( !strcmp(data,"startref") ) {
      fscanf(in,"%d",&info->nsref);
      info->sref = (int*)calloc(info->nsref+1,sizeof(double));
      assert(info->sref);

      for(k=1; k<=info->nsref; k++) {
        fscanf(in,"%d ",&info->sref[k]);
      }
    }
  }
  fclose(in);

  return(1);
}

/* Write statistics about the supplied mesh(es) */
static void stats(pMesh mesh1,pMesh mesh2) {
  fprintf(stdout,"     NUMBER OF GIVEN VERTICES   %8d  %8d\n",mesh1->np,mesh2->np);
  if ( mesh1->dim == 2 )
    fprintf(stdout,"     NUMBER OF GIVEN ELEMENTS   %8d  %8d\n",mesh1->nt,mesh2->na);
  else
    fprintf(stdout,"     NUMBER OF GIVEN ELEMENTS   %8d  %8d\n",mesh1->ne,mesh2->nt);
}


/* Set function pointers */
int setfunc(Info info,int dim) {
  if ( dim == 2 ) {
    newBucket = newBucket_2d;

    if ( info.option == 1 ) {
      inidist = inidist_2d;
      inidistpcloud = inidistpcloud_2d;
    }
    else if ( info.option == 3 ){
      inireftrias = inireftrias_2d;
      iniencdomain = iniencdomain_2d;
    }
    else
      iniredist = iniredist_2d;

    ppgdist = ppgdist_2d;
    ppgdistfmm = ppgdistfmm_2d;
  }
  
  /* Three-dimensional case */
  else {
    newBucket = newBucket_3d;

    /* 3d surface mesh */
    if ( info.dsurf ) {
      iniencdomain = iniencdomain_s;
      iniredist = iniredist_s;
      ppgdist    = ppgdistfmm_s;
      ppgdistfmm = ppgdistfmm_s;
    }
    
    /* 3d volume mesh */
    else {
      if ( info.option == 1 ) {
        inidist = inidist_3d;
        inidistpcloud = inidistpcloud_3d;
      }
      else if ( info.option == 3 ) {
        inireftrias = inireftrias_3d;
        iniencdomain = iniencdomain_3d;
      }
      else
        iniredist = iniredist_3d;

      ppgdist = ppgdist_3d;
      ppgdistfmm = ppgdistfmm_3d;
    }
  }

  return(1);
}

/* Generate distance function according to the selected mode */
static int mshdis1(Info info,pMesh mesh1,pMesh mesh2,pSol sol1) {
  pBucket  bucket;
  int      ier;

  if(info.specdist){
    printf("GENERATING HOLES \n");
    //genHolesPCB_2d(mesh1,sol1);
    //genHolesRadia_2d(mesh1,sol1);
    //anafuncbuddha(mesh1,sol1);
    //anafuncsq(mesh1,sol1);
    gen2Holes_2d(mesh1,sol1);
    //holeCl_3d(mesh1,sol1);
    //anafunchelix(mesh1,sol1);
    //holeClCrane(mesh1,sol1);
    //holeClGrip(mesh1,sol1);
    //holeGripIni(mesh1,sol1);
    //gen1Hole_2d(mesh1,sol1);
    //anafuncspherelag(mesh1,sol1);
    //holeCraneIni(mesh1,sol1);
    //genHolesMast_3d(mesh1,sol1);
    //holeClMast_3d(mesh1,sol1);
    //holeClStarfish_3d(mesh1,sol1);
    //genHolesStarfish_3d(mesh1,sol1);
    //genHolesBridge_3d(mesh1,sol1);
    //holeClBridge_3d(mesh1,sol1);
    //holeClBridge2_3d(mesh1,sol1);
    //genHolesCanti_2d(mesh1,sol1);
    //genHolesMast_2d(mesh1,sol1);
    //genHolesCanti_3d(mesh1,sol1);
    //holeLBeamIni(mesh1,sol1);
    //holeLBBBeamIni(mesh1,sol1);
    //holeCl_Chair(mesh1,sol1);
    //holeChairIni(mesh1,sol1);
    //genHolesCantia_3d(mesh1,sol1);
    //genholes_3d(mesh1,sol1,3,3,3);
    //genHolesSCantia_3d(mesh1,sol1);
    return(1);
  }

  /* Distance initialization, depending on the option */
  /* Signed distance generation from the contour of mesh2 */
  if ( info.option == 1 ) {
    /* bucket sort */
    bucket = newBucket(mesh1,BUCKSIZ);
    if ( !bucket )  return(0);
    
    if ( info.imprim )  fprintf(stdout,"  ** Initialization\n");
    if ( info.pcloud )
      ier = inidistpcloud(mesh1,mesh2,sol1,bucket);
    else {
      ier = inidist(info,mesh1,mesh2,sol1,bucket);
    }
    if ( info.imprim )  fprintf(stdout,"  ** Sign identification\n");
    
    /* Put a sign to the initial distance field */
    if ( !info.pcloud ) {
      ier = mesh1->dim == 2 ? sgndist_2d(info,mesh1,mesh2,sol1,bucket) : sgndist_3d(info,mesh1,mesh2,sol1,bucket);
      if ( !ier )  return(0);
    }
  }

  /* Signed distance generation from entities enclosed in mesh1 */
  else if ( info.option == 3 ){
    if( info.startref ) {
      if ( info.imprim )  fprintf(stdout,"  ** Generation of signed distance function from surface edg/tria with prescribed refs\n");
      ier = inireftrias(info,mesh1,sol1);
    }
  	else {
      if ( info.imprim )  fprintf(stdout,"  ** Generation of signed distance function from tria/tets with ref 3\n");
      ier = iniencdomain(info,mesh1,sol1);
    }
    assert(ier);
  }

  /* Redistancing */
  else {
    if ( info.imprim )  fprintf(stdout,"  ** Redistancing\n");
    ier = iniredist(info,mesh1,sol1);
  }
    
  /* Signed distance propagation */
  if ( ier > 0 ) {
    chrono(ON,&info.ctim[4]);
    if ( info.imprim )  fprintf(stdout,"  ** Propagation [%d cpu]\n",info.ncpu);
    
    if ( info.fmm )
      ier = ppgdistfmm(mesh1,sol1);
    else
      ier = ppgdist(info,mesh1,sol1);
    chrono(OFF,&info.ctim[4]);
  }
  else if ( ier < 0 )
    fprintf(stdout,"  ## Problem in sign function\n");

  if ( info.option == 1 )  freeBucket(bucket);

  return(ier);
}


/* Main program */
int main(int argc,char **argv) {
  Info     info;
  Mesh     mesh1,mesh2;
  Sol      sol1;
  int      k,ier,*perm;
  char     stim[32];

  fprintf(stdout,"  -- MSHDIST, Release %s (%s) \n",D_VER,D_REL);
  fprintf(stdout,"     %s\n",D_CPY);
  fprintf(stdout,"    %s\n",COMPIL);

  /* Trap exceptions */
  signal(SIGABRT,excfun);
  signal(SIGFPE,excfun);
  signal(SIGILL,excfun);
  signal(SIGSEGV,excfun);
  signal(SIGTERM,excfun);
  signal(SIGINT,excfun);

  tminit(info.ctim,TIMEMAX);
  chrono(ON,&info.ctim[0]);

  /* Default values */
  memset(&mesh1,0,sizeof(Mesh));
  memset(&mesh2,0,sizeof(Mesh));
  memset(&sol1,0,sizeof(Sol));
  memset(&info,0,sizeof(Info));

  info.imprim = -99;
  info.option = 1;
  info.ddebug = 0;
  info.ncpu   = 1;
  info.res    = EPS;
  info.nexp   = 0;
  info.dt     = 0.001;
  info.maxit  = 1000;
  info.size   = SIZE;
  info.nintel = -1;
  
  /* Parse command line arguments */
  if ( !parsar(argc,argv,&info,&mesh1,&sol1,&mesh2) )  return(1);
  
  /* Load data */
  if ( info.imprim )   fprintf(stdout,"\n  -- INPUT DATA\n");
  chrono(ON,&info.ctim[1]);
  
  if ( !loadMesh(&info,&mesh1,&mesh2) )  return(1);
    
  /* Load solution or allocate memory */
  if ( info.option == 2 ) {
    if ( !loadSol(&sol1) )  return(1);
  }
  else if ( ( info.option == 1 ) || ( info.option == 3 ) ) {
    sol1.bin = mesh1.bin;
    sol1.ver = GmfDouble;
    sol1.size    = 1;
    sol1.type[0] = 1;
    sol1.typtab[0][0] = GmfSca;
    sol1.dim = mesh1.dim;
    sol1.np  = mesh1.np;
    sol1.val = (double*)malloc((sol1.np+1)*sizeof(double));
    assert(sol1.val);
  }

  /* Set pointers */
  if ( !setfunc(info,mesh1.dim) )  return(1);

  chrono(OFF,&info.ctim[1]);
  if ( info.imprim )
    stats(&mesh1,&mesh2);
  
  /* Read file DEFAULT.mshdist, if any */
  parsop(&info,&mesh1);
  
  /* Default value for the interior domain, if none supplied (used in -dom option only) */
  if ( info.nintel < 0 ) {
    info.nintel = 1;
    info.intel = (int*)calloc(1,sizeof(int));
    info.intel[0] = REFINT;
  }
  
  /* Default value for the starting point if none supplied (used in generating signed distance only) */
  if ( !info.nexp ) {
    info.nexp = -1;
    info.exp  = (double*)calloc(mesh1.dim,sizeof(double));
    for (k=0; k<mesh1.dim; k++)
      info.exp[k] = 0.01;
  }

  /* Te be deleted, ultimately */
  if ( info.startref && info.nsref < 1 ) {
    printf("    *** No starting ref for triangles found.\n ");
    exit(0);
  }
  
  /* Compress mesh in the surface case if tetrahedra are supplied */
  if ( info.zip ) {
    perm = (int*)calloc(mesh1.np+1,sizeof(int));
    if ( !pack_s(&mesh1,&sol1,perm) ) {
      printf("    *** Impossible to pack mesh; abort.\n ");
      exit(0);
    }
  }

  printim(info.ctim[1].gdif,stim);
  fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);

  /* Analysis */
  fprintf(stdout,"\n  %s\n   MODULE MSHDIST-LJLL : %s (%s)\n  %s\n",D_STR,D_VER,D_REL,D_STR);
  chrono(ON,&info.ctim[2]);
  if ( info.imprim )   fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  
  if ( !info.noscale || !info.specdist )
    if ( !scaleMesh(&info,&mesh1,&mesh2,&sol1) )  return(1);

  info.nexp = info.nexp == -1 ? 1 : info.nexp;
  
  /* Create adjacencies */
  ier = ( mesh1.dim == 2 || ( mesh1.dim ==3 && info.dsurf ) )? hashelt_2d(&mesh1) : hashelt_3d(&mesh1);
  if ( !ier )  return(1);

  chrono(OFF,&info.ctim[2]);
  if ( info.imprim ) {
    printim(info.ctim[2].gdif,stim);
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);
  }
  
  /* Distance calculation */
  if ( info.imprim )   fprintf(stdout,"\n  -- PHASE 2 : DISTANCING\n");
  chrono(ON,&info.ctim[3]);

  if ( !mshdis1(info,&mesh1,&mesh2,&sol1) )  return(1);
  /* free memory */
  if ( info.nintel )  free(info.intel);
  if ( info.nst )  free(info.st);
  if ( info.nsa )  free(info.sa);
  if ( info.nsp )  free(info.sp);
  if ( info.nexp )  free(info.exp);

  chrono(OFF,&info.ctim[3]);
  if ( info.imprim ) {
    printim(info.ctim[3].gdif,stim);
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
  }
  
  fprintf(stdout,"\n  %s\n   END OF MODULE MSHDIST \n  %s\n",D_STR,D_STR);

  /* Save file */
  if ( info.imprim )  fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",sol1.name);
  chrono(ON,&info.ctim[1]);
  if ( !info.noscale || !info.specdist )
    if ( !unscaleSol(info,&sol1) )  return(1);
  
  /* Unpacking */
  if ( info.zip ) {
    if ( !unpack_s(&mesh1,&sol1,perm) ) {
      printf("    *** Impossible to pack mesh; abort.\n ");
      exit(0);
    }
    free(perm);
  }

  if ( !saveSol(&sol1) )     return(1);
  chrono(OFF,&info.ctim[1]);
  if ( info.imprim )  fprintf(stdout,"  -- WRITING COMPLETED\n");

  free(mesh1.point);
  free(mesh1.adja);
  if ( mesh2.point )  free(mesh2.point);
  if ( mesh1.dim == 3 ) {
    free(mesh1.tetra);
    if ( mesh2.tria )  free(mesh2.tria);
  }
  else {
    free(mesh1.tria);
    if ( mesh2.edge )  free(mesh2.edge);
  }
  free(sol1.val);

  return(0);
}
