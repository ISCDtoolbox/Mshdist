#include "mshdist.h"

/* Read mesh data */
int loadMesh(Info *info,pMesh mesh1,pMesh mesh2) {
  pPoint       ppt;
  pTetra       pt;
  pTria        pt1;
  pEdge        pe;
  float        fp1,fp2,fp3;
  int          inm,i,k,ref,nh,nq;
  char         *ptr,*name,data[128];

  name = mesh1->name;
  strcpy(data,name);
  
  /* Open mesh file (or .meshb) */
  ptr = strstr(data,".mesh");
  mesh1->bin = 0;
  if ( !ptr ) {
    strcat(data,".meshb");
    if (!(inm = GmfOpenMesh(data,GmfRead,&mesh1->ver,&mesh1->dim)) ) {
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if (!(inm = GmfOpenMesh(data,GmfRead,&mesh1->ver,&mesh1->dim)) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        return(0);
      }
      mesh1->bin = 0;
    }
  }
  else if (!(inm = GmfOpenMesh(data,GmfRead,&mesh1->ver,&mesh1->dim)) ) {
    fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  if ( abs(info->imprim) > 3 )
    fprintf(stdout,"  -- READING DATA FILE %s\n",data);

  /* Read numbers of tetrahedra and vertices in the case of a full (volumic) 3d mesh */
  if ( mesh1->dim == 3 && !info->dsurf ) {
    mesh1->np = GmfStatKwd(inm,GmfVertices);
    mesh1->ne = GmfStatKwd(inm,GmfTetrahedra);
    nh = GmfStatKwd(inm,GmfHexahedra);
    if ( nh ) {
      fprintf(stdout,"  ## Non simplicial mesh\n");
      return(0);
    }
    if( info->option == 3 && info->startref ) mesh1->nt = GmfStatKwd(inm,GmfTriangles);
  }
  
  /* Read numbers of Triangles and vertices in 2d mesh or 3d surface triangulation */
  else {
    mesh1->np = GmfStatKwd(inm,GmfVertices);
    /* Need to compress the surface mesh if tetras are also supplied */
    if ( mesh1->dim == 3 && info->dsurf )
      if ( GmfStatKwd(inm,GmfTetrahedra) ) info->zip = 1;
        
    if ( (info->option == 3) || (info->bbbc) || (info->hausdorff))
      mesh1->na = GmfStatKwd(inm,GmfEdges);
    
    mesh1->nt = GmfStatKwd(inm,GmfTriangles);
    nq = GmfStatKwd(inm,GmfQuadrilaterals);
    if ( nq ) {
      fprintf(stdout,"  ## Non simplicial mesh\n");
      return(0);
    }
  }
  
  /* Invalid mesh */
  if ( !mesh1->np || (mesh1->ne+mesh1->nt == 0)  ) {
    fprintf(stdout,"  ** MISSING DATA (expecting elements)\n");
    //return(0);
  }

  /* Memory allocation for vertices, triangles and tetrahedra */
  mesh1->point = (pPoint)calloc(mesh1->np+1,sizeof(Point));
  assert(mesh1->point);
  
  if ( mesh1->dim == 3 && !info->dsurf ) {
    mesh1->tetra = (pTetra)calloc(mesh1->ne+1,sizeof(Tetra));
    assert(mesh1->tetra);
    
    if( info->option == 3 && info->startref ) {
      mesh1->tria  = (pTria)calloc(mesh1->nt+1,sizeof(Tria));
      assert(mesh1->tria);
    }
  }
  else {
    mesh1->tria  = (pTria)calloc(mesh1->nt+1,sizeof(Tria));
    assert(mesh1->tria);
    
    if ((( mesh1->na )&&(info->option == 3))||(( mesh1->na )&&(info->bbbc))||((info->hausdorff)&&(mesh1->na))) {
	    mesh1->edge  = (pEdge)calloc(mesh1->na+1,sizeof(Edge));
	    assert(mesh1->edge);
  	}
  }
  
  /* Memory allocation for adjacencies */
  if ( mesh1->dim == 2 || ( mesh1->dim ==3 && info->dsurf) ) {
    mesh1->adja = (int*)calloc(3*mesh1->nt+5,sizeof(int));
    assert(mesh1->adja);
  }
  else {
    mesh1->adja = (int*)calloc(4*mesh1->ne+5,sizeof(int));
    assert(mesh1->adja);
  }

  /* Read mesh vertices */
  GmfGotoKwd(inm,GmfVertices);
  if ( mesh1->dim == 3 ) {
    for (k=1; k<=mesh1->np; k++) {
      ppt = &mesh1->point[k];
      ppt->flag = 0;
      if ( mesh1->ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
        ppt->c[2] = fp3;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->c[2],&ref);
    }
  }
  else {
    for (k=1; k<=mesh1->np; k++) {
      ppt = &mesh1->point[k];
      ppt->flag = 0;
      if ( mesh1->ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ref);
    }
  }
  
  /* Read mesh edges */
  if ( ((mesh1->dim == 2)&&(info->option == 3)&&(mesh1->na)) || (((mesh1->dim == 2)&&(info->bbbc)&&(mesh1->na))) || ((info->hausdorff)&&(mesh1->na))){
	GmfGotoKwd(inm,GmfEdges);
    for (k=1; k<=mesh1->na; k++) {
      pe = &mesh1->edge[k];
      GmfGetLin(inm,GmfEdges,&pe->v[0],&pe->v[1],&pe->ref);
    }
  }
  
  /* Read mesh elements in the case of a 3d mesh */
  if ( mesh1->dim == 3 && !info->dsurf ) {
    GmfGotoKwd(inm,GmfTetrahedra);
    for (k=1; k<=mesh1->ne; k++) {
      pt = &mesh1->tetra[k];
      if (info->option ==3)
        GmfGetLin(inm,GmfTetrahedra,&pt->v[0],&pt->v[1],&pt->v[2],&pt->v[3],&pt->ref);
      else
        GmfGetLin(inm,GmfTetrahedra,&pt->v[0],&pt->v[1],&pt->v[2],&pt->v[3],&ref);
      
      /* Attribute the number of one tetra containing p0 */
      for (i=0; i<4; i++) {
        ppt = &mesh1->point[pt->v[i]];
        if ( !ppt->s )  ppt->s = k;
      }
    }
    
    /* In case of option sref, read also mesh triangles */
    if( info->option == 3 && info->startref ) {
      GmfGotoKwd(inm,GmfTriangles);
      for (k=1; k<=mesh1->nt; k++) {
        pt1 = &mesh1->tria[k];
        if (info->option ==3)
          GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
        else
          GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&ref);
        
        /* Attribute the number of one tria containing p0 */
        for (i=0; i<3; i++) {
          ppt = &mesh1->point[pt1->v[i]];
          if ( !ppt->s )  ppt->s = k;
        }
      }
    }    
  }
  /* Read mesh elements in the case of a 2d mesh or a 3d surface triangulation */
  else {
    GmfGotoKwd(inm,GmfTriangles);
    for (k=1; k<=mesh1->nt; k++) {
      pt1 = &mesh1->tria[k];
      if ( info->option == 3 )
        GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
      else
        GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&ref);
      
	  for (i=0; i<3; i++) {    
        ppt = &mesh1->point[pt1->v[i]];
        if ( !ppt->s )  ppt->s = k;
      }
    }
  }
  GmfCloseMesh(inm);

  if ( ( info->option == 2 ) || ( info->option == 3 ) )  return(1);
	
  /* Load mesh 2 */
  name = mesh2->name;
  strcpy(data,name);
  ptr = strstr(data,".mesh");
  mesh2->bin = 0;
  if ( !ptr ) {
    strcat(data,".meshb");
    if (!(inm = GmfOpenMesh(data,GmfRead,&mesh2->ver,&mesh2->dim)) ) {
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if (!(inm = GmfOpenMesh(data,GmfRead,&mesh2->ver,&mesh2->dim)) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        return(0);
      }
      mesh2->bin = 0;
    }
  }
  else if (!(inm = GmfOpenMesh(data,GmfRead,&mesh2->ver,&mesh2->dim)) ) {
    fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  if ( abs(info->imprim) > 3 )
    fprintf(stdout,"  -- READING DATA FILE %s\n",data);

  if ( mesh2->dim != mesh1->dim ) {
		fprintf(stdout,"  %%%% Wrong dimension %d %d\n",mesh1->dim,mesh2->dim);
		return(0);
	}
	
  if ( mesh2->dim == 3 ) {
    mesh2->np = GmfStatKwd(inm,GmfVertices);
    mesh2->nt = GmfStatKwd(inm,GmfTriangles);
  }
  else {
    mesh2->np = GmfStatKwd(inm,GmfVertices);
    mesh2->na = GmfStatKwd(inm,GmfEdges);
  }
  
  if ( (!info->pcloud && ( !mesh2->np || (mesh2->nt+mesh2->na == 0) )) || ( info->pcloud && !mesh2->np ) ) {
    fprintf(stdout,"  ** MISSING DATA (expecting elements)\n");
    return(0);
  }

  /* Memory allocation */
  mesh2->point = (pPoint)calloc(mesh2->np+1,sizeof(Point));
  assert(mesh2->point);
  
  /* Read mesh vertices */
  GmfGotoKwd(inm,GmfVertices);
  if ( mesh2->dim == 3 ) {
    for (k=1; k<=mesh2->np; k++) {
      ppt = &mesh2->point[k];
      if ( mesh2->ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
        ppt->c[2] = fp3;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->c[2],&ref);
    }
  }
  else {
    for (k=1; k<=mesh2->np; k++) {
      ppt = &mesh2->point[k];
      if ( mesh2->ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
      }
      else {
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ref);
      }
    }
  }

  /* Allocate memory and read mesh elements */
  if ( !info->pcloud ) {
    if ( mesh2->dim == 3 ) {
      mesh2->tria = (pTria)calloc(mesh2->nt+1,sizeof(Tria));
      assert(mesh2->tria);
    }
    else {
      mesh2->edge  = (pEdge)calloc(mesh2->na+1,sizeof(Edge));
      assert(mesh2->edge);
    }
    
    if ( mesh2->dim == 3 ) {
      GmfGotoKwd(inm,GmfTriangles);
      for (k=1; k<=mesh2->nt; k++) {
        pt1 = &mesh2->tria[k];
        GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
      }
    }
    else {
      GmfGotoKwd(inm,GmfEdges);
      for (k=1; k<=mesh2->na; k++) {
        pe = &mesh2->edge[k];
        GmfGetLin(inm,GmfEdges,&pe->v[0],&pe->v[1],&pe->ref);
      }
    }
  }

  GmfCloseMesh(inm);
  return(1);
}


/* save (part of) mesh to disk */
int saveMesh(pMesh mesh,char *fileout) {
	pPoint     ppt;
	pTria      pt;
	pTetra     ptt;
	float      fp1,fp2,fp3;
	int        outm,k,ver,ref;
	
	ver = GmfFloat;
	if ( !(outm = GmfOpenMesh(fileout,GmfWrite,ver,mesh->dim)) ) {
		fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",fileout);
		return(0);
	}
	fprintf(stdout,"  %%%% %s OPENED\n",fileout);
	
	GmfSetKwd(outm,GmfVertices,mesh->np);
	for (k=1; k<=mesh->np; k++) {
		ppt = &mesh->point[k];
		if ( mesh->dim == 2 ) {
			fp1 = ppt->c[0];
			fp2 = ppt->c[1];
			ref = 0;
			GmfSetLin(outm,GmfVertices,fp1,fp2,ref);
		}
		else {
			fp1 = ppt->c[0];
			fp2 = ppt->c[1];
			fp3 = ppt->c[2];
			ref = 0;
			GmfSetLin(outm,GmfVertices,fp1,fp2,fp3,ref);
		}
	}
	
	/* write triangles */
	GmfSetKwd(outm,GmfTriangles,mesh->nt);
	for (k=1; k<=mesh->nt; k++) {
		pt = &mesh->tria[k];
		if ( !pt->v[0] )  continue;
		ref = 0;
		GmfSetLin(outm,GmfTriangles,pt->v[0],pt->v[1],pt->v[2],ref);
	}
	
	/* write tetrahedra */
	GmfSetKwd(outm,GmfTetrahedra,mesh->ne);
	for (k=1; k<=mesh->ne; k++) {
		ptt = &mesh->tetra[k];
		ref = 0;
		GmfSetLin(outm,GmfTetrahedra,ptt->v[0],ptt->v[1],ptt->v[2],ptt->v[3],ref);
	}
	
		fprintf(stdout,"     TOTAL NUMBER OF VERTICES   %8d\n",mesh->np);
		fprintf(stdout,"     TOTAL NUMBER OF TRIANGLES  %8d\n",mesh->nt);
		fprintf(stdout,"     TOTAL NUMBER OF TETRA      %8d\n",mesh->ne);
	
	GmfCloseMesh(outm);
	return(1);
}


/* Load function defined at the verticles of the mesh */
int loadSol(pSol sol) {
  double       dbuf[ GmfMaxTyp ];
  float        fbuf[ GmfMaxTyp ];
  int          k,inm;
  char        *ptr,data[128];
	
  strcpy(data,sol->name);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  strcat(data,".solb");
  if (!(inm = GmfOpenMesh(data,GmfRead,&sol->ver,&sol->dim)) ) {
	ptr  = strstr(data,".solb");
	*ptr = '\0';
	strcat(data,".sol");
	if (!(inm = GmfOpenMesh(data,GmfRead,&sol->ver,&sol->dim)) ) {
	  fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
	  return(0);
	}
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);
  fprintf(stdout,"  -- READING DATA FILE %s\n",data);
	
  sol->np = GmfStatKwd(inm,GmfSolAtVertices,sol->type,&sol->size,sol->typtab);
  if ( !sol->np ) {
	fprintf(stdout,"  ** MISSING DATA.\n");
	return(0);
  }
  if ( sol->size > 1 ) {
	fprintf(stdout,"  ** DATA NOT COMPATIBLE\n");
	return(0);
  }

  /* mem alloc */
  sol->val  = (double*)calloc(sol->np+1,sizeof(double));
  assert(sol->val);

  /* read sol */
  GmfGotoKwd(inm,GmfSolAtVertices);
  for (k=1; k<=sol->np; k++) {
	if ( sol->ver == GmfFloat ) {
	  GmfGetLin(inm,GmfSolAtVertices,fbuf);
	  sol->val[k] = fbuf[0];
	}
	else {
	  GmfGetLin(inm,GmfSolAtVertices,dbuf);
	  sol->val[k] = dbuf[0];
	}
  }
	
  GmfCloseMesh(inm);
  return(1);
}

/* Save solution, defined at the vertices of the mesh */
int saveSol(pSol sol) {
  float        fbuf;
  int          k,inm;
  char        *ptr,data[128];
  double       dbuf;	

  strcpy(data,sol->name);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
	strcat(data,sol->bin ? ".solb" : ".sol");
  if (!(inm = GmfOpenMesh(data,GmfWrite,sol->ver,sol->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s\n",data);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /* write sol */
  GmfSetKwd(inm,GmfSolAtVertices,sol->np,sol->type[0],sol->typtab[0]);
  if ( sol->ver == GmfFloat ) {
	  for (k=1; k<=sol->np; k++) {
	    fbuf = sol->val[k];
	    GmfSetLin(inm,GmfSolAtVertices,&fbuf);
	  }
  }
  else {
	  for (k=1; k<=sol->np; k++) {
      GmfSetLin(inm,GmfSolAtVertices,&sol->val[k]);
    }
  }

  GmfCloseMesh(inm);
  return(1);
}

