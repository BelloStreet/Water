/*#####################################*/
/*#Este es para la version petsc-3.7.5#*/
/*#Para la version petsc-3.2 xchem    #*/
/*#####################################*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscsys.h>
#include <petscviewer.h>
#include <petscts.h>
int main(int argc,char **argv){
  Vec y;
  PetscScalar *array;
  PetscViewer viewen;
  FILE *file;
  int i,nvec,size,nlocal,rank,root=0;
  _Complex double *wloc,*wtot;
  char solution[50]="bound.m";
  PetscInitialize(&argc,&argv,0,0);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  if(rank==root){
    printf("Numbers of Processors=%d\n",size);
    fflush(stdout);
  }
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,solution,FILE_MODE_READ,&viewen);
  VecCreate(PETSC_COMM_WORLD,&y);
  VecSetType(y,VECMPI);
  VecLoad(y,viewen);
  PetscViewerDestroy(&viewen);
  VecGetSize(y,&nvec);
  VecGetLocalSize(y,&nlocal);
  wloc=(_Complex double*)calloc(nlocal,sizeof(_Complex double));
  wtot=(_Complex double*)calloc(nvec,sizeof(_Complex double));
  VecGetArray(y,&array);
  for(i=0;i<nlocal;i++){
    wloc[i]=array[i];
  }
  VecRestoreArray(y,&array);
  MPI_Gather(wloc,nlocal,MPI_DOUBLE_COMPLEX,wtot,nlocal,MPI_DOUBLE_COMPLEX,root,PETSC_COMM_WORLD);
  if(rank==root){
    file=fopen("bound_ascii.dat","w");
    fprintf(file,"%d\n",nvec);
    fprintf(file,"%d\n",2);
    for(i=0;i<nvec;i++){
      fprintf(file,"%d % .16le\n",i+1,creal(wtot[i]));
    }
    fclose(file);
  }
  free(wtot);
  free(wloc);
  VecDestroy(&y);
  PetscFinalize();
  return 0;
}

