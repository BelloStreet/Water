#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_eigen.h>
#include "utils.h"
#include "sasha.h"
#include "class_grid.h"
#include "orb.h"
int main(int argc,char **argv){
  FILE *file;
  int i,j,n,l,i1,i2,j1,j2,tag=99;
  int id=0,numprocs,root=0,nbas,lamax,lmax,io,*lo=NULL,*mo=NULL,lmorb;
  int *la,*ma;
  double z1,z2,z3,z4,qa=0.0,qb=60.0,**co,*C3j,*C3i;
  cmplx *Xi,*Wi,*TXX,*TIXX,*H1,*cilm=NULL,*K,*J;
  char command[150],fichero[50]="he_no.molden";//"lih_states3.02.molden2015";
  //  char symm[5]="C2v",irrep[5]="a1",str[2]="c";
  char symm[5]="D2h",irrep[5]="ag",str[2]="c";
  orbitals od[10],ord;  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&id);
  printf("id=%d numprocs=%d\n",id,numprocs);

  if(id==root){
    strcpy(command,"./Hatom");
    system(command);
  }
  file=fopen("quadrature.dat","r");
  fscanf(file,"%d",&lamax);
  fscanf(file,"%d",&nbas);
  Xi=(cmplx*)calloc(nbas,sizeof(cmplx));
  Wi=(cmplx*)calloc(nbas,sizeof(cmplx));
  for(i=0;i<nbas;i++){
    fscanf(file,"%le%le%le%le",&z1,&z2,&z3,&z4);
    Xi[i]=z1+z2*I;
    Wi[i]=z3+z4*I;
  }
  fclose(file);

  lmax=lamax/2;
  TXX=(cmplx*)calloc(nbas*nbas*lamax,sizeof(cmplx));
  TIXX=(cmplx*)calloc(nbas*nbas*lamax,sizeof(cmplx));
  file=fopen("kinetic.dat","r");
  for(i=0;i<nbas*nbas*lamax;i++){
    fscanf(file,"%le%le",&z1,&z2);
    TXX[i]=z1+z2*I;
  }
  for(i=0;i<nbas*nbas*lamax;i++){
    fscanf(file,"%le%le",&z1,&z2);
    TIXX[i]=z1+z2*I;
  }
  fclose(file);  

  GRID CORR(nbas,lamax,Xi,Wi,TXX,TIXX);

  H1=(cmplx*)calloc(nbas*nbas,sizeof(cmplx));
  //  CORR._HH1(id,1.0,H1);

  //read orbitales from molpro, put them on the grid  
  //  for(io=0;io<1;io++){
  io=0;
  if(id==root){
    CORR._GETORB(fichero,symm,irrep,io,str,&lmorb,&lo,&mo,&cilm);
  }
  MPI_Bcast(&lmorb,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&str,strlen(str)+1,MPI_CHAR,root,MPI_COMM_WORLD);
  MPI_Bcast(&symm,strlen(symm)+1,MPI_CHAR,root,MPI_COMM_WORLD);
  MPI_Bcast(&irrep,strlen(irrep)+1,MPI_CHAR,root,MPI_COMM_WORLD);
  if(id==root){
    printf("Done Bcasting\n");
    fflush(stdout);
  }
  if(id!=root){
    lo=(int*)calloc(lmorb,sizeof(int));
    mo=(int*)calloc(lmorb,sizeof(int));
    cilm=(cmplx*)calloc(nbas*lmorb,sizeof(cmplx));
  }
  MPI_Bcast(lo,lmorb,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(mo,lmorb,MPI_INT,root,MPI_COMM_WORLD);
  if(id==root){
    printf("Done Bcasting\n");
    fflush(stdout);
  }
  MPI_Bcast(cilm,nbas*lmorb,MPI_DOUBLE_COMPLEX,root,MPI_COMM_WORLD);
  if(id==root){
    printf("Done Bcasting\n");
    fflush(stdout);
  }
  n=(lamax-1)*(lamax+1)+1;
  C3j=(double*)calloc(lmorb*lmorb*n,sizeof(double));
  C3i=(double*)calloc(lmorb*lmorb*n,sizeof(double));
  _GETC3J(lamax,str,str,lmorb,lmorb,lo,mo,lo,mo,C3j,C3i);
  //     _INITORB(symm,irrep,str,lmorb,lo,mo,cilm,&ord);

  K=(cmplx*)calloc(nbas*nbas*lmorb*lmorb,sizeof(cmplx));
  J=(cmplx*)calloc(nbas*nbas*lmorb*lmorb,sizeof(cmplx));
  CORR._KO(id,numprocs,str,lmorb,lo,mo,cilm,C3j,C3i,K);
  CORR._JO(id,numprocs,str,lmorb,lo,mo,cilm,C3j,C3i,J);
  CORR._JMK(id,numprocs,lmorb,cilm,K,J);
  
  CORR._HH1(2.0,id,numprocs,lmorb,lo,mo,K,J,H1);
  
  free(lo);
  free(mo);
  free(cilm);
  free(C3i);
  free(C3j);
  //  }  
  
  //two electron hamiltonian

  free(Xi);
  free(Wi);
  free(TXX);
  free(TIXX);
  free(H1);

  MPI_Finalize();
  return 0;
}
