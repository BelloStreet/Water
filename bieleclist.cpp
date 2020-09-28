#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void _bieleclist_0(int n,int *ijkl){ //Wilson 257 pseudo-code
  int ii,jj,i,j,k,l,lmax,kount,m;
  kount=0;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      for(k=0;k<n;k++){
  	for(l=0;l<n;l++){
  	  ijkl[0+4*kount]=i;
  	  ijkl[1+4*kount]=j;
  	  ijkl[2+4*kount]=k;
  	  ijkl[3+4*kount]=l;
  	  kount=kount+1;
  	}
      }
    }
  }
}


void _bieleclist_2(int n,int *ijkl){ //Wilson 257 pseudo-code
  int ii,jj,i,j,k,l,lmax,kount,m;
  kount=0;
  for(ii=0;ii<n*n;ii++){
    i=floor(ii/n);
    j=ii-n*i;
    for(jj=ii;jj<n*n;jj++){
      k=floor(jj/n);
      l=jj-n*k;
      ijkl[0+4*kount]=i;
      ijkl[1+4*kount]=j;
      ijkl[2+4*kount]=k;
      ijkl[3+4*kount]=l;
      kount=kount+1;
    }
  }
}


void _bieleclist_4(int n,int *ijkl){ //Wilson 257 pseudo-code
  int ii,jj,i,j,k,l,lmax,kount,m;
  kount=0;
  for(i=0;i<n;i++){
    for(j=i;j<n;j++){
      for(k=0;k<n;k++){
  	for(l=k;l<n;l++){
  	  ijkl[0+4*kount]=i;
  	  ijkl[1+4*kount]=j;
  	  ijkl[2+4*kount]=k;
  	  ijkl[3+4*kount]=l;
  	  kount=kount+1;
  	}
      }
    }
  }  
}


void _bieleclist_8(int n,int *ijkl){ //Wilson 257 pseudo-code
  int i,j,k,l,lmax,kount,m;
  kount=0;
  for(i=0;i<n;i++){
    for(j=0;j<=i;j++){
      for(k=0;k<=i;k++){
  	if(i==k){
	  lmax=j;
  	}
  	else{
  	  lmax=k;
  	}
  	for(l=0;l<=lmax;l++){
	  ijkl[0+4*kount]=i;
	  ijkl[1+4*kount]=j;
	  ijkl[2+4*kount]=k;
	  ijkl[3+4*kount]=l;
  	  kount=kount+1;
  	}
      }
    }    
  }
}


void _mblock(int nf,int nb,int *nmkl){
  int i,j,k=0,*b,*n1,*n2;
  b=(int*)calloc(nb,sizeof(int));
  n1=(int*)calloc(nb,sizeof(int));
  n2=(int*)calloc(nb,sizeof(int));
  for(i=0;i<nb;i++){
    b[i]=nf*nf;
  }
  for(i=0;i<nb;i++){
    k=k+b[i];
    n2[i]=k;
    n1[i]=n2[i]-b[i]+1;
  }
  k=0;
  for(i=0;i<nb;i++){
    for(j=i;j<nb;j++){
      nmkl[0+4*k]=n1[i];
      nmkl[1+4*k]=n2[i];
      nmkl[2+4*k]=n1[j];
      nmkl[3+4*k]=n2[j];
      k=k+1;
    }
  }
  free(b);
  free(n1);
  free(n2);
}

int _dpr(char *symm,int il, int jk){
  int i,j,k;
  if(strcmp(symm,"C2v")==0){
    //a1=0,b1=1,b2=2,a2=3
    int dpr[4][4]={0,1,2,3,
		   1,0,3,2,
		   2,3,0,1,
		   3,2,1,0};
    k=dpr[il][jk];
  }
  if(strcmp(symm,"D2h")==0){
    //ag=0,b3u=1,b2u=2,b1g=3,b1u=4,b2g=5,b3g=6,au=7
    int dpr[8][8]={0,1,2,3,4,5,6,7,
		   1,0,3,2,5,4,7,6,
		   2,3,0,1,6,7,4,5,
		   3,2,1,0,7,6,5,4,
		   4,5,6,7,0,1,2,3,
		   5,4,7,6,1,0,3,2,
		   6,7,4,5,2,3,0,1,
		   7,6,5,4,3,2,1,0};
    k=dpr[il][jk];
  }
  return k;
}



  
