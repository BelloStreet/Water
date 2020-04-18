program Hatom_fem
  use numbers
  use dvrecs
  implicit none
#include "petsc/finclude/petsc.h"
#include "slepc/finclude/slepc.h"
#include "petsc/finclude/petscvec.h90"
  Mat AA,VVA
  Vec x,y,z,vxx,phisc
  PetscInt d_nz,o_nz
  EPS eps
  PetscScalar vAA,xx(1),Er,Ei
  PetscScalar,pointer::xx_v(:),ll_v(:)
  PetscReal nrm,errest,relerr,tol
  PetscViewer viewen
  integer::i,j,k,n=600,lmax=12,v,l,lda,ldb,info,li,lf,ni,nf,nloc,id,numprocs
  integer::jj,ii,ierr,root=0,its,nconv,lwork,ncount
  integer,allocatable::ipiv(:)
  real(kind=8)::Cval,Ebs,e
  real(kind=8),allocatable::rwork(:),yy(:),wav(:)
  complex*16,allocatable::HH(:,:,:),VV(:,:),zImat(:,:),vr(:,:),vl(:,:)
  complex*16,allocatable::ee(:,:),c_nl(:,:),c_ni(:,:),w(:),work(:),A(:,:)
  character:: jobz='V',uplo='L'
  call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,numprocs,ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,id,ierr)
  if(numprocs/=lmax) then
     call abort()
  end if
  call setupfls(90)
  call get_grid_param()
  call build_grids()
  call build_KEmats(id)
  n=nel*(nq-1)-1 !!!este funciona en el lobatto
  !!  n=nbas !! este deberia funcionar en todos
  allocate(VV(n,n))
  allocate(HH(n,n,lmax))
  HH(:,:,id+1)=TXX(:,:,1)
  do i=1,n
     do j=1,n
        if(i==j) then
           VV(i,j)=xz(i)
        end if
     end do
  end do  
  if(flag==0) then
     allocate(A(n,n))
     allocate(ee(n,n))
     allocate(c_nl(n,n))
     allocate(c_ni(n,n))
     allocate(vr(n,n))
     li=id+1
     do i=1,n
        do j=1,n
           A(i,j)=HH(i,j,li)
        end do
     end do
     allocate(w(n))
     allocate(work(2))
     allocate(rwork(2*n))
     lwork=-1
     lda=n
     ldb=n
     call zgeev('N','V',lda,A,lda,w,vl,lda,vr,lda,work,lwork,rwork,info)
     lwork=work(1)
     deallocate(work)
     allocate(work(lwork))
     call zgeev('N','V',lda,A,lda,w,vl,lda,vr,lda,work,lwork,rwork,info)
     ee=0.d0
     do i=1,lda
        if(dble(w(i))<50.d0) then
           ee(i,i)=w(i)
        else
           ee(i,i)=50.d0
        end if
     end do
     deallocate(rwork)
     deallocate(work)
     deallocate(w)
     do i=1,lda
        do j=1,lda
           c_nl(j,i)=vr(j,i)
        end do
     end do
     allocate(zImat(n,n))
     zImat=0.d0
     do i=1,lda
        zImat(i,i)=1.d0
     end do
     allocate(ipiv(lda))
     call zgesv(lda,lda,c_nl,lda,ipiv,zImat,lda,info)
     c_ni=zImat
     do i=1,lda
        do j=1,lda
           c_nl(j,i)=vr(j,i)
        end do
     end do
     deallocate(ipiv)
     zImat=0.d0
     A=0.d0
     zImat=matmul(ee,c_ni)
     A=matmul(c_nl,zImat)
     do i=1,n
        do j=1,n
           HH(i,j,li)=A(i,j)
        end do
     end do
     deallocate(vr)
     deallocate(zImat)
     deallocate(ee)
     deallocate(A)
     deallocate(c_nl)
     deallocate(c_ni)
  end if
  nloc=n
  lda=nloc*lmax
  d_nz=nloc
  o_nz=nloc*(nloc-1)
  call MatCreateAIJ(PETSC_COMM_WORLD,nloc,nloc,lda,lda,d_nz,PETSC_NULL_INTEGER,o_nz,PETSC_NULL_INTEGER,AA,ierr)
  call MatSetFromOptions(AA,ierr)
  li=id+1
  do i=1,nloc
     do j=1,nloc
        do lf=1,lmax
           ii=i+nloc*(li-1)
           jj=j+nloc*(lf-1)
           if(lf==li) then
              vAA=HH(i,j,lf)
              call MatSetValues(AA,1,ii-1,1,jj-1,vAA,INSERT_VALUES,ierr)
           end if
        end do
     end do
  end do
  call MatAssemblyBegin(AA,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(AA,MAT_FINAL_ASSEMBLY,ierr)  

  call EPSCreate(PETSC_COMM_WORLD,eps,ierr)
  call EPSSetOperators(eps,AA,PETSC_NULL_DOUBLE,ierr)
  call EPSSetProblemType(eps,EPS_HEP,ierr)
  call EPSSetFromOptions(eps,ierr)
  call EPSSolve(eps,ierr)
  call EPSGetIterationNumber(eps,its,ierr)
  call EPSGetConverged(eps,nconv,ierr)
  call EPSGetErrorEstimate(eps,0,errest,ierr)
  call EPSComputeError(eps,0,EPS_ERROR_RELATIVE,relerr,ierr)
  call VecCreateMPI(PETSC_COMM_WORLD,nloc,lda,x,ierr)
  call VecDuplicate(x,y,ierr)
  if(id==root)then
     write(6, *)
     write(6, *) 'Eigenvalues Computed'
  end if
  do i=0,nconv-1
     call EPSGetEigenpair(eps,i,Er,Ei,x,y,ierr)
     if(id==root) then
        write(6,fmt='(20e16.8)') Er
     end if
  end do
  call EPSGetEigenpair(eps,0,Er,Ei,x,y,ierr)
  Ebs=Er
  call VecNorm(x,NORM_2,nrm,ierr)
  write(*,fmt='(20e16.8)') nrm
  call VecGetArrayF90(x,xx_v,ierr)
  allocate(wav(lda))
  allocate(yy(nloc))
  do i=1,nloc
     yy(i)=xx_v(i)
  end do
  call VecRestoreArrayF90(x,xx_v,ierr)
  call MPI_GATHER(yy,nloc,MPI_DOUBLE_PRECISION,wav,nloc,MPI_DOUBLE_PRECISION,root,PETSC_COMM_WORLD,ierr)

  call VecCreateMPI(PETSC_COMM_WORLD,nloc,lda,z,ierr)
  if(id==0)then
     nrm=sum(yy(:)**2)     
     yy=yy/dsqrt(nrm)
     open(601,file='bound_ascii.dat',status='replace')
     write(601,*) nloc,lmax
     do i=1,nloc
        write(601,fmt='(20e16.8)') yy(i)
        vAA=yy(i)
        call VecSetValue(z,i-1,vAA,ADD_VALUES,ierr);
     end do
     close(601)
  end if
  call VecAssemblyBegin(z,ierr)
  call VecAssemblyEnd(z,ierr);
  deallocate(wav)
  deallocate(yy)

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,"ground.m",FILE_MODE_WRITE,viewen,ierr)
  call VecView(x,viewen,ierr)
  call PetscViewerDestroy(viewen,ierr)
  call VecDestroy(x,ierr)
  call VecDestroy(y,ierr)
  call EPSDestroy(eps,ierr)

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,"HH.m",FILE_MODE_WRITE,viewen,ierr);
  call MatView(AA,viewen,ierr)
  call PetscViewerDestroy(viewen,ierr)
  call MatDestroy(AA,ierr)     

  call MatCreateAIJ(PETSC_COMM_WORLD,nloc,nloc,lda,lda,d_nz,PETSC_NULL_INTEGER,o_nz,PETSC_NULL_INTEGER,VVA,ierr)
  call MatSetFromOptions(VVA,ierr)
  li=id+1
  do i=1,nloc
     do j=1,nloc
        do lf=1,lmax
           ii=i+nloc*(li-1)
           jj=j+nloc*(lf-1)
           call Gcoeff(li-1,0,1,0,lf-1,0,Cval)
           vAA=VV(i,j)*Cval*dsqrt(4.d0*pi/3.d0)
           call MatSetValues(VVA,1,ii-1,1,jj-1,vAA,INSERT_VALUES,ierr)
        end do
     end do
  end do
  call MatAssemblyBegin(VVA,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(VVA,MAT_FINAL_ASSEMBLY,ierr)
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,"VV.m",FILE_MODE_WRITE,viewen,ierr);
  call MatView(VVA,viewen,ierr)
  call PetscViewerDestroy(viewen,ierr)
  
  call VecDuplicate(z,x,ierr)
  call VecZeroEntries(x,ierr)
  call MatMult(VVA,z,x)
  call VecGetArrayF90(x,ll_v,ierr)
  allocate(wav(lda))
  allocate(yy(nloc))
  do i=1,nloc
     yy(i)=ll_v(i)
  end do
  call VecRestoreArrayF90(x,ll_v,ierr)
  call MPI_GATHER(yy,nloc,MPI_DOUBLE_PRECISION,wav,nloc,MPI_DOUBLE_PRECISION,root,PETSC_COMM_WORLD,ierr)
  if(id==0)then
     open(601,file='vector_ascii.dat',status='replace') 
     write(601,*) lda
     do i=1,lda
        write(601,fmt='(20e16.8)') wav(i)
     end do
     close(601)
  end if

  deallocate(HH)
  deallocate(VV)
  deallocate(wav)
  deallocate(yy)
  call VecDestroy(z,ierr)     
  call VecDestroy(x,ierr)     
  call MatDestroy(VVA,ierr)     
  call SlepcFinalize(ierr)  

end program Hatom_fem
