program Hatom_fem
  use numbers
  use dvrecs
  implicit none
  integer::i,j,k,n,v,l,lda,ldb,info,li,lf,ni,nf,nloc,id,numprocs
  integer::jj,ii,ierr,root=0,its,nconv,lwork,ncount
  integer,allocatable::ipiv(:)
  real(kind=8)::Cval,Ebs,e
  real(kind=8),allocatable::rwork(:),yy(:),wav(:)
  complex*16,allocatable::HH(:,:,:),VV(:,:),zImat(:,:),vr(:,:),vl(:,:)
  complex*16,allocatable::ee(:,:),c_nl(:,:),c_ni(:,:),w(:),work(:),A(:,:)
  character:: jobz='V',uplo='L'
  call setupfls(90)
  call get_grid_param()
  call build_grids()
  call build_KEmats()
  open(601,file="quadrature.dat",status='replace')
  write(601,*) lamax+1 
  write(601,*) nbas
  do i=1,nbas
     write(601,fmt='(20e16.8)') xz(i),wz(i)
  end do
  close(601)
  open(601,file="kinetic.dat",status='replace')
  do l=0,lamax
     do i=1,nbas
        do j=1,nbas
           write(601,fmt='(20e16.8)') TXX(j,i,l)
        end do
     end do
  end do
  do l=0,lamax
     do i=1,nbas
        do j=1,nbas
           write(601,fmt='(20e16.8)') TIXX(j,i,l)
        end do
     end do
  end do
  close(601)
end program Hatom_fem
