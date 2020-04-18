!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module dvrecs

  implicit none
  
  integer,parameter::unit_stdout=6
  integer nq, mq, nel, rel, nbas, nbas2, rbas, ndvr,flag
  integer lamax
  real*8 theta, R0, Zcharge
  complex*16 eit, eit_conj, eit2

  !element boundaries
  real*8, allocatable:: ar(:)
  complex*16, allocatable:: az(:)
  
  !one element DVR
  real*8, allocatable:: xq_l(:), wq_l(:), xq_r(:), wq_r(:)
  complex*16, allocatable:: xqz(:)
  complex*16, allocatable:: xs(:), ws(:)
  
  !full grid
  complex*16, allocatable:: xz(:), wz(:), siwz(:)
  
  !2 x kinetic energy and inverse
  complex*16, allocatable:: TXX(:, :, :), TIXX(:, :, :)

contains

subroutine get_grid_param()

  use numbers
  implicit none
  
  integer i
  
  open(unit=501,file='indvr.main')
  read(501, *)
  read(501, *)
  read(501, *) nq, mq !Labatto order, Radau order
  read(501, *)
  read(501, *) nel, rel
  read(501, *)
  allocate(ar(0:nel), az(0:nel))
  do i=0, nel
     read(501, *) ar(i)
  end do
  read(501, *)
  read(501, *) theta,flag
  
  R0=ar(rel)
  eit=exp(eye*(pi*theta/180.d0))
  az(0:rel)=ar(0:rel)
  az(rel+1:nel)=R0+eit*(ar(rel+1:nel)-R0)
  
  rbas=rel*(nq-1)-1
  if(rel==nel)then
     nbas=nel*(nq-1)-1
  else
     ndvr=(nel-1)*(nq-1)+mq !Generalized for Radua with different order mq
     nbas=ndvr-2 ! Only taking off basis function at the origin
  end if
  nbas2=nbas**2
  
  
end subroutine get_grid_param

end module dvrecs
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
