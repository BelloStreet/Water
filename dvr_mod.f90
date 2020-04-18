!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module dvrecs

  implicit none
  
  integer nq, nel, rel, nbas, nbas2, rbas, lamax,lmax,flag
  real*8 theta, R0

  !element boundaries
  real*8, allocatable:: ar(:),ar1(:)
  complex*16, allocatable:: az(:)
  
  !one element DVR
  real*8, allocatable:: xq(:), wq(:)
  complex*16, allocatable:: xqz(:)
  
  !full grid
  complex*16, allocatable:: xz(:), wz(:)
  
  !2 x kinetic energy and inverse
  complex*16, allocatable:: TXX(:, :, :), TIXX(:, :, :)

  contains
    !----------------------------------!
    subroutine get_grid_param()
      
      use numbers
      implicit none
      
      integer i
      complex*16 eit
      
      open(unit=501,file='indvr.main')
      read(501, *)
      read(501, *)
      read(501, *) nq
      read(501, *)
      read(501, *) nel, rel
      read(501, *)      
      allocate(ar(0:nel), az(0:nel))
      do i=0, nel
         read(501, *) ar(i)
      end do
      read(501, *)
      read(501, *) theta,flag
      read(501, *) lmax
      lamax=2*lmax
      close(501)

      R0=ar(rel) 
      eit=exp(eye*(pi*theta/180.d0))
      az(0:rel)=ar(0:rel)
      az(rel+1:nel)=R0+eit*(ar(rel+1:nel)-R0)

      rbas=rel*(nq-1)-1
      nbas=nel*(nq-1)-1  !!Numero de FEM X (Numero de DVR-1)-1
      nbas2=nbas**2


    end subroutine get_grid_param
    !----------------------------------!
    

end module dvrecs
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
