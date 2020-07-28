!----------------------------------!
subroutine get_grid_param_Radau()

  use numbers, only: eye, pi
  use dvrecsRad
  implicit none

  integer i, kel, bel, ilm
  real*8 rval
  character*1 keep

  ! Aug 2018 - allows for following possible grid calculations (keep is flag):
  !   DEFAULT - pure Lobatto with ECS (nel=kel, theta as read in)
  !   keep=b  - pure Lobatto on a subset of the grid as read (nel=bel,theta=0)
  !   keep=t  - pure Lobatto on a subset of the grid as read (nel=rel,theta=0)
  !   keep=L  - one Radau-Laguerre tail appended after rel (nel=rel+1)
  read(501, *)
  read(501, *)
  read(501, *) nq, mq ! Lobatto order, and Radau order for last FEM (if used)
  read(501, *)
  read(501, *) theta
  read(501, *)
  read(501, *) keep, kel, rel, bel

  select case (keep)
     case('b')
        ! boundstate solution (Lobatto only)
        rel=abs(bel)
        nel=rel
        theta=0.d0
        nbas=nel*(nq-1)-1
        mq=0
     case('t')
        ! Lobatto-only real grid for time propagation (depricated Aug.2018)
        nel=rel
        theta=0.d0
        nbas=nel*(nq-1)-1
        mq=0
     case('L')
        ! keep a REAL (non-rotated) grid consisting of 
        ! all real Lobatto FEMs being read in and add one FEM with Laguerre tail
        nel=rel+1
        nbas=(nel-1)*(nq-1)-1+mq ! adding Radau tail, keeping last element
     case default
        ! default case is to keep a Lobatto-only ECS rotated grid
        nel=kel
        !rel is kept as read in
        nbas=nel*(nq-1)-1
        mq=0
  end select

  read(501, *)      
  allocate(ar(0:nel), az(0:nel))
  do i=0, kel
     read(501, *) rval
     if(i.le.nel) ar(i)=rval
  end do

  R0=ar(rel) 
  eit=exp(eye*(pi*theta/180.d0))
  az(0:rel)=ar(0:rel)
  az(rel+1:nel)=R0+eit*(ar(rel+1:nel)-R0)
  !!! *** NOTE az(nel) affects Poisson solve integrals ... boundary conditions
  !!! *** for Laguerre tail, az(nel) and ar(nel) will be rewritten in setgridRad

  rbas=rel*(nq-1)-1
  nbas2=nbas**2

  ! alphaRad is needed for the complex Radau tail (ignored otherwise) 
  !           0.3d0 by Scrinzi used for Coulomb problems
  read(501, *) 
  read(501, *) alphaRad  

end subroutine get_grid_param_Radau
!----------------------------------!
