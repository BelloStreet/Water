!!! these routines permit the last FEM appended to be Gauss Radau to infinity
!!! updated August 2018
!-----------------------------------!
subroutine gl_quad(n, t, w)

  implicit none
  real*8 t(n), w(n), scr(n)
  real*8 alpha, beta, endpts(2)
  integer n, kind, kpts
  
 !generate Gauss-Lobatto quadrature points and weights
  alpha=0.d0
  beta=0.d0
  endpts(1)=-1.d0
  endpts(2)=1.d0
  kind=1
  kpts=2
  call gaussq(kind, n, alpha, beta, kpts, endpts, scr, t, w)
  t(1)=-1.d0
  t(n)=1.d0

end subroutine gl_quad
!-----------------------------------!
!-----------------------------------!
subroutine gleg_quad(n, t, w)

  implicit none
  real*8 t(n), w(n), scr(n)
  real*8 alpha, beta, endpts(2)
  integer n, kind, kpts
  
 !generate Gauss-Lobatto quadrature points and weights
  alpha=0.d0
  beta=0.d0
  endpts(1)=-1.d0
  endpts(2)=1.d0
  kind=1
  kpts=0
  call gaussq(kind, n, alpha, beta, kpts, endpts, scr, t, w)
  !t(1)=-1.d0
  !t(n)=1.d0

end subroutine gleg_quad
!-----------------------------------!
!-----------------------------------!
subroutine gr_quad(n, t, w)

  implicit none
  real*8 t(n), w(n), scr(n)
  real*8 alpha, beta, endpts(2)
  integer n, kind, kpts

 !generate Gauss-Radau quadrature points and weights
  alpha=0.d0
  beta=0.d0
  endpts(1) = 0.d0
  endpts(2) = 1.d0 
  kind=6 !Laguerre boundaries
  kpts=1 !Fixes left endpoint, i.e. Radau quadrature
  call gaussq(kind, n, alpha, beta, kpts, endpts, scr, t, w)
  t(1)=0.d0

end subroutine gr_quad
!-----------------------------------!
!-----------------------------------!
subroutine build_grids_Radau()

  use dvrecsRad
  implicit none
  complex*16 xs(nq), ws(nq), radscale
  integer i, j, ii

  allocate(xz(nbas), wz(nbas))
  !create nth order DVR on [-1, 1]
  allocate(xq(nq), wq(nq))
  call gl_quad(nq, xq, wq)
  allocate(xqz(nq))
  xqz=xq

  ii=0
  !do all of the FEMs up to the last one defined 
  do i=1, nel-1  
     xs=az(i-1)+0.5d0*(xq+1)*(az(i)-az(i-1))
     ws=0.5*wq*(az(i)-az(i-1))
     if(ii.ne.0) wz(ii)=wz(ii)+ws(1)
     do j=2, nq
        ii=ii+1
        xz(ii)=xs(j)  
        wz(ii)=ws(j)
     end do
  end do
  ! ii now sits at ending boundary before last FEM
  ! last point written at xz(ii) is the penultimate FEM right endpoint

  if(mq.eq.0) then
     ! last finite element is Lobatto, handled just as before
     xs=az(nel-1)+0.5d0*(xq+1)*(az(nel)-az(nel-1))
     ws=0.5*wq*(az(nel)-az(nel-1))
     wz(ii)=wz(ii)+ws(1)
     do j=2, nq-1  ! throwing away last DVR to enforce zero at boundary
        ii=ii+1
        xz(ii)=xs(j)
        wz(ii)=ws(j)
     end do
     !! *** end of grid, az(nel), already set in get_grid_param_Radau() 
  else
     ! last finite element is done with Radau-Lagurerre quadrature
     ! incorporating Radau tail as (rel+1) FEM

     ! prepare mth order Radau quadrature on [0,infty)
     allocate(xqR(mq), wqR(mq), xqzR(mq), wqzR(mq))
     call gr_quad(mq, xqR, wqR)

     radscale=eit/(2.d0*alphaRad)
     xqzR=radscale*xqR+R0  ! complex-scaled Radau points used in KE matrix
     wqzR=radscale*wqR

     wz(ii)=wz(ii)+wqzR(1)  ! adjust connection to last FEM-DVR at R0
     do j=2, mq ! keeping last Radau point on the grid (mq complex DVRs)
        ii=ii+1
        xz(ii)=xqzR(j)
        wz(ii)=wqzR(j)
    end do
    do i = 1, nbas
      print *, xz(ii), wz(ii)
    end do
    !! *** Set the end of the last FEM to the last quadrature points
    !!  *** this is now done here, and not in get_grid_param_Radau
    az(nel)=xqzR(mq)
    ! map last complex point back to real r (same ECS formula as always)
    ar(nel)=(az(nel)-R0)*conjg(eit)+R0

 end if

end subroutine build_grids_Radau
!-----------------------------------!
!-----------------------------------!
subroutine build_dq_nw(n, d)
  
  ! NOTE ... output matrix d is REAL (it's Real-valued along contour in Lobatto)
  ! xqz, defined in build_grids, has nq-th order complex Lobatto pts)
  implicit none
  integer n, i, j, k
  real*8 d(n, n), x(n), w(n), pr
  
  call gl_quad(n, x, w)
  do i=1, n
     d(i, i)=0.d0
     do k=1, n
        if(i.ne.k) d(i, i)=d(i, i)+1.d0/(x(i)-x(k))
     end do
     do j=1, n
        if(i.ne.j)then
           pr=1.d0/(x(j)-x(i))
           do k=1, n
              if((k.ne.i).and.(k.ne.j)) pr=pr*(x(i)-x(k))/(x(j)-x(k))
           end do
           d(i, j)=pr ! check the ordering ... may not be symmetrically defined
        end if
     end do
  end do

end subroutine build_dq_nw
!-----------------------------------!
!-----------------------------------!
subroutine build_dq_Radau(n, d, flag)
  
  ! Produce the derivative matrix in COMPLEX array d 
  ! (assume Radau quadrature already rotated into complex FEM (xqzR, wqzR)
  ! This version scales by the sqrt(wi/wj) when flag=1 (not usually done)

  use dvrecsRad, only : mq, xqzR, wqzR
  implicit none
  integer n, i, j, k, flag
  complex*16 d(n, n), pr 
  
  if (n.ne.mq) stop ' mismatch in Radau Quadrature order on building KE'
  do i=1, n
     d(i, i)=0.d0
     do k=1, n
        if(i.ne.k) d(i, i)=d(i, i)+1.d0/(xqzR(i)-xqzR(k))
     end do
     do j=1, n
        if(i.ne.j)then
           pr=1.d0/(xqzR(j)-xqzR(i))
           do k=1, n
              if((k.ne.i).and.(k.ne.j))pr=pr*(xqzR(i)-xqzR(k))/(xqzR(j)-xqzR(k))
           end do
           d(i, j)=pr ! check the ordering ... may not be symmetrically defined
           if(flag.eq.1) d(i, j)=pr*sqrt(wqzR(i)/wqzR(j))
        end if
     end do
  end do

end subroutine build_dq_Radau
!-----------------------------------!
!+++++++++++++++++++++++++++++++++++!
function xi(i, x, xqz, n)

  implicit none
  integer i, n, j
  complex*16 xi, pr, x, xqz(n), xm
  
  !xm=xqz(i)
  !xi=product(((x-xqz)/(xm-xqz)), mask= xqz.ne.xm)

  pr=1.d0
  do j=1, i-1
     pr=pr*(x-xqz(j))/(xqz(i)-xqz(j))
  end do
 do j=i+1, n
     pr=pr*(x-xqz(j))/(xqz(i)-xqz(j))
  end do
  xi=pr

  return
  
end function xi
!+++++++++++++++++++++++++++++++++++!
!-----------------------------------!
subroutine build_KEmats_Radau()
  
  use dvrecsRad
  implicit none

  integer i, j, l, ii, jj
  integer i1, i2, o, m, lam, info
  complex*16 de

  integer ipiv(nbas)
  real*8 Dq(nq, nq)
  complex*16 DqRad(mq, mq), D2(nbas, nbas), eitm2
  complex*16 zImat(nbas, nbas), zR1(nbas, nbas)
  logical appendRad

  call build_dq_nw(nq, Dq)

  eitm2=1.d0/(eit**2)  ! same as e^{-2i*theta}

  if (mq.eq.0) then !decide if pure Lobatto or append Radau tail
     lam=nel
     appendRad=.false.
  else
     lam=nel-1
     appendRad=.true.
     call build_dq_Radau(mq, DqRad, 0)
  end if

  D2=0.d0
  do l=1, lam 
     i1=1
     i2=nq
     if(l==1) i1=2          ! omit first DVR point
     if(l==nel) i2=nq-1     ! omit last DVR point on purely real FEM Lobatto 
     o=(nq-1)*(l-1)
     do i=i1, i2
        ii=(i-1)+o
        do j=i1, i2
           jj=(j-1)+o
           de=0.d0
           do m=1, nq
              de=de+dq(m, i)*dq(m, j)*wq(m)
           end do
           D2(ii, jj)=D2(ii, jj)+de*2.d0/((az(l)-az(l-1))*sqrt(wz(ii)*wz(jj)))
        end do
     end do
  end do

  ! now compute and append the last Radau FEM (if needed)
  if(appendRad) then
     i1=1
     i2=mq
     o=(nq-1)*(nel-1)
     do i=i1, i2
        ii=(i-1)+o
        do j=i1, i2
           jj=(j-1)+o
           de=0.d0
           do m=1, mq
              de=de+dqRad(m, i)*dqRad(m, j)*wqzR(m)
           end do
           if(i==j) then
              de=de-alphaRad*conjg(eit)&
                   *(wqzR(j)*dqRad(i, j)+wqzR(i)*dqRad(j, i))&
                   +(alphaRad**2)*eitm2*wqzR(i)
              D2(ii,jj)=D2(ii,jj)+de/sqrt(wz(ii)*wz(jj))
           else
              de=de-alphaRad*conjg(eit)&
                   *(wqzR(j)*dqRad(j, i)+wqzR(i)*dqRad(i, j))
              D2(ii,jj)=D2(ii,jj)+de/sqrt(wz(ii)*wz(jj))
           end if
        end do
     end do
  end if

  allocate(TXX(nbas, nbas, 0:lamax))
  TXX=0.d0
  do lam=0, lamax
     TXX(:, :, lam)=D2
     do i=1, nbas
        TXX(i, i, lam)=TXX(i, i, lam)+dble(lam*(lam+1))/(xz(i)**2)
     end do
  end do
  zImat=0.d0
  do i=1, nbas
     zImat(i, i)=1.d0
  end do
  allocate(TIXX(nbas, nbas, 0:lamax))
  do lam=0, lamax
     D2=TXX(:, :, lam)
     zR1=zImat
     call zgesv(nbas, nbas, D2, nbas, ipiv, zR1, nbas, info)
     TIXX(:, :, lam)=zR1
  end do
  
end subroutine Build_KEmats_Radau
!-----------------------------------!
