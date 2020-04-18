!-----------------------------------!
subroutine build_grids()

  use dvrecs
  implicit none
  integer,parameter:: debug=0 !For debugging
  complex*16 xs_l(nq), ws_l(nq), xs_r(mq), ws_r(mq)!Temps 
  real*8 lag_alpha
  complex*16 lscal
  integer i, j, ii, pp

  allocate(xz(nbas), wz(nbas))

  !create nth order DVR on [-1, 1]
  allocate(xq_l(nq), wq_l(nq))!Raw points and weights
  call gl_quad(nq, xq_l, wq_l)
  allocate(xqz(nq))
  xqz=xq_l
  !Write raw points and weights if desired
  if(debug/=0) then
      write(unit_stdout,*) 'Raw DVR points and weights on -1,1'
      write(unit_stdout,'(6f11.6)')(xq_l(i),i=1,nq)
      write(unit_stdout,*)
      write(unit_stdout,'(6f11.6)')(wq_l(i),i=1,nq)
      write(unit_stdout,*)
    endif

  !create nth order Radau DVR on [0, Inf]
  allocate(xq_r(mq), wq_r(mq))!Raw points and weights
  call gr_quad(mq, xq_r, wq_r)
  !Write raw points and weights if desired
  if(debug/=0) then
      write(unit_stdout,*) 'Raw DVR points and weights on 0,inf'
      write(unit_stdout,'(6f11.6)')(xq_r(i),i=1,mq)
      write(unit_stdout,*)
      write(unit_stdout,'(6f11.6)')(wq_r(i),i=1,mq)
      write(unit_stdout,*)
    endif

    lag_alpha=0.3d0 !Scrinzi says this was a good choice.
    lscal=eit/(2*lag_alpha) !Transformation to Laguerre boundaries on a complex ray

!*******************************************
!Creating Grid from raw points--------------
!*******************************************
if(rel==nel)then
    pp=rel-1
else
    pp=rel
endif
ii=0
do i=1, pp !Labatto body 
        xs_l=az(i-1)+0.5d0*(xq_l+1)*(az(i)-az(i-1))
        ws_l=0.5*wq_l*(az(i)-az(i-1))
        if(ii.ne.0) wz(ii)=wz(ii)+ws_l(1)
        do j=2, nq
            ii=ii+1
            xz(ii)=xs_l(j)  
            wz(ii)=ws_l(j)
        end do
end do

if (rel==nel) then
    xs_l=az(i-1)+0.5d0*(xq_l+1)*(az(i)-az(i-1))
    ws_l=0.5*wq_l*(az(i)-az(i-1))
    if(ii.ne.0) wz(ii)=wz(ii)+ws_l(1)
    do j=2, nq-1
        ii=ii+1
        xz(ii)=xs_l(j)  
        wz(ii)=ws_l(j)
    end do
else !Radau tail
    xs_r=lscal*xq_r+R0
    ws_r=lscal*wq_r
    wz(ii)=wz(ii)+ws_r(1)
    do j=2, mq-1
        ii=ii+1
        xz(ii)=xs_r(j)
        wz(ii)=ws_r(j)
    end do
endif

end subroutine build_grids
!-----------------------------------!
!-----------------------------------!
subroutine build_dq_nw(n, d, quad)
  
  use dvrecs
  implicit none
  integer n, i, j, k, quad
  real*8 d(n, n), x(n), w(n), pr
  
  if (quad==1) then 
    call gl_quad(n, x, w) !Evaluate Labatto derivative
  else
    call gr_quad(n, x, w) !Evaluate Radau derivative
  end if

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
           d(i, j)=pr
        end if
    end do
  end do

end subroutine build_dq_nw
!-----------------------------------!
!-----------------------------------!
subroutine build_dq(n, d, quad)
  
  implicit none
  integer n, i, j, k, quad
  real*8 d(n, n), x(n), w(n), pr
  
  if (quad==1) then
    call gl_quad(n, x, w) !Evaluate Labatto derivative
  else
    call gr_quad(n, x, w) !Evaluate Radau derivative
  end if

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
           d(i, j)=pr*sqrt(w(i)/w(j))
        end if
     end do
  end do

end subroutine build_dq
!-----------------------------------!
!-----------------------------------!
subroutine derivmatr(nel_b, Dmat)

  use dvrecs
  implicit none
  
  integer i, j, l, ii, jj, o, i1, i2
  integer nel_b
  real*8 Dmat(nel_b*(nq-1)-1, nel_b*(nq-1)-1)
  real*8 Dq(nq, nq)
  

  call build_dq_nw(nq, Dq, 1)

  Dmat=0.d0
  do l=1, nel_b
     i1=1
     i2=nq
     if(l==1) i1=2
     if(l==nel_b) i2=nq-1
     o=(nq-1)*(l-1)
     do i=i1, i2
        ii=(i-1)+o
        do j=i1, i2
           jj=(j-1)+o
           Dmat(ii, jj)=Dmat(ii, jj)+&
                Dq(i, j)*sqrt(wz(ii)/wz(jj))*2.d0/(az(l)-az(l-1))
        end do
     end do
  end do
  do l=1, nel_b-1
     ii=(nq-1)*l
     Dmat(ii, :)=0.5d0*Dmat(ii, :)
  end do

end subroutine derivmatr
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
subroutine build_KEmats(id)
  
  use dvrecs
  use numbers, only: pi, eye
  implicit none

  integer i, j, l, ii, jj
  integer i1, i2, o, m, lam, info
  complex*16 de

  integer ipiv(nbas),id
  real*8 dq_l(nq, nq), dq_r(mq, mq), lag_alpha
  complex*16 D2(nbas, nbas), scal_w
  complex*16 zImat(nbas, nbas), zR1(nbas, nbas)

  lag_alpha=0.3d0 !Scrinzi
  eit2=exp(eye*(pi*theta*(-2)/180.d0))
  eit_conj=conjg(eit)

  call build_dq_nw(nq, dq_l, 1)
  call build_dq_nw(mq, dq_r, 2)
  
  D2=0.d0
    do l=1, nel
        i1=1
        i2=nq
        o=(nq-1)*(l-1)
        if(l==1) i1=2
        if(rel==nel) then
            if(l==nel) i2 = nq-1
        else
            if(l==nel) i2 = mq-1 
        endif
        do i=i1, i2
            ii=(i-1)+o
            do j=i1, i2
                jj=(j-1)+o
                de=0.d0
                if (l<=rel) then !Labatto body
                    scal_w=2.d0/(az(l)-az(l-1))/sqrt(wz(ii)*wz(jj))
                    do m=1, nq
                        de=de+dq_l(m, i)*dq_l(m, j)*scal_w*wq_l(m)
                    end do
                    D2(ii, jj)=D2(ii, jj)+de
                else !Radau tail
                    scal_w=2*lag_alpha/(eit*sqrt(wz(ii)*wz(jj))) 
                    do m=1, mq
                        de=de+dq_r(m, i)*dq_r(m, j)*wq_r(m)*scal_w
                    end do
                    if (i==j) then
                        D2(ii,jj)=D2(ii,jj)+de- &
                        lag_alpha*eit_conj*scal_w*(wz(jj)*dq_r(j,i)+wz(ii)*dq_r(i,j))+ & 
                        lag_alpha**2*eit2
                    else
                        D2(ii,jj)=D2(ii,jj)+de- &
                        lag_alpha*eit_conj*scal_w*(wz(jj)*dq_r(j,i)+wz(ii)*dq_r(i,j))
                    end if
                end if
            end do
        end do
    end do

  ! allocate(TXX(nbas, nbas, 0:lamax))
  ! TXX=0.d0
  ! do lam=0, lamax
  !    TXX(:, :, lam)=D2
  !    do i=1, nbas
  !       TXX(i, i, lam)=TXX(i, i, lam)+dble(lam*(lam+1))/(xz(i)**2)
  !    end do
  ! end do
  ! zImat=0.d0
  ! do i=1, nbas
  !    zImat(i, i)=1.d0
  ! end do
  ! allocate(TIXX(nbas, nbas, 0:lamax))
  ! do lam=0, lamax
  !    D2=TXX(:, :, lam)
  !    zR1=zImat
  !    call zgesv(nbas, nbas, D2, nbas, ipiv, zR1, nbas, info)
  !    TIXX(:, :, lam)=zR1
  ! end do

  allocate(TXX(nbas, nbas, 1))
  TXX=0.d0
  TXX(:, :, 1)=0.5d0*D2
  do i=1, nbas
     TXX(i, i, 1)=TXX(i, i, 1)+0.5d0*(dble(id*(id+1))/(xz(i)**2))-(1.d0/xz(i)) !!evidentemente kinetic energy plus centrifuga
  end do
  
end subroutine Build_KEmats
