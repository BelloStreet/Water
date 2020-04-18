!-----------------------------------!
subroutine build_grids()
  
  use dvrecs
  implicit none
  complex*16 xs(nq), ws(nq) 
  integer i, j, ii

  allocate(xz(nbas), wz(nbas))
  !create nth order DVR on [-1, 1]
  allocate(xq(nq), wq(nq))
  call gl_quad(nq, xq, wq) !!!!llama a la cuadratura de gauss-lobbato
  allocate(xqz(nq))
  xqz=xq
  
  !!basicamente calcula los ptos de cuadratura entre las fonteras del FEM
  ii=0
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
  xs=az(nel-1)+0.5d0*(xq+1)*(az(nel)-az(nel-1))
  ws=0.5*wq*(az(nel)-az(nel-1))
  wz(ii)=wz(ii)+ws(1)
  do j=2, nq-1
     ii=ii+1
     xz(ii)=xs(j)
     wz(ii)=ws(j)
  end do
  
end subroutine build_grids
!-----------------------------------!
!-----------------------------------!
subroutine build_dq_nw(n, d)
  
  implicit none
  integer n, i, j, k
  real*8 d(n, n), x(n), w(n), pr
  
  call gl_quad(n, x, w) !!llama a la cuadratura y luego calcula las derivadas de los lagrange
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
subroutine build_dq(n, d)
  
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
  

  call build_dq_nw(nq, Dq)

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
!--------------------------------------!
subroutine build_KEmats()
  
  use dvrecs
  use molecule
  implicit none

  integer i, j, l, ii, jj
  integer i1, i2, o, m, lam, info
  complex*16 de

  integer ipiv(nbas),lwork,lda,ldb,ncount
  real*8 Dq(nq, nq),alpha
  complex*16 D2(nbas, nbas)
  complex*16, allocatable::c_nl(:,:),c_ni(:,:),wener(:),work(:),vr(:,:),vl(:,:),ee(:,:)
  complex*16 zImat(nbas, nbas),zR1(nbas,nbas)
  real(kind=8),allocatable::rwork(:)
  character:: jobz='V',uplo='L'
  call build_dq_nw(nq, Dq) !!calcula las derivadas del Lagrange y las guarda en Dq
  D2=0.d0
  do l=1, nel
     i1=1
     i2=nq
     if(l==1) i1=2
     if(l==nel) i2=nq-1
     o=(nq-1)*(l-1)
     do i=i1, i2
        ii=(i-1)+o
        do j=i1, i2
           jj=(j-1)+o
           de=0.d0
           do m=1, nq
              de=de-dq(m, i)*dq(m, j)*2.d0/(az(l)-az(l-1))*(wq(m)/sqrt(wz(ii)*wz(jj)))  !!supongo que generaliza las derivadas
           end do
           D2(ii, jj)=D2(ii, jj)-de
        end do
     end do
  end do
  allocate(TXX(nbas, nbas, 0:lamax))
  TXX=0.d0
  do lam=0, lamax
     TXX(:, :, lam)=D2
     do i=1, nbas
        TXX(i, i, lam)=TXX(i, i, lam)+dble(lam*(lam+1))/(xz(i)**2) !!evidentemente kinetic energy plus centrifuga
     end do                                                        !!equation 64 del topical
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
     TIXX(:, :, lam)=zR1                                          !!se invierte la matriz para calcularse el electron-electron 
  end do

  do lam=0, lamax
     D2=TXX(:, :, lam)
     allocate(ee(nbas,nbas))
     allocate(c_nl(nbas,nbas))
     allocate(c_ni(nbas,nbas))
     allocate(wener(nbas))
     allocate(work(2))
     allocate(vr(nbas,nbas))
     allocate(rwork(2*lda))
     lwork=-1
     lda=nbas
     ldb=nbas
     call zgeev('N','V',lda,D2,lda,wener,vl,lda,vr,lda,work,lwork,rwork,info)
     lwork=work(1)
     deallocate(work)
     allocate(work(lwork))
     call zgeev('N','V',lda,D2,lda,wener,vl,lda,vr,lda,work,lwork,rwork,info)
     deallocate(rwork)
     ee=0.d0
     do i=1,lda
        ! ee(i,i)=wener(i)
        if(dble(wener(i))<150.0d0) then
           ee(i,i)=wener(i)
        else
           ee(i,i)=150.0d0
        end if
     end do
     deallocate(work)
     deallocate(wener)
     do i=1,lda
        do j=1,lda
           c_nl(j,i)=vr(j,i)
        end do
     end do
     zImat=0.d0
     do i=1,lda
        zImat(i,i)=1.d0
     end do
     call zgesv(lda,lda,c_nl,lda,ipiv,zImat,lda,info)
     c_ni=zImat
     do i=1,lda
        do j=1,lda
           c_nl(j,i)=vr(j,i)
        end do
     end do
     zImat=matmul(ee,c_ni)
     D2=matmul(c_nl,zImat)
     TXX(:, :, lam)=D2

     deallocate(vr)
     deallocate(ee)
     deallocate(c_nl)
     deallocate(c_ni)
  end do
  
end subroutine Build_KEmats
!--------------------------------------!
