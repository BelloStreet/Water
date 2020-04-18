!++++++++++++++++++++++++++++++++++++++++++++++!
function Ylm(l, m, theta, phi)

  use numbers
  implicit none
  integer l, m, absm
  real*8 theta, phi
  real*8 rnrm, cost, Pm(0:l+1, 0:l+1)
  complex*16 Ylm


  absm=abs(m)
  Ylm=0.d0
  if(absm>l) goto 86
  rnrm=dble(2*l+1)*rfl(l-absm)/rfl(l+absm)
  rnrm=sqrt(0.25d0*dble(rnrm)/pi)
  if(m<0) rnrm=rnrm*((-1.d0)**m)
  cost=cos(theta)
  call Plm(l+1, absm, cost, Pm)
  Ylm=rnrm*Pm(absm, l)*exp(eye*dble(m)*phi)

86 continue

  return

end function Ylm
!++++++++++++++++++++++++++++++++++++++++++++++!
!++++++++++++++++++++++++++++++++++++++++++++++!
function Ylmc(l, m, cost, phi)

  use numbers
  implicit none
  integer l, m, absm
  real*8 phi
  real*8 rnrm, cost, Pm(0:l+1, 0:l+1)
  complex*16 Ylmc

  Ylmc=0.d0
  absm=abs(m)
  if(absm>l) goto 86
  rnrm=dble(2*l+1)*rfl(l-absm)/rfl(l+absm)
  rnrm=sqrt(0.25d0*dble(rnrm)/pi)
  if(m<0) rnrm=rnrm*((-1.d0)**m)
  call Plm(l+1, absm, cost, Pm)
  Ylmc=rnrm*Pm(absm, l)*exp(eye*dble(m)*phi)

  86 continue

  return

end function Ylmc
!++++++++++++++++++++++++++++++++++++++++++++++!
!----------------------------------------------!
subroutine Plm(l, am, x, Pm)

  implicit none
  integer l, am, i, j
  real*8 x, Pm(0:l, 0:l)
  real*8 xq, xs, ls

  Pm(0:l, 0:l)=0.d0
  Pm(0, 0)=1.d0
  if(abs(x)==1.d0)then
     do i=1, l
        Pm(0, i)=x**i
     end do
  else
     ls=1.d0
     if (abs(x)>1.0d0) ls=-1.d0
     xq=sqrt(ls*(1.0d0-x**2))
     xs=ls*(1.0d0-x**2)
     do i=1,am
        pm(i,i)=-ls*(2.0d0*i-1.0d0)*xq*pm(i-1,i-1)
     end do
     do i=0,am
        pm(i,i+1)=(2.0d0*i+1.0d0)*x*pm(i,i)
     end do
     do i=0,am   
       do j=i+2,l
           pm(i,j)=((2.0d0*j-1.0d0)*x*pm(i,j-1)-(i+j-1.0d0)*pm(i,j-2))/(j-i)
        end do
     end do
  end if

end subroutine Plm
!----------------------------------------------!
