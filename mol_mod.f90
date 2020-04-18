!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module molecule
  
  implicit none

  integer nnuc
  real*8, allocatable:: cent(:, :), rnuc(:, :)
  
  integer nlmps, ellmax
  integer, allocatable:: lmps(:, :)


contains
  !----------------------------------!
  subroutine get_nucs()
    
    use numbers
    implicit none
    
    integer inuc

    open(501, file='inrun.inp')
    
    read(501, *)
    read(501, *) nnuc !!!numero de nucleos
    read(501, *)
    allocate(cent(4, nnuc))
    do inuc=1, nnuc
       read(501, *) cent(1, inuc), cent(2, inuc),  &  !!!Posicion de los nucleos y carga
            cent(3, inuc), cent(4, inuc)
    end do

    !rnuc is polar coordinates (r, theta, phi) of the nuclei
    allocate(rnuc(3, nnuc)) 
    do inuc=1, nnuc
       rnuc(1, inuc)=dsqrt(cent(1, inuc)**2+cent(2, inuc)**2+cent(3, inuc)**2)
       if(rnuc(1, inuc).ne.0.d0)then
          rnuc(2, inuc)=acos(cent(3, inuc)/rnuc(1, inuc))
          if(cent(1, inuc)==0.d0)then 
             if(cent(2, inuc)==0.d0)then 
                rnuc(3, inuc)=0.d0 !x=0, y=0
             else 
                rnuc(3, inuc)=pi*(1.d0-sign(0.5d0, cent(2, inuc))) !x=0, y.ne.0
             end if
          else
             if(cent(2, inuc)==0.d0)then
                rnuc(3, inuc)=0.5d0*pi*(1.d0-sign(1.d0, cent(2, inuc))) !x.ne.0, y=0
             else
                rnuc(3, inuc)=atan(cent(2, inuc)/cent(1, inuc)) !x.ne.0, y.ne.0
             end if
          end if
       else
          rnuc(:, inuc)=0.d0
       end if
    end do
    

  end subroutine get_nucs
  !----------------------------------!
  !----------------------------------!
  subroutine get_lmps()
    
    use dvrecs, only: lamax
    
    implicit none
    
    integer i, j
    
    read(501, *) 
    read(501, *)
    read(501, *) nlmps !numero de configuariones l1,m1,l2,m2
    read(501, *) 

    allocate(lmps(4, nlmps))! guarda las configuaraciones l1, m1, l2, m2, n
    do i=1, nlmps
       read(501, *) lmps(1, i), lmps(2, i), & 
            lmps(3, i), lmps(4, i)
    end do

    j=0
    do i=1, nlmps
       if((lmps(1, i)<0).or.(lmps(3, i)<0))then  !!!!condiciones que l1,l2, deben ser mayores que 0
          write(6, *) 'INPUT ERROR (l < 0)!', i  
          j=1
       end if
       if((abs(lmps(2, i))>lmps(1, i)).or.(abs(lmps(4, i))>lmps(3, i)))then !!!!condiciones que |m1,m2|>|l1,l2|
          write(6, *) 'INPUT ERROR (|m| > l)!', i
          j=1
       end if
    end do
    if(j.ne.0)stop
    
    ellmax=maxval(lmps)
    lamax=2*ellmax

   
  end subroutine get_lmps
!----------------------------------!

end module molecule
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
