subroutine inith 
  !     ------------------------------------------------                  
  USE header
  ! Initialize h using thermal wind balance
  ! Assume the level of no motion at z=-D  (bottom)
  ! Any density distribution is fine. 
  ! findzall has been called, so zf can be used. sigma not yet called, so wz 
  ! not usable
  implicit none
  integer i,j,k
  double precision fu(0:NJ),const,hsum,hmean,dz,drho

  const=DL/R0
  !fu is actually fu / g and is at cell faces
  do i=1,NI
     do j=1,NJ-1
        fu(j)= 0.d0
        do k=1,NK
           !fu is at cell faces,is actually fu/g
           dz= zf(i,j,k)-zf(i,j,k-1)
           drho= (rho(i,j+1,k)- rho(i,j,k))*vy(i,j)/LEN
           fu(j)= fu(j) + const*drho*dz
        end do
     end do
     fu(0)= fu(1)
     fu(NJ)= fu(NJ-1)

! at k=NK, fu = g*hy
     h(i,0)= 0.d0
     do j=1,NJ+1
        h(i,j)= h(i,j-1) - (LEN/vy(i,j))*fu(j-1)
     end do
  end do
  do j=0,NJ+1
     h(0,j)= h(NI,j)
     h(NI+1,j)= h(1,j)
  end do

  hsum=0.d0 
  do  i=1,NI 
     do  j=1,NJ 
        hsum= hsum +h(i,j)
     end do
  end do
  hmean= hsum/dble(NI*NJ) 
  h=h-hmean
  h=h/HL  
  return 
END subroutine inith
