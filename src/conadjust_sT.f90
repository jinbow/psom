subroutine conadjust(n) 
  !----------------------------------------------------                   
  USE header
  !     T(i,j,k) is set to 0, so Treal is T0                              
  !     Performs convective adjustment                                    
  implicit logical (a-z) 
  integer i,j,k,n,it,step 
  double precision rh,rhtop,zt,zb,dz,zinv 
  double precision dzupper 
  double precision sreal,Treal,rhreal,pbar 
  ! assume that evalrho has been called
  !                                                                       
  do j=1,NJ 
     do i=1,NI 
        do k=NK-1,1,-1 
           conv(i,j,k)= 0 
                                                                        
           rhtop= rho(i,j,k+1) 
           !     top                                                               
           dzupper = zf(i,j,k+1) -zf(i,j,k) 
                                                                        
           !     this level                                                        
           dz = zf(i,j,k) -zf(i,j,k-1) 
           rh =  rho(i,j,k) 
                                                                        
           if (rh.lt.rhtop) then 
!     mix                                                               
              conv(i,j,k)= 1 
              zinv= 1.d0/(dzupper +dz) 
              s(i,j,k,n)=(dzupper*s(i,j,k+1,n) +dz*s(i,j,k,n))*zinv 
              s(i,j,k+1,n)= s(i,j,k,n) 
              T(i,j,k,n)=(dzupper*T(i,j,k+1,n) +dz*T(i,j,k,n))*zinv 
              T(i,j,k+1,n)= T(i,j,k,n) 
              rho(i,j,k)= potdens(s(i,j,k,n),T(i,j,k,n)) 
              rho(i,j,k+1) = rho(i,j,k)
              do it=1,ntr 
                 Tr(it,i,j,k,n)=(dzupper*Tr(it,i,j,k+1,n)             &
                      +dz*Tr(it,i,j,k,n))*zinv                       
                 Tr(it,i,j,k+1,n)= Tr(it,i,j,k,n) 
              end do
           end if
                                                                        
        end do
     end do
  end do
                                                                        
  if (mod((step-1),100).eq. 0) then 
     con100= 0
  else 
     con100= con100 +conv
  end if
!     re-evaluate density (if needed)                                   
!      call evalrho(rho,n)                                              
                                                                        
!-      call sTbc_periodicew(n)                                         
  return 
END subroutine conadjust
