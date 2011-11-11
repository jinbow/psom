subroutine calcn2 
  !     ---------------------------                                       
  USE header
  !     Calculate N2  (units: per second)                                 
                                                                        
      implicit none 
      integer i,j,k 
      double precision DLinv,rdz,depth,dz 
      double precision avezf(NK),aveN2(NK) 
                                                                        
!     evalrho was called from rpevalgrad                                
      DLinv= 1.d0/DL 
      do k=1,NK 
         do j=1,NJ 
            do i=1,NI 
                                                                        
               if (k.eq.NK) then 
                  rdz= (rho(i,j,k) -rho(i,j,k-1))*wz(i,j,k)*DLinv 
               else if (k.eq.1) then 
                  rdz= (rho(i,j,k+1) -rho(i,j,k))*wz(i,j,k)*DLinv 
               else 
                  rdz= 0.5*(rho(i,j,k+1) -rho(i,j,k-1))*wz(i,j,k)*DLinv 
               end if 
                                                                        
               freqN2(i,j,k)=(-gpr*10.d0/R0)*rdz 
                                                                        
            end do 
         end do 
      end do 
                                                                        
      return 
      END                                           
