subroutine energy(step) 
!---------------------------------------------------                    
  USE header
      implicit logical (a-z) 
!                                                                       
      integer i,j,k,step 
      double precision const,enkin,enin,enout,enbal 
!                                                                       
      const= EPS*EPS*delta*delta 
      enkin= 0.d0 
      do k=1,NK 
         do j=1,NJ 
            do i=1,NI 
               enkin= Jac(i,j,k)*(u(i,j,k,0)*u(i,j,k,0) +v(i,j,k,0)*    &
     &              v(i,j,k,0) + const*w(i,j,k,0)*w(i,j,k,0))           &
     &              +enkin                                              
            end do 
         end do 
      end do 
!                                                                       
!c     flow in                                                          
!      enin= 0.d0                                                       
!      do k=1,NK                                                        
!         do i=ient1,ient2                                              
!            enin= enin + vfbcs(i,k)*dtf*(uent(i,k)*uent(i,k) +         
!     &           vent(i,k)*vent(i,k) )                                 
!                                                                       
!         end do                                                        
!      end do                                                           
!c                                                                      
!c     flow out                                                         
!      enout = 0.d0                                                     
!      do k=1,NK                                                        
!         do j=jex1,jex2                                                
!            enout= enout +ufbce(j,k)*dtf*(u(NI,j,k,0)*u(NI,j,k,0) +    
!     &           v(NI,j,k,0)*v(NI,j,k,0) +const*w(NI,j,k,0)*           
!     &           w(NI,j,k,0) )                                         
!         end do                                                        
!      end do                                                           
!c                                                                      
!      enbal= enin + enout + enkin                                      
!                                                                       
!      write(901,11) step,enin,enout,enkin,enbal                        
      write(91,*) step,enkin 
   11 format(I5,5x,4(e10.4,4x)) 
      return 
      END                                           
