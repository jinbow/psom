subroutine inith 
!     ------------------------------------------------                  
  USE header
!     Initialize h for the sb front from the density                    
!     assuming geostrophy and assuing the geos. velocity at z=h         
!     to be the same as that at z=0.                                    
!     This is for the case of the simple vertical front: for more compli
!     density distributions, this routine will have to be refined.      
!                                                                       
      implicit double precision (a-h,o-z) 
                                                                        
!      dh= 0.016d0 ! surf vel of 1.2 m/s                                
                                ! for dx=1km                            
      dh= 0.008d0 
!      dh= 0.005d0             !for 125m                                
!c      dh= 0.0005d0 largely barotropic front for 125 m reslution       
      do i=0,NI+1 
! perturb                                                               
!         dscl=dyfront*(0.5*(yc(NJ)+yc(NJ+1))-0.5*(yc(0)+yc(1)))        
         dscl=dyfront* sclwidth 
!==         dscl=(1.d0 + 1.d-4*ran3(iseed))*dscl                        
         rdscl=1.d0/dscl 
         dscl=dscl*1.d-2 
                                                                        
         hmean= 0.d0 
                                                                        
!     old profile : (1-exp y)/(1 +exp y) works better                   
                                                                        
         goto 111 
!     new method (- 1+ exp(-y)) or 1 - exp(-y)                          
         do j=0,NJ+1 
            yarg = yc(j) - yfront 
                                                                        
            if(yarg.lt.0.d0) then 
               thy = (-1. + dexp(-rdscl*dabs(yarg))) 
            else 
               thy = -(-1. + dexp(-rdscl*dabs(yarg))) 
            end if 
                                                                        
            h(i,j)= -dh*thy 
            hmean = hmean + h(i,j) 
         end do 
  111    continue 
                                                                        
!         goto 112                                                      
!     old method - Steve's curve (1-exp y)/(1 + exp y)                  
         do j=0,NJ+1 
!            yarg=yc(j)-  0.5*(yc(nj/2)+yc(nj/2+1))                     
!            yarg = yc(j) - yfront                                      
!            ex2y=dexp(-rdscl*dabs(yarg))                               
!            thy=(1.d0-ex2y)/(1.d0+ex2y)                                
!=                                                                      
            thy = tanh(tightness*(yc(j)-yfront)*PI) 
!=same curve as in stprofile                                            
!            write(6,*) j,thy,s(10,j,24,0)                              
                                                                        
!           sechy2=1.d0-thy*thy                                         
            h(i,j)= -dh*thy 
            hmean = hmean + h(i,j) 
         end do 
  112    continue 
                                                                        
         hmean = hmean/dble(NJ+2) 
         do j=0,NJ+1 
            h(i,j)= (h(i,j) - hmean) 
!            if (i.eq.1) write(6,*) h(i,j)                              
         end do 
      end do 
                                                                        
                                                                        
      return 
      END                                           
