subroutine coriolis(n) 
!----------------------------------------------------                   
  USE header
  USE rpgrads
!     Evaluate the Coriolis terms si,sj                                 
      implicit logical (a-z) 
      integer i,j,k,n 
      double precision pxi,peta,psig,px,py,fac,fac2,ainv,betainv 
                                                                        
      fac= EPS*delta 
      fac2= EPS*lambda 
      ainv= 1.d0/apr 
      betainv= 1.0/beta 
                                                                        
!     We are using the values at the ghost points.                      
      do 10 k=1,NK 
         do 20 j=1,NJ 
            do 30 i=1,NI 
!=               pxi= 0.5d0*(p(i+1,j,k)-p(i-1,j,k))                     
!=               peta= 0.5d0*(p(i,j+1,k)-p(i,j-1,k))                    
!=               psig= 0.5d0*(p(i,j,k+1)-p(i,j,k-1))                    
!=               px= ux(i,j)*pxi +vx(i,j)*peta +wx(i,j,k)*psig          
!=               py= uy(i,j)*pxi +vy(i,j)*peta +wy(i,j,k)*psig          
                                                                        
               si(i,j,k)= -ffc(i,j)*v(i,j,k,n)                          &
     &              +fac*bbc(i,j)*w(i,j,k,n)                            &
     &              + drpx(i,j,k)                                       &
     &              - EPS*uvis(i,j,k)                                   
!     &             + qpr* px  now put in uvchy (so that si can comprise
               sj(i,j,k)= ffc(i,j)*u(i,j,k,n)                           &
     &              + drpy(i,j,k)                                       &
     &              - EPS* vvis(i,j,k)                                  
!     &             + qpr*py  this is now put in uvchy                  
               sk(i,j,k)= fnhhy*(-bbc(i,j)*u(i,j,k,n) -fac2*            &
     &              (u(i,j,k,n)*u(i,j,k,n) +                            &
     &              v(i,j,k,n)*v(i,j,k,n))*ainv)                        &
     &              - betainv*wvis(i,j,k)                               
   30       continue 
   20    continue 
   10 continue 
!                                                                       
!      do k=1,NK                                                        
!         write(100,*) 'k=', k                                          
!         do j=1,NJ                                                     
!            write(100,*) 'j=', j                                       
!      k=NK                                                             
!      j=18                                                             
!      write(100,*) 'j=',j                                              
!            do i= 1,NI                                                 
!               write(100,*) si(i,j,k),uvis(i,j,k),drpx(i,j,k)          
!            end do                                                     
!         end do                                                        
!      end do                                                           
!      stop                                                             
!                                                                       
      return 
      END                                           
