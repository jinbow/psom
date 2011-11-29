subroutine uvchy(dtime) 
!----------------------------------------------------                   
  USE header
  USE rpgrads
!     compute the tilde velocities at cell centers and store them in n=1
      implicit logical (a-z) 
      integer i,j,k 
      double precision dte,dtime,hxi,heta,hx,hy,kaph1,                  &
     &     Ub,Vb,Wb,d2,c1,c2,epsinv,temp1,temp2,wfk                     
!     operation:  u(n)= u(0) + dtime* f(u(m))                           
!                                                                       
      double precision pxi,peta,psig,px,py 
!                                                                       
      dte= dtime/EPS 
      epsinv= 1.d0/EPS 
      kaph1= 1.d0 - kappah 
!                                                                       
      do 10 j=1,NJ 
         do 20 i=1,NI 
            hxi= 0.5d0*( h(i+1,j)-h(i-1,j) ) 
            heta= 0.5d0*( h(i,j+1)-h(i,j-1) ) 
            hx= ux(i,j)*hxi +vx(i,j)*heta 
            hy= uy(i,j)*hxi +vy(i,j)*heta 
            do 30 k=1,NK 
                                                                        
               pxi= 0.5d0*(p(i+1,j,k)-p(i-1,j,k)) 
!     IF STATEMENTS ADDED jUN 8,2005 (doesn't make much diff)           
               if (j.eq.1) then 
                  peta= p(i,j+1,k)-p(i,j,k) 
               else if (j.eq.NJ) then 
                  peta= p(i,j,k)-p(i,j-1,k) 
               else 
                  peta= 0.5d0*(p(i,j+1,k)-p(i,j-1,k)) 
               endif 
               psig= 0.5d0*(p(i,j,k+1)-p(i,j,k-1)) 
               px= ux(i,j)*pxi +vx(i,j)*peta +wx(i,j,k)*psig 
               py= uy(i,j)*pxi +vy(i,j)*peta +wy(i,j,k)*psig 
                                                                        
!     cx and cy contain the convective terms.                           
               cx(i,j,k)=  cx(i,j,k) -dte*( gpr*(kappah*hx              &
     &              + kaph1*gradhn(i,j,1)) + si(i,j,k)                  &
     &              + qpr *px )                                         
               cy(i,j,k)=  cy(i,j,k) -dte*( gpr*(kappah*hy              &
     &              + kaph1*gradhn(i,j,2)) + sj(i,j,k)                  &
     &              + qpr *py )                                         
!               cz(i,j,k)=  cz(i,j,k)                                   
   30       continue 
            gradhn(i,j,1)= hx 
            gradhn(i,j,2)= hy 
   20    continue 
   10 continue 
!                                                                       
!     czf is already computed                                           
      do 11 j=1,NJ 
         do 21 i=1,NI 
            do 31 k=1,NK 
               wfk= 0.5d0*(czf(i,j,k) +czf(i,j,k-1)) 
               cz(i,j,k)= (wfk/Jac(i,j,k) -cx(i,j,k)*wx(i,j,k)          &
     &              -cy(i,j,k)*wy(i,j,k) )/(EPS*wz(i,j,k))              
   31       continue 
   21    continue 
   11 continue 
      return 
                                                                        
!     Boundaries (these are not probably needed)                        
      do j=1,NJ 
         do i=1,NI 
            cx(i,j,0)= cx(i,j,1) 
            cx(i,j,NK+1)= cx(i,j,NK) 
            cy(i,j,0)= cy(i,j,1) 
            cy(i,j,NK+1)= cy(i,j,NK) 
         end do 
      end do 
      do k=0,NK+1 
         do i=1,NI 
            cx(i,0,k)= cx(i,1,k) 
            cx(i,NJ+1,k)= cx(i,NJ,k) 
            cy(i,0,k)= cy(i,1,k) 
            cy(i,NJ+1,k)= cy(i,NJ,k) 
         end do 
      end do 
      do k=0,NK+1 
         do j=0,NJ+1 
            cx(0,j,k)= cx(NI,j,k) 
            cx(NI+1,j,k)= cx(1,j,k) 
            cy(0,j,k)= cy(NI,j,k) 
            cy(NI+1,j,k)= cy(1,j,k) 
         end do 
      end do 
!                                                                       
!c     Boundaries                                                       
!      i=NI+1                                                           
!      do 35 j=jex1,jex2                                                
!         hxi= h(i,j)-h(i-1,j)                                          
!         if (j.eq.jex1) then                                           
!            heta= h(i,j+1)-h(i,j)                                      
!         else if (j.eq.jex2) then                                      
!            heta= h(i,j)-h(i,j-1)                                      
!         else                                                          
!            heta= 0.5*(h(i,j+1)-h(i,j-1))                              
!         endif                                                         
!         hx= ux(i,j)*hxi +vx(i,j)*heta                                 
!         hy= uy(i,j)*hxi +vy(i,j)*heta                                 
!         do 36  k=1,NK                                                 
!            cx(i,j,k)=  cx(i,j,k) -dte*( gpr*(kappah*hx                
!     &           + kaph1*gradhn(i,j,1)) + si(i,j,k) )                  
!            cy(i,j,k)=  cy(i,j,k) -dte*( gpr*(kappah*hy                
!     &           + kaph1*gradhn(i,j,2)) + sj(i,j,k) )                  
!c            cz(i,j,k)=  cz(i,j,k)                                     
! 36      continue                                                      
!         gradhn(i,j,1)= hx                                             
!         gradhn(i,j,2)= hy                                             
! 35   continue                                                         
!                                                                       
!c     Need not fill the ghost points. The ghost points must not be used
!c     in newcor.                                                       
!c     *******************************                                  
!c                                                                      
!c      do k=1,NK                                                       
!c         do j=1,NJ                                                    
!c            if (ufbcw(j,k).lt.0.d0) then                              
!c               cx(0,j,k)= uwest(j,k)                                  
!c               cy(0,j,k)= cy(1,j,k)                                   
!c               cz(0,j,k)= cz(1,j,k)                                   
!c            else                                                      
!cc     incoming at west boundary                                       
!c               cx(0,j,k)= uwest(j,k)                                  
!c               cy(0,j,k)= vwest(j,k)                                  
!c               cz(0,j,k)= wwest(j,k)                                  
!c            end if                                                    
!c            if (ufbce(j,k).gt.0.d0) then                              
!c               cx(NI+1,j,k)= ueast(j,k)                               
!c               cy(NI+1,j,k)= cy(NI,j,k)                               
!c               cz(NI+1,j,k)= cz(NI,j,k)                               
!c            else                                                      
!cc     incoming at west boundary                                       
!c               cx(NI+1,j,k)= ueast(j,k)                               
!c               cy(NI+1,j,k)= veast(j,k)                               
!c               cz(NI+1,j,k)= weast(j,k)                               
!c            end if                                                    
!c         end do                                                       
!c      end do                                                          
!c      do k=1,NK                                                       
!c         do i=1,NI                                                    
!c            if (vfbcs(i,k).lt.0.d0) then                              
!c               cy(i,0,k)= vsouth(i,k)                                 
!c               cx(i,0,k)= cx(i,1,k)                                   
!c               cz(i,0,k)= cz(i,1,k)                                   
!c            else                                                      
!cc     incoming at south boundary                                      
!c               cy(i,0,k)= vsouth(i,k)                                 
!c               cx(i,0,k)= usouth(i,k)                                 
!c               cz(i,0,k)= wsouth(i,k)                                 
!c            end if                                                    
!c            if (vfbcn(i,k).gt.0.d0) then                              
!c               cy(i,NJ+1,k)= vnorth(i,k)                              
!c               cx(i,NJ+1,k)= cx(i,NJ,k)                               
!c               cz(i,NJ+1,k)= cz(i,NJ,k)                               
!c            else                                                      
!cc     incoming at north boundary                                      
!c               cy(i,NJ+1,k)= vnorth(i,k)                              
!c               cx(i,NJ+1,k)= unorth(i,k)                              
!c               cz(i,NJ+1,k)= wnorth(i,k)                              
!c            end if                                                    
!c         end do                                                       
!c      end do                                                          
!cc     vertical edges  (use east and west boundary arrays)             
!c      do k=1,NK                                                       
!c         cx(0,0,k)= uwest(0,k)                                        
!c         cx(NI+1,0,k)= ueast(0,k)                                     
!c         cx(0,NJ+1,k)= uwest(NJ+1,k)                                  
!c         cx(NI+1,NJ+1,k)= ueast(NJ+1,k)                               
!c         cy(0,0,k)= vwest(0,k)                                        
!c         cy(NI+1,0,k)= veast(0,k)                                     
!c         cy(0,NJ+1,k)= vwest(NJ+1,k)                                  
!c         cy(NI+1,NJ+1,k)= veast(NJ+1,k)                               
!c         cz(0,0,k)= wwest(0,k)                                        
!c         cz(NI+1,0,k)= weast(0,k)                                     
!c         cz(0,NJ+1,k)= wwest(NJ+1,k)                                  
!c         cz(NI+1,NJ+1,k)= weast(NJ+1,k)                               
!c      end do                                                          
!c      do j=0,NJ+1                                                     
!c         do i=0,NI+1                                                  
!c            cx(i,j,NK+1)= cx(i,j,NK)                                  
!c            cy(i,j,NK+1)= cy(i,j,NK)                                  
!c            cz(i,j,NK+1)= cz(i,j,NK)                                  
!c            cx(i,j,0)= cx(i,j,1)                                      
!c            cy(i,j,0)= cy(i,j,1)                                      
!c            cz(i,j,0)= -cz(i,j,1)                                     
!c         end do                                                       
!c      end do                                                          
!                                                                       
                                                                        
                                                                        
                                                                        
                                                                        
      return 
      END                                           
