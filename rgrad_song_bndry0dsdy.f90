subroutine rpevalgrad(n) 
!     ---------------------------------------------                     
  USE header
  USE rpgrads
!     modified for periodicew bc                                        
!     ---------------------------                                       
!     MODIFIED FOR THE SHELFBREAK (SLOPE in Y) TOPOGRAPHY               
      implicit logical (a-z) 
      integer i,j,k,n 
      double precision const,dz,Jx(NI,NK),Jy(0:NJ,0:NK),Jy2(NJ,0:NK),   &
     &     dzeta(NJ,0:NK),dzsig(NJ,0:NK),Jsum                           
                                                                        
!     Evaluate Jacobian at the edge-centers                             
!     multiply press and press gradient by const                        
!     const= G*gpr/P1                                                   
      const= 10.d0*gpr/P1 
!     const = (approx)  1.d-3                                           
                                                                        
      call evalrho(rho,n) 
!      do k=1,NK                                                        
!         do j=1,NJ                                                     
!            do i=1,NI                                                  
!               rho(i,j,k)= s(i,j,k,n)                                  
!            end do                                                     
!         end do                                                        
!      end do                                                           
                                                                        
!     x-direction                                                       
!     -----------                                                       
!     Assumes sigma surfaces are horizontal.                            
                                                                        
      do j=1,NJ 
!     periodic-ew boundaries                                            
         do k=1,NK 
            Jx(1,k)= 0.5*(rho(2,j,k) -rho(NI,j,k))*                     &
     &           (zf(1,j,k) -zf(1,j,k-1))*DL                            
            do i=2,NI-1 
               dz= (zf(i,j,k) -zf(i,j,k-1))*DL 
               Jx(i,k)= 0.5*(rho(i+1,j,k) -rho(i-1,j,k))*dz 
            end do 
            Jx(NI,k)= 0.5*(rho(1,j,k) -rho(NI-1,j,k))*                  &
     &           (zf(NI,j,k)-zf(NI,j,k-1))*DL                           
         end do 
         do i=1,NI 
                                                                        
            drpx(i,j,NK)= 0.5*const*ux(i,j)*Jx(i,NK) 
            do k=NK-1,1,-1 
               drpx(i,j,k)= drpx(i,j,k+1) +                             &
     &              0.5d0*(Jx(i,k+1) +Jx(i,k)) *const*ux(i,j)           
            end do 
         end do 
!         if (j.eq.40) then                                             
!            i = 25                                                     
!            do k=NK-1,1,-1                                             
!               write(6,*) 'rgrad',k,drpx(i,j,k+1),Jx(i,k),Jx(i,k+1),   
!     &              0.5d0*(Jx(i,k+1) +Jx(i,k)) *const*ux(i,j)          
!            end do                                                     
!         end if                                                        
      end do 
!                                                                       
!                                                                       
!     grpifc  (at cell faces)                                           
!     -----------------------                                           
      do j=1,NJ 
!     periodic-ew boundaries                                            
         do k=1,NK 
            do i=1,NI-1 
               dz= 0.5d0*(zf(i,j,k) +zf(i+1,j,k) -                      &
     &              (zf(i,j,k-1) +zf(i+1,j,k-1)) )*DL                   
               Jx(i,k)= (rho(i+1,j,k) -rho(i,j,k))*dz*gi(i,j,k,1) 
            end do 
            dz= 0.5d0*(zf(1,j,k) +zf(NI,j,k) -                          &
     &           (zf(1,j,k-1)+zf(NI,j,k-1)))*DL                         
            Jx(NI,k)= (rho(1,j,k) -rho(NI,j,k))*dz*gi(NI,j,k,1) 
         end do 
                                                                        
         do i=1,NI 
            grpifc(i,j,NK)= 0.5*const*Jx(i,NK) 
            do k=NK-1,1,-1 
               grpifc(i,j,k)= grpifc(i,j,k+1) +                         &
     &              0.5d0*const*(Jx(i,k+1) +Jx(i,k))                    
            end do 
         end do 
         do k=1,NK 
            grpifc(0,j,k)= grpifc(NI,j,k) 
         end do 
      end do 
                                                                        
                                                                        
!     y-direction                                                       
!     -----------                                                       
!     sloping sigma surfaces                                            
!     solid boundaries                                                  
                                                                        
      do i=1,NI 
         do j=1,NJ-1 
            k=0 
            Jy(j,k)= rho(i,j+1,k+1) - rho(i,j,k+1) 
            dzeta(j,k)= (zc(i,j+1,k+1) -zc(i,j,k+1) )*DL 
            do k=1,NK-1 
               Jy(j,k)= 0.5*(rho(i,j+1,k)+rho(i,j+1,k+1) -              &
     &              (rho(i,j,k) +rho(i,j,k+1) ))                        
               dzeta(j,k)= 0.5*(zc(i,j+1,k) +zc(i,j+1,k+1) -            &
     &              (zc(i,j,k) +zc(i,j,k+1) ) )*DL                      
            end do 
            k= NK 
            Jy(j,k)= rho(i,j+1,k) -rho(i,j,k) 
            dzeta(j,k)= (zc(i,j+1,k) -zc(i,j,k))*DL 
         end do 
                                                                        
         do j=1,NJ-1 
            k=0 
            Jy2(j,0)= 0.5*(rho(i,j,k+2) +rho(i,j+1,k+2) -               &
     &              (rho(i,j,k+1) + rho(i,j+1,k+1)) )                   
            dzsig(j,1)= 0.5*(zc(i,j,2) +zc(i,j+1,2) -                   &
     &           (zc(i,j,1) + zc(i,j+1,1) ))*DL                         
            do k=1,NK-1 
               Jy2(j,k)= 0.5*(rho(i,j,k+1) +rho(i,j+1,k+1) -            &
     &              (rho(i,j,k) + rho(i,j+1,k)) )                       
               dzsig(j,k)= 0.5d0*(zc(i,j,k+1) +zc(i,j+1,k+1)            &
     &              - (zc(i,j,k) +zc(i,j+1,k)) )*DL                     
            end do 
            Jy2(j,NK)= 0.5*(rho(i,j,NK) +rho(i,j+1,NK) -                &
     &           (rho(i,j,NK-1) + rho(i,j+1,NK-1) ) )                   
            dzsig(j,NK)= 0.5*(zc(i,j,NK) +zc(i,j+1,NK) -                &
     &           (zc(i,j,NK-1) + zc(i,j+1,NK-1) ) )*DL                  
         end do 
                                                                        
         do k=0,NK 
            do j=1,NJ-1 
               Jy(j,k)= Jy(j,k)*dzsig(j,k) 
               Jy2(j,k)= Jy2(j,k)*dzeta(j,k) 
               Jy(j,k)= Jy(j,k) -Jy2(j,k) 
            end do 
         end do 
! Solid boundaries                                                      
         do k=0,NK 
            Jy(0,k)= 0.0 
            Jy(NJ,k)= 0.0 
         end do 
                                                                        
         do j=1,NJ-1 
            k= NK 
            Jsum= Jy(j,k)*0.5 
            grpjfc(i,j,NK)= Jsum*const*gj(i,j,k,2) 
            do k=NK-1,1,-1 
               Jsum= Jsum + Jy(j,k) 
               grpjfc(i,j,k)= Jsum*const*gj(i,j,k,2) 
            end do 
         end do 
         do k=1,NK 
            grpjfc(i,0,k)= 0.d0 
            grpjfc(i,NJ,k)= 0.d0 
         end do 
         do j=1,NJ 
            k= NK 
            drpy(i,j,k)= const*vy(i,j)*0.5                              &
     &           *( Jy(j,NK) + Jy(j-1,NK))*0.5                          
            do k=NK-1,1,-1 
               drpy(i,j,k)= drpy(i,j,k+1) +const*vy(i,j)*0.5            &
     &           *( Jy(j,k) + Jy(j-1,k))                                
            end do 
         end do 
      end do 
                                                                        
      return 
                                                                        
!      goto 202                                                         
      do k=1,NK 
         do j=1,NJ 
            do i=1,NI 
               if( drpy(i,j,k).ne. 0.0) then 
                  write(6,*) 'drpy' 
                  stop 
               end if 
               if (drpx(i,j,k).ne. 0.0) then 
                  write(6,*) 'drpx' 
                  stop 
               end if 
            end do 
         end do 
         do j=0,NJ 
            do i=1,NI 
               if (grpjfc(i,j,k).ne.0.d0) then 
                  write(6,*) 'grpjfc' 
                  stop 
               end if 
            end do 
         end do 
         do j=1,NJ 
            do i=0,NI 
               if (grpifc(i,j,k).ne.0.d0) then 
                  write(6,*) 'grpifc' 
                  stop 
               end if 
            end do 
         end do 
      end do 
                                                                        
                                                                        
      return 
!     ***********                                                       
                                                                        
!     grpjfc  (at cell faces)                                           
!     -----------------------                                           
      do i=1,NI 
!     Solid boundaries at j=0 and NJ                                    
!     grpjfc(0,k)= grpjfc(NJ,k)= 0                                      
         do k=1,NK 
            do j=1,NJ-1 
!               dz= 0.5d0*(zf(i,j,k) +zf(i,j+1,k) -                     
!     &              (zf(i,j,k-1) +zf(i,j+1,k-1)) )*DL                  
               Jy(j,k)= (rho(i,j+1,k) -rho(i,j,k)) 
               dzeta(j,k)= (zc(i,j+1,k) -zc(i,j,k))*DL 
!                                                                       
            end do 
         end do 
                                                                        
         do j=1,NJ-1 
            k=1 
            Jy2(j,k)= 0.5d0*(rho(i,j+1,k+1)+rho(i,j,k+1)                &
     &           -(rho(i,j+1,k) +rho(i,j,k)))*dzeta(j,k)                
            dzsig(j,k)= 0.5d0*(zc(i,j+1,k+1)+zc(i,j,k+1)                &
     &           -(zc(i,j+1,k) +zc(i,j,k)))*DL                          
            do k=2,NK-1 
               Jy2(j,k)= 0.25d0*(rho(i,j+1,k+1)+rho(i,j,k+1)            &
     &              -(rho(i,j+1,k-1) +rho(i,j,k-1)))*dzeta(j,k)         
               dzsig(j,k)= 0.25d0*(zc(i,j+1,k+1)+zc(i,j,k+1)            &
     &              -(zc(i,j+1,k-1) +zc(i,j,k-1)))*DL                   
            end do 
            k=NK 
            Jy2(j,k)= 0.5d0*(rho(i,j+1,k)+rho(i,j,k)                    &
     &           -(rho(i,j+1,k-1) +rho(i,j,k-1)))*dzeta(j,k)            
            dzsig(j,k)= 0.5d0*(zc(i,j+1,k)+zc(i,j,k)                    &
     &           -(zc(i,j+1,k-1) +zc(i,j,k-1)))*DL                      
         end do 
         do k=1,NK 
            do j=1,NJ-1 
               Jy(j,k)= (Jy(j,k)*dzsig(j,k) -Jy2(j,k))*gj(i,j,k,2) 
!               if (i.eq.1) write(60,*) Jy(j,k),Jy2(j,k)                
            end do 
         end do 
                                                                        
         do j=1,NJ-1 
            grpjfc(i,j,NK)= 0.5*const*Jy(j,NK) 
            do k=NK-1,1,-1 
               grpjfc(i,j,k)= grpjfc(i,j,k+1) +                         &
     &              0.5d0*const*(Jy(j,k+1) +Jy(j,k))                    
            end do 
         end do 
         do k=1,NK 
            grpjfc(i,0,k)= 0.d0 
            grpjfc(i,NJ,k)= 0.d0 
         end do 
      end do 
                                                                        
                                                                        
      return 
!     *************                                                     
  202 do k=1,NK 
         write(200,*) 'k= ',k 
         do j=1,NJ 
            i=8 
            write(200,11) j,drpy(i,j,k),grpjfc(i,j,k) 
   11       format(I4,2(2x,E16.8)) 
         end do 
      end do 
                                                                        
      do k=1,NK 
         write(300,*) 'k= ',k 
         do i=1,NI 
            do j=1,NJ 
               if (drpx(i,j,k).ne.0.d0) then 
                  write(320,*) i,j,k,drpx(i,j,k) 
               endif 
               if (grpifc(i,j,k).ne.0.d0) then 
                  write(330,*) i,j,k,grpifc(i,j,k) 
               endif 
            end do 
            j=10 
            write(300,*) i,drpx(i,j,k),grpifc(i,j,k) 
         end do 
      end do 
                                                                        
      write(6,*) 'stopping in rpevalgrad_spline' 
      stop 
!     **************                                                    
      return 
      END                                           
