subroutine newsrc 
!----------------------------------------------------                   
  uSE header
!     We interpolate the source terms onto the cell faces               
      implicit logical (a-z) 
      integer i,j,k 
      double precision be2,uxi,uyi,vxj,vyj,wzsk,ainv,                   &
     &     wxsk,wysk,fac2,Jack,con                                      
!     We are using the values at the ghost points.                      
      fac2= EPS*lambda 
      ainv= 1.d0/apr 
      be2= beta*EPS*EPS 
      con= qpr/(1.d0-qpr) 
!                                                                       
      do k=1,NK 
         do j=1,NJ 
            do i=0,NI 
               sifc(i,j,k) = 0.0 
            end do 
         end do 
      end do 
                                                                        
      do k=1,NK 
         do j=0,NJ 
            do i=1,NI 
               sjfc(i,j,k)= 0.0 
            end do 
         end do 
      end do 
                                                                        
      do k=0,NK 
         do j=1,NJ 
            do i=1,NI 
               skfc(i,j,k)= 0.0 
            end do 
         end do 
      end do 
                                                                        
      return 
                                                                        
                                                                        
                                                                        
      do 10 i=0,NI 
         do 20 j=1,NJ 
            uxi= 0.5*(ux(i+1,j)+ux(i,j)) 
            uyi= 0.5*(uy(i+1,j)+uy(i,j)) 
            do 30 k=1,NK 
!     vint= 0.5*(cy(i+1,j,k)-v(i+1,j,k) +cy(i,j,k)-v(i,j,k))            
!     uint= 0.5*(cx(i+1,j,k)-u(i+1,j,k) +cx(i,j,k)-u(i,j,k))            
!     wint= 0.5*(cz(i+1,j,k) +cz(i,j,k))                                
!     sifc(i,j,k)= (uxi*(-ffi(i,j)*vint +fac*bbi(i,j)*wint) +           
!     &              uyi*ffi(i,j)*uint)*Jifc(i,j)                       
!-               sifc(i,j,k)= 0.5d0*(uxi*(si(i+1,j,k)+si(i,j,k))        
!-     &              +  uyi*(sj(i+1,j,k)+sj(i,j,k)) )*Jifc(i,j,k)      
               sifc(i,j,k)= con*sifc(i,j,k) 
   30       continue 
   20    continue 
   10 continue 
!                                                                       
      do 40 i=1,NI 
         do 50 j=0,NJ 
            vxj= 0.5*(vx(i,j+1)+vx(i,j)) 
            vyj= 0.5*(vy(i,j+1)+vy(i,j)) 
            do 60 k=1,NK 
!     vint= 0.5*(cy(i,j+1,k)-v(i,j+1,k) +cy(i,j,k)-v(i,j,k))            
!     uint= 0.5*(cx(i,j+1,k)-u(i,j+1,k) +cx(i,j,k)-u(i,j,k))            
!     wint= 0.5*(cz(i,j+1,k) +cz(i,j,k))                                
!     sjfc(i,j,k)= (-vxj*ffj(i,j)*vint +                                
!     &              vyj*ffj(i,j)*uint )*Jjfc(i,j)                      
!-               sjfc(i,j,k)= 0.5d0*( vxj*(si(i,j+1,k)+si(i,j,k))       
!-     &              + vyj*(sj(i,j+1,k)+sj(i,j,k)) )*Jjfc(i,j,k)       
               sjfc(i,j,k)= con*sjfc(i,j,k) 
   60       continue 
   50    continue 
   40 continue 
!                                                                       
!--      goto 99                                                        
!     Since the boundary conditions on cx,cy are the same as on u,v     
!     sifc,sjfc are zero at the boundaries                              
!     comment the following out because of periodic-ew boundaries       
!-      do k=1,NK                                                       
!-         do j=1,NJ                                                    
!-            do i=0,NI,NI                                              
!-               sifc(i,j,k)= 0.d0                                      
!-            end do                                                    
!-         end do                                                       
!-      end do                                                          
      do k=1,NK 
         do j=0,NJ,NJ 
            do i=1,NI 
               sjfc(i,j,k)= 0.d0 
            end do 
         end do 
      end do 
   99 continue 
!                                                                       
!     Notice that skfc goes from i=0,NI+1 and j=0,NJ+1. This is         
!     because in mgpfill we wish to use skfc at i=0,NI+1,j=0,NJ+1       
!     at k=0 and NK levels to fill p along the horizontal edges         
!                                                                       
      do 70 i=1,NI 
         do 80 j=1,NJ 
            do 90 k=1,NK-1 
               wxsk= 0.25d0*(wx(i,j,k+1) +wx(i,j,k))*(si(i,j,k+1)+      &
     &              si(i,j,k) )                                         
               wysk= 0.25d0*(wy(i,j,k+1) +wy(i,j,k))*(sj(i,j,k+1)+      &
     &              sj(i,j,k) )                                         
               wzsk= 0.25d0*(wz(i,j,k+1) +wz(i,j,k))*(sk(i,j,k+1)+      &
     &              sk(i,j,k) )                                         
               Jack= 0.5d0*(Jac(i,j,k+1) + Jac(i,j,k)) 
               skfc(i,j,k)= ( be2*wzsk +wxsk +wysk )*Jack 
   90       continue 
            k=0 
            wxsk= wx(i,j,k+1)*si(i,j,k+1) 
            wysk= wy(i,j,k+1)*sj(i,j,k+1) 
            wzsk= wzk(i,j,k)*sk(i,j,k+1) 
            Jack= Jac(i,j,k+1) 
            skfc(i,j,k)= ( be2*wzsk +wxsk +wysk )*Jack 
            k= NK 
            wxsk= 0.5d0*(wx(i,j,k+1) +wx(i,j,k))*si(i,j,k) 
            wysk= 0.5d0*(wy(i,j,k+1) +wy(i,j,k))*sj(i,j,k) 
            wzsk= wzk(i,j,k)*sk(i,j,k) 
            Jack= Jac(i,j,k) 
            skfc(i,j,k)= ( be2*wzsk +wxsk +wysk )*Jack 
                                                                        
! july 11, 2001;  compute skfc more accurately at top boundary          
!c            skfc(i,j,0)= skfc(i,j,1)                                  
!c            skfc(i,j,NK)= skfc(i,j,NK-1)                              
   80    continue 
   70 continue 
!                                                                       
      return 
      END                                           
