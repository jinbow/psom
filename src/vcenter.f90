      subroutine vcenter(pf,dtime,n) 
!----------------------------------------------------                   
      USE header
!     compute the final vel (n+1 th step) at the cell centers           
      implicit logical (a-z) 
      integer i,j,k,n 
      double precision dtbeta,pxi,peta,psig,c1,c2,dte,                  &
     &     dtime,Ub,Vb,Wb,d2,Wr,Wtemp(0:NI+1,0:NJ+1),                   &
     &     px,py,pz,fac,pf(0:NI+1,0:NJ+1,0:NK+1),wfk                    
!     operation:  u(n)= u(0) + dtime* f(u(m))                           
      double precision z,seval,bu0(NK),cu0(NK),du0(NK),                 &
     &     bv0(NK),cv0(NK),dv0(NK),bw0(NK),cw0(NK),dw0(NK),             &
     &     dep(NK),uvert(NK),vvert(NK),wvert(NK)                        
!                                                                       
      dte= dtime/EPS 
      dtbeta= dtime*beta 
!                                                                       
      do 10 j=1,NJ 
         do 20 i=1,NI 
            do 30 k=1,NK 
               pxi= 0.5d0*(pf(i+1,j,k)-pf(i-1,j,k)) 
               peta= 0.5d0*(pf(i,j+1,k)-pf(i,j-1,k)) 
               psig= 0.5d0*(pf(i,j,k+1)-pf(i,j,k-1)) 
               px= ux(i,j)*pxi +vx(i,j)*peta +wx(i,j,k)*psig 
               py= uy(i,j)*pxi +vy(i,j)*peta +wy(i,j,k)*psig 
               pz= wz(i,j,k)*psig 
               u(i,j,k,n)=  cx(i,j,k) -dte*(qpr*px +si(i,j,k) ) 
               v(i,j,k,n)=  cy(i,j,k) -dte*(qpr*py +sj(i,j,k) ) 
!c                                                                      
!=               w(i,j,k,n)=  cz(i,j,k) -dtbeta*( pz +sk(i,j,k) )       
   30       continue 
   20    continue 
   10 continue 
!                                                                       
!=      return                                                          
!=                                                                      
!     vface is already called - compute w from wf                       
      do 11 j=1,NJ 
         do 21 i=1,NI 
            do 31 k=1,NK 
               wfk= 0.5d0*(wf(i,j,k) +wf(i,j,k-1)) 
               w(i,j,k,n)= (wfk/Jac(i,j,k) -u(i,j,k,n)*wx(i,j,k)        &
     &              -v(i,j,k,n)*wy(i,j,k) )/(EPS*wz(i,j,k))             
   31       continue 
   21    continue 
   11 continue 
!                                                                       
!==      call velbc_periodicew(n)                                       
!      j=17                                                             
!      i=15                                                             
!      do k=1,NK                                                        
!         write(6,*) u(i,j,k,n),v(i,j,k,n),w(i,j,k,n)                   
!      end do                                                           
                                                                        
      return 
                                                                        
                                                                        
!     Boundaries                                                        
!     k=0,NK                                                            
      k=0 
      do 70 i=1,NI 
         do 75 j=1,NJ 
            u(i,j,k,n)= u(i,j,1,n) 
            v(i,j,k,n)= v(i,j,1,n) 
            w(i,j,k,n)= -w(i,j,1,n) 
   75    continue 
   70 continue 
                                                                        
                                                                        
!-      do k=0,NK                                                       
!-         do i=1,NI                                                    
!-            u(i,0,k,n)= u(i,1,k,n)                                    
!-            v(i,0,k,n)= v(i,1,k,n)                                    
!-            w(i,0,k,n)= w(i,1,k,n)                                    
!-            u(i,NJ+1,k,n)= u(i,NJ,k,n)                                
!-            v(i,NJ+1,k,n)= v(i,NJ,k,n)                                
!-            w(i,NJ+1,k,n)= w(i,NJ,k,n)                                
!-         end do                                                       
!-      end do                                                          
!                                                                       
!     FIT SPLINES to EXTRAPOLATE                                        
!     ----------------------------                                      
      do i=1,NI 
         do k=1,NK 
            dep(k)= zc(i,1,k) 
            uvert(k)= u(i,1,k,n) 
            vvert(k)= v(i,1,k,n) 
            wvert(k)= w(i,1,k,n) 
         end do 
         call spline(NK,dep,uvert,bu0,cu0,du0) 
         call spline(NK,dep,vvert,bv0,cv0,dv0) 
         call spline(NK,dep,vvert,bw0,cw0,dw0) 
!                                                                       
         j=0 
         z= zc(i,j,k) 
         u(i,j,k,n)= seval(NK,z,dep,uvert,bu0,cu0,du0) 
         v(i,j,k,n)= seval(NK,z,dep,vvert,bv0,cv0,dv0) 
         w(i,j,k,n)= seval(NK,z,dep,wvert,bw0,cw0,dw0) 
!                                                                       
         do k=1,NK 
            dep(k)= zc(i,NJ,k) 
            uvert(k)= u(i,NJ,k,n) 
            vvert(k)= v(i,NJ,k,n) 
            wvert(k)= w(i,NJ,k,n) 
         end do 
         call spline(NK,dep,uvert,bu0,cu0,du0) 
         call spline(NK,dep,vvert,bv0,cv0,dv0) 
         call spline(NK,dep,vvert,bw0,cw0,dw0) 
!                                                                       
         j=NJ+1 
         z= zc(i,j,k) 
         u(i,j,k,n)= seval(NK,z,dep,uvert,bu0,cu0,du0) 
         v(i,j,k,n)= seval(NK,z,dep,vvert,bv0,cv0,dv0) 
         w(i,j,k,n)= seval(NK,z,dep,wvert,bw0,cw0,dw0) 
      end do 
!                                                                       
                                                                        
                                                                        
!                                                                       
!     vertical edges                                                    
      do k=1,NK 
         u(0,0,k,n)= u(1,0,k,n) 
         v(0,0,k,n)= v(1,0,k,n) 
         w(0,0,k,n)= w(1,0,k,n) 
         u(0,NJ+1,k,n)= u(1,NJ+1,k,n) 
         v(0,NJ+1,k,n)= v(1,NJ+1,k,n) 
         w(0,NJ+1,k,n)= w(1,NJ+1,k,n) 
         u(NI+1,NJ+1,k,n)= u(NI,NJ+1,k,n) 
         v(NI+1,NJ+1,k,n)= v(NI,NJ+1,k,n) 
         w(NI+1,NJ+1,k,n)= w(NI,NJ+1,k,n) 
         u(NI+1,0,k,n)= u(NI,0,k,n) 
         v(NI+1,0,k,n)= v(NI,0,k,n) 
         w(NI+1,0,k,n)= w(NI,0,k,n) 
      end do 
!                                                                       
!     edges                                                             
      do 71 j=0,NJ+1 
         u(0,j,0,n)=  0.d0 
         v(0,j,0,n)=  0.d0 
         w(0,j,0,n)=  0.d0 
         u(NI+1,j,0,n)=  0.d0 
         v(NI+1,j,0,n)=  0.d0 
         w(NI+1,j,0,n)=  0.d0 
   71 continue 
      do 72 i=1,NI 
         u(i,0,0,n)= 0.d0 
         v(i,0,0,n)= 0.d0 
         w(i,0,0,n)= 0.d0 
         u(i,NJ+1,0,n)= 0.d0 
         v(i,NJ+1,0,n)= 0.d0 
         w(i,NJ+1,0,n)= 0.d0 
   72 continue 
                                                                        
      do 78 j=1,NJ 
         do 78 i=1,NI 
!            Wtemp(i,j)= wfbct(i,j)                                     
            Wtemp(i,j)= wf(i,j,NK) 
   78 continue 
      do 79 j=1,NJ 
         Wtemp(0,j)= 2.d0*wf(1,j,NK) -wf(2,j,NK) 
         Wtemp(NI+1,j)= 2.d0*wf(NI,j,NK) -wf(NI-1,j,NK) 
   79 continue 
      do 81 i=1,NI 
         Wtemp(i,0)= 2.d0*wf(i,1,NK) -wf(i,2,NK) 
         Wtemp(i,NJ+1)= 2.d0*wf(i,NJ,NK) -wf(i,NJ-1,NK) 
   81 continue 
      Wtemp(0,0)= 0.5*(Wtemp(1,0)+Wtemp(0,1)) 
      Wtemp(NI+1,0)= 0.5*(Wtemp(NI,0)+Wtemp(NI+1,1)) 
      Wtemp(0,NJ+1)= 0.5*(Wtemp(0,NJ)+Wtemp(1,NJ+1)) 
      Wtemp(NI+1,NJ+1)= 0.5*(Wtemp(NI,NJ+1)+Wtemp(NI+1,NJ)) 
!                                                                       
      k=NK+1 
      do 80 i=0,NI+1 
         do 85 j=0,NJ+1 
            u(i,j,k,n)= u(i,j,NK,n) 
            v(i,j,k,n)= u(i,j,NK,n) 
            w(i,j,k,n)= u(i,j,NK,n) 
   85    continue 
   80 continue 
!                                                                       
      k=NK+1 
      do 89 j=1,NJ 
         w(0,j,k,n)= w(1,j,k,n) 
         w(NI+1,j,k,n)=w(NI,j,k,n) 
   89 continue 
      do 91 i=1,NI 
         w(i,0,k,n)= w(i,1,k,n) 
         w(i,NJ+1,k,n)=w(i,NJ,k,n) 
   91 continue 
      w(0,0,k,n)= 0.5*(w(1,0,k,n)+w(0,1,k,n)) 
      w(NI+1,0,k,n)= 0.5*(w(NI,0,k,n)+w(NI+1,1,k,n)) 
      w(0,NJ+1,k,n)= 0.5*(w(0,NJ,k,n)+w(1,NJ+1,k,n)) 
      w(NI+1,NJ+1,k,n)= 0.5*(w(NI,NJ+1,k,n)+w(NI+1,NJ,k,n)) 
!                                                                       
                                                                        
!     periodicew                                                        
      do 35 k=0,NK+1 
         do 36  j=0,NJ+1 
            u(0,j,k,n)= u(NI,j,k,n) 
            v(0,j,k,n)= v(NI,j,k,n) 
            w(0,j,k,n)= w(NI,j,k,n) 
!                                                                       
            u(NI+1,j,k,n)= u(1,j,k,n) 
            v(NI+1,j,k,n)= v(1,j,k,n) 
            w(NI+1,j,k,n)= w(1,j,k,n) 
   36    continue 
   35 continue 
                                                                        
      return 
      END                                           
