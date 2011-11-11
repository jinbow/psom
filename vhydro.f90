subroutine vhydro(dtime) 
!     -----------------------------------------------                   
  USE header
  USE rpgrads
!     modified for periodicew                                           
!     computes, u,v,w -tilde and writes these over cxf,cyf,czf.         
!     czf is unchanged.                                                 
      implicit logical (a-z) 
      integer i,j,k 
!      double precision v1(0:NI+1,0:NJ+1,0:NK+1),                       
!     &     v2(0:NI+1,0:NJ+1,0:NK+1),v3(0:NI+1,0:NJ+1,0:NK+1),          
!     &     v4(0:NI+1,0:NJ+1,0:NK+1),v5(0:NI+1,0:NJ+1,0:NK+1)           
      double precision hxi,heta,dte,dtime,gradh,kaph1,temp,pz,          &
     &     sumuf(0:NI,NJ),sumvf(NI,0:NJ),hx,hy,Jack,wfbct,ufdu,vfdv     
      dte= dtime/EPS 
      kaph1= 1.d0 -kappah 
!                                                                       
      do 10 j=1,NJ 
         do 20 i=1,NI 
            hxi= h(i+1,j)- h(i,j) 
            heta= 0.25*(h(i+1,j+1) +h(i,j+1) -h(i+1,j-1)                &
     &            -h(i,j-1))                                            
!=            sumuf(i,j)= 0.d0                                          
            do 30 k=1,NK 
               hx= gi(i,j,k,1)*hxi +gi(i,j,k,2)*heta 
               gradh= gpr*(kappah*hx + kaph1*hxn(i,j,k)) 
               cxf(i,j,k)= cxf(i,j,k) -dte*(gradh +sifc(i,j,k)) 
!=               sumuf(i,j)= sumuf(i,j) + kappah*cxf(i,j,k)             
!=     &              + kaph1*uf(i,j,k)                                 
               hxn(i,j,k)= hx 
   30       continue 
   20    continue 
!                                                                       
!=         sumuf(0,j)= sumuf(NI,j)                                      
         do k=1,NK 
            cxf(0,j,k)= cxf(NI,j,k) 
            hxn(0,j,k)= hxn(NI,j,k) 
         end do 
   10 continue 
!     ************added on Mar 21,94 ***********                        
!      do 15 j=1,NJ                                                     
!         sumuf(0,j)= 0.d0                                              
!         do  16 k=1,NK                                                 
!            cxf(0,j,k)= ufbcw(j,k)                                     
!            sumuf(0,j)= sumuf(0,j) + kappah*cxf(0,j,k)                 
!     &           + kaph1*uf(0,j,k)                                     
! 16      continue                                                      
!         sumuf(NI,j)= 0.d0                                             
!         do 17 k=1,NK                                                  
!            cxf(NI,j,k)= ufbce(j,k)                                    
!            sumuf(NI,j)= sumuf(NI,j) + kappah*cxf(NI,j,k)              
!     &           + kaph1*uf(NI,j,k)                                    
! 17      continue                                                      
! 15   continue                                                         
!     **************************************                            
!                                                                       
      do 40 i=1,NI 
         do 50 j=0,NJ 
!            if (i.eq.16) write(500,*) 'j = ',j                         
            hxi= 0.25*(h(i+1,j+1) +h(i+1,j) -h(i-1,j+1) -h(i-1,j)) 
            heta= h(i,j+1) -h(i,j) 
!=            sumvf(i,j)= 0.d0                                          
            do 60 k=1,NK 
               hy= gj(i,j,k,1)*hxi +gj(i,j,k,2)*heta 
               gradh= gpr*(kappah*hy + kaph1*hyn(i,j,k)) 
               cyf(i,j,k)= cyf(i,j,k) -dte*(gradh +sjfc(i,j,k)) 
!=               sumvf(i,j)= sumvf(i,j) + kappah*cyf(i,j,k)             
!=     &              + kaph1*vf(i,j,k)                                 
               hyn(i,j,k)= hy 
   60       continue 
   50    continue 
   40 continue 
!                                                                       
!     ************added on Mar 21,94 ***********                        
!     June 6,2005                                                       
!     There's a problem with the bc on p at j=0,NJ+1 (solid boundaries),
!     which traces back to the bc on pcorr in mgpfill. If cyf is set    
!     equal to vfbcs at this stage, then dp/d(eta) =0 at the N,S        
!     boundaries. This make p discountinuous and creates oscillations   
!     at the boundary.                                                  
!     So we are going to stop setting cyf to vfbcs in vhydro. The       
!     bc on p, will ensure that the final vf= vfbcs                     
!      goto 399  - this just doesn't work, so I'm giving in to setting c
      do 25 i=1,NI 
!=         sumvf(i,NJ)= 0.d0                                            
         do  26 k=1,NK 
            cyf(i,NJ,k)= vfbcn(i,k) 
!=            sumvf(i,NJ)= sumvf(i,NJ) + kappah*cyf(i,NJ,k)             
!=     &           + kaph1*vf(i,NJ,k)                                   
   26    continue 
!=         sumvf(i,0)= 0.d0                                             
         do 27 k=1,NK 
            cyf(i,0,k)= vfbcs(i,k) 
!=            sumvf(i,0)= sumvf(i,0) + kappah*cyf(i,0,k)                
!=     &           + kaph1*vf(i,0,k)                                    
   27    continue 
   25 continue 
!     **** added July 1, 2002                                           
  399 goto 401 
!     **************************************                            
!     Evaluate b.c. on wf(at k=NK) (excludes wt). Bc on wf(k=0)= 0      
!     This would be the Bc if the residual of the h-eqn were zero       
!     For the edges (i=0,NI+1,j=0,NJ+1) hdt is computed in hsolve.      
      temp= dtime/HDL 
      do 310 j=1,NJ 
         do 310 i=1,NI 
            wfbct= -(sumuf(i,j) -sumuf(i-1,j) +sumvf(i,j)               &
     &           -sumvf(i,j-1)) +wfbcb(i,j)                             
            hdt(i,j)= wfbct/J2d(i,j) 
!     hdt is actually (dh'/dt')*(H/D)                                   
!     coorect h(i,j), so that the residual does not accumulate          
!==            h(i,j)= oldh(i,j) + temp*hdt(i,j)                        
  310 continue 
!     added 12/99                                                       
!==      call hfill(dtime,h)                                            
!                                                                       
      do 330 i=1,NI 
         hdt(i,0)= hdt(i,1) 
         hdt(i,NJ+1)= hdt(i,NJ) 
  330 continue 
      do 320 j=0,NJ+1 
         hdt(0,j)= hdt(NI,j) 
         hdt(NI+1,j)= hdt(1,j) 
  320 continue 
!                                                                       
!-      call findztop (no longer use this - now move the whole grid)    
!===                                                                    
!     no longer call this here - remain explicit in wz,wx,wy,Jac -      
!     these are called only before the next time step in momentum.      
!===      call findzall                                                 
!===      call sigma                                                    
!                                                                       
!     Compute czf                                                       
!     If HY, then pz=0, skfc=0                                          
!                                                                       
  401 do j=1,NJ 
         do i=1,NI 
! czf must be 0 only in intpol, not here.  czf(i,j,0)= wfbcb(i,j)       
! czf will become 0 if HY. Else, it may be non-zero at this stage,but 0 
!            if (i.eq.16)          write(200,*) 'j = ',j                
            do k=0,NK 
               pz= (p(i,j,k+1) -p(i,j,k))*gqk(i,j,k,3) +                &
     &              0.25*(p(i+1,j,k+1)                                  &
     &              +p(i+1,j,k)-p(i-1,j,k+1)-p(i-1,j,k))*gqk(i,j,k,1)   &
     &              +0.25*(p(i,j+1,k+1)+p(i,j+1,k)-p(i,j-1,k+1)         &
     &              -p(i,j-1,k))*gqk(i,j,k,2)                           
!==                                                                     
!               if ((k.eq.NK).and.(j.eq.40)) then                       
!                  pz= (p(i,j,k+1) -p(i,j,k))*gqk(i,j,k,3) +            
!     &                 0.25*(p(i+1,j,k+1)+p(i+1,j,k)                   
!     &                 -p(i-1,j,k+1)-p(i-1,j,k))*gqk(i,j,k,1)          
!     &                 +0.25*(p(i,j+1,k+1)+p(i,j+1,k)-p(i,j-1,k+1)     
!     &                 -p(i,j-1,k))*gqk(i,j,k,2)                       
!                                                                       
!                  write(6,*)                                           
!     &              i,czf(i,j,k),pz,skfc(i,j,k),                       
!     &              czf(i,j,k)-dte*(pz+skfc(i,j,k))                    
!               end if                                                  
!==                                                                     
               czf(i,j,k)= czf(i,j,k) -dte*(pz +skfc(i,j,k)) 
!               ufdu= (cxf(i,j,k)-cxf(i-1,j,k))                         
!               vfdv= (cyf(i,j,k)-cyf(i,j-1,k))                         
!               czf(i,j,k)= czf(i,j,k-1) - ufdu -vfdv                   
!               if (i.eq.16)  then                                      
!                  write(200,*) pz,skfc(i,j,k),czf(i,j,k)               
!               endif                                                   
            end do 
         end do 
      end do 
!      write(6,*) 'writing done in vhydro'                              
                                                                        
!     call uvchy after czf is computed                                  
      call uvchy(dtime) 
                                                                        
!8**************                                                        
!      do k=1,NK                                                        
!         do j=1,NJ                                                     
!            do i=1,NI                                                  
!               v1(i,j,k)= cxf(i,j,k) - cxf(i-1,j,k)                    
!               v2(i,j,k)= cyf(i,j,k)- cyf(i,j-1,k)                     
!               v3(i,j,k)= czf(i,j,k)                                   
!               v4(i,j,k)= hxn(i,j,k)                                   
!               v5(i,j,k)= hyn(i,j,k)- hyn(i,j-1,k)                     
!            end do                                                     
!         end do                                                        
!      end do                                                           
!                                                                       
!      call outarray(v1,v2,v3,v4,v5)                                    
!      write(6,*) 'stopping in vhydro'                                  
!      stop                                                             
!********************                                                   
      return 
  101 continue 
                                                                        
!     U-tilde and V-tilde are computed from the equation in vhydro.f    
!     W-tilde is found from interpolation (and stored in czf)           
      do 210 k=1,NK-1 
         do 220 j=1,NJ 
            do 230 i=1,NI 
               Jack= 0.5d0*(Jac(i,j,k) +Jac(i,j,k+1)) 
               czf(i,j,k)= 0.5*Jack*( wx(i,j,k)*cx(i,j,k) +             &
     &              wx(i,j,k+1)*cx(i,j,k+1) +                           &
     &              wy(i,j,k)*cy(i,j,k) +wy(i,j,k+1)*cy(i,j,k+1) +      &
     &              EPS*(wz(i,j,k)*cz(i,j,k)                            &
     &              +wz(i,j,k+1)*cz(i,j,k+1)) )                         
  230       continue 
  220    continue 
  210 continue 
!                                                                       
!     boundaries                                                        
      k= NK 
      do 240 j=1,NJ 
         do 250 i=1,NI 
            czf(i,j,0)= wfbcb(i,j) 
!*            czf(i,j,NK)= wfbct(i,j)                                   
!     Use cx,cy,cz from k=NK, but use mean of wx,wy,wz,Jac in case      
!     grid is sretched.                                                 
               Jack= 0.5d0*(Jac(i,j,k) +Jac(i,j,k+1)) 
               czf(i,j,k)= 0.5*Jack*( (wx(i,j,k) +                      &
     &              wx(i,j,k+1))*cx(i,j,k) +                            &
     &              (wy(i,j,k) +wy(i,j,k+1))*cy(i,j,k) +                &
     &              EPS*((wz(i,j,k) +wz(i,j,k+1))*cz(i,j,k) ) )         
!-            czf(i,j,NK)= Jac(i,j,NK)*( wx(i,j,NK)*cx(i,j,NK) +        
!-     &           wy(i,j,NK)*cy(i,j,NK) +                              
!-     &           EPS*wz(i,j,NK)*cz(i,j,NK) )                          
  250    continue 
  240 continue 
!                                                                       
      return 
      END                                           
