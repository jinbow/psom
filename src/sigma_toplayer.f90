subroutine sigma 
!     ------------------------                                          
  USE header
!     ONLY FOR THE TOPMOST GRID CELL                                    
!     modified a for periodicew bc                                      
!     This subroutine updates all those quantities that are a function  
!     of time and  (time and depth). It is called  every time step.     
!     The surface bc  wfbct is updated here.                            
      implicit double precision (a-h,o-z) 
      integer i,j,k 
      double precision hpd,hpdinv,hu,hv,hx,hy,hxpdx,hypdy,z,temp,       &
     &     be2,wxk,wyk,d2,dnk,dnkm1,sig,zpd,                            &
     &     g13(0:NI+1,0:NJ+1,0:NK+1),g23(0:NI+1,0:NJ+1,0:NK+1)          
!     Ddx and Ddy are known from init                                   
!                                                                       
      dnkm1= dble(NK-1) 
      dnk= dble(NK) 
!                                                                       
      be2= beta*EPS*EPS 
!                                                                       
      do 10 j=0,NJ+1 
         do 10 i=0,NI+1 
!     All these variables are functions of time                         
            hpd= h(i,j)*HDL +dztop 
            hpdinv= 1.d0/hpd 
            if (i.eq.0) then 
               hu= HDL*(h(i+1,j)-h(NI-1,j)) 
            else if  (i.eq.NI+1) then 
               hu= HDL*(h(2,j)-h(i-1,j)) 
            else 
               hu= 0.5d0*HDL*(h(i+1,j)-h(i-1,j)) 
            endif 
            if (j.eq.0) then 
               hv= HDL*(h(i,j+1)-h(i,j)) 
            else if (j.eq.NJ+1) then 
               hv= HDL*(h(i,j)-h(i,j-1)) 
            else 
               hv= 0.5d0*HDL*(h(i,j+1)-h(i,j-1)) 
            endif 
            hx= hu*ux(i,j) + hv*vx(i,j) 
            hy= hu*uy(i,j) + hv*vy(i,j) 
!     wz is not a function of depth when the sigma lines are equally spa
!     Hence wz is wz(i,j,time). For a stretched grid wz would be w(i,j,k
!     Then Jac which is now Jac(i,j,time) would  become  Jac(i,j,k,time)
!                                                                       
!            write(100,*) J2d(i,j)                                      
!            write(200,*) Jac(i,j)                                      
!                                                                       
            do 20 k=NK,NK+1 
               wz(i,j,k)= hpdinv 
               Jac(i,j,k)= J2d(i,j)/wz(i,j,k) 
!     All these variables are functions of time and depth               
!     z is  dimensionaless and is equal to  z*/DL                       
!     now hdt computed in vhydro already contains HDL                   
               sig= dfloat(k)-0.5d0 
               z= (sig -dnkm1)*hpd -dztop 
               wx(i,j,k)= (dnkm1-sig)*hx*hpdinv 
               wy(i,j,k)= (dnkm1-sig)*hy*hpdinv 
               g13(i,j,k)= ux(i,j)*wx(i,j,k) +uy(i,j)*wy(i,j,k) 
               g23(i,j,k)= vx(i,j)*wx(i,j,k) +vy(i,j)*wy(i,j,k) 
!                                                                       
!               g33(i,j,k)= wx(i,j,k)*wx(i,j,k) +wy(i,j,k)*wy(i,j,k) +  
!     &              be2*wz(i,j)*wz(i,j)                                
   20       continue 
!     wt is defined at cell faces                                       
            do k=NK-1,NK 
               sig= dfloat(k) 
               wt(i,j,k)= (dnkm1-sig)*hdt(i,j)*hpdinv 
            end do 
            wzk(i,j,NK)= hpdinv 
            wzk(i,j,NK-1)= 0.5*(wz(i,j,NK) + wz(i,j,NK-1)) 
!                                                                       
!     We evaluate gk at i=0,NI+1,j=0,NJ+1 only because these are used   
!     at k=0,NK to fill the horizontal edges in mgpfill or hnbc.        
            do 15 k=NK-1,NK 
               sig= dfloat(k) 
               wxk= (dnkm1-sig)*hx*hpdinv 
               wyk= (dnkm1-sig)*hy*hpdinv 
               gqk(i,j,k,1)= qpr*Jac(i,j,k)*(ux(i,j)*wxk +uy(i,j)*wyk) 
               gqk(i,j,k,2)= qpr*Jac(i,j,k)*(vx(i,j)*wxk +vy(i,j)*wyk) 
               gqk(i,j,k,3)= Jac(i,j,k)*(qpr*(wxk*wxk +wyk*wyk) +       &
     &              be2*wz(i,j,k)*wz(i,j,k))                            
   15       continue 
   10 continue 
!      do j=0,NJ+1                                                      
!         do i=0,NI+1                                                   
!            write(120,*) gk(i,j,0,1), gk(i,j,3,1)                      
!            write(130,*) gk(i,j,0,2), gk(i,j,3,2)                      
!            write(140,*) gk(i,j,0,3), gk(i,j,3,3)                      
!         end do                                                        
!      end do                                                           
!      stop                                                             
!                                                                       
!                                                                       
      k= NK 
!     ======                                                            
      do 21 i=0,NI 
         do 31 j=1,NJ 
            Jifc(i,j,k)= 0.5d0*(Jac(i,j,k)+ Jac(i+1,j,k)) 
            gi(i,j,k,1)= 0.5d0*(g11(i,j) +g11(i+1,j))*Jifc(i,j,k) 
            gi(i,j,k,2)= 0.5d0*(g12(i,j) +g12(i+1,j))*Jifc(i,j,k) 
            gqi(i,j,k,1)= qpr*gi(i,j,k,1) 
            gqi(i,j,k,2)= qpr*gi(i,j,k,2) 
   31    continue 
   21 continue 
                                                                        
      do 22 i=1,NI 
         do 32 j=0,NJ 
            Jjfc(i,j,k)= 0.5d0*(Jac(i,j,k)+ Jac(i,j+1,k)) 
            gj(i,j,k,1)= 0.5d0*(g12(i,j) +g12(i,j+1))*Jjfc(i,j,k) 
            gj(i,j,k,2)= 0.5d0*(g22(i,j) +g22(i,j+1))*Jjfc(i,j,k) 
            gqj(i,j,k,1)= qpr*gj(i,j,k,1) 
            gqj(i,j,k,2)= qpr*gj(i,j,k,2) 
   32    continue 
   22 continue 
!                                                                       
      do 30 j=1,NJ 
         do 30 i=0,NI 
               gi3(i,j,k)= 0.5d0*(g13(i,j,k) +g13(i+1,j,k))*Jifc(i,j,k) 
               gqi3(i,j,k)= qpr*gi3(i,j,k) 
   30 continue 
      do 40 j=0,NJ 
         do 40 i=1,NI 
               gj3(i,j,k)= 0.5d0*(g23(i,j,k) +g23(i,j+1,k))*Jjfc(i,j,k) 
               gqj3(i,j,k)= qpr*gj3(i,j,k) 
   40 continue 
!                                                                       
      return 
      END                                           
