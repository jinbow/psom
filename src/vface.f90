      subroutine vface(pf,dtime) 
!----------------------------------------------------                   
      USE header
!     modified for periodicew bc                                        
!     --------------------------                                        
!     compute the final vel (n+1 th step) at the cell faces             
!     applies the boundary conditions                                   
      implicit logical (a-z) 
      integer i,j,k,nexit,iflag 
      double precision dte,px,py,pz,dtime,                              &
     &     kaph1,flowin,flowout,const,                                  &
     &     pf(0:NI+1,0:NJ+1,0:NK+1)                                     
!     &     oldh(0:NI+1,0:NJ+1),                                        
!      double precision dte,px,py,pz,dtime                              
!                                                                       
      dte= dtime/EPS 
      kaph1= 1.d0 -kappah 
!      k=12                                                             
!      j=40                                                             
!      do i=0,NI+1                                                      
!         write(6,*) 'pcorr',i,pf(i,j,k)                                
!      end do                                                           
!      write(6,*) 'stopping in vface'                                   
!      stop                                                             
!                                                                       
!     uf compute here is used as the flux in convecn. In intpol it is   
!     used just at the outflow as the value of Uf at the prev time step.
!     In hbc and pbc, hnbc and pnbc, ufbce,ufbcw,vfbcn,vfbcs are used fo
      do 10 j=1,NJ 
!c         do 15 i=0,NI                                                 
         do 15 i=1,NI 
            do 16 k=1,NK 
!               if (i.eq.0) then                                        
!                  uf(i,j,k)= ufbcw(j,k)                                
!               else if (i.eq.NI) then                                  
!                     uf(i,j,k)= ufbce(j,k)                             
!               else                                                    
               px= (pf(i+1,j,k) -pf(i,j,k))*gqi(i,j,k,1) +0.25*         &
     &           (pf(i+1,j+1,k)+pf(i,j+1,k)                             &
     &              -pf(i+1,j-1,k)-pf(i,j-1,k))                         &
     &            *gqi(i,j,k,2)+0.25*(pf(i+1,j,k+1)+pf(i,j,k+1)-        &
     &            pf(i+1,j,k-1)-pf(i,j,k-1))*gqi3(i,j,k)                
               uf(i,j,k)= cxf(i,j,k) -dte*( px +sifc(i,j,k)) 
!               endif                                                   
   16       continue 
   15    continue 
   10 continue 
      do k=1,NK 
         do j=1,NJ 
            uf(0,j,k)=uf(NI,j,k) 
         end do 
      end do 
!                                                                       
!c     periodicity check                                                
!      iflag=0                                                          
!      do k=1,NK                                                        
!         do j=1,NJ                                                     
!            if (uf(NI,j,k).ne.uf(0,j,k)) then                          
!               write(6,*) uf(NI,j,k),uf(0,j,k),j,k                     
!               write(6,*) cxf(NI,j,k),cxf(0,j,k),                      
!     &              sifc(NI,j,k),sifc(0,j,k)                           
!               iflag=1                                                 
!            endif                                                      
!         end do                                                        
!      end do                                                           
!      if (iflag.eq.1) stop                                             
                                                                        
!     y-direction                                                       
!     -----------                                                       
      do 20 i=1,NI 
         do 25 j=0,NJ 
            do 26 k=1,NK 
               if (j.eq.0) then 
                  vf(i,j,k)= vfbcs(i,k) 
               else if (j.eq.NJ) then 
                  vf(i,j,k)= vfbcn(i,k) 
               else 
               py= (pf(i,j+1,k)                                         &
     &                 -pf(i,j,k))*gqj(i,j,k,2) +0.25*(pf(i+1,j+1,k)    &
     &          +pf(i+1,j,k)-pf(i-1,j+1,k)-pf(i-1,j,k))*gqj(i,j,k,1)    &
     &          +0.25*(pf(i,j+1,k+1)+pf(i,j,k+1)-pf(i,j+1,k-1)          &
     &          -pf(i,j,k-1))*gqj3(i,j,k)                               
               vf(i,j,k)= cyf(i,j,k) -dte*( py +sjfc(i,j,k)) 
               endif 
   26       continue 
   25    continue 
   20 continue 
!                                                                       
      do 30 j=1,NJ 
         do 30 i=1,NI 
            wf(i,j,0)= wfbcb(i,j) 
!     wf(i,j,NK)= wfbct(i,j)                                            
!     do 35 k=1,NK-1                                                    
!+    const= dtime/(J2d(i,j)*HDL)                                       
            do 35 k=1,NK 
               pz= (pf(i,j,k+1) -pf(i,j,k))*gqk(i,j,k,3) +              &
     &           0.25*(pf(i+1,j,k+1)                                    &
     &           +pf(i+1,j,k)-pf(i-1,j,k+1)-pf(i-1,j,k))*gqk(i,j,k,1)   &
     &           +0.25*(pf(i,j+1,k+1)+pf(i,j+1,k)-pf(i,j-1,k+1)         &
     &           -pf(i,j-1,k))*gqk(i,j,k,2)                             
               wf(i,j,k)= czf(i,j,k) -dte*(pz +skfc(i,j,k)) 
   35       continue 
   30 continue 
!                                                                       
!      j=17                                                             
!      i=15                                                             
!      do k=1,NK                                                        
!         write(6,*) uf(i,j,k),vf(i,j,k)                                
!      end do                                                           
!      write(6,*) 'czf'                                                 
!      do k=0,NK                                                        
!         write(6,*) czf(i,j,k),wf(i,j,k)                               
!      end do                                                           
!      stop                                                             
      return 
      END                                           
