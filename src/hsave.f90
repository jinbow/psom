subroutine hsave 
!----------------------------------------------------                   
  USE header
!     saves the hydrostatic pressure gradients at the initial time step 
      implicit logical (a-z) 
      integer i,j,k 
      double precision kap1,hxi,heta 
!                                                                       
      kap1= 1.d0 -kappah 
!     hxn,hyn,hzn need to be evaluated at the boundaries too, because   
!     they are used in mghfill.                                         
!                                                                       
!     x-direction                                                       
!     -----------                                                       
      do 10 j=1,NJ 
         do 10 i=0,NI 
            hxi= h(i+1,j) -h(i,j) 
            heta= 0.25*(h(i+1,j+1)+h(i,j+1)-h(i+1,j-1)-h(i,j-1)) 
            do 15 k=1,NK 
               hxn(i,j,k)= hxi*gi(i,j,k,1) +heta*gi(i,j,k,2) 
   15       continue 
   10 continue 
!                                                                       
!     y-direction                                                       
!     -----------                                                       
      do 20 j=0,NJ 
         do 20 i=1,NI 
            heta= h(i,j+1) -h(i,j) 
            hxi= 0.25*(h(i+1,j+1)+h(i+1,j)-h(i-1,j+1)-h(i-1,j)) 
            do 25 k=1,NK 
               hyn(i,j,k)= heta*gj(i,j,k,2) + hxi*gj(i,j,k,1) 
!               if ((i.eq.1).and.(k.eq.1).and.((j.eq.0).or.(j.eq.NJ)))  
!     &             write(6,*) 'hsave',i,j,k,hyn(i,j,k),heta,gj(i,j,k,2)
   25       continue 
   20 continue 
!                                                                       
!                                                                       
!     Now save grad h as is done in uvchy.f for later time steps.       
      do 40 j=1,NJ 
         do 40 i=1,NI 
            hxi= 0.5d0*(h(i+1,j)-h(i-1,j)) 
            heta= 0.5d0*(h(i,j+1)-h(i,j-1)) 
            gradhn(i,j,1)= ux(i,j)*hxi +vx(i,j)*heta 
            gradhn(i,j,2)= uy(i,j)*hxi +vy(i,j)*heta 
   40 continue 
!                                                                       
      return 
      END                                           
