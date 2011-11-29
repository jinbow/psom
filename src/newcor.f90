subroutine newcor(dtime,n) 
!----------------------------------------------------                   
  USE header
!     modified for periodicew bcs                                       
!     Evaluate the Coriolis terms si,sj                                 
      implicit logical (a-z) 
      integer i,j,k,n 
      double precision fac2,ainv,dtime,con 
!     We are using the values at the ghost points.                      
      fac2= EPS*lambda 
      ainv= 1.d0/apr 
!                                                                       
      do 10 k=1,NK 
         do 20 j=1,NJ 
            do 30 i=1,NI 
!      do 10 k=0,NK+1                                                   
!         do 20 j=0,NJ+1                                                
!            do 30 i=0,NI+1                                             
!     since w has not changed we cannot correct fac*b*w                 
!-               si(i,j,k)= -ffc(i,j)*(cy(i,j,k) -v(i,j,k,n))           
!-               sj(i,j,k)= ffc(i,j)*(cx(i,j,k) -u(i,j,k,n))            
!     Doesn't work with si,sj=0. Oscillations produced.                 
               si(i,j,k)= 0.d0 
               sj(i,j,k)= 0.d0 
               sk(i,j,k)= 0.d0 
!c               si(i,j,k)= con*si(i,j,k)                               
!c               sj(i,j,k)= con*sj(i,j,k)                               
!c               sk(i,j,k)= -bbc(i,j)*u(i,j,k,n) -fac2*                 
!c     &              (u(i,j,k,n)*u(i,j,k,n) +                          
!c     &              v(i,j,k,n)*v(i,j,k,n))*ainv                       
   30       continue 
   20    continue 
   10 continue 
!                                                                       
      goto 999 
!     Since velocities outside the boundary are not correct we          
!     set si,sj,sk outside equal to values inside. (Feb 2)              
      do 170 j=1,NJ 
         do 180 i=1,NI 
            si(i,j,0)= si(i,j,1) 
            sj(i,j,0)= sj(i,j,1) 
            sk(i,j,0)= sk(i,j,1) 
!                                                                       
            si(i,j,NK+1)= si(i,j,NK) 
            sj(i,j,NK+1)= sj(i,j,NK) 
            sk(i,j,NK+1)= sk(i,j,NK) 
  180    continue 
  170 continue 
!                                                                       
      do 150 k=0,NK+1 
         do 160 i=1,NI 
            si(i,0,k)= si(i,1,k) 
            sj(i,0,k)= sj(i,1,k) 
            sk(i,0,k)= sk(i,1,k) 
!                                                                       
            si(i,NJ+1,k)= si(i,NJ,k) 
            sj(i,NJ+1,k)= sj(i,NJ,k) 
            sk(i,NJ+1,k)= sk(i,NJ,k) 
  160    continue 
  150 continue 
                                                                        
      do 130 k=0,NK+1 
         do 140 j=0,NJ+1 
            si(0,j,k)= si(NI,j,k) 
            sj(0,j,k)= sj(NI,j,k) 
            sk(0,j,k)= sk(NI,j,k) 
!                                                                       
            si(NI+1,j,k)= si(1,j,k) 
            sj(NI+1,j,k)= sj(1,j,k) 
            sk(NI+1,j,k)= sk(1,j,k) 
  140    continue 
  130 continue 
!                                                                       
  999 return 
      END                                           
