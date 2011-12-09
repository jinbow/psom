      subroutine efill(nxm,nym,nzm,p) 
!     --------------------------------------------------------          
!     call efill(nx(m),ny(m),nz(m),p(loco(m)))                          
!     fills in the values at ficticious points outside the domain for   
!     grid levels other than 1                                          
      integer i,j,k,nxm,nym,nzm 
      double precision p(0:nxm+1,0:nym+1,0:nzm+1) 
!                                                                       
!                                                                       
      do 10 k=1,nzm 
         do 10 i=1,nxm 
            p(i,0,k)= p(i,1,k) 
            p(i,nym+1,k)= p(i,nym,k) 
   10 continue 
      do 20 k=1,nzm 
         do 20 j=1,nym 
            p(0,j,k)= p(1,j,k) 
            p(nxm+1,j,k)= p(nxm,j,k) 
   20 continue 
      do 30 j=1,nym 
         do 30 i=1,nxm 
            p(i,j,0)= p(i,j,1) 
!            p(i,j,nzm+1)= p(i,j,nzm)                                   
   30 continue 
!                                                                       
      do 35 j=0,nym+1 
         do 35 i=0,nxm+1 
            p(i,j,nzm+1)= -p(i,j,nzm) 
   35 continue 
!                                                                       
      do 40 k=1,nzm 
         p(0,0,k)= p(1,1,k) 
         p(nxm+1,0,k)= p(nxm,1,k) 
         p(0,nym+1,k)= p(1,nym,k) 
         p(nxm+1,nym+1,k)= p(nxm,nym,k) 
   40 continue 
      do 50 i=1,nxm 
         p(i,0,0)= p(i,1,1) 
!         p(i,0,nzm+1)= p(i,1,nzm)                                      
         p(i,nym+1,0)= p(i,nym,1) 
!         p(i,nym+1,nzm+1)= p(i,nym,nzm)                                
   50 continue 
      do 60 j=1,nym 
         p(0,j,0)= p(1,j,1) 
!         p(0,j,nzm+1)= p(1,j,nzm)                                      
         p(nxm+1,j,0)= p(nxm,j,1) 
!         p(nxm+1,j,nzm+1)= p(nxm,j,nzm)                                
   60 continue 
!                                                                       
      p(0,0,0)= p(1,1,1) 
!      p(0,0,nzm+1)= p(1,1,nzm)                                         
      p(0,nym+1,0)= p(1,nym,1) 
      p(nxm+1,0,0)= p(nxm,1,1) 
      p(nxm+1,nym+1,0)= p(nxm,nym,1) 
!      p(nxm+1,0,nzm+1)= p(nxm,1,nzm)                                   
!      p(0,nym+1,nzm+1)= p(1,nym,nzm)                                   
!      p(nxm+1,nym+1,nzm+1)= p(nxm,nym,nzm)                             
!                                                                       
      return 
      END                                           
