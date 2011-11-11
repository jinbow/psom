subroutine cfdiv(maxdiv) 
  !     ------------------                                                
  USE header
  !     checks divergence of cxf,cyf,czf                                  
  !     czf does not contain  J*wt and is not exactly the contravariant ve
      integer i,j,k,imax,jmax,kmax 
      double precision cxfdx,cyfdy,czfdz,div,maxdiv 
!     edd is a scaling factor that makes this divergence                
!     the same as the residual from the pressure-Poisson equation.      
!     This is a check on the code.                                      
!                                                                       
      maxdiv= 0.d0 
      do 10 k=1,NK 
         do 20 j=1,NJ 
            do 30 i=1,NI 
               cxfdx= (cxf(i,j,k)-cxf(i-1,j,k)) 
               cyfdy= (cyf(i,j,k)-cyf(i,j-1,k)) 
               czfdz= (czf(i,j,k)-czf(i,j,k-1)) 
                                                                        
               div= dabs(cxfdx+ cyfdy + czfdz) 
               if (div.gt.maxdiv) then 
                  maxdiv=div 
                  imax=i 
                  jmax=j 
                  kmax=k 
               end if 
   30       continue 
   20    continue 
   10 continue 
!      write(6,*) 'in facediv, i,j,k',imax,jmax,kmax                    
      return 
      END                                           
