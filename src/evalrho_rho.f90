subroutine evalrho(rhonew,n) 
!     ---------------------------------------------                     
  USE header
!     s is the pot. density, so it is just copied into rho              
      implicit logical (a-z) 
      double precision rhonew(0:NI+1,0:NJ+1,0:NK+1) 
      integer i,j,k,n 
      do k=0,NK+1 
         do j=0,NJ+1 
            do i=0,NI+1 
               rhonew(i,j,k)= s(i,j,k,n) 
            end do 
         end do 
      end do 
      return 
      END                                           
