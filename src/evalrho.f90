subroutine evalrho(rhonew,n) 
!     ---------------------------------------------                     
  USE header
!  s is salinity, T is temp, rhonew is pot density
  implicit none
  double precision rhonew(0:NI+1,0:NJ+1,0:NK+1),potdens
  integer i,j,k,n 
  do k=0,NK+1 
     do j=0,NJ+1 
        do i=0,NI+1 
           rhonew(i,j,k)= potdens(s(i,j,k,n),T(i,j,k,n))
        end do
     end do
  end do
  return 
END subroutine evalrho
