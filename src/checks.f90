subroutine checks 
  !     ------------------                                                
  USE header
  integer i,j,k 
  !     Checks and writes a few things                                    
                                                                        
                                                                        
      do j=1,NJ 
         do i=0,NI 
            if (dabs(D(i+1,j)-D(i,j)).gt.1.d-4) then 
               write(6,*) i,j,D(i+1,j),D(i,j) 
               write(6,*) 'need to modify rpevalgrad (rgradients) for   &
     &              slope in x direcn: grpifc needs to be modified'     
               write(6,*) 'need to modify velbc_periodicew' 
               stop 
            end if 
         end do 
      end do 
                                                                        
      do j=0,NJ+1 
         do i=0,NI+1 
            if ((uy(i,j).ne.0.d0).or.(vx(i,j).ne.0.d0)) then 
               write(6,*) 'Need to modify biharmonic for curvi grid' 
               stop 
            end if 
         end do 
      end do 
                                                                        
      if (rect.eqv.(.false.)) then 
!        need to modify grpifc,grpjfc to contain cross diff terms       
         write(6,*) 'modify grpifc,grpjfc, stop in checks' 
         stop 
      end if 
!                                                                       
!      open(unit=88,file='umax.out',status='new')                       
!      open(unit=89,file='umin.out',status='new')                       
!      close(88 89)                                                     
!      open(unit=79,file='vtimeser.out',status='new')                   
!      open(unit=78,file='utimeser.out',status='new')                   
!      close(78 79)                                                     
                                                                        
      return 
      END                                           
