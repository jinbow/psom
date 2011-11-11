subroutine findzall 
  !     -------------------------------------------------------           
  USE header
  !     finds the value of z (non-dim by DL) given the value of sigma.    
  !     ht and dep are non-dim by DL.                                     
  !     At every (i,j), the column is divided into NK equal-depth cells.  
  integer i,j,k 
  double precision sigma,dep,ht,epm1,dnkm1,em1inv,hpd,         &
       dnkm1inv,xfac                                                
  !    pfac is the stretching in z. higher pfac gives more points near surf.

  pfac= 2.0d0 
  !      pfac= 5.d0  !c=NK32,dztop=0.5m                                   
  !==      pfac= 4.d0  !c=NK32,dztop=0.5m USED FOR SEVERAL MLI runs       
  !=      pfac=3.d0  !c used with NK=32 for first set of runs             
  dnkm1= dfloat(NK-1) 
  dnkm1inv= 1.d0/dnkm1 
  epm1= dexp(pfac) -1.d0 
  epm1inv= 1.d0/(dexp(pfac) -1.d0) 
  do j=0,NJ+1 
     do i=0,NI+1 
        !     In the surface layer                                              
        hpd= h(i,j)*HDL +dztop 
        do k=NK,NK+1 
           sigma= dble(k)-0.5 
           zc(i,j,k)= (sigma -dnkm1)*hpd -dztop 
        end do
        do k=NK-1,NK+1 
           sigma= dble(k) 
           zf(i,j,k)= (sigma -dnkm1)*hpd -dztop 
        end do
        !                                                                       
        !     Below the surface layer                                           
        do k=0,NK-1 
           sigma= dble(k)-0.5 
           xfac = (dnkm1 -sigma)*dnkm1inv 
           zc(i,j,k)= (dexp(pfac*xfac)-1.d0)*epm1inv*               &
                (D(i,j)+dztop) -dztop                               
        end do
        do k=-1,NK-2 
           sigma= dble(k) 
           xfac = (dnkm1 -sigma)*dnkm1inv 
           zf(i,j,k)= (dexp(pfac*xfac)-1.d0)*epm1inv*               &
                (D(i,j)+dztop) -dztop                               
        end do
     end do
  end do
  
  return 
end subroutine findzall
