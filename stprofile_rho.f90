subroutine stprofile 
  !     --------------------                                              
  USE header
!     initializes s as pot.density with a single vertical profile       
!     Across-front use the TANH profile with tightness as a measure of f
!     spread                                                            
!     Larger the factor, tighter the front and larger b_xx.             
!     tightness=4, gives the needed b_xx to make b_xx H/ (2f^2) < 1.    
      implicit logical (a-z) 
        integer  i,j,k,it,iseed,npm1 
      double precision rhovert(np),Tvert(np),dep(np),depoff(np) 
      double precision bs1(np),cs1(np),ds1(np),bT1(np),cT1(np),dT1(np), &
     &     z,seval,zmm,sbkgrnd,z1,z2                                    
      double precision slfac,dum,dscl,rdscl,yarg,ex2y,thy,ran3,         &
     &     perturb,slfacnew,zoffset,dz,bfsqbkgrnd                       
!     tightness = 10 represents a very tight front, =1 loose front      
!=      parameter (zoffset= 200.)  also used zoffset=50 ! if zoffset=0 m
!=       parameter (zoffset= 200.d0, tightness=0.03)                    
!=       parameter (zoffset= -10.d0)                                    
!       parameter (zoffset= -20.d0) 
       parameter (zoffset= 0.d0) 
      data dep /-1100.,-800.,-600.,-450.,-350.,-305.,-275.,             &
     &    -250.,-240.,-230.,-220.,-210.,-200.,-190.,-180.,              &
     &    -170.,-160.,-150.,-140.,-130.,-120.,-110.,-100.,              &
     &    -90.,-80.,-70.,-60.,-50.,-40.,-30.,-20.,-10., 0./             
      data rhovert /   1025.9728, 1025.9728, 1025.9689,                 &
     &     1025.9612, 1025.9496, 1025.9342, 1025.9142,                  &
     &     1025.8949, 1025.8748, 1025.8556, 1025.8324,                  &
     &     1025.8062, 1025.7754, 1025.7414, 1025.7029,                  &
     &     1025.6643, 1025.6258, 1025.5872, 1025.5487,                  &
     &     1025.5102, 1025.4716, 1025.4331, 1025.3945,                  &
     &     1025.3560, 1025.3174, 1025.2789, 1025.2481,                  &
     &     1025.2249, 1025.2095, 1025.2018, 1025.2018,                  &
     &     1025.2018, 1025.2018 /                                       
!     The z values must be in STRICTLY INCREASING ORDER and hence should
!     be written from ocean bottom up.                                  
                                                                        
                          !for larger domain                            
      tightness= 0.03d0 
      bfsqbkgrnd = 1.d-6 
!      bfsqbkgrnd = 1.d-5                                               
      write(6,*) 'ML background stratific N2 = ', bfsqbkgrnd 
      do k=1,np 
         depoff(k)= dep(k)- zoffset 
      end do 
                                                                        
      call spline (np,depoff,rhovert,bs1,cs1,ds1) 
                                                                        
!     sclwidth is the scale width of the channel (i.e. orig width = 48km
!     yfront is the distance of the front from the coast                
!     dyfront is the scale width of the front as a ratio of the domain w
      sclwidth = 48.0 
      yfront = 0.5*(yc(NJ+1) +yc(0)) 
!=      yfront = (yc(NJ+1) +yc(0))/3.d0                                 
      dyfront = 1./24.0 
!     dyfront = 1./48.0 works well for this set up                      
!     dscl=dyfront* sclwidth  this sets the width of the front (exp cuto
                                                                        
!     ADD a PERTURBATION to the width of the front                      
      iseed= 44294 
      dum = ran3(iseed) 
                                                                        
!     z1 is the depth of the ml (diff rho on both sides above this depth
!     z2 is the vertical extent of the density anamoly (it is gradually 
!        linearly anihillated with depth).                              
!     Orig vals. z1= 50. z2=250.                                        
      z1= 50. + zoffset 
      z2= 250. + zoffset 
!      slfac= -0.15d0                                                   
!      slfac= -0.12d0 This was used for SB front
!      slfac= 0.12d0      ! Change sign for Arctic shelf (if dense on shelf)
      slfac=0.d0
!     0.12 in pot dens, 0.15 in salinity                                
      do i=0,NI+1 
! perturb                                                               
!         dscl=dyfront* sclwidth                                        
!         dscl=(1.d0 + 1.d-4*ran3(iseed))*dscl                          
!         rdscl=1.d0/dscl                                               
!         dscl=dscl*1.d-2                                               
                                                                        
         do j=0,NJ+1 
            thy = tanh(tightness*(yc(j)-yfront)*PI) 
            do k=0,NK+1 
               z= DL*zc(i,j,k) 
                  if (z.ge.-z1) then 
                     sbkgrnd =                                          &
     &                    seval(np,-1.*z1,depoff,rhovert,bs1,cs1,ds1)   &
     &                    - (z+z1)*bfsqbkgrnd*R0/(gpr*10.)              
                  else 
                     sbkgrnd =                                          &
     &                    seval(np,z,depoff,rhovert,bs1,cs1,ds1)        
                  end if 
                  if (z.ge.-z1) then 
                     slfacnew = slfac 
                  else if (z.ge.-z2) then 
                     slfacnew = slfac*(z+z2)/(z2-z1) 
                  else 
                     slfacnew = 0. 
                  end if 
                  s(i,j,k,0)= -slfacnew*thy + sbkgrnd 
                                                                        
!-               else                                                   
!-                  s(i,j,k,0)=                                         
!-     &                 seval(np,z,dep,svert,bs1,cs1,ds1) -S0          
!-               end if                                                 
            end do 
         end do 
      end do 
                                                                        
      return 
      END                                           
