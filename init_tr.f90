      subroutine init(pcorr) 
!      subroutine init(p,vfent,ufex)                                    
!     ------------------------------------------------                  
!     For the curvilinear grid.                                         
!     sigma levels are evenly spaced.  Hence wz is a function of        
!     time but not of z.   Here u,v,w refer to the xi,eta,sigma direcito
!     Those metric quantities that do not change with time are evaluated
!     The rest are evaluated in "sigma.f" which will be called at every 
!     step.                                                             
!     Also it is absolutely essential that the integral of the flux     
!     over the entrance section is equal to that over the exit section. 
!     -No longer necessary with the free surface.                       
!     This may require some special attention when the free surface is  
!     moving at the exit. At the entrance vf is held fixed.             
!                                                                       
      USE header 
      implicit double precision (a-h,o-z) 
!      double precision lat(0:NI+1,0:NJ+1),lon(0:NI+1,0:NJ+1)           
!      double precision lat(0:NJ+1)                                     
      double precision xdu(0:NI+1,0:NJ+1),ydu(0:NI+1,0:NJ+1),           &
     &     xdv(0:NI+1,0:NJ+1),ydv(0:NI+1,0:NJ+1),C,fconst,              &
     &     phi0,cosphi0,sinphi0,sumuf,sumvf,dthet,dtheta,dphi,          &
     &     gin,c1,c2,temp,temp2,dnkinv,tmplon,tmplat,cjent,             &
     &     phi0deg,pcorr(maxout),dep,y,DLinv                            
!                                                                       
!     FPAR is the value used for scaling Coriolis parameters f and b    
      integer i,j,k,jmid,it,iseed 
                                                                        
!      open(unit=40, file='xyD.dat')                                    
!      open(unit=45, file='du.metrics.dat')                             
!      open(unit=50, file='dv.metrics.dat')                             
!      open(unit=55, file='partialsD.dat')                              
!      open(unit=110,file='ufopen.in')                                  
      write(6,*) 'in init' 
      DLinv= 1.d0/DL 
      C= PI/180.d0 
      gin= 1.d0/gpr 
      dnkinv= 1.d0/dble(NK) 
!                                                                       
!     dx,dy are the grid dims (in m) read in in main.f                  
!     phi0 is the central lat, dphi,dtheta grid spacing(angle)          
!=      phi0deg= 25.d0                                                  
      phi0deg=35.d0 
!      phi0deg= 3.d0                                                    
      phi0= phi0deg*C 
      dphi = dy/(apr*AL) 
      write(6,*) 'dphi = ',dphi, 'in deg ',dphi/C 
!     The domain should be within the latitude limits  in stinit.f.     
!     lat0 is the center latitude about which we linearize              
      jmid= NJ/2 
      write(6,*) 'jmid=',jmid 
      do 10 j=0,NJ+1 
         latrad(j)= phi0 + dble(j-jmid)*dphi 
!     f const                                                           
!         latrad(j)= phi0                                               
   10 continue 
!      do 10 j=0,NJ+1                                                   
!         do 20 i=0,NI+1                                                
!            read(40,*) lon(i,j),lat(i,j),D(i,j)                        
!            lat(i,j)= lat(i,j)*C                                       
!     rescale lat and lon so that the grid spacing is about 1/8 degree  
!            lon(i,j)= lon(i,j)*C                                       
! 20      continue                                                      
! 10   continue                                                         
!      close(40)                                                        
!c      cosphi0= dCos(phi0)                                             
!c      sinphi0= dSin(phi0)                                             
!c      write(6,*) 'cosphi0=',cosphi0                                   
!      lon0= 260.d0*C                                                   
!      lat0= 20.d0*C                                                    
!      dtheta= C/8.d0                                                   
!      dphi= C/8.d0                                                     
!      phi0= lat0 + 0.5d0*(dphi*dfloat(NJ))                             
!c     for the rectilinear case, dx,dy,dz are constant                  
!      dx= AL*apr*cosphi0*dtheta /LEN                                   
!      dy= AL*apr*dphi/LEN                                              
!      dz= 1000.d0/(dfloat(NK)*DL)                                      
!      write(6,*) 'dx,dy,dz',dx,dy,dz                                   
!      do 10 i=-1,NI+1                                                  
!         do 10 j=-1,NJ+1                                               
!            lat(i,j)= lat0 + dphi*dfloat(j)                            
!            lon(i,j)= lon0 + dtheta*dfloat(i)                          
! 10   continue                                                         
!     These values must be read in for the curvilinear grid.            
!c      c1= apr*AL*cosphi0/LEN                                          
!c      c2= apr*AL/LEN                                                  
      do 30 j=0,NJ+1 
         do 30 i=0,NI+1 
!            xdu(i,j)= c1*dtheta                                        
            xdu(i,j)= dx/LEN 
            ydu(i,j)= 0.d0 
            xdv(i,j)= 0.d0 
            ydv(i,j)= dy/LEN 
   30 continue 
!                                                                       
      yc(0) = -0.5*dy*1.d-3 
      do j=1,NJ+1 
         yc(j)= yc(j-1) + dy*1.d-3 
      end do 
      xc(0) = -0.5*dx*1.d-3 
      do i=1,NI+1 
         xc(i)= xc(i-1) + dx*1.d-3 
      end do 
!                                                                       
      do 50 j=0,NJ+1 
         do 60 i=0,NI+1 
            J2d(i,j)= xdu(i,j)*ydv(i,j) -xdv(i,j)*ydu(i,j) 
!            write(300,13) xdu(i,j),ydv(i,j),xdv(i,j),ydu(i,j),J2d(i,j) 
! 13         format(5(e8.2,2x))                                         
            ux(i,j)= ydv(i,j)/J2d(i,j) 
            vx(i,j)= -ydu(i,j)/J2d(i,j) 
            uy(i,j)= -xdv(i,j)/J2d(i,j) 
            vy(i,j)= xdu(i,j)/J2d(i,j) 
            g11(i,j)= ux(i,j)*ux(i,j) +uy(i,j)*uy(i,j) 
            g12(i,j)= ux(i,j)*vx(i,j) +uy(i,j)*vy(i,j) 
            g22(i,j)= vx(i,j)*vx(i,j) +vy(i,j)*vy(i,j) 
   60    continue 
   50 continue 
!                                                                       
!     D(i,j) is -ve and  non-dim by DL                                  
      !for MLI runs                           
!===      dep= 500.d0 
      dep= 1000.d0 
      D= -dep*DLinv
!     Compute partial  dD/dx,dD/dy                                      
      call smooth 
      write(6,*) 'smmoth called' 
                                                                        
! 101  call findzall                                                    
  101 continue 
!      close(55)                                                        
!c      open (unit= 31, file='xyD.dat')                                 
!c      do j=0,NJ+1                                                     
!c         tmplat= phi0deg + dble(j-jmid)*dthet                         
!c         yc(j)= tmplat                                                
!c         do i=0,NI+1                                                  
!c            tmplon= 300.d0 +dthet*dble(i)                             
!c            write(31,*) tmplon,tmplat,D(i,j)                          
!c            if (j.eq.0) xc(i)= tmplon                                 
!c         end do                                                       
!c      end do                                                          
!c      close(31)                                                       
!     In smooth we smooth D and then use central differencing to evaluat
!     Ddx and Ddy. The values obtained are much lower than those obtaine
!     by smoothing Ddx, Ddy evaluated from the bspline - probably bec   
!     we had a bug with a factor of 3 which is now corrected. (values   
!     differ even after the correction and we use the numerial values.) 
!c      call smooth                                                     
!*    flat bottom                                                       
! 999  do j=0,NJ+1                                                      
!         do i=0,NI+1                                                   
!            Ddx(i,j)= 0.d0                                             
!            Ddy(i,j)= 0.d0                                             
!            D(i,j)= 1.d0                                               
!            write(150,*) Ddx(i,j)                                      
!            write(250,*) Ddy(i,j)                                      
!            write(350,*) D(i,j)                                        
!         end do                                                        
!      end do                                                           
!      stop                                                             
!                                                                       
      do 55 j=0,NJ+1 
         do 65 i=0,NI+1 
            do 75 k=0,NK+1 
               conv(i,j,k)= 0 
               con100(i,j,k)=0 
               u(i,j,k,0)= 0.d0 
               v(i,j,k,0)= 0.d0 
               w(i,j,k,0)= 0.d0 
               s(i,j,k,0)= 0.d0 
               T(i,j,k,0)= 0.d0 
               do it=1,ntr 
                  Tr(it,i,j,k,0)= 0.d0 
               end do 
!     p is initialized in geostroph                                     
!              p(i,j,k)= 0.d0                                           
   75       continue 
   65    continue 
   55 continue 
      do k=1,maxout 
         pcorr(k)= 0.d0 
      end do 
!                                                                       
      do k=1,NK 
         do j=1,NJ 
            do i=1,NI 
               uvis(i,j,k)= 0.d0 
               vvis(i,j,k)= 0.d0 
               wvis(i,j,k)= 0.d0 
            end do 
         end do 
      end do 
!                                                                       
!     specify the initial fields                                        
!     read in the cell-centered pressure and u,v velocities             
      do 81 j=0,NJ+1 
         do 82 i=0,NI+1 
!            hdt(i,j)= 0.d0                                             
            h(i,j)= 0.d0 
   82    continue 
   81 continue 
!                                                                       
      do 202 k=1,NK 
         do 202 j=1,NJ 
            do 202 i=0,NI 
               uf(i,j,k)= 0.d0 
  202 continue 
!                                                                       
      do 203 k=1,NK 
         do 203 j=0,NJ 
            do 203 i=1,NI 
               vf(i,j,k)= 0.d0 
  203 continue 
!                                                                       
      do 105 k=0,NK 
         do 105 j=1,NJ 
            do 105 i=1,NI 
               wf(i,j,k)= 0.d0 
  105 continue 
!      close(110)                                                       
!                                                                       
!     fill in the vfbc arrays                                           
      do 130 k=1,NK 
         do 131 j=1,NJ 
            ufbce(j,k)= uf(NI,j,k) 
            ufbcw(j,k)= uf(0,j,k) 
  131    continue 
  130 continue 
      do 140 k=1,NK 
         do 141 i=1,NI 
            vfbcn(i,k)= vf(i,NJ,k) 
            vfbcs(i,k)= vf(i,0,k) 
  141    continue 
  140 continue 
!     at the sea bed, wf=0                                              
      do 150 j=1,NJ 
         do 151 i=1,NI 
!            wfbct(i,j)= 0.d0                                           
            wfbcb(i,j)= 0.d0 
  151    continue 
  150 continue 
!                                                                       
!     Evaluate the Coriolis parameter f' and fill the arrays ffi,ffj,ffc
!     f= f0+By= 2Omega(Sin(phi0) +(phi-phi0)Cos(phi0))                  
!     b= b0+By= 2Omega(Cos(phi0) -(phi-phi0)Sin(phi0))                  
      fconst= 2.d0*OMEGA/FPAR 
      do 230 j=1,NJ 
         ffi(0,j)= fconst*dSin(latrad(j)) 
         bbi(0,j)= fnhhy*fconst*dCos(latrad(j)) 
!         do 230 i=0,NI                                                 
         do 230 i=1,NI 
         ffi(i,j)= ffi(0,j) 
         bbi(i,j)= bbi(0,j) 
!         ffi(i,j)= fconst*(sinphi0 +(latrad(i,j)-phi0)*cosphi0)        
!         bbi(i,j)= fconst*(cosphi0 -(latrad(i,j)-phi0)*sinphi0)        
  230 continue 
      do 240 j=0,NJ 
            ffj(1,j)= fconst*dSin(0.5d0*(latrad(j+1)+latrad(j))) 
            bbj(1,j)= fnhhy*fconst*dCos(0.5d0*(latrad(j+1)+latrad(j))) 
!         do 240 i=1,NI                                                 
         do 240 i=2,NI 
            ffj(i,j)= ffj(1,j) 
            bbj(i,j)= bbj(1,j) 
  240 continue 
      do 250 j=0,NJ+1 
         ffc(0,j)= fconst*dSin(latrad(j)) 
         bbc(0,j)= fnhhy*fconst*dCos(latrad(j)) 
!         do 250 i=0,NI+1                                               
         do 250 i=1,NI+1 
         ffc(i,j)= ffc(0,j) 
         bbc(i,j)= bbc(0,j) 
  250 continue 
!                                                                       
!     Initialize s,T                                                    
      call findzall 
      call stprofile 
      call evalrho(rho,0) 
      call inith 
      call findzall 
!     write out z-grid                                                  
!     -----------------                                                 
      open (unit=60,file='zgrid.out') 
      write(60,*) '#vertical grid' 
      do k=0,NK+1 
         write(60,*) k,zc(10,10,k)*1000. 
      end do 
      write(60,*) '#face values' 
      do k=-1,NK+1 
         write(60,*) k,zf(10,10,k)*1000. 
      end do 
      close(60) 
!     -----------                                                       
                                                                        
      do m=1,3 
         advecpv(m)= 0.d0 
         friction(m)= 0.d0 
         diabatic(m)= 0.d0 
      end do 
                                                                        
      return 
                                                                        
      drho= rho(NI/2,5*NJ/6,NK)-rho(NI/2,NJ/6,NK) 
      write(6,*) 'drho = ',drho 
      hlfac = 0.5*drho*100./DL/HL 
!     100 is the depth of the front : upper ocean  depth                
                                                                        
!     ADD a PERTURBATION to the width of the front                      
      iseed= 44294 
      dum = ran3(iseed) 
                                                                        
!cc      hlfac=.1d0/HL  (set above)                                     
      do i=0,NI+1 
! perturb                                                               
         dscl=0.1d0*(0.5*(yc(NJ)+yc(NJ+1))-0.5*(yc(0)+yc(1))) 
         dscl=(1.d0 + 1.d-4*ran3(iseed))*dscl 
         rdscl=1.d0/dscl 
         dscl=dscl*1.d-2 
         if (i.eq.10) open(65) 
!         write(6,*) 'dscl',dscl                                        
         do j=0,NJ+1 
!     ulfac=gpr*hlfac/(dscl*ffc(0,j))                                   
            yarg=yc(j)-  0.5*(yc(nj/2)+yc(nj/2+1)) 
            ex2y=exp(-2.d0*rdscl*abs(yarg)) 
            thy=(1.d0-ex2y)/(1.d0+ex2y) 
            if(yarg.lt.0.d0) thy=-thy 
            sechy2=1.d0-thy*thy 
                                                                        
            hdt(i,j)= 0.d0 
            h(i,j)= -hlfac*thy 
         end do 
         if (i.eq.10) write(65,*) j,h(10,j),u(10,j,15,0) 
      end do 
      close(65) 
                                                                        
      return 
                                                                        
  201 hlfac=.1d0/HL 
      dscl=0.1d0*(0.5*(yc(NJ)+yc(NJ+1))-0.5*(yc(0)+yc(1))) 
!cmy      dscl=0.25d0*(0.5*(yc(NJ)+yc(NJ+1))-0.5*(yc(0)+yc(1)))         
      rdscl=1.d0/dscl 
      dscl=dscl*1.d-2 
      open(65) 
      do j=0,NJ+1 
!         ulfac=gpr*hlfac/(dscl*ffc(0,j))                               
         yarg=yc(j)-  0.5*(yc(nj/2)+yc(nj/2+1)) 
         ex2y=exp(-2.d0*rdscl*abs(yarg)) 
         thy=(1.d0-ex2y)/(1.d0+ex2y) 
         if(yarg.lt.0.d0) thy=-thy 
         sechy2=1.d0-thy*thy 
         do i=0,NI+1 
            hdt(i,j)= 0.d0 
            h(i,j)= -hlfac*thy 
!            do k=0,nk+1                                                
!               u(i,j,k,0)=ulfac*sechy2                                 
!            enddo                                                      
         end do 
         write(65,*) j,h(10,j),u(10,j,15,0) 
      end do 
      close(65) 
                                                                        
      return 
                                                                        
!     Perturb the front boundary                                        
!     for a quicker onset of instability                                
!     ----------------------------------------------------              
      iseed= 44294 
      dum = ran3(iseed) 
!                                                                       
      do i=0,NI+1 
         perturb=  1.d-5*ran3(iseed) 
         do j=1,NJ 
            h(i,j)= h(i,j) +perturb 
         end do 
      end do 
!                                                                       
!     call evalrho(0)                                                   
!      write(6,*) 'in initsq, evalrho called'                           
!      write(6,*) 'writing out layer depths'                            
!      do k=0,NK                                                        
!         temp= findz(dble(k),dztop,D(10,10),0.d0,NK)                   
!         write(6,*) 'k,z',k,temp                                       
!      end do                                                           
      write(6,*) 'to leave init' 
      return 
!                                                                       
      return 
      END                                           
