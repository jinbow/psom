program main 
!---------------------------------------------------                    
!     test model                                                        
      USE header
#ifdef runtracmass
      USE mod_vel
      USE interface
#endif
      implicit logical (a-z) 
!                                                                       
!      integer NI,NJ,NK,lexp,np                                         
!     parameters:                                                       
!     NI,NJ,NK are the no. of points in each of the 3 co-ordintate direc
!     L, DL are the characteristic length and depth scales.             
!     R is the radius of the earth = 6371 km.                           
!     NEQ is the number of equations and the number or unknowns         
!     EPS is the Rosby number                                           
!     AH and AV are the eddy diffusivities in the horizontal and vertica
!     lexp is a constant in the equations denoting the length scale     
!     LEN= 10**(4+lexp),  0<=lexp<=1                                    
!     dtime is the time step                                            
!     S0 and T0 are the mean salinity and temp. sreal= S0+s,Treal=T0+T  
!     np is the number of points which makes up the salinity,temp,dens v
!     depth data curve                                                  
!     OMEGA is the angular velocity of rotation of the earth            
!     dtf is the(dimensionless) time step                               
!     delta is the ratio of horizontal to vertical length scale = L/D   
!     delinv= 1/delta                                                   
!     qpr is the ratio of non-hyd to hyd pressure = Q/P                 
!     qprinv = 1/qpr                                                    
!     lambda is the ration of hor length scale to earth's rad = L/A     
!     kappah is the implcitness parameter in the Crank-Nicolson scheme f
!     the h eqn. kaphinv is its inverse.                                
!     set kappah=1/2.  kaphinv=1/kappah                                 
      integer step,nsteps,n,i,j,k,nblock,ik,inew,iold,frame_int,        &
     &     ngraph,m,budget_int                                          
!      double precision fdiv,p(maxout),coslat(0:NJ+1),                  
!      double precision fdiv,p(maxout),colati(0:NJ+1),                  
      double precision fdiv,pcorr(maxout),                              &
     &     ctrdiv,hsum,hmean                                            
      real dtime, tarray(2), tim 
      CHARACTER(len=10) :: stepchar

      open(unit=91,file='energy.out') 
      open(unit=31, file='gulf.in',status='old') 
      read(31,*) 
      read(31,*) nsteps,dtf,ngraph,frame_int 
      read(31,*) 
      read(31,*) fnhhy 
      read(31,*) 
      read(31,*) dx,dy 
      close(31) 
      write(6,*) 'nsteps= ',nsteps 
                                                                        
#ifdef runtracmass
      call initialize_tracmass
      call seed_particles
#endif

      budget_int= 10 
!     define the parameters                                             
!=      frame_int= 100  read in gulf.in                                 
!      EPS= 0.1d0 defined in header                                     
      delta= DL/LEN 
      delinv= LEN/DL 
!     set qpr to zero for hydrostatic, delta for nonhydrostatic         
!      qpr= 0.d0                                                        
      qpr= fnhhy*delta 
!      qprinv= 1.d0/qpr not to be used                                  
      kappah= 0.65d0 
      kaphinv= 1.d0/kappah 
!=      kappah= 0.5d0                                                   
!=      kaphinv= 2.d0                                                   
      lambda= LEN/AL 
!     beta= 1.d0/(EPS*EPS*delta)  - orignal value                       
      beta= 1.d0/(EPS*EPS*delta) 
!      beta=100.d0                                                      
      UL= FPAR *LEN*EPS 
!     UL= F*LEN*EPS                                                     
      P1= R0*UL*UL*EPS**(-1) 
!     P1= RU^2/EPS                                                      
      HL= P1/(R0*10.d0) 
!     R0*G*HL= P1                                                       
!     HDL= 1.d-4 or 1.d-3                                               
      HDL= HL/DL 
      TL= LEN/UL 
      WL= EPS*delta*UL 
      write(6,*) 'parameters' 
      write(6,*) 'fnhhy=',fnhhy,'eps=',EPS,'delta=',delta,'qpr',qpr,    &
     &     'lambda=',lambda,'beta=',beta,'dtime=',dtf,'WL',WL           
      write(6,*) 'R0, UL, P1, HL, HL/DL', R0,UL,P1,HL,HDL 
!      call potfunc                                                     
!      call init(p,vfent,ufex)                                          
      call init(pcorr) 
!==      call stprofile (now called from init_tr.f)                     
      call meanh(NI,NJ,h,hmean) 
      write(6,*) 'initial  hmean',hmean 
!-      call stinit    -- for variation with latitude                   
      call sigma 
      call staticsigma 
!-      call staticsigma (move the whole grid)                          
!     important to call sigma before staticsigma and stprofile          
!=      call stprofile                                                  
!=      call bioinit                                                    
                                                                        
      call tracerinit(0) 
!      call setbcuf                                                     
      call hsave 
!     hsave must be called before geostroph                             
!      do n=1,2                                                         
!         write(6,*) 'n= ',n,'geostroph'                                
      call geostroph(0) 
!      end do                                                           
!     sigma needs to be called in momentum after each time advancement  
!     check the divergence of the initial vel field                     
      call facediv(EPS,fdiv) 
      call cdiv(EPS,ctrdiv,0) 
      write(6,*) ' fdiv  ',fdiv,' cdiv  ',ctrdiv 
      write(6,*) 'init called' 
!      call output(p,0,0)                                               
!      call finalout(p,0,0)                                             
!     write out initial values                                          
      call vort(0) 
!     evalrho is called in outcdf                                       
                                                                        
      call evalrho(rho,0)
!     call outcdf(xc,yc,zc,p,h,consump,s,T,rho,Tr,u,v,w,vor,pv,uf,vf,wf,Kz,  &
!    &     con100,conv,0)                                               
!      call writestrain(0)
      step= 0 
      call calcn2 
      call n2budget(step) 
      ksurf= NK 
!     call writeksurf(frame_int,step,ksurf,h,consump,Tr,s,T,rho, &
!    &     u,v,w,p,vor,strain,freqN2,xc,yc,zc,DL,LEN,Jac,dtf*TL)     
!      call writeisopycnal(frame_int,sigrelease,stressx,T,rho,vor,       &
!     &     strain,shear,freqN2,xc,yc,zc,zf,DL,step,TL*dtf)              
!     call writeslice(frame_int,step,consump,Tr,s,T,rho,u,v,w,vor,freqN2,    &
!          shear,fricb,stressx,yc,zc,dtf*TL)  
!       call  writeslicenut(frame_int,step,prod,divreyn,divmean,dcdt,T,   &
!            rho,u,v,w,vor,freqN2,stressx,yc,zc,dtfsec)          
!=      call writeyslice(frame_int,step,consump,T,rho,u,v,w,vor,xc,zc)  
!                                                                       
!     nblock is the number of time steps we can advance with each semtne
!     data file (since each contains 10 time frames, it is undatefreq*10
!      nblock= 10*updatefreq                                            
      write(6,*) 'to start momentum' 
      tim= dtime(tarray) 
!      ik= 1                                                            
!     BCs set just once in the beginning                                
      call setbc(step) 
      call correctbc 
                                                                        
      call checks 

   50 step=step+1 
     WRITE(stepchar,'(I10.10)') step
      !FOR tracmass      
#ifdef runtracmass
      call set_tracmass_fluxes(1)
#endif
      call momentum(pcorr,step) 
#ifdef runtracmass
      call set_tracmass_fluxes(2)
      call tracmass
#endif

!===
!     calcn2 must be called before n2budget. Called at end of momentum  
      call n2budget(step) 
!     compute the KE                                                    
      call energy(step) 
      call meanh(NI,NJ,h,hmean) 
                                                                        
                                            !budget_int=10              
!      if (mod(step,budget_int).eq.0) then 
!         call writen2budget(budget_int,step,dtf,TL,stressmax,mldn2,     &
!     &     zbtop,zbbot,advecPV,friction,diabatic)                       
!      end if 
                                                                        
      if (mod(step,frame_int).eq.0) then 
                                                                        
         call vort(0) 
         call evalrho(rho,0) 
           CALL save3d(NI+2,NJ+2,NK+2,u(:,:,:,1),'op.u.'//stepchar//'.bin')
           CALL save3d(NI+2,NJ+2,NK+2,v(:,:,:,1),'op.v.'//stepchar//'.bin')
           CALL save3d(NI+2,NJ+2,NK+2,w(:,:,:,1),'op.w.'//stepchar//'.bin')
           CALL save3d(NI+2,NJ+2,NK+2,rho,'op.rho.'//stepchar//'.bin')
!         call writeframe(frame_int,step,h,s,T,rho,u,v,w,vor,xc,yc,zc)  
!        call writeksurf(frame_int,step,ksurf,h,consump,Tr,rho,u,v,w,    &
!    &        p,vor,strain,freqN2,xc,yc,zc,DL,LEN,Jac,dtf*TL)           
!         call writeisopycnal(frame_int,sigrelease,stressx,T,rho,vor,    &
!     &        strain,shear,freqN2,xc,yc,zc,zf,DL,step,TL*dtf)           
!        call writeslice(frame_int,step,consump,Tr,s,T,rho,u,v,w,vor,freqN2, &
!    &        shear,fricb,stressx,yc,zc,dtf*TL) 
!       call  writeslicenut(frame_int,step,prod,divreyn,divmean,dcdt,T,   &
!            rho,u,v,w,vor,freqN2,stressx,yc,zc,dtfsec)          
!=         call writeyslice(frame_int,step,consump,T,rho,u,v,w,vor,xc,zc
      endif 
      write(6,*) 'steps',step, ' hmean',hmean 
      if (mod(step,100).eq.0) then 
         tim= dtime(tarray) 
         write(6,*) 'cpu time for 100 steps=', tim 
      endif 
      if (step.eq.nsteps) goto 100 
      if (mod(step,ngraph).eq.0) goto 100 
      goto 50 
!                                                                       
  100 n=0 
      call vort(n) 
!     call outcdf(xc,yc,zc,p,h,consump,s,T,rho,Tr,u,v,w,vor,pv,uf,vf,wf,Kz, &
!    &     conv,con100,step)                                            
      if (step.lt.nsteps) goto 50 
      close(91) 
      stop 
      END                                           
