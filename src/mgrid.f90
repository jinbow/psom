      subroutine mgrid(p,dtime,edt,cfcdiv) 
      implicit logical (a-z) 
!     --------------                                                    
!     multigrid solver for the elliptic equation Del2 p= f              
!     Run pgrogram "preprocess" to get the value of maxout              
!     input :                                                           
!     -----                                                             
!     call the subroutine mgrid with p_initial on fine grid             
!     maxout : max dimension on the one-dim array - includes outer point
!     maxint: max dimension of the one-dim array for variables at interi
!              points of grid                                           
!     int1 :  number of internal points on the  fine grid               
!     ngrid : number of grid levels                                     
!     nx(m),ny(m),nz(m) = number of int grid points at level m. m=1...ng
!     m=1 is the finest grid, m=ngrid refers to the coarsest grid.      
!     ntout(m) is the total number of grid points (storage locations) at
!               grid level m - including the outer ficticious points    
!     ntint(m) is the total number of interior grid points (storage loca
!               at grid level m                                         
!     loco(m) is the storage location of a scalar variable in a one-dim 
!             of dimension  Summation(m=1,ngrid) {ntout(m)}             
!     loci(m) is the storage location of a scalar specified only at inte
!             grid points                                               
!     loccp(m) is the storage location for cp on the m-th level         
!     res(i)  temporary stores the residual before it interpolated onto 
!             the coarse grid                                           
!     dxf,dyf,dzf  are the grid spacings on the fine grid               
      integer ngrid,m,maxout,maxint,int1,l,ncycle,ll,                   &
     &     nu1,nu2,noc                                                  
!c     48x32x16 grid                                                    
!c      parameter(ngrid=4,maxout=36312,maxint=28080,int1=24576)         
!     48x24x32 grid                                                     
!      parameter(ngrid=4,maxout=13272,maxint=9360,int1=8192)            
!-      parameter(ngrid=5,maxout=46416,maxint=37448,int1=32768)         
!1k     96x96x24                                                        
!      parameter(ngrid=4,maxout=291092,maxint=252720,int1=221184)       
!.5k  96x192x32                                                         
      parameter(ngrid=5,maxout=750240,maxint=674064,int1=589824) 
!.5k  96x192x64 - increased vertical resol
!!      parameter(ngrid=6,maxout=1449264,maxint=1348164,int1=1179648) 
!.big  192x384x32                                                       
!      parameter(ngrid=5,maxout=2946528,maxint=2696256,int1=2359296)    
!.big  96x384x32                                                        
!      parameter(ngrid=5,maxout=1491264,maxint=1348128,int1=1179648)    
!       192x192x32                                                      
!      parameter(ngrid=5,maxout=1482336,maxint=1348128,int1=1179648)    
!4k   12x24x247884                                                      
!      parameter(ngrid=3,maxout=11352,maxint=7884,int1=6912)            
!2k   24x48x24                                                          
!      parameter(ngrid=4,maxout=39992,maxint=31590,int1=127648)         
!1k   48x96x24                                                          
!      parameter(ngrid=4,maxout=149072,maxint=126360,int1=110592)       
!1k   48x96x32                                                          
!      parameter(ngrid=5,maxout=194472,maxint=168516,int1=147456)       
!.25k  192x384x24  (insuff. memory on vayu)                             
!      parameter(ngrid=4,maxout=2258852,maxint=2021760,int1=1769472)    
!.5km  96x192x24                                                        
!=      parameter(ngrid=4,maxout=575132,maxint=505440,int1=442368)      
!test125m  16x96x24                                                     
!      parameter(ngrid=4,maxout=54392,maxint=42120,int1=36864)          
!.5k  192x96x16                                                         
!      parameter(ngrid=4,maxout=400472,maxint=336960,int1=294912)       
!-     48x48x24                                                         
!-      parameter(ngrid=4,maxout=76352,maxint=63180,int1=55296)         
!     32x32x16                                                          
!      parameter(ngrid=4,maxout=24792,maxint=18720,int1=16384)          
!=32x64x16      parameter(ngrid=4,maxout=47832,maxint=37440,int1=32768) 
!=      parameter(ngrid=4,maxout=92312,maxint=74880,int1=65536)         
!      parameter(ngrid=4,maxout=92312,maxint=74880,int1=65536)          
      integer ntout(ngrid),ntint(ngrid),nx(ngrid),ny(ngrid),nz(ngrid),  &
     &     loco(ngrid),loci(ngrid),loccp(ngrid)                         
      double precision cp(19*maxint),rhs(maxint),p(maxout),             &
     &     res(int1),dtime,tol,maxres,edt,ratio,oldres,cfcdiv           
!     edt= EPS/dtime. If the tolerance is on (u_x +v_y +eps*w_z)        
!     then  the tolerance on the residual r= edt*(u_x+ v_y + eps*w_z)   
!     is equal to edt*(the specified tolerance)                         
!                                                                       
      open (unit=90, file='mg.in') 
      read(90,*) nx(1),ny(1),nz(1),noc,nu1,nu2,tol 
      close(90) 
                                                                        
      if (cfcdiv.le.tol) then 
         do m=1,maxout 
            p(m)= 0.d0 
         end do 
         maxres = edt*cfcdiv 
         ncycle = 0 
         goto 101 
      end if 
                                                                        
!     redefine the tolerance as tol*edt                                 
      tol= edt*tol 
      ntint(1)= nx(1)*ny(1)*nz(1) 
      ntout(1)= (nx(1)+2)*(ny(1)+2)*(nz(1)+2) 
!                                                                       
      loco(1)= 1 
      loci(1)= 1 
      loccp(1)= 1 
      do m=2,ngrid 
         if (mod(nx(m-1),2).ne.0) goto 20 
         if (mod(ny(m-1),2).ne.0) goto 20 
         if (mod(nz(m-1),2).ne.0) goto 20 
         nx(m)= nx(m-1)/2 
         ny(m)= ny(m-1)/2 
         nz(m)= nz(m-1)/2 
         ntint(m)= nx(m)*ny(m)*nz(m) 
         ntout(m)= (nx(m)+2)*(ny(m)+2)*(nz(m)+2) 
         loco(m)= loco(m-1) + ntout(m-1) 
         loci(m)= loci(m-1) + ntint(m-1) 
!     for the 19 point stencil                                          
         loccp(m)= loccp(m-1) + 19*ntint(m-1) 
      end do 
!                                                                       
!     compute fine grid coefficients and rhs on fine grid               
      call cpfine(dtime,cp,rhs) 
!      call checkcp(1,nx(1),ny(1),nz(1),cp(loccp(1)))                   
!      write(6,*) 'cpfine checked'                                      
!     **********************************                                
      do m=2,ngrid 
         call cpcors(nx(m),ny(m),nz(m),cp(loccp(m-1)),cp(loccp(m))) 
!         call checkcp(m,nx(m),ny(m),nz(m),cp(loccp(m)))                
      end do 
!                                                                       
      do 1000 ncycle=1,noc 
!+         if (mod(ncycle,3).eq.0) then                                 
!+c     then we update the rhs                                          
!+            call intermq(dtime)                                       
!+            call remcor                                               
!+            call newsrc                                               
!+            call exitp                                                
!+            call newqfn(dtime,rhs)                                    
!+         endif                                                        
!w         write(6,*) 'cycle=',ncycle                                   
!     initialize the coarse grid values of p                            
         do l=loco(2),maxout 
            p(l)= 0.d0 
         end do 
!     DOWN CYCLE                                                        
!     for m=1                                                           
      m=1 
!     -----------------------                                           
!wc     write out starting resid                                        
!w      call resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),                    
!w     &        p(loco(m)),rhs(loci(m)),res,maxres)                     
!w      write(6,*) 'starting resid=',maxres                             
!     ---------------------------------                                 
      do 11 l=1,nu1 
         call linerelax(nx(m),ny(m),nz(m),cp(loccp(m)),p(loco(m)),      &
     &        rhs(loci(m)))                                             
!     *******try and put mgpfill here ******                            
         call mgpfill(dtime,p) 
!         call mgpfill(dtime,p)                                         
   11 continue 
!c      call mgpfill(dtime,p)                                           
      call resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),                      &
     &        p(loco(m)),rhs(loci(m)),res,maxres)                       
!w      if (ncycle.ne.1) then                                           
!w         ratio = maxres/oldres                                        
!w         write(6,*) ncycle,'maxres',maxres,'conv ratio=',ratio        
!w      endif                                                           
!w      oldres= maxres                                                  
!w      write(6,*) 'm= ',m, 'maxres=',maxres                            
      if (maxres.lt.tol) goto 101 
      call restrict(nx(m),ny(m),nz(m),res,rhs(loci(m+1))) 
!      do 100 m=1,ngrid-1                                               
      do 100 m=2,ngrid-1 
         do 21 l=1,nu1 
            call linerelax(nx(m),ny(m),nz(m),cp(loccp(m)),p(loco(m)),   &
     &           rhs(loci(m)))                                          
   21    continue 
         call resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),                   &
     &        p(loco(m)),rhs(loci(m)),res,maxres)                       
!w         write(6,*) 'm= ',m, 'maxres=',maxres                         
         call restrict(nx(m),ny(m),nz(m),res,rhs(loci(m+1))) 
  100 continue 
!                                                                       
!     UP CYCLE                                                          
      do 200 m= ngrid,2,-1 
         if (m.eq.ngrid) then 
   25       call sor(nx(m),ny(m),nz(m),cp(loccp(m)),p(loco(m)),         &
     &           rhs(loci(m)))                                          
         else 
            do 41 l=1,nu2 
               call linerelax(nx(m),ny(m),nz(m),cp(loccp(m)),           &
     &              p(loco(m)),rhs(loci(m)))                            
   41       continue 
         endif 
!     this call to resid is not necessary - temporary for res check     
!w         call resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),                 
!w     &        p(loco(m)),rhs(loci(m)),res,maxres)                     
!w         write(6,*) 'm= ',m, 'maxres=',maxres                         
         call efill(nx(m),ny(m),nz(m),p(loco(m))) 
         call prolong(nx(m),ny(m),nz(m),p(loco(m)),p(loco(m-1)) ) 
  200 continue 
!      do l=1,3                                                         
!         call  mgpfill(dtime,p)                                        
!      end do                                                           
 1000 continue 
!*      write(6,*) '100 mg V cycles, not converged'                     
!     added Nov 17,1993                                                 
      m=1 
      do 111 l=1,nu1 
         call linerelax(nx(m),ny(m),nz(m),cp(loccp(m)),p(loco(m)),      &
     &        rhs(loci(m)))                                             
!     *******try and put mgpfill here ******                            
         call  mgpfill(dtime,p) 
  111 continue 
!                                                                       
      m=1 
      call resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),                      &
     &     p(loco(m)),rhs(loci(m)),res,maxres)                          
!w      write(6,*) 'max-res=',maxres                                    
!     normliz can be called only for the case of dirichlet bcs.         
  101 continue 
! 101  call normliz(nx(1),ny(1),nz(1),p)                                
!      do l=1,3                                                         
!           call  mgpfill(dtime,p)                                      
!      end do                                                           
!     maxres/edt is the value of (u_x +v_y +ep*w_Z)                     
      maxres= maxres/edt 
      write(6,*) 'mg-cycle ', ncycle,'  maxres=',maxres 
!     this call to resid is not necessary - temporary for res check     
!      call resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),                     
!     &     p(loco(m)),rhs(loci(m)),res,maxres)                         
      goto 201 
   20 write(6,*) 'cannot coarsen this grid as specified' 
                                                                        
  201 return 
      END                                           
