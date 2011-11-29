subroutine momentum(pcorr,step) 
!----------------------------------------------------                   
  USE header
      implicit logical (a-z) 
      integer step,i,j,k 
      double precision dtim,fdiv,cfcdiv,ctrdiv,edt,pcorr(maxout),temp 
      real dtime, tarray(2),tim 
!                                                                       
!     --------------------------------------------------------          
!     initialized with u(t=0) and uf(t=0), ie. u(n=0),uf(n=0)           
!     ------------------------------------------------------            
!     get ready for the 1st step                                        
!     evalrp and srcface were callled at the finish or the previous step
!     or in init                                                        
!                                                                       
!     save old value of h  in array oldh
      oldh= h

!     1st R-K step                                                      
!     -----------------                                                 
      dtim= dtf/3.d0 
!     compute new zc,zf and wz,wx,wy and Jac (since grid has moved)     
      call findzall 
      call sigma 
!      call outflowbc(0,dtim)                                           
!     advecn(m,n,dtim) means use level m and write s,T to level n       
      call advecn(0,1,dtim,step) 
!     call intpol with the time level at which uf is known              
!     besides interpolating the convective terms, intpol computes u*    
      call intpol 
!     we call the following 3 after advecn, in which case               
!     evalrp would use the new value of s,T(n+1), but old value of h(n) 
      call rpevalgrad(1) 
!      call rpflatbtm(1)                                                
!c      call evalrp(1)                                                  
!     evalrp must be called before srcface or coriolis                  
!      call viscous(0)                                                  
!c      call gradrp                                                     
!-      call vertmix(0)                                                 
      call coriolis(0) 
      call srcface(0,step) 
!                                                                       
!      open(unit=70,file='hbef.out')                                    
!      do j=0,NJ+1                                                      
!         write(70,*) h(NI/2,j)                                         
!      end do                                                           
!      close(70)                                                        
      call hsolve(h,oldh,hdt,dtim) 
      call calcskfc 
      call vhydro(dtim) 
      call cfdiv(cfcdiv) 
!==                                                                     
!      write(6,*) 'mom, frame_int',frame_int                            
!      call writeslice(1,2,consump,T,rho,cx,cy,cz,vor,freqN2,           
!     &     fricb,yc,zc,dtf*TL)                                         
!      write(6,*) 'from mom, call outcdf and stop'                      
!      call outcdf(xc,yc,zc,p,h,consump,s,T,cx,cy,cz,vor,pv,cxf,cyf,czf,
!     &     conv,con100,1)                                              
!      stop                                                             
! v,w, present beacuse cannot bal both v and vf in geostrophic perfectly
!     in calling newcor, n should be the same as in coriolis            
!c      call evalrp(1)                                                  
      call newcor(dtim,0) 
      call newsrc 
!     compute q (correction to NH press p)                              
!     cpuf must be called before mgrid which calls vface and overwrites 
!c      call exitp(p)                                                   
      edt= EPS/dtim 
      call mgrid(pcorr,dtim,edt,cfcdiv) 
!      call pfilter(pcorr)                                              
!c      call linesolve(p,dtim,step,edt)                                 
!      call psolve(p,dtim,step,edt)                                     
!c      call pnsolve(p,dtim,step,edt)                                   
!      call outcdf(xc,yc,zc,pcorr,h,s,T,cx,cy,cz,hdt,cxf,cyf,czf,Kz,conv
!      stop                                                             
      call vface(pcorr,dtim) 
!c      call vnface(p,dtim)                                             
!     fill f(n=1), use f(n=0)                                           
!     call vcenter with the value of n that we wish to fill             
      call vcenter(pcorr,dtim,1) 
!     computes the final vel (at n+1 -th time step) at the cell centers 
!     Now correct the NH pressure                                       
      if (fnhhy.ne.0) call pcorrect(pcorr) 
      call facediv(dtim,fdiv) 
      call cdiv(dtim,ctrdiv,1) 
      write(6,*)  'rk-1 cfcdiv ',cfcdiv,' fdiv ',fdiv,' cdiv ',ctrdiv 
!      call outcdf(xc,yc,zc,pcorr,h,s,T,u(0,0,0,1),v(0,0,0,1),          
!     &     w(0,0,0,1),vor,uf,vf,wf,Kz,conv,con100,1)                   
!      write(6,*) 'stopping in momentum.f'                              
!      stop                                                             
!                                                                       
!     2nd R-K step                                                      
!     -----------------                                                 
      dtim= 0.5d0*dtf 
!     compute new zc,zf and wz,wx,wy and Jac (since grid has moved)     
      call findzall 
      call sigma 
!      call outflowbc(1,dtim)                                           
!     advecn(m,n,dtim) means use level m and write s,T to level n       
      call advecn(1,1,dtim,step) 
!     call intpol with the time level at which uf is known              
!     besides interpolating the convective terms, intpol computes u*    
      call intpol 
!     we call the following 3 after advecn, in which case               
!     evalrp would use the new value of s,T(n+1), but old value of h(n) 
!      tim= dtime(tarray)                                               
      call rpevalgrad(1) 
!      call rpflatbtm(1)                                                
!      tim= dtime(tarray)                                               
!      write(6,*) 'time for rpevalgrad',tim                             
!c      call evalrp(1)                                                  
!     evalrp must be called before srcface or coriolis                  
!      call viscous(1)                                                  
!c      call gradrp                                                     
!-      call vertmix(1)                                                 
      call coriolis(1) 
      call srcface(1,step) 
!cc      call bcmake(dtim)                                              
      call hsolve(h,oldh,hdt,dtim) 
      call calcskfc 
!c      call hnsolve(h,oldh,hdt,dtim)                                   
      call vhydro(dtim) 
      call cfdiv(cfcdiv) 
!      call outcdf(xc,yc,zc,p,h,s,T,cx,cy,cz,hdt,cxf,cyf,czf,Kz,conv,con
!      stop                                                             
!      call ufcorrector                                                 
!     in calling newcor, n should be the same as in coriolis.           
!c      call evalrp(1)                                                  
      call newcor(dtim,1) 
      call newsrc 
!c      call cpsi                                                       
!     compute q                                                         
!      call cpuf                                                        
!c      call exitp(p)                                                   
      edt= EPS/dtim 
      call mgrid(pcorr,dtim,edt,cfcdiv) 
!      call pfilter(pcorr)                                              
!c      call linesolve(p,dtim,step,edt)                                 
!-      call psolve(p,dtim,step,edt)                                    
!c      call pnsolve(p,dtim,step,edt)                                   
!     ********************                                              
      call vface(pcorr,dtim) 
!c      call vnface(p,dtim)                                             
!     computes the vel at the cell faces                                
!     fill f(n=1), use f(n=0)                                           
!     call vcenter with the value of n that we wish to fill             
      call vcenter(pcorr,dtim,1) 
!     computes the final vel (at n+1 -th time step) at the cell centers 
!     Now correct the NH pressure                                       
      if (fnhhy.ne.0) call pcorrect(pcorr) 
      call facediv(dtim,fdiv) 
!      if (step.eq.27) write(6,*) 'rk-2, step=27 to call cdiv,fdiv',fdiv
      call cdiv(dtim,ctrdiv,1) 
      write(6,*)  'rk-2 cfcdiv ',cfcdiv,' fdiv ',fdiv,' cdiv ',ctrdiv 
!      write(6,*)  'rk-2  fdiv  ',fdiv,' cdiv  ',ctrdiv                 
!                                                                       
!     3rd R-K step                                                      
!     -----------------                                                 
      dtim= dtf 
!     compute new zc,zf and wz,wx,wy and Jac (since grid has moved)     
      call findzall 
      call sigma 
!      call outflowbc(1,dtim)                                           
!     advecn(m,n,dtim) means use level m and write s,T to level n       
      call advecn(1,0,dtim,step) 
!=      call ecol(1,0,dtim)                                             
!c      call correctbc                                                  
!     call intpol with the time level at which uf is known              
!     besides interpolating the convective terms, intpol computes u*    
      call intpol 
      call rpevalgrad(0) 
!      call rpflatbtm(0)                                                
!c      call evalrp(0)                                                  
!     evalrp must be called before srcface or coriolis                  
!      call viscous(1)                                                  
!c      call gradrp                                                     
!-      call vertmix(1)                                                 
      call coriolis(1) 
      call srcface(1,step) 
!cc      call bcmake(dtim)                                              
      call hsolve(h,oldh,hdt,dtim) 
      call calcskfc 
!c      call hnsolve(h,oldh,hdt,dtim)                                   
      call vhydro(dtim) 
      call cfdiv(cfcdiv) 
!      call ufcorrector                                                 
!     in calling newcor, n should be the same as in coriolis.           
!c      call evalrp(0)                                                  
      call newcor(dtim,1) 
      call newsrc 
!c      call cpsi                                                       
!     compute q                                                         
!      call cpuf                                                        
!c      call exitp(p)                                                   
      edt= EPS/dtim 
      call mgrid(pcorr,dtim,edt,cfcdiv) 
!      call pfilter(pcorr)                                              
!c      call linesolve(p,dtim,step,edt)                                 
!-      call psolve(p,dtim,step,edt)                                    
!c      call pnsolve(p,dtim,step,edt)                                   
!     *********************                                             
      call vface(pcorr,dtim) 
!c      call vnface(p,dtim)                                             
!     computes the vel at the cell faces                                
!     fill f(n=1), use f(n=0)                                           
!     call vcenter with the value of n that we wish to fill             
      call vcenter(pcorr,dtim,0) 
!     computes the final vel (at n+1 -th time step) at the cell centers 
!     Now correct the NH pressure (remains 0 if hy)                     
      if (fnhhy.ne.0) call pcorrect(pcorr) 
      call facediv(dtim,fdiv) 
      call cdiv(dtim,ctrdiv,0) 
      write(6,*)  'rk-3 cfcdiv ',cfcdiv,' fdiv ',fdiv,' cdiv ',ctrdiv 
!      write(6,*)  'rk-3  fdiv  ',fdiv,' cdiv  ',ctrdiv                 
!                                                                       
!     Convective adjustment carried out every 9 time steps.             
!      if (mod(step,9).eq.0) then                                       
      call evalrho(rho,0) 
      call conadjust(0) 
!      end if                                                           
      call calcn2 
!     This was eliminated for non reactive tracer - ONR expts.          
      call tracersource(0,dtim)   
!     call eddydivergence(dtim) - needs cbar
!                                                                       
!=      call resettr(step)                                              
!                                                                       
!     compute and write out extreme values                              
!=      call extremes(step) !open statement in extremes problematic with
                                                                        
      return 
      END                                           
