      subroutine facediv(dtime,maxdiv) 
!     ------------------                                                
!     checks divergence                                                 
!     wf does not contain  J*wt and is not exactly the contravariant vel
      USE header
      integer i,j,k,imax,jmax,kmax 
      double precision ufdx,vfdy,wfdz,div,maxdiv,dtime 
!     edd is a scaling factor that makes this divergence                
!     the same as the residual from the pressure-Poisson equation.      
!     This is a check on the code.                                      
!                                                                       
!      double precision dtheta,dphi,phi0,cosphi0,sinphi0,C,dx,dy,dz,    
!     & lon0,lat0                                                       
!      C= PI/180.d0                                                     
!c                                                                      
!      lon0= 260.d0*C                                                   
!      lat0= 20.d0*C                                                    
!      dtheta= C/8.d0                                                   
!      dphi= C/8.d0                                                     
!      phi0= lat0 + 0.5d0*(dphi*dfloat(NJ))                             
!      cosphi0= dCos(phi0)                                              
!      sinphi0= dSin(phi0)                                              
!      write(6,*) dphi,dtheta,phi0,cosphi0,sinphi0                      
!c     for the rectilinear case, dx,dy,dz are constant                  
!      dx= R*cosphi0*dtheta /LEN                                        
!      dy= R*dphi/LEN                                                   
!      dz= 1000.d0/(dfloat(NK)*DL)                                      
!      write(6,*) 'dx dy',dx,dy                                         
!                                                                       
!      edd= EPS/(dtime*delta*kappa)                                     
!     We have now changed the defn of tol so that it is just (u_x +v_y +
!     We do not need edd any more.                                      
!      edd= EPS*delinv/dtime                                            
      maxdiv= 0.d0 
      do 10 k=1,NK 
         do 20 j=1,NJ 
            do 30 i=1,NI 
               ufdx= (uf(i,j,k)-uf(i-1,j,k)) 
               vfdy= (vf(i,j,k)-vf(i,j-1,k)) 
               wfdz= (wf(i,j,k)-wf(i,j,k-1)) 
!               div= dabs((ufdx+ vfdy + wfdz)*edd)                      
               div= dabs(ufdx+ vfdy + wfdz) 
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
