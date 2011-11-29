subroutine stprofile 
  !     --------------------                                              
!  USE header
  ! Initializes s,T from Andrey's section gridded for model grid
  ! s_init.dat and T_init.dat  written in ascii - NK+2 x NJx2
!  implicit none
!  integer n,i,j,k

  !     --------------------                                              
  USE header, ONLY: NI,NJ,NK,rho,s,T,pi,zc,DL,LEN,yc,dy
  implicit none
  INTEGER n,i,j,k,Nu,k0
  REAL*8 :: A,B,G,C,potdens,y,z,tmp,yfront,mldepth,width,slfac
  INTEGER,PARAMETER :: seed = 86456  
  CALL RANDOM_SEED()

  ! assign the density or temperature and salinity by either analytic functions or 
  ! any particular hydrographic section.

  ! t=A/B exp(A*z)erf(B*y) + G*(1+z)
  ! z=-1:0, y = -0.5,0.5, 

  ! ===========================
  ! cases 
  ! 0: jet
  ! 1: idealized surface jet
  ! 2: krushio
  SELECT CASE (1)
  CASE(1)
     A = 4.0d0
     B = 10d3 !jet has 2B width
     C = 1d0
     G = 10d0
     !Nu =NK-9
     mldepth = 70d0
     s(:,:,:,0) = 34d0
     DO k = NK, 0, -1
        DO j = 0, NJ+1
           DO i = 0, NI+1
              CALL RANDOM_NUMBER(tmp)
              !y = -(REAL(j)-REAL(NJ+1)/2d0)/REAL(NJ+1)
              y = -(dble(j)-dble(NJ+1)/2d0)*dy
              !IF (k>Nu) THEN
              IF (zc(1,1,k)*DL>-mldepth) THEN
                 z=0d0
                 Nu=k
              ELSE
                 z = dble(k-Nu)/dble(Nu)
                 z = zc(1,1,k)*DL+mldepth
              ENDIF
              T(i,j,k,0) = C * EXP(z/200d0) * erf(y/B) + G*(1+z/1000d0) !+ 0.02*tmp ! + 1d-2*SIN(DBLE(i)/DBLE(NI+1)*2d0*pi*1d0)
              !T(i,j,k,0) = C * EXP(A*z) * erf(B*y) + G*(1+z) + 1d-3*SIN(DBLE(i)/DBLE(NI+1)*2d0*pi*3d0)
              rho(i,j,k) = potdens(s(0,j,k,0),T(0,j,k,0))
              if (k == 0) then
                  T(i,j,k,0) = C * EXP(z/200d0)  + G*(1+z/1000d0) !+ 0.02*tmp ! + 1d-2*SIN(DBLE(i)/DBLE(NI+1)*2d0*pi*1d0)
                  !T(i,j,k,0) = C * EXP(A*z) * erf(B*y) + G*(1+z) + 1d-3*SIN(DBLE(i)/DBLE(NI+1)*2d0*pi*3d0)
                  rho(i,j,k) = potdens(s(0,j,k,0),T(0,j,k,0))
              endif
           ENDDO

        ENDDO
     ENDDO
     T(:,0,:,0)= T(:,1,:,0)
     T(:,NJ+1,:,0)= T(:,NJ,:,0)
     s(:,0,:,0)= s(:,1,:,0)
     s(:,NJ+1,:,0)= s(:,NJ,:,0)
     rho(:,0,:)= rho(:,1,:)
     rho(:,NJ+1,:)= rho(:,NJ,:)
     T(:,:,NK+1,:)=T(:,:,NK,:)
     s(:,:,NK+1,:)=s(:,:,NK,:)
     rho(:,:,NK+1)=rho(:,:,NK)
     
     DO k = 1, NK
        rho(:,:,k) = 0.25d0*(rho(:,:,k-1)+2d0*rho(:,:,k)+rho(:,:,k+1))
        s(:,:,k,0) = 0.25d0*(s(:,:,k-1,0)+2d0*s(:,:,k,0)+s(:,:,k+1,0))
     ENDDO
     call save3d(NI+2,NJ+2,NK+2,rho(:,:,:),'initial-rho.bin')
  END SELECT

END SUBROUTINE stprofile
