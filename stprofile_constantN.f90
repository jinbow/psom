subroutine stprofile 
  !     --------------------                                              
  USE header, ONLY: NI,NJ,NK,rho,s,T,pi,zc,DL,LEN,yc,dy
  implicit none
  INTEGER n,i,j,k,Nu,k0
  REAL*8 :: A,B,G,C,potdens,y,z,tmp,yfront,mldepth,width,slfac
  INTEGER,PARAMETER :: seed = 86456  

  SELECT CASE (1)
  CASE(1)
     s(:,:,:,0) = 34d0
     DO k = NK+1, 0, -1
     z = zc(1,1,k)*DL
              T(:,:,k,0) = 29d0 - z/800d0 * 10d0
     ENDDO
     do i = 0, NI+1
        do j = 0, NJ+1
           do k = 0, NK+1
              rho(i,j,k) = potdens(s(i,j,k,0),T(i,j,k,0))
           enddo
        enddo
     enddo
     call save3d(NI+2,NJ+2,NK+2,rho(:,:,:),'initial-rho.bin')
  END SELECT

END SUBROUTINE stprofile
