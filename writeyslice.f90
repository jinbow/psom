subroutine writeyslice(frame_int,step,consump,T,rho,u,v,w,vor,    &
     xc,zc)                                                       
  !     ------------------------------------------------------------      
  USE dimensions
  !     reads in a netcdf file and interpolates the data from the sigma gr
  !     onto a z level                                                    
                                                                        
      INCLUDE '/usr/local/netcdf-3.6.2/include/netcdf.inc' 
      integer NI2,NJ2,NK2,i,j,k,step,frame_int 
      parameter ( NI2=NI+2,NJ2=NJ+2,NK2=NK+2) 
      double precision DL 
      parameter (DL= 1000.d0) 
      double precision svert(NK),Tvert(NK),rvert(NK),uvert(NK),         &
     &     vvert(NK),wvert(NK),zvert(NK),seval                          
                                                                        
      double precision                                                  &
     &     rho(0:NI+1,0:NJ+1,0:NK+1),u(0:NI+1,0:NJ+1,0:NK+1),           &
     &     v(0:NI+1,0:NJ+1,0:NK+1),w(0:NI+1,0:NJ+1,0:NK+1),             &
     &     consump(0:NI+1,0:NJ+1,0:NK+1,ntr),                           &
     &     T(ntr,0:NI+1,0:NJ+1,0:NK+1),                                 &
     &     vor(0:NI+1,0:NJ+1,0:NK+1),                                   &
     &     zc(0:NI+1,0:NJ+1,0:NK+1),                                    &
     &     xc(0:NJ+1),bs1(NK),cs1(NK),ds1(NK)                           
      real*4 xslice(NI),zslice(NI,0:NK+1),sslice(ntr,NI,0:NK+1),        &
     &     Tslice(ntr,NI,0:NK+1),rslice(NI,0:NK+1),uslice(NI,0:NK+1),   &
     &     vslice(NI,0:NK+1),wslice(NI,0:NK+1),vorslice(NI,0:NK+1)      
!                                                                       
                                                                        
      character*80 outname 
                                                                        
      integer start(3),count(3),dims(3),dims4(4),start2d(2),count4(4),  &
     &     start4(4)                                                    
                                                                        
                                                                        
      j= NJ/2 +1 
      do k=1,NK 
         do i=1,NI 
            zslice(i,k)= zc(i,j,k)*DL 
            xslice(i)= xc(i) 
                                                                        
!=            sslice(i,k)= s(i,j,k) +S0                                 
!=            Tslice(i,k)= T(i,j,k) +T0                                 
            do itr=1,ntr 
               sslice(itr,i,k)= consump(i,j,k,itr) 
               Tslice(itr,i,k)= T(itr,i,j,k) 
            end do 
            rslice(i,k)= rho(i,j,k) 
            uslice(i,k)= u(i,j,k) 
            vslice(i,k)= v(i,j,k) 
            wslice(i,k)= w(i,j,k) 
            vorslice(i,k)= vor(i,j,k) 
                                                                        
         end do 
      end do 
!     FOR THE SAKE OF BETTER PLOTS                                      
      do i=1,NI 
         zslice(i,NK+1)= 0.d0 
         zslice(i,0)= 0.5*DL*(zc(i,j,0)+zc(i,j,1)) 
         rslice(i,NK+1)= rslice(i,NK) 
         rslice(i,0)= rslice(i,1) 
         uslice(i,NK+1)= uslice(i,NK) 
         uslice(i,0)= 2.0*uslice(i,1) -uslice(i,2) 
         vslice(i,NK+1)= 2.0*vslice(i,NK) -vslice(i,NK-1) 
         vslice(i,0)= 2.0*vslice(i,1) -vslice(i,2) 
         wslice(i,NK+1)= 2.0*wslice(i,NK) -wslice(i,NK-1) 
         wslice(i,0)= 0.0 
         vorslice(i,NK+1)= 2.0*vorslice(i,NK) -vorslice(i,NK-1) 
         vorslice(i,0)= 2.0*vorslice(i,1) -vorslice(i,2) 
         do itr=1,ntr 
            sslice(itr,i,NK+1)= sslice(itr,i,NK) 
            sslice(itr,i,0)= sslice(itr,i,1) 
            Tslice(itr,i,NK+1)= Tslice(itr,i,NK) 
            Tslice(itr,i,0)= Tslice(itr,i,1) 
         end do 
      end do 
                                                                        
!     write to a netcdf file                                            
                                                                        
      outname= 'yslicemovie.cdf' 
                                                                        
      if (step.eq.0) then 
         idDatFile =  nccre(outname,NCCLOB,rcode) 
                                                                        
         dims(1) = ncddef(idDatFile,'x',NI,rcode) 
         dims(2) = ncddef(idDatFile,'z',NK+2,rcode) 
         dims(3) = ncddef(idDatFile,'time',NCUNLIM,rcode) 
                                                                        
         dims4(1) = dims(1) 
         dims4(2) = ncddef(idDatFile,'ntr',ntr,rcode) 
         dims4(3) = dims(2) 
         dims4(4) = dims(3) 
                                                                        
         idvy = ncvdef(idDatFile,'xc',NCFLOAT,1,dims(1),rcode) 
         idvz = ncvdef(idDatFile,'zc',NCFLOAT,2,dims,rcode) 
         idvs = ncvdef(idDatFile,'consump',NCFLOAT,4,dims4,rcode) 
!=         idvt = ncvdef(idDatFile,'temp',NCFLOAT,3,dims,rcode)         
         idvt = ncvdef(idDatFile,'tr',NCFLOAT,4,dims4,rcode) 
         idvrho = ncvdef(idDatFile,'rho',NCFLOAT,3,dims,rcode) 
         idvu = ncvdef(idDatFile,'u',NCFLOAT,3,dims,rcode) 
         idvv = ncvdef(idDatFile,'v',NCFLOAT,3,dims,rcode) 
         idvw = ncvdef(idDatFile,'w',NCFLOAT,3,dims,rcode) 
         idvvor = ncvdef(idDatFile,'vor',NCFLOAT,3,dims,rcode) 
                                                                        
         CALL ncendf(idDatFile,rcode) 
                                                                        
      else 
         idDatFile = ncopn(outname, NCWRITE,rcode) 
!         call ncredf(idDatFile,rcode)                                  
         idvs = NCVID(idDatFile, 'consump', RCODE) 
!=         idvt = NCVID(idDatFile, 'temp', RCODE)                       
         idvt = NCVID(idDatFile, 'tr', RCODE) 
         idvrho = NCVID(idDatFile, 'rho', RCODE) 
         idvu = NCVID(idDatFile, 'u', RCODE) 
         idvv = NCVID(idDatFile, 'v', RCODE) 
         idvw = NCVID(idDatFile, 'w', RCODE) 
         idvvor = NCVID(idDatFile, 'vor', RCODE) 
      endif 
                                                                        
      count(1)= NI 
      count(2)= NK+2 
      count(3)= 1 
                                                                        
      count4(1)= NI 
      count4(2)= ntr 
      count4(3)= NK+2 
      count4(4)= 1 
                                                                        
      start2d(1)= 1 
      start2d(2)= 1 
                                                                        
      start(1)= 1 
      start(2)= 1 
      start(3)= step/frame_int +1 
                                                                        
      start4(1)= 1 
      start4(2)= 1 
      start4(3)= 1 
      start4(4)= step/frame_int +1 
                                                                        
      if (step.eq.0) then 
         CALL ncvpt(idDatFile,idvy, start(1), count(1), xslice, rcode) 
         CALL ncvpt(idDatFile,idvz, start, count, zslice, rcode) 
      endif 
      CALL ncvpt(idDatFile,idvs, start4, count4, sslice, rcode) 
!=      CALL ncvpt(idDatFile,idvT, start, count, Tslice, rcode)         
      CALL ncvpt(idDatFile,idvT, start4, count4, Tslice, rcode) 
      CALL ncvpt(idDatFile,idvrho, start, count, rslice, rcode) 
      CALL ncvpt(idDatFile,idvu, start, count, uslice, rcode) 
      CALL ncvpt(idDatFile,idvv, start, count, vslice, rcode) 
      CALL ncvpt(idDatFile,idvw, start, count, wslice, rcode) 
      CALL ncvpt(idDatFile,idvvor, start, count, vorslice, rcode) 
                                                                        
      CALL ncclos(idDatFile,rcode) 
                                                                        
      END                                           
