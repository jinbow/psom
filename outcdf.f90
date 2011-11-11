subroutine outcdf(xc,yc,zc,p,h,consump,s,temp,rho,Tr,u,v,w,          &
     vor,pv,uf,vf,wf,Kz,conv,con100,step)                         
  !----------------------------------------------------------             
  !     outputs the solution as a netCDF file                             
  USE dimensions
  integer i,j,k,step,nstp,idigit,ipos,it 
  integer conv(0:NI+1,0:NJ+1,0:NK+1),con100(0:NI+1,0:NJ+1,0:NK+1) 
  double precision rho(0:NI+1,0:NJ+1,0:NK+1),                       &
       p(0:NI+1,0:NJ+1,0:NK+1),u(0:NI+1,0:NJ+1,0:NK+1),             &
       v(0:NI+1,0:NJ+1,0:NK+1),w(0:NI+1,0:NJ+1,0:NK+1),             &
       vor(0:NI+1,0:NJ+1,0:NK+1),pv(0:NI+1,0:NJ+1,0:NK+1),          &
       s(0:NI+1,0:NJ+1,0:NK+1),temp(0:NI+1,0:NJ+1,0:NK+1),          &
       Tr(ntr,0:NI+1,0:NJ+1,0:NK+1),                                &
       Trwrite(0:NI+1,0:NJ+1,0:NK+1,ntr),                           &
       h(0:NI+1,0:NJ+1),uf(0:NI,NJ,NK),                             &
       vf(NI,0:NJ,NK),wf(NI,NJ,0:NK),Kz(NI,NJ,0:NK),                &
       consump(0:NI+1,0:NJ+1,0:NK+1,nconsume)                       
  double precision xc(0:NI+1),yc(0:NJ+1),                           &
       zc(0:NI+1,0:NJ+1,0:NK+1)                                     
  real*8 z(0:NI+1,0:NJ+1,0:NK+1) 
  real*8 smax,Tmax,umax,vmax,wmax,smin,Tmin,umin,vmin,wmin 
  character  outname*12, facename*12 
                                                                        
  INCLUDE '/usr/local/netcdf-3.6.2/include/netcdf.inc' 
!      INCLUDE '/sw/lib/netcdf-gfortran/include/netcdf.inc'             
  
  integer start(3), count(3), start2d(2), count2d(2),count2dtr(3),  &
       countuf(3), countvf(3), countwf(3),count4(4),start4(4),      &
       start2dtr(3),count4consump(4)                                
  integer dims(3), dims2d(2),dims2dtr(3),dimuf(3), dimvf(3),        &
       dimwf(3),dims4(4),dimsconsump(4)                             
                                                                        
                                                                        
  DATA start /1, 1, 1/ 
  DATA start4 /1, 1, 1, 1/ 
  DATA start2d /1, 1/ 
  DATA start2dtr /1, 1, 1/ 
                                                                        
                                                                        
!      do j=0,NJ+1                                                      
!         do i=0,NI+1                                                   
!            x(i,j)= xc(i)                                              
!            y(i,j)= yc(j)                                              
!         end do                                                        
!      end do                                                           
  do k=0,NK+1 
     do j=0,NJ+1 
        do i=0,NI+1 
           z(i,j,k)= zc(i,j,k) 
           do it=1,ntr 
              Trwrite(i,j,k,it)= Tr(it,i,j,k) 
           end do
        end do
     end do
  end do
                                                                        
!     FOR THE SAKE OF BETTER PLOTS                                      
  do j=0,NJ+1 
     do i=0,NI+1 
        z(i,j,NK+1)= 0.d0 
        z(i,j,0)= 0.5*(z(i,j,0)+z(i,j,1)) 
        s(i,j,NK+1)= s(i,j,NK) 
        s(i,j,0)= s(i,j,1) 
        temp(i,j,NK+1)= temp(i,j,NK) 
        temp(i,j,0)= temp(i,j,1) 
        rho(i,j,NK+1)= rho(i,j,NK) 
        rho(i,j,0)= rho(i,j,1) 
        u(i,j,NK+1)= u(i,j,NK) 
        u(i,j,0)= u(i,j,1) 
        v(i,j,NK+1)= v(i,j,NK) 
        v(i,j,0)= v(i,j,1) 
        w(i,j,NK+1)= w(i,j,NK) 
        w(i,j,0)= w(i,j,1) 
        vor(i,j,NK+1)= vor(i,j,NK) 
        vor(i,j,0)= vor(i,j,1) 
        pv(i,j,NK+1)= pv(i,j,NK) 
        pv(i,j,0)= pv(i,j,1) 
        do it=1,ntr 
           Trwrite(i,j,NK+1,it)= Trwrite(i,j,NK,it) 
           Trwrite(i,j,0,it)= Trwrite(i,j,1,it) 
        end do
     end do
  end do
                                                                        
      count(1)= NI+2 
      count(2)= NJ+2 
      count(3)= NK+2 
                                                                        
      count4(1)= NI+2 
      count4(2)= NJ+2 
      count4(3)= NK+2 
      count4(4)= ntr 
                                                                        
      count4consump(1)= NI+2 
      count4consump(2)= NJ+2 
      count4consump(3)= NK+2 
      count4consump(4)= nconsume 
                                                                        
      count2d(1)= NI+2 
      count2d(2)= NJ+2 
                                                                        
      count2dtr(1)= NI+2 
      count2dtr(2)= NJ+2 
      count2dtr(3)= ntr 
                                                                        
      countuf(1)= NI+1 
      countuf(2)= NJ 
      countuf(3)= NK 
                                                                        
      countvf(1)= NI 
      countvf(2)= NJ+1 
      countvf(3)= NK 
                                                                        
      countwf(1)= NI 
      countwf(2)= NJ 
      countwf(3)= NK+1 
                                                                        
!     call evalrho(rho,0)  should be called from outside
                                                                        
!     name the output file                                              
!                                                                       
      outname = 'outxxxxx.cdf' 
!      facename = 'vfcxxxxx.cdf'                                        
      facename = 'vfc.cdf' 
      nstp = step 
      do  ipos = 8, 4, -1
	 idigit = mod(nstp,10)
	 outname(ipos:ipos) = char(ichar('0') + idigit)
!	 facename(ipos:ipos) = char(ichar('0') + idigit)
	 nstp = nstp / 10
      end do

!      write(6,*) 'in outcdf, outname',outname
!      write(6,*) 'in outcdf, facename',facename

      idDatFile =  nccre(outname,NCCLOB,rcode)      

      dims(1) = ncddef(idDatFile,'x',NI+2,rcode)
      dims(2) = ncddef(idDatFile,'y',NJ+2,rcode)
      dims(3) = ncddef(idDatFile,'sigma',NK+2,rcode)

      dims4(1) = dims(1)
      dims4(2) = dims(2)
      dims4(3) = dims(3)
      dims4(4) = ncddef(idDatFile,'ntr',ntr,rcode)

      dimsconsump(1) = dims(1)
      dimsconsump(2) = dims(2)
      dimsconsump(3) = dims(3)
      dimsconsump(4) = ncddef(idDatFile,'ntrcon',nconsume,rcode)

      dims2d(1)= dims(1)
      dims2d(2)= dims(2)

      dims2dtr(1)= dims(1)
      dims2dtr(2)= dims(2)
      dims2dtr(3)= dims4(4)

      idvx = ncvdef(idDatFile,'xc',NCDOUBLE,1,dims(1),rcode)
      idvy = ncvdef(idDatFile,'yc',NCDOUBLE,1,dims(2),rcode)
      idvz = ncvdef(idDatFile,'zc',NCDOUBLE,3,dims,rcode)
      idvh = ncvdef(idDatFile,'h',NCDOUBLE,2,dims2d,rcode)
      idvc = ncvdef(idDatFile,'consump',NCDOUBLE,4,dimsconsump,rcode)
      idvtr = ncvdef(idDatFile,'tr',NCDOUBLE,4,dims4,rcode)
      idvs = ncvdef(idDatFile,'s',NCDOUBLE,3,dims,rcode)
      idvt = ncvdef(idDatFile,'temp',NCDOUBLE,3,dims,rcode)
      idvrho = ncvdef(idDatFile,'rho',NCDOUBLE,3,dims,rcode)
      idvp = ncvdef(idDatFile,'p',NCDOUBLE,3,dims,rcode)
      idvu = ncvdef(idDatFile,'u',NCDOUBLE,3,dims,rcode)
      idvv = ncvdef(idDatFile,'v',NCDOUBLE,3,dims,rcode)
      idvw = ncvdef(idDatFile,'w',NCDOUBLE,3,dims,rcode)
      idvz3 = ncvdef(idDatFile,'vor',NCDOUBLE,3,dims,rcode)
      idvpv = ncvdef(idDatFile,'pv',NCDOUBLE,3,dims,rcode)
      idvcon = ncvdef(idDatFile,'conv',NCshort,3,dims,rcode)
      idvcon100 = ncvdef(idDatFile,'con100',NCshort,3,dims,rcode)

      CALL ncendf(idDatFile,rcode)
!      call calc_dvdy(dvdy,0)


      CALL ncvpt(idDatFile,idvx, start(1), count(1), xc, rcode)
      CALL ncvpt(idDatFile,idvy, start(2), count(2), yc, rcode)
      CALL ncvpt(idDatFile,idvz, start, count, z, rcode)
      CALL ncvpt(idDatFile,idvh, start2d, count2d, h, rcode)
      CALL ncvpt(idDatFile,idvc,start4,count4consump,consump,rcode)
      CALL ncvpt(idDatFile,idvtr, start4, count4, Trwrite, rcode)
      CALL ncvpt(idDatFile,idvs, start, count, s, rcode)
      CALL ncvpt(idDatFile,idvt, start, count, temp, rcode)
      CALL ncvpt(idDatFile,idvrho, start, count, rho, rcode)
      CALL ncvpt(idDatFile,idvp, start, count, p, rcode)
      CALL ncvpt(idDatFile,idvu, start, count, u, rcode)
      CALL ncvpt(idDatFile,idvv, start, count, v, rcode)
      CALL ncvpt(idDatFile,idvw, start, count, w, rcode)
      CALL ncvpt(idDatFile,idvz3, start, count, vor, rcode)
      CALL ncvpt(idDatFile,idvpv, start, count, pv, rcode)
      CALL ncvpt(idDatFile,idvcon, start, count, conv, rcode)
      CALL ncvpt(idDatFile,idvcon100, start, count, con100, rcode)

      CALL ncclos(idDatFile,rcode)

!     ----------------------------------------------------------------
!     write face velocities, uf,vf,wf
!     -------------------------------

      idFaceFile =  nccre(facename,NCCLOB,rcode)      

      dimuf(1) = ncddef(idFaceFile,'xi-ew',NI+1,rcode)
      dimuf(2) = ncddef(idFaceFile,'eta-ew',NJ,rcode)
      dimuf(3) = ncddef(idFaceFile,'sigma-ew',NK,rcode)

      dimvf(1) = ncddef(idFaceFile,'xi-ns',NI,rcode)
      dimvf(2) = ncddef(idFaceFile,'eta-ns',NJ+1,rcode)
      dimvf(3) = ncddef(idFaceFile,'sigma-ns',NK,rcode)

      dimwf(1) = ncddef(idFaceFile,'xi-tb',NI,rcode)
      dimwf(2) = ncddef(idFaceFile,'eta-tb',NJ,rcode)
      dimwf(3) = ncddef(idFaceFile,'sigma-tb',NK+1,rcode)

      idvuf = ncvdef(idFaceFile,'uf',NCDOUBLE,3,dimuf,rcode)
      idvvf = ncvdef(idFaceFile,'vf',NCDOUBLE,3,dimvf,rcode)
      idvwf = ncvdef(idFaceFile,'wf',NCDOUBLE,3,dimwf,rcode)
      idvKf = ncvdef(idFaceFile,'Kz',NCDOUBLE,3,dimwf,rcode)

      CALL ncendf(idFaceFile,rcode)

      CALL ncvpt(idFaceFile,idvuf, start, countuf, uf, rcode)
      CALL ncvpt(idFaceFile,idvvf, start, countvf, vf, rcode)
      CALL ncvpt(idFaceFile,idvwf, start, countwf, wf, rcode)
      CALL ncvpt(idFaceFile,idvKf, start, countwf, Kz, rcode)

      CALL ncclos(idFaceFile,rcode)

      return
      end
      
