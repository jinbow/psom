subroutine save3d(Nx,Ny,Nz,var,filename)
   integer :: Nx,Ny,Nz
   character(len=*) :: filename
   real*8 :: var(Nx,Ny,Nz)
   open(3333,file=trim(filename),form='unformatted',access='direct',recl=Nx*Ny*Nz,status='replace')
   write(3333,rec=1) real(var,4)
   close(3333)
end subroutine 

subroutine save2d(Nx,Ny,var,filename)
   !Nx,Ny are not necessary in zonal and meridional directions.
   integer :: Nx,Ny
   character(len=*) :: filename
   real*8 :: var(Nx,Ny)
   open(3333,file=trim(filename),form='unformatted',access='direct',recl=Nx*Ny,status='replace')
   write(3333,rec=1) real(var,4)
   close(3333)
end subroutine 

