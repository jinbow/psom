subroutine stprofile 
  !     --------------------                                              
  USE header
!     initializes s as pot.density with a single vertical profile       
!     Across-front use the TANH profile with tightness as a measure of f
!     spread                                                            
!     Larger the factor, tighter the front and larger b_xx.             
!     tightness=4, gives the needed b_xx to make b_xx H/ (2f^2) < 1.    
  implicit none
  INTEGER, PARAMETER                 :: np=100
!  REAL*8 :: bs(np),cs(np),ds(np), bt(np),ct(np),dt(np),                 &
!       xs(np),ys(np),xt(np),yt(np)
  integer  i,j,k,it,iseed,npm1,n
  double precision rhovert(np),Tvert(np),dep(np),depoff(np) 
  double precision bs1(np),cs1(np),ds1(np),bT1(np),cT1(np),dT1(np), &
       sal(np),temp(np),z,seval,zmm,sbkgrnd,z1,z2,potdens
  double precision slfac,dum,dscl,rdscl,yarg,ex2y,thy,ran3,         &
       perturb,slfacnew,zoffset,dz,bfsqbkgrnd                       
!     tightness = 10 represents a very tight front, =1 loose front      
!=      parameter (zoffset= 200.)  also used zoffset=50 ! if zoffset=0 m
!=       parameter (zoffset= 200.d0, tightness=0.03)                    
!=       parameter (zoffset= -10.d0)                                    
!       parameter (zoffset= -20.d0) 
  parameter (zoffset= 0.d0) 
  open(unit=40,file='profile_sTarctic.in')
  do i=np,1,-1
     read(40,*) dep(i),temp(i),sal(i)
     dep(i)= -dep(i)
  end do
  close(40)
  !     The z values must be in STRICTLY INCREASING ORDER and hence should
  !     be written from ocean bottom up.                                  
                                                                        
  n=0
  call spline (np,dep,sal,bs1,cs1,ds1) 
  call spline (np,dep,temp,bt1,ct1,dt1) 
  do k=0,NK+1 
     do j=0,NJ+1 
        i=1
        z= DL*zc(i,j,k) 
        s(i,j,k,n)= seval(np,z,dep,sal,bs1,cs1,ds1)   
        T(i,j,k,n)= seval(np,z,dep,temp,bt1,ct1,dt1) 
        rho(i,j,k)= potdens(s(i,j,k,n),T(i,j,k,n))
        do i=2,NI+1
           s(i,j,k,n)=s(1,j,k,n)
           T(i,j,k,n)=T(1,j,k,n)
           rho(i,j,k)= rho(1,j,k)
        end do
     end do
  end do
!  call evalrho(rho,n)

  return 
END subroutine stprofile
