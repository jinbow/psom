MODULE mod_param
 INTEGER                                   :: IMT, JMT, KM
! INTEGER, POINTER                          :: IMT, JMT, KM
 INTEGER                                   :: JMAX, LBT
 INTEGER, PARAMETER                        :: MR=1001
 INTEGER                                   :: NEND
 INTEGER, PARAMETER                        :: NST=2,NNRJ=8,NTRJ=7
 INTEGER, PARAMETER                        :: NTRACmax= 1 
 REAL*4, ALLOCATABLE, DIMENSION(:,:)       :: trj
 INTEGER, ALLOCATABLE, DIMENSION(:,:)       :: nrj

 INTEGER                                   :: ncoor,kriva,iter,ngcm
 REAL*8                                    :: tseas,tday,tyear,dtmin,voltr
 REAL*8                                    :: tstep,dstep,tss,partQuant
 REAL*8, PARAMETER                         :: UNDEF=1.d20 
 REAL*8, PARAMETER                         :: PI =3.14159265358979323846d0
ENDMODULE mod_param

MODULE mod_grid
  REAL*4, ALLOCATABLE, DIMENSION(:,:)       :: dxv, dyu
  REAL*8, ALLOCATABLE, DIMENSION(:)         :: dz
  REAL*8, ALLOCATABLE, DIMENSION(:,:)       :: dxdy
  REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: dvol
  REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: dzt
  
  REAL*8                                    :: rmin ,tmin ,smin
  REAL*8                                    :: dr ,dtemp ,dsalt
  REAL*8                                    :: arcscale
  INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: kmt
  INTEGER                                   :: subGrid     ,subGridID
  INTEGER                                   :: subGridImin ,subGridImax
  INTEGER                                   :: subGridJmin ,subGridJmax
  CHARACTER(LEN=200)                        :: SubGridFile 
ENDMODULE mod_grid

MODULE mod_vel
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:)    :: uflux ,vflux
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)    :: wflux
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)      :: hs  
  REAL*8                                     :: ff
ENDMODULE mod_vel
