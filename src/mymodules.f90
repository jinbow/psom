MODULE dimensions
!  INTEGER,PARAMETER                  :: NI=96,NJ=192,NK=64,ntr=1,nconsume=1  
!  INTEGER, PARAMETER                 :: maxout=1449264
  INTEGER,PARAMETER                  :: NI=96,NJ=192,NK=32,ntr=1,nconsume=1  
  INTEGER, PARAMETER                 :: maxout=750240  ! for NK=32
!  INTEGER,PARAMETER                  :: IMT=96,JMT=192,KM=32
end MODULE dimensions

MODULE rpgrads
  USE dimensions
  REAL*8 :: drpx(NI,NJ,NK),drpy(NI,NJ,NK),                   &
       grpifc(0:NI,NJ,NK),grpjfc(NI,0:NJ,NK)                        
end MODULE rpgrads

MODULE common2tracmass
  USE dimensions
  REAL*8, PARAMETER                  :: LEN= 1.d5, DL=1.d3 
  REAL*8 ux(0:NI+1,0:NJ+1),vy(0:NI+1,0:NJ+1),                        &
       J2d(0:NI+1,0:NJ+1),wz(0:NI+1,0:NJ+1,0:NK+1),                         &
       Jac(0:NI+1,0:NJ+1,0:NK+1),                                           &
       uf(0:NI,NJ,NK), vf(NI,0:NJ,NK),wf(NI,NJ,0:NK) 
end MODULE common2tracmass


MODULE header
  USE dimensions
  USE common2tracmass
  REAL*8, PARAMETER                  :: S0=35.7d0, T0=15.d0 ,R0=1027.d0
  REAL*8, PARAMETER                  :: EPS= 0.1d0,&
       AL=1.d7, FPAR=1.d-4, PI=3.14159265358979323846, OMEGA=7.272d-5, &
       gpr= 0.981,apr=0.6371,dztop=1.d-3,z0= 0.2d0,zm= 0.1d0
  REAL*8                             :: fnhhy,dx,dy,yfront,dyfront,sclwidth, &
       tightness,mldepth,stressmax,pfac,sigrelease(ntr)
  !     sigrelease is the isopycnal of tracer release
  !     dztop is the depth of the top layer non-dim by DL
  !     pfac is the grid stretching factor in z, used in findzall and sigma
  LOGICAL, parameter                :: rect = .true., periodicew = .true.
  INTEGER                           :: conv(0:NI+1,0:NJ+1,0:NK+1), &
       &     con100(0:NI+1,0:NJ+1,0:NK+1)
  REAL*8                            ::  xc(0:NI+1),yc(0:NJ+1),stressx(NJ), &
       &     zc(0:NI+1,0:NJ+1,0:NK+1),zf(0:NI+1,0:NJ+1,-1:NK+1)
  REAL*8                            :: uy(0:NI+1,0:NJ+1), &
       &     vx(0:NI+1,0:NJ+1),                           &
       &     wt(0:NI+1,0:NJ+1,0:NK),wx(0:NI+1,0:NJ+1,0:NK+1),               &
       &     wy(0:NI+1,0:NJ+1,0:NK+1),                                      &
       &     wzk(0:NI+1,0:NJ+1,0:NK),                                       &
       &     ffi(0:NI,NJ),bbi(0:NI,NJ),bbj(NI,0:NJ),                        &
       &     ffj(NI,0:NJ),ffc(0:NI+1,0:NJ+1),bbc(0:NI+1,0:NJ+1),            &
       &     latrad(0:NJ+1)                                                 
  REAL*8                          ::  u(0:NI+1,0:NJ+1,0:NK+1,0:1),          &
       &     v(0:NI+1,0:NJ+1,0:NK+1,0:1),                   &
       &     w(0:NI+1,0:NJ+1,0:NK+1,0:1),                   &
       &     p(0:NI+1,0:NJ+1,0:NK+1),strain(0:NI+1,0:NJ+1,0:NK+1),          &
       &     shear(0:NI+1,0:NJ+1,0:NK+1),                                   &
       &     rho(0:NI+1,0:NJ+1,0:NK+1),vor(0:NI+1,0:NJ+1,0:NK+1),           &
       &     freqN2(0:NI+1,0:NJ+1,0:NK+1),pv(0:NI+1,0:NJ+1,0:NK+1),         &
       &     gradhn(0:NI+1,0:NJ+1,2),hxn(0:NI,NJ,NK),hyn(NI,0:NJ,NK),       &
       &     si(0:NI+1,0:NJ+1,0:NK+1),sj(0:NI+1,0:NJ+1,0:NK+1),             &
       &     sk(0:NI+1,0:NJ+1,0:NK+1),oldh(0:NI+1,0:NJ+1),                  &
       &     sifc(0:NI,NJ,NK),sjfc(NI,0:NJ,NK),skfc(0:NI+1,0:NJ+1,0:NK)     
  REAL*8    ::          consump(0:NI+1,0:NJ+1,0:NK+1,nconsume),             &
       &     s(0:NI+1,0:NJ+1,0:NK+1,0:1),                                   &
       &     T(0:NI+1,0:NJ+1,0:NK+1,0:1),                                   &
       &     Tr(ntr,0:NI+1,0:NJ+1,0:NK+1,0:1),rp(0:NI+1,0:NJ+1,0:NK+1),     &
       &     wtr(ntr,0:NI+1,0:NJ+1,0:NK+1),                                 &
       &     cx(0:NI+1,0:NJ+1,0:NK+1),cy(0:NI+1,0:NJ+1,0:NK+1),             &
       &     cz(0:NI+1,0:NJ+1,0:NK+1),cxf(0:NI,NJ,NK),cyf(NI,0:NJ,NK),      &
       &     czf(NI,NJ,0:NK),                                               &
       &     h(0:NI+1,0:NJ+1),hdt(0:NI+1,0:NJ+1)
  REAL*8    ::   ufbce(NJ,NK),ufbcw(NJ,NK),vfbcn(NI,NK),vfbcs(NI,NK),       &
       &     wfbcb(NI,NJ),ueast(0:NJ+1,NK),uwest(0:NJ+1,NK),                &
       &     vnorth(NI,NK),vsouth(NI,NK),ssouth(NI,NK),sbackgrnd
  REAL*8     trinit(NJ,NK)
  REAL*8 :: D(0:NI+1,0:NJ+1),Ddx(0:NI+1,0:NJ+1),                            &
       &     Ddy(0:NI+1,0:NJ+1),                                            &
       &     g11(0:NI+1,0:NJ+1),g22(0:NI+1,0:NJ+1),                         &
       &     g12(0:NI+1,0:NJ+1),gi(0:NI,NJ,NK,2),gj(NI,0:NJ,NK,2),          &
       &     gi3(0:NI,NJ,NK),gj3(NI,0:NJ,NK),gqi(0:NI,NJ,NK,2),             &
       &     gqj(NI,0:NJ,NK,2),gqi3(0:NI,NJ,NK),                            &
       &     gqj3(NI,0:NJ,NK),gqk(0:NI+1,0:NJ+1,0:NK,3),                    &
       &     Jifc(0:NI,NJ,NK),Jjfc(NI,0:NJ,NK)                             
  REAL*8 ::   beta,dtf,P1,delta,delinv,qpr,lambda,                          &
       &     UL,WL,TL,HL,HDL,kappah,kaphinv
  REAL*8 :: uvis(NI,NJ,NK),vvis(NI,NJ,NK),wvis(NI,NJ,NK),                  &
       &     Kz(NI,NJ,0:NK)
  REAL*8 :: divreyn(NJ,NK),divmean(NJ,NK),dcdt(NJ,NK),prod(NJ,NK)
  REAL*8 :: fricu(NI,NJ,NK),fricv(NI,NJ,NK),fricw(NI,NJ,NK),               &
       &     fricb(NI,NJ,NK),rhoadv(NI,NJ,NK),rhoprev(NI,NJ,NK),           &
       &     mldn2(3),mldn2init(3),zbtop(3),zbbot(3),frictop(3),diatop(3), &
       &     zbtopinit(3),zbbotinit(3),advecpv(3),friction(3),             &
       &     diabatic(3),diabot(3),fricbot(3)
  integer ::  ktest(3)

end MODULE header

