
CASE         = gulf


#F90COMPILER  = "ifort"
F90COMPILER  = "gfortran"

#================================================================

NC_LIB           = $(shell nc-config --flibs)
NC_INC           = $(shell nc-config --fflags)
#TRACMASS_FLAGS    = -Dzgrid3d -Druntracmass -Dfull_wflux
#TRACMASS_FLAGS    = -Dzgrid3d  -Dfull_wflux

#================================================================

CASE_FLAG         = -DCASE_NAME=\'$(CASE)\'
ARG_FLAGS         = -DINARG_1=KALLE

ifeq ($(F90COMPILER),"ifort")
        FFLAGS = -fast -fpp
        CFLAGS = -O3 -c
        FC = ifort
endif
ifeq ($(F90COMPILER),"gfortran")
	CFLAGS = -c -x f95-cpp-input -fconvert=big-endian -gdwarf-2   
	FFLAGS =-fno-underscoring  
	FC = gfortran #$(LIB_DIR) $(INC_DIR) $(F90_FLAGS) $(ORM_FLAGS)
endif

VPATH = src

%.o : %.f90
	$(FC) $(FFLAGS) $(CFLAGS) $(CASE_FLAG) $(ARG_FLAGS) $(TRACMASS_FLAGS) $< -o $@

 NHO  = mymodules.o  main.o init_tr.o  momentum.o potdens.o facediv.o intpol_3rdorder.o srcface_3rdorder.o vcenter.o vface.o advect_quick.o cdiv.o cpcors.o cpfine.o efill.o mgrid.o prolong.o resid.o restrict.o mgpfill.o pbc.o sor.o spline.o seval.o sigma_toplayer.o staticsigma.o hsolve.o hfill.o chfine.o hbc.o vhydro.o evalrho.o psolve.o uvchy.o newcor.o newsrc.o coriolis.o hsave.o smooth.o setbc.o linerelax_periodicew.o dgtsl.o energy.o findz_topmoves.o meanh.o rpflattopog.o correctbc.o surfaceflux.o conadjust_sT.o diffusion_wind.o extremes_tr.o checks.o biharmonic.o  density.o geostroph.o ran3.o tracerinit_nut.o mprove.o dgtsv.o pcorrect.o stprofile_deepfront_jbw.o inith_thermalwind.o vort.o   calcskfc.o cfdiv.o pfilter.o viscous.o prepvisc.o srcface_nopy.o tracersource_nut.o windstress.o calc_dvdy.o streamfunction.o n2budget_topbot.o calcn2.o pvcalc.o   analytic_eval.o sub_io.o

 NHTRMASSO  = mymodules.o tracmass_modules.o tracmass_interface.o tracmass.o seed_particles.o time_subs.o main.o init_tr.o outcdf.o momentum.o potdens.o facediv.o intpol_3rdorder.o srcface_3rdorder.o vcenter.o vface.o advect_quick.o cdiv.o cpcors.o cpfine.o efill.o mgrid.o prolong.o resid.o restrict.o mgpfill.o pbc.o sor.o spline.o seval.o sigma_toplayer.o staticsigma.o hsolve.o hfill.o chfine.o hbc.o vhydro.o evalrho_rho.o psolve.o uvchy.o newcor.o newsrc.o coriolis.o hsave.o smooth.o setbc.o linerelax_periodicew.o dgtsl.o energy.o findz_topmoves.o meanh.o rgrad_song_bndry0dsdy.o correctbc.o surfaceflux.o conadjust_rho.o diffusion_wind.o extremes_tr.o checks.o biharmonic.o writeslice.o density.o geostroph.o ran3.o tracerinit_dots.o mprove.o dgtsv.o pcorrect.o stprofile_lightintrusion.o topog.o inith_fixdh.o vort.o writeksurf.o writeyslice.o calcskfc.o cfdiv.o pfilter.o viscous.o prepvisc.o srcface_nopy.o tracersource_isopyc.o windstress.o calc_dvdy.o streamfunction.o n2budget_topbot.o calcn2.o pvcalc.o writen2budget.o outpv.o writeisopycnal.o analytic_eval.o tracerrelease_dots.o

 NHCNO  =  maincn.o initcn_tr.o outcdf.o momentum.o potdens.o facediv.o intpol_3rdorder.o srcface_3rdorder.o vcenter.o vface.o advect_quick.o cdiv.o cpcors.o cpfine.o efill.o mgrid.o prolong.o resid.o restrict.o mgpfill.o pbc.o sor.o spline.o seval.o sigma_toplayer.o staticsigma.o hsolve.o hfill.o chfine.o hbc.o vhydro.o evalrho_rho.o psolve.o uvchy.o newcor.o newsrc.o coriolis.o hsave.o smooth.o setbc.o linerelax_periodicew.o dgtsl.o energy.o findz_topmoves.o meanh.o rgrad_song_bndry0dsdy.o correctbc.o surfaceflux.o conadjust_tr.o diffusion_wind.o extremes_tr.o checks.o biharmonic.o writeframe.o writeslice.o density.o geostroph.o ran3.o tracerinit.o mprove.o dgtsv.o pcorrect.o topog.o inith_fixdh.o vort.o writeksurf.o writeyslice.o calcskfc.o cfdiv.o pfilter.o viscous.o prepvisc.o readcdfcn_tr.o tracersource.o srcface_nopy.o windstress.o calc_dvdy.o streamfunction.o n2budget.o calcn2.o pvcalc.o writen2budget.o

 VORO  =  vort_writecdf.o initcn_tr.o seval.o sigma_toplayer.o staticsigma.o findzall_toplayer_moves.o density.o readcdfcn_vor.o evalrho_tr.o smooth.o ran3.o meanh.o potdens.o

 VORXTO  =  vorticity_xt.o initcn_tr.o seval.o sigma.o findzall.o density.o readcdfcn_vor.o evalrho_tr.o smooth.o ran3.o meanh.o potdens.o

 VELXTO  =  vel_xt.o initcn_tr.o seval.o sigma.o findzall.o density.o readcdfcn_vor.o evalrho_tr.o smooth.o ran3.o meanh.o potdens.o

nh: $(NHO)
	$(FC) $(FFLAGS)  -o nh $(NHO) $(NC_LIB) $(NC_INC) 

nhcn: $(NHCNO)
	$(FC) $(FFLAGS)  -o nhcn $(NHCNO) -L/usr/local/lib -lnetcdf

writevor: $(VORO)
	gfortran $(FFLAGS)  -o writevor $(VORO) 
vorxt: $(VORXTO)
	gfortran $(FFLAGS)  -o vorxt $(VORXTO)
velxt: $(VELXTO)
	gfortran $(FFLAGS)  -o velxt $(VELXTO) -L/usr/local/lib -lnetcdf

clean :
	rm *.o *.mod


