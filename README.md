
PSOM: Process Study Ocean Model
===============================

PSOM, pronounced "soam" (the nectar derived from the churning of the oceans in Indian mythology), stands for Process Study Ocean Model. It is a versatile, three-dimensional, non-hydrostatic, computational fluid dynamical model for oceanographic (as well as other) applications (Mahadevan et al., 1996a,b). The model uses the finite volume method on a structured grid with boundary fitted coordinates (topography conforming sigma grid in the vertical, and boundary conforming in the horizontal). The model has a free-surface. It can be used for large- and small-scale phenomena and can be run in hydrostatic or non-hydrostatic mode (Mahadevan, 2006). It uses a highly efficient solution procedure that exploits the skewness arising from the small geometrical aspect (depth to length scale) ratio of the ocean to speed up the solution of the non-hydrostatic pressure, which is solved by the multigrid method. The model has been used for a number of process studies, including investigation of the vertical transport of nutrients for phytoplankton production (Mahadevan and Archer, 2000) and the dynamics of submesoscale processes (Mahadevan and Tandon, 2006; Mahadevan, Tandon and Ferrari, 2010). Since the non-hydrostatic model is well-posed with open boundaries, it can be used as a nested high-resolution model with time-varying boundary conditions applied from a coarser resolution general circulation model (Mahadevan and Archer, 1998). The model is thus ideally suited for high-resolution, limited-region modeling studies. 



This code was created for Takeyoshi
stprofile_rho was replaced by stprofile_Kuroshio

to test Jinbo's initialization, replace stprofile_kuroshio with stprofile_deepfront_jbw.f90