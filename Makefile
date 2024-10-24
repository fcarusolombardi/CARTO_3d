#=====================================================
# Note: use -C=all to catch array out of bounds errors
#============================================================================ 
# Compiler options
#============================================================================
# server='tesla'
# server='m100'
# server='dgx'
# server='uge'
# server='cpuuge'
# server='astra'
server='lngs'

ifeq ($(server), 'tesla')
# GPU mode Tesla
FC=mpif90 -r8 -pgf90libs
FC+=-Mnovect
FC+=-Mcuda=cuda11.3,cc70,ptxinfo -DUSE_NVTX -DUSE_CUDA
FC+=-Mcudalib=cufft
FC+=-mp
endif

ifeq ($(server), 'm100')
# GPU mode M100
# FC=mpipgifort -r8 -pgf90libs
# FC+=-Mnovect
# FC+=-Mcuda=cuda11.0,cc70,ptxinfo -DUSE_NVTX -DUSE_CUDA
# FC+=-Mcudalib=cufft,nccl
# FC+=-mp

FC=mpifort -r8
FC+=-Mnovect
FC+=-gpu=cuda11.0,cc70,ptxinfo -DUSE_NVTX -DUSE_CUDA
FC+=-cudalib=cufft,nccl
#,nvtx                                                                                                                                                          
FC+=-cuda
FC+=-mp
FC+=-DCUDECOMPAVAIL

# only if mkl is available (for now everywhere except m100)
FC+=-DMKLNOTAVAIL
# only if nccl is available (for now only on m100)
FC+=-DNCCLAVAIL
endif

ifeq ($(server), 'dgx')
# GPU mode Tesla                                                                                                                                                                  
FC=mpif90 -r8 -pgf90libs
FC+=-Mnovect
FC+=-Mcuda=cuda11.0,cc80,ptxinfo -DUSE_NVTX -DUSE_CUDA
FC+=-Mcudalib=cufft,nccl
FC+=-mp

FC+=-DMKLNOTAVAIL
FC+=-DNCCLAVAIL
endif



ifeq ($(server), 'uge')
FC=mpif90 -r8 -pgf90libs
FC+=-Mnovect
FC+=-gpu=cuda11.7,cc80,ptxinfo -DUSE_NVTX -DUSE_CUDA
FC+=-cudalib=cufft,nccl
FC+=-mp
FC+=-cuda
FC+=-DMKLNOTAVAIL
FC+=-DNCCLAVAIL
FC+=-DCUDECOMPAVAIL
endif


ifeq ($(server), 'astra')
FC=mpif90 -O0 -r8 -pgf90libs 
FC+=-Mnovect
FC+=-gpu=cuda11.8,cc80,ptxinfo -DUSE_NVTX -DUSE_CUDA
FC+=-cudalib=cufft,nccl
FC+=-ffpe-trap=overflow
FC+=-mp
FC+=-cuda
FC+=-DMKLNOTAVAIL
FC+=-DNCCLAVAIL
FC+=-DCUDECOMPAVAIL
#FC+=-g 
# #Debugging flags
FC+=-Kieee
# FC+=-Ktrap=fp
# FC+=-Ktrap=divz,inv,ovf,inexact,unf

# PGI_TERM=trace
# export PGI_TERM
# FC=mpif90 -r8 -pgf90libs -C -traceback
# FC+=-mp

endif

ifeq ($(server), 'debtesla')
# CPU degub mode
PGI_TERM=trace
export PGI_TERM
FC=mpif90 -r8 -pgf90libs -C -traceback

# CPU mode
# FC=mpif90 -r8 -pgf90libs
# FC+=-fast #-Mnovect
# FC+=-mp 
endif

ifeq ($(server), 'cpuuge')
# CPU mode
FC=mpif90 -r8 -pgf90libs
FC+=-fast #-Mnovect
FC+=-mp 
endif


ifeq ($(server), 'lngs')
# GPU mode          
FC=mpif90 -O0 -r8 -pgf90libs
FC+=-Mnovect
FC+=-Mcuda=cuda11.8,cc80,ptxinfo -DUSE_NVTX -DUSE_CUDA
FC+=-Mcudalib=cufft,nccl
FC+=-mp

FC+=-DMKLNOTAVAIL
FC+=-DNCCLAVAIL
#FC+=-DCUDECOMPAVAIL
endif


#->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
# FSI COUPLING: only one should be active
FC+=-DLOOSE
# FC+=-DSTRONG
# FC+=-DSTRONGRV
# FC+=-DEXPLICIT
# FC+=-DSCALINGTEST
# FC+=-DLUINVERSE

# CONTACT MODEL ACTIVE
# FC+=-DCONTACT
FC+=-DCONTACTV

# BODY OFFSET: can be either active or not
# FC+=-DOFFSETBODY

# ELECTRO BIDOMAIN MODEL
FC+=-DELECTRO
FC+=-DELEGEO0
# FC+=-DHYPERELASTIC
# FC+=-DBIDOMAIN

FC+=-DFIBERS
FC+=-DCARTO
FC+=-DEGM
# FC+=-DCORONARYVEINS

# Geometry type (Specify if the geometry has been created from ansa)
FC+=-DGEO_ANSA

# FSEI part
#FC+=-DFLUID_STRUCTURE_SOLVER

# -------------------------------------------
# Cellular models
# ventricles
# FC+=-DTP06
FC+=-DMINIMAL_MODEL
#FC+=-DMITCHELL_SCHAEFFER
# atria
#FC+=-DCOURTEMANCHE
#Purkinje
FC+=-DMV_PURK
#FC+=-DSTEWART09
# -------------------------------------------

# Fibrillation starting procedure
# FC+=-S1S2-3D

FC+=-DSOLOUNO

#->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
# flow config
#FC+=-DMOVIE
# FC+=-DDEBUG

# calculate temperature and scalars
#FC+=-DCALTEMP
#FC+=-DCALSCAL
#->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->


ifeq ($(server), 'tesla')
# Library Tesla
LINKS=-L/home/francesco/Programs/hdf5-1.10.5/hdf5/lib -lhdf5_fortran -lhdf5
LINKS+=-lfftw3 #-lfftw3_omp
LINKS+=-L/opt/intel/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
LINKS+=-L/usr/local/cuda/lib64 -lnvToolsExt
INCS=-I/home/francesco/Programs/hdf5-1.10.5/hdf5/include
INCS+=-I/opt/intel/compilers_and_libraries_2019.5.281/linux/mkl/include
endif

ifeq ($(server), 'm100')
# LibraryM100
# LINKS=-L/cineca/prod/opt/libraries/hdf5/1.12.0--spectrum_mpi--10.4.0/hpc-sdk--2021--binary/lib -lhdf5_fortran -lhdf5 
LINKS=-L/m100/home/userexternal/fviola00/Programs/hdf5/lib -lhdf5_fortran -lhdf5 
LINKS+=-L/cineca/prod/opt/libraries/fftw/3.3.8/spectrum_mpi--10.4.0--binary/lib -lfftw3
LINKS+=-L/cineca/prod/opt/libraries/lapack/3.9.0/pgi--19.10--binary -lblas -llapack
LINKS+=-L/cineca/prod/opt/compilers/cuda/11.0/none/lib64 -lnvToolsExt
LINKS+=-L/m100/home/userexternal/fviola00/Programs/cuDecomp/build/lib -lcudecomp -lcudecomp_fort

# INCS=-I/cineca/prod/opt/libraries/hdf5/1.12.0--spectrum_mpi--10.4.0/hpc-sdk--2021--binary/include
INCS=-I/m100/home/userexternal/fviola00/Programs/hdf5/include
INCS+=-I/cineca/prod/opt/libraries/fftw/3.3.8/spectrum_mpi--10.4.0--binary/include
INCS+=-I/cineca/prod/opt/libraries/lapack/3.9.0/pgi--19.10--binary
INCS+=-I/cineca/prod/opt/libraries/essl/6.2.1-0/binary/include/
INCS+=-I/m100/home/userexternal/fviola00/Programs/cuDecomp/build/include



# LINKS=-L${HDF5_HOME}/lib -lhdf5_fortran -lhdf5
# INCS=-I${HDF5_HOME}/include

# LINKS+=-L${NVHPC_HOME}/compilers/cuda/lib64 -lnvToolsExt

# LINKS+=-L${CUDECOMP_HOME}/lib -lcudecomp -lcudecomp_fort
# INCS+=-I${CUDECOMP_HOME}/include

# LINKS+=-L/cineca/prod/opt/libraries/fftw/3.3.8/spectrum_mpi--10.4.0--binary/lib -lfftw3
# LINKS+=-L/cineca/prod/opt/libraries/lapack/3.9.0/pgi--19.10--binary -lblas -llapack
endif

ifeq ($(server), 'dgx')
LINKS=-L/dgx/home/userexternal/fviola00/Programs/hdf5FV-1.12.1/lib -lhdf5_fortran -lhdf5
INCS=-I/dgx/home/userexternal/fviola00/Programs/hdf5FV-1.12.1/include
LINKS+=-L/cineca/prod/opt/libraries/lapack/3.9.0/pgi--19.10--binary/lib -lblas -llapack
INCS+=-I/cineca/prod/opt/libraries/lapack/3.9.0/pgi--19.10--binary/include
LINKS+=-L/dgx/home/userexternal/fviola00/Programs/fftwFV-3.3.10/lib -lfftw3
INCS+=-I/dgx/home/userexternal/fviola00/Programs/fftwFV-3.3.10/include
LINKS+=-L/cineca/prod/opt/compilers/nvhpc/2021/binary/Linux_x86_64/2021/cuda/lib64 -lnvToolsExt
endif

ifeq ($(server), 'uge')
# link both lapack and mkl, then the compilation flag will choose which one to be used
LINKS=-L/home/viola/Programs/hdf5-1.12.2/hdf5/lib -lhdf5_fortran -lhdf5
LINKS+=-lfftw3 #-lfftw3_omp
LINKS+=-L/usr/local/cuda-11.7/targets/x86_64-linux/lib -lblas -llapack
LINKS+=-L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/lib64 -lnvToolsExt
LINKS+=-L/opt/intel/oneapi/mkl/2022.1.0/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
LINKS+=-L/home/viola/Programs/cuDecomp/build/lib -lcudecomp -lcudecomp_fort

INCS=-I/home/viola/Programs/hdf5-1.12.2/hdf5/include
INCS+=-I/usr/local/cuda-11.7/targets/x86_64-linux/include
INCS+=-I/opt/intel/oneapi/mkl/2022.1.0/include
INCS+=-I/home/viola/Programs/cuDecomp/build/include
endif

ifeq ($(server), 'cpuuge')
# link both lapack and mkl, then the compilation flag will choose which one to be used
LINKS=-L/home/viola/Programs/hdf5-1.12.2/hdf5/lib -lhdf5_fortran -lhdf5
LINKS+=-lfftw3 #-lfftw3_omp
LINKS+=-L/usr/local/cuda-11.7/targets/x86_64-linux/lib -lblas -llapack
LINKS+=-L/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/lib64 -lnvToolsExt
LINKS+=-L/opt/intel/oneapi/mkl/2022.1.0/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

INCS=-I/home/viola/Programs/hdf5-1.12.2/hdf5/include
INCS+=-I/usr/local/cuda-11.7/targets/x86_64-linux/include
INCS+=-I/opt/intel/oneapi/mkl/2022.1.0/include
endif

ifeq ($(server), 'astra')
# link both lapack and mkl, then the compilation flag will choose which one to be used
LINKS=-L/home/franvio/Programs/hdf5-1.14.0_loc/lib -lhdf5_fortran -lhdf5
LINKS+=-lfftw3 #-lfftw3_omp
LINKS+=-L/usr/local/cuda-11.8/targets/x86_64-linux/lib -lblas -llapack
LINKS+=-L/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/cuda/lib64  -lnvToolsExt
LINKS+=-L/opt/intel/oneapi/mkl/2023.0.0/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
LINKS+=-L/home/franvio/Programs/cuDecomp/build/lib -lcudecomp -lcudecomp_fort

INCS=-I/home/franvio/Programs/hdf5-1.14.0_loc/include
INCS+=-I/usr/local/cuda-11.8/targets/x86_64-linux/lib/include
INCS+=-I/opt/intel/oneapi/mkl/2023.0.0/include
INCS+=-I/home/franvio/Programs/cuDecomp/build/include
endif


ifeq ($(server), 'debtesla')
# Library Tesla
LINKS=-L/home/francesco/Programs/hdf5-1.10.5/hdf5/lib -lhdf5_fortran -lhdf5
LINKS+=-lfftw3 #-lfftw3_omp
LINKS+=-L/opt/intel/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
LINKS+=-L/usr/local/cuda/lib64 -lnvToolsExt
INCS=-I/home/francesco/Programs/hdf5-1.10.5/hdf5/include
INCS+=-I/opt/intel/compilers_and_libraries_2019.5.281/linux/mkl/include
endif


ifeq ($(server), 'lngs')
PROGRAMS_HOME=/projects/gssi-fluids/Programs
LINKS=-L${PROGRAMS_HOME}/loc_hdf5-1.14.0/lib -lhdf5_fortran -lhdf5
INCS=-I${PROGRAMS_HOME}/loc_hdf5-1.14.0/include
LINKS+=-L/usr/local/cuda/targets/x86_64-linux/lib/ -lblas -llapack
INCS+=-I/usr/local/cuda/targets/x86_64-linux/include/
LINKS+=-L${PROGRAMS_HOME}/loc_hdf5-1.14.0/lib -lfftw3
INCS+=-I${PROGRAMS_HOME}/loc_hdf5-1.14.0/include
LINKS+=-L/software/prod/spack/opt/spack/linux-rocky8-zen/gcc-8.5.0/nvhpc-23.3-ftx6zua52fjmyibsnyeiut4ljzhnqaty/Linux_x86_64/23.3/cuda/11.8/lib64 -lnvToolsExt

LINKS+=-L${PROGRAMS_HOME}/cuDecomp/build/lib -lcudecomp -lcudecomp_fort
INCS+=-I${PROGRAMS_HOME}/cuDecomp/build/include


endif


#============================================================================ 
# make PROGRAM   
#============================================================================

PROGRAM=ibmheart_gpu.e

OBJECTS=continua_structures.o indson.o nvtx.o trans.o tridiag.o allocate_mls_local.o allocate_trigeo.o findindices.o \
	  TilingRoutines.o  mlsForceTiled3Comp.o velforce3Comp.o mlsForceTiledOff3Comp.o findindicesTiled.o \
         mlsStruc3Comp.o  \
	TriAuxRoutines.o MLSAuxRoutines.o write_to_vtk.o\
	  AllBoundCond.o um_ib.o GeoAuxRoutines.o internalforce.o internalforce_3d.o \
        balance.o cfl.o cflr.o coetar.o cordin.o continua.o densbo.o densmc.o divg.o divgck.o \
        divgloc.o dsalbo.o scar_tagging.o gcurv.o hdnlq1.o hdnlq2.o hdnlq3.o hdnlte.o hdnlsa.o \
        indic.o inirea.o initia.o inqpr.o invtrq1.o invtrq2.o \
        invtrq3.o explicitNS.o invtrte.o invtrsa.o matrix_transpose.o meshes.o \
        mpi_routines.o openfi.o param.o papero.o phcalc.o phini.o \
        prcalc.o solq1k.o solq2k.o solq3k.o soltek.o solsak.o solxi.o solxj.o solxri.o \
        solxrj.o stst.o preambolotsch.o tsch.o updvp.o vmaxv.o \
        coetarr.o indicr.o slab_rcd.o mkmov_dsal.o velbc.o toteng.o \
        mgrd_idc.o mgrd_mem.o mgrd_velitp.o mgrd_dsalc.o meanprofs.o \
        iniitp_grid.o iniitp_dens.o iniitp_velx.o iniitp_vely.o iniitp_velz.o iniitp_dsal.o \
	tripvmy_line.o electroEF_1d.o electroEF_2d.o electroEF_3d.o  deallocate_stuff.o postpro.o \
        mlsInterp3Comp.o 


MODULES=param.o

#============================================================================ 
# Linking    
#============================================================================

$(PROGRAM): $(MODULES) $(OBJECTS)
	$(FC) $(OBJECTS) $(LINKS) $(INCS) -o $@

#============================================================================
#  Dependencies
#============================================================================

param.o: param.F90
	$(FC) $(INCS) -c param.F90

%.o: %.F $(MODULES)
	$(FC) $(INCS) -c $<
%.o: %.F90 $(MODULES)
	$(FC) $(INCS) -c $<

#============================================================================
#  Clean up
#============================================================================

clean :
	rm *.o 
	rm *.mod
	rm ibmheart_gpu.e
