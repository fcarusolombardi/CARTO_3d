!***********************************************************************
!  this subroutine perform the calculation of trigz for temperton fft
!
      subroutine fftqua
      use param
      integer :: n2mh,n2mp,j,i,n1mh,n1mp
      n1mh=n1m/2+1
      n1mp=n1mh+1
      n2mh=n2m/2+1
      n2mp=n2mh+1
!
!     wave number definition
!
      do i=1,n1mh
        ao(i)=(i-1)*2.d0*pi
      enddo
      do i=n1mp,n1m
        ao(i)=-(n1m-i+1)*2.d0*pi
      enddo
      do i=1,n1m
        ak1(i)=2.d0*(1.d0-dcos(ao(i)/n1m))*(dble(n1m)/rext1)**2
      enddo

      do j=1,n2mh
        ap(j)=(j-1)*2.d0*pi
      enddo
      do j=n2mp,n2m
        ap(j)=-(n2m-j+1)*2.d0*pi
      enddo
      do j=1,n2m
        ak2(j)=2.d0*(1.d0-dcos(ap(j)/n2m))*(dble(n2m)/rext2)**2
      enddo

      return
      end
!=====================================================
      subroutine phini
      use param
      use mpi_param
      use mpih
#ifdef USE_CUDA
      use cudafor
      use cufft
#endif
      implicit none
      integer FFTW_EXHAUSTIVE
      parameter(FFTW_EXHAUSTIVE=64)
      real(DP), dimension(m2m,m1m) :: xr
      complex(DP), dimension(m2mh,m1m) :: xa
#ifdef _OPENMP
      integer :: nt,fftw_info
#endif
#ifdef USE_CUDA
      integer :: istat
      integer(int_ptr_kind()) :: worksize
      integer :: rank(2)
      integer :: inembed(2), onembed(2)
#endif
    
!   Initialize tridiag matrices
      call tridiag_matrices   

!    Initialize FFTW
      call fftqua

#ifndef MKLNOTAVAIL      
#ifdef _OPENMP
!$OMP PARALLEL &
!$OMP REDUCTION(+:nt)
      nt = nt + 1
!$OMP  END PARALLEL
      call dfftw_init_threads(fftw_info)
      if(fftw_info.eq.0) write(*,*) "ERROR: FFTW THREAD INIT FAIL"
      call dfftw_plan_with_nthreads(nt)
#endif
#endif


      call dfftw_plan_dft_r2c_2d(fwd_plan,m2m,m1m,xr,xa,FFTW_EXHAUSTIVE)
      call dfftw_plan_dft_c2r_2d(bck_plan,m2m,m1m,xa,xr,FFTW_EXHAUSTIVE)

#ifdef USE_CUDA
      rank(1) = m1m
      rank(2) = m2m
      inembed(1) = m1m
      inembed(2) = m2m
      onembed(1) = m1m
      onembed(2) = m2mh

      istat = cufftPlanMany(cufft_fwd_plan, 2, rank, inembed, 1,  &
      m1m*m2m, onembed, 1, m1m*m2mh, CUFFT_D2Z, kend-kstart+1)
     
      istat = cufftPlanMany(cufft_bck_plan, 2, rank, onembed, 1,  &
      m1m*m2mh, inembed, 1, m1m*m2m, CUFFT_Z2D, kend-kstart+1)
#endif

      return
      end
      
!=======================================================================
      subroutine tridiag_matrices
      use param
      implicit none
      integer  :: kc,km,kp
      real(DP) :: ugmmm,a33icc,a33icp
!
!
!   tridiagonal matrix coefficients at each k and i
!   x1 and x3 cartesian coordinates
!
      do kc=1,n3m
        km=kmv(kc)
        kp=kpv(kc)
        a33icc=kmc(kc)*dx3q/g3rc(kc)
        a33icp=kpc(kc)*dx3q/g3rc(kp)
        ugmmm=1.0d0/g3rm(kc)
        amphk(kc)=a33icc*ugmmm
        apphk(kc)=a33icp*ugmmm
        acphk(kc)=-(amphk(kc)+apphk(kc))
      enddo

      end subroutine tridiag_matrices
!==================================================
      subroutine phend
      use param
      implicit none

      call dfftw_destroy_plan(fwd_plan)
      call dfftw_destroy_plan(bck_plan)

#ifndef MKLNOTAVAIL
#ifdef _OPENMP
      call dfftw_cleanup_threads()
#endif
#endif

      return
      end subroutine phend
