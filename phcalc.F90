!***********************************************************************
!  this subroutine perform the calculation of dph , periodic direction
!  along x3 and x1to use the real(DP) fourier transform
!
      subroutine phcalc
      use constants
      use param
      use local_arrays, only: dph, dpht
      ! JR NOTE: fft plans are setup using m and not n variables, meaning xr=>rhs_t only works
      ! since n and m vars are equal. We should remove all usage of the m variables in favor of
      ! dynamic sized allocations everywhere.
      !use local_arrays, only: xr=>rhs_t, xa=>cwork1_3D, rwork1      
      use local_arrays, only: xr=>rhs_t, rwork2, rwork1
      use mpi_param
      use mpih
      use nvtx
#ifdef USE_CUDA
      use cufft
      use trans
      use tridiag
#endif
!@cuf      use cudafor
      implicit none
      integer :: i,j,k,info
      real(DP) :: coefnorm,pmed
!#ifndef USE_CUDA
      real(DP) :: amphT(m3m-1),acphT(m3m),apphT(m3m-1),drhs(m3m),acphT_b(m3m)
!#endif
      integer :: n2mh,jmh
      real(DP), pointer, contiguous ::  xa(:,:,:)      
#ifdef USE_CUDA
      attributes(device) :: pmed,xa
#endif
!@cuf integer :: istat

      if (allocated(rwork2) .and. size(rwork2) < 2*m2mh*m1m*(kend-kstart+1)) deallocate(rwork2)
      if (.not. allocated(rwork2)) allocate(rwork2(2*m2mh*m1m*(kend-kstart+1)))
      xa(1:2*m2mh, 1:m1m, kstart:kend) => rwork2(:)

      n2mh=n2m/2+1

      coefnorm = 1.d0/(dble(n1m)*dble(n2m))

!   fft applied to the x2 direction to the
!   complex(DP) coeff. from cos fft
!   from physical to wave number space

      !dph_d = dph
      !call trans_xy(dph_d, xr_d, n1m, n2m, kend-kstart+1) 
     
#ifdef USE_CUDA
      !$cuf kernel do (3)
#endif
      do k=kstart,kend
#ifndef USE_CUDA
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(j,i)
#endif
        do j=1,n2m
          do i=1,n1m
           xr(j,i,k)=dph(i,j,k)
          enddo
        enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
     enddo
     
#ifdef USE_CUDA
      istat = cufftExecD2Z(cufft_fwd_plan, xr, xa)
#else
      do k=kstart,kend
        call dfftw_execute_dft_r2c(fwd_plan,xr(:,:,k),xa(:,:,k))
      enddo
#endif
!@cuf istat = cudaDeviceSynchronize() ! JDR TMP

!=================================================

      call nvtxStartRange("PackZ_UnpackRP", 1)
      call PackZ_complex_UnpackRP(xa,dpht)
      !call PackZ_UnpackRP(dpho,dpht)
      call nvtxEndRange

!
!m==================================================================
!     inversion of the matrix in the x2 and x3 directions (real part)
!

#ifdef USE_CUDA
      if (allocated(rwork1) .and. size(rwork1) < size(dpht)) deallocate(rwork1)
      if (.not. allocated(rwork1)) allocate(rwork1(size(dpht)))
      call tepDgtsv_pres_nopivot(n3m,n1m*(jendp-jstartp+1),n1m,(jendp-jstartp+1),&
                                 amphk,acphk,apphk,ak1,ak2,jmhv,dpht,&
                                 n3m,rwork1,jstartp)
#else
      call mkl_set_num_threads(1)
      do j=jstartp,jendp
        jmh=jmhv(j)
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(i,k,drhs,acphT_b,info) &
!$OMP PRIVATE(amphT,acphT,apphT)
        do i=1,n1m
          do k = 1,n3m
            acphT_b(k)=1.d0/(acphk(k)-ak2(jmh)-ak1(i))
            drhs(k)=dpht(k,i,j)*acphT_b(k)
          enddo

          amphT(1:n3m-1)=amphk(2:n3m)*acphT_b(2:n3m)
          apphT(1:n3m-1)=apphk(1:n3m-1)*acphT_b(1:n3m-1)
          acphT(1:n3m)=1.d0

          call ddttrfb(n3m,amphT,acphT,apphT,info)
          call ddttrsb('N',n3m,1,amphT,acphT,apphT,drhs,n3m,info)

          do k=1,n3m
            dpht(k,i,j) = drhs(k)
          enddo
        enddo
!$OMP  END PARALLEL DO
      enddo

      call mkl_set_num_threads(numthreads)
#endif
!@cuf istat = cudaDeviceSynchronize() ! JDR TMP

!=============================================

call nvtxStartRange("PackR_UnpackZP", 2)
call PackR_UnpackZP_complex(dpht,xa)
call nvtxEndRange
!
!================================================================
!   inverse fft applied to the phi x1 direction
!   from wave number space to physical space
!
#ifdef USE_CUDA
      istat = cufftExecZ2D(cufft_bck_plan, xa, xr)
#else
      do k=kstart,kend
        call dfftw_execute_dft_c2r(bck_plan,xa(:,:,k),xr(:,:,k))
      end do
#endif

#ifdef USE_CUDA
      !$cuf kernel do (3)
#endif
      do k=kstart,kend
#ifndef USE_CUDA
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(j,i)
#endif
       do j=1,n2m
         do i=1,n1m
           dph(i,j,k)=xr(j,i,k)
         enddo
       end do
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      end do

!@cuf istat = cudaDeviceSynchronize()

#ifdef USE_CUDA
      if (myid .eq. 0) istat = cudaMemcpy(pmed, dph(1,1,1), 1, cudaMemcpyDeviceToDevice)
#else
      if (myid .eq. 0) pmed = dph(1,1,1)
#endif

      call MPI_Bcast(pmed, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)

#ifdef USE_CUDA
      !$cuf kernel do(3)
#endif
      do k = kstart, kend
        do j = 1, n2m
          do i = 1, n1m
            dph(i,j,k) = dph(i,j,k) - pmed
          end do
        end do
      end do

!@cuf istat = cudaDeviceSynchronize() ! JDR TMP
      return
      end
