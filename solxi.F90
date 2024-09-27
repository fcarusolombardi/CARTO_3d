      subroutine solxi(betadx,ind)
!   Solves tridiagonal system in i direction
      use param
      use tridiag
      use local_arrays, only : rwork1, rwork2, rwork3
      use local_arrays, only : rhs, rhs_t, forclo
      use local_arrays, only : amil=>aml, apil=>apl, acil=>acl
      use mpi_param, only: kstart,kend
      use ibm_param
!@cuf use cudafor
      implicit none
      integer :: jc,kc,ic,n
      real(DP),intent(in) :: betadx
      integer,intent(in) :: ind
      real(DP) :: ackl_b, offdiag
!@cuf integer :: istat

      ! Set common tridiagonal system (with all nonzeros)
      ackl_b = 1.0/(1.0+2.0*betadx)
      offdiag=-betadx*ackl_b
#ifdef USE_CUDA
      !$cuf kernel do(1)
#endif
      do ic=1,n1m
        apil(ic)=offdiag
        acil(ic)=1.0d0
        amil(ic)=offdiag
      enddo

      ! Modify rhs
#ifdef USE_CUDA
      !$cuf kernel do(3)
#endif
      do kc=kstart,kend
#ifndef USE_CUDA
!$OMP PARALLEL DO &
!$OMP DEFAULT(NONE) &
!$OMP SHARED(n2m,n1m,kc,betadx,rhs,forclo,ind) &
!$OMP PRIVATE(jc,ic,ackl_b)
#endif
        do jc=1,n2m
          do ic=1,n1m
            ackl_b = 1.0/(1.0+2.0*betadx*forclo(ic,jc,kc,ind))
            rhs(ic,jc,kc)=rhs(ic,jc,kc)*ackl_b
          enddo
        enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      enddo

      if (allocated(rwork1) .and. size(rwork1) < size(rhs)) deallocate(rwork1)
      if (.not. allocated(rwork1)) allocate(rwork1(size(rhs)))
      if (allocated(rwork2) .and. size(rwork2) < size(rhs)) deallocate(rwork2)
      if (.not. allocated(rwork2)) allocate(rwork2(size(rhs)))
      if (allocated(rwork3) .and. size(rwork3) < size(rhs)) deallocate(rwork3)
      if (.not. allocated(rwork3)) allocate(rwork3(size(rhs)))
#ifdef USE_CUDA
      call trisolve_periodic_3D_gpu(n1m,n2m,kend-kstart+1,amil,acil,apil,rhs(1,1,kstart),forclo(1,1,kstart,ind),rhs_t(1,1,kstart),rwork1,rwork2,rwork3)
#else
      call trisolve_periodic_cpu(n1m,n2m*(kend-kstart+1),amil,acil,apil,rhs,n1m,forclo(:,:,:,ind),rhs_t,rwork1,rwork2)
#endif

!@cuf istat = cudaDeviceSynchronize() ! JDR TMP
 
      return
      end
