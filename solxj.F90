      subroutine solxj(betadx,ind)
!   Solves tridiagonal system in j direction
      use param
      use tridiag
      use local_arrays, only : rwork1, rwork2, rwork3
      use local_arrays, only : rhs, rhs_t
      use local_arrays, only : forclo, forclo_t
      use local_arrays, only : amjl=>aml, apjl=>apl, acjl=>acl
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
      do jc=1,n2m
        apjl(jc)=offdiag
        acjl(jc)=1.0d0
        amjl(jc)=offdiag
      enddo

#ifdef USE_CUDA
      ! Modify rhs without transpose
      !$cuf kernel do(3)
      do kc=kstart,kend
        do jc=1,n2m
          do ic=1,n1m
            ackl_b = 1.0/(1.0+2.0*betadx*forclo(ic,jc,kc,ind))
            rhs(ic,jc,kc)=rhs(ic,jc,kc)*ackl_b
          enddo
        enddo
      enddo
#else
      ! Modify rhs with transpose
      do kc=kstart,kend
!$OMP PARALLEL DO &
!$OMP DEFAULT(NONE) &
!$OMP SHARED(kc,forclo,forclo_t,n1m,n2m,rhs,rhs_t,betadx,ind) &
!$OMP PRIVATE(jc,ic,ackl_b)
        do ic=1,n1m
          do jc=1,n2m
            forclo_t(jc,ic,kc)=forclo(ic,jc,kc,ind)
            ackl_b = 1.0/(1.0+2.0*betadx*forclo_t(jc,ic,kc))
            rhs_t(jc,ic,kc)=rhs(ic,jc,kc)*ackl_b
          enddo
        enddo
!$OMP  END PARALLEL DO
      enddo
#endif

      if (allocated(rwork1) .and. size(rwork1) < size(rhs)) deallocate(rwork1)
      if (.not. allocated(rwork1)) allocate(rwork1(size(rhs)))
      if (allocated(rwork2) .and. size(rwork2) < size(rhs)) deallocate(rwork2)
      if (.not. allocated(rwork2)) allocate(rwork2(size(rhs)))
      if (allocated(rwork3) .and. size(rwork3) < size(rhs)) deallocate(rwork3)
      if (.not. allocated(rwork3)) allocate(rwork3(size(rhs)))
#ifdef USE_CUDA
      call trisolve_periodic_3D_notrans_gpu(n2m,n1m,(kend-kstart+1),amjl,acjl,apjl,rhs(1,1,kstart),forclo(1,1,kstart,ind),rwork1,rwork2,rwork3)
#else
      call trisolve_periodic_cpu(n2m,n1m*(kend-kstart+1),amjl,acjl,apjl,rhs_t,n2m,forclo_t,rhs,rwork1,rwork2)
#endif

#ifndef USE_CUDA
      ! Transpose rhs back to original configuration
      do kc=kstart,kend
!$OMP PARALLEL DO &
!$OMP DEFAULT(NONE) &
!$OMP SHARED(kc,n1m,n2m,rhs,rhs_t) &
!$OMP PRIVATE(jc,ic)
        do jc=1,n2m
          do ic=1,n1m
            rhs(ic,jc,kc)=rhs_t(jc,ic,kc)
          enddo
        enddo
!$OMP  END PARALLEL DO
      enddo
#endif

!@cuf istat = cudaDeviceSynchronize() ! JDR TMP
      return
      end
