!
!***********************************************************************
!   this subroutine performs the inversion of the q2 momentum equation
!   by a factored implicit scheme, only the derivatives 11,22,33 of q2
!   are treated implicitly
!   in the first part the rhs is calculated
!        direction x3
!
      subroutine solq1k(ind)
      use param
      use tridiag
      use local_arrays, only : rhs, q1, forclo
      use local_arrays, only : rhst, forclot
      use local_arrays, only : rwork1, rwork2
      use local_arrays, only : amkl=>aml, apkl=>apl, ackl=>acl
      use mpi_param
      use mpih
      use outflow_vars
      use ibm_param
!@cuf use cudafor
      implicit none
      integer :: jc,kc,info,ipkv(m3m),ic,n
      integer, intent(in) :: ind
      real(DP) :: betadx
      real(DP) :: ugkks, ugkkn
      real(DP) :: ackl_b
!@cuf integer :: istat

      betadx=beta*al

      ugkks = 1./(g3rm(1)*g3rc(1))*dx3q*2.*betadx
      ugkkn = 1./(g3rm(n3m)*g3rc(n3))*dx3q*2.*betadx

      call PackZ_UnpackR(rhs(:,:,kstart:kend),rhst(:,:,jstart:jend))
      call PackZ_UnpackR(forclo(:,:,kstart:kend,ind),forclot(:,:,jstart:jend))

!
!  ************ compute dq2 sweeping along the x3 direction
!               periodic
!

      ! Set common tridiagonal system (with all nonzeros)
#ifdef USE_CUDA
      !$cuf kernel do(1)
#endif
      do kc=1,n3m
        ackl_b=1.0d0/(1.0d0-ac3sk(kc)*betadx)
        amkl(kc)=-am3sk(kc)*betadx*ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap3sk(kc)*betadx*ackl_b
      enddo

      ! Modify rhs
#ifdef USE_CUDA
      !$cuf kernel do(3)
#endif
      do jc=jstart,jend
#ifndef USE_CUDA
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(kc,ic,ackl_b)
#endif
        do ic=1,n1m
          do kc=1,n3m
            ackl_b=1.0d0/(1.0d0-ac3sk(kc)*betadx*forclot(kc,ic,jc))
            if (kc .eq. 1) then
              rhst(kc,ic,jc)=rhst(kc,ic,jc)*ackl_b + dqb1s(ic,jc)*ugkks
            elseif (kc .eq. n3m) then
              rhst(kc,ic,jc)=rhst(kc,ic,jc)*ackl_b + dqb1n(ic,jc)*ugkkn
            else
              rhst(kc,ic,jc)=rhst(kc,ic,jc)*ackl_b
            endif
          enddo
        enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      enddo

      if (allocated(rwork1) .and. size(rwork1) < size(rhst)) deallocate(rwork1)
      if (.not. allocated(rwork1)) allocate(rwork1(size(rhst)))
      if (allocated(rwork2) .and. size(rwork2) < size(rhst)) deallocate(rwork2)
      if (.not. allocated(rwork2)) allocate(rwork2(size(rhst)))

#ifdef USE_CUDA
      call trisolve_gpu(n3m,n1m*(jend-jstart+1),amkl,ackl,apkl,rhst,n3m,forclot,rwork1,rwork2)
#else
      call trisolve_cpu(n3m,n1m*(jend-jstart+1),amkl,ackl,apkl,rhst,n3m,forclot,rwork1)
#endif
!@cuf istat = cudaDeviceSynchronize() ! JDR TMP

      call PackR_UnpackZ(rhst(:,:,jstart:jend),rhs(:,:,kstart:kend))

#ifdef USE_CUDA
      !$cuf kernel do(3)
#endif
      do kc=kstart,kend
#ifndef USE_CUDA
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ic,jc)
#endif
        do jc=1,n2m
          do ic=1,n1m
            q1(ic,jc,kc) = q1(ic,jc,kc) + rhs(ic,jc,kc)
          enddo
        enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      enddo

!@cuf istat = cudaDeviceSynchronize() ! JDR TMP

      return
      end
