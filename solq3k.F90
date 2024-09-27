!***********************************************************************
!   this subroutine performs the inversion of the q3 momentum equation
!   by a factored implicit scheme, only the derivatives 11,22,33 of q3
!   are treated implicitly
!       direction x3
!
      subroutine solq3k(ind)
      use param
      use tridiag
      use local_arrays, only : rhs, q3, forclo
      use local_arrays, only : rwork1, rwork2
      use local_arrays, only : rhst, forclot
      use local_arrays, only : amkl=>aml, apkl=>apl, ackl=>acl
      use mpi_param
      use mpih
      use outflow_vars
      use ibm_param
!@cuf use cudafor
      implicit none
      integer :: jc,kc,info,ipkv(m3),ic,n,ke
      integer, intent(in) :: ind
      real(DP) :: betadx
      real(DP) :: ackl_b
!@cuf integer :: istat

      betadx=beta*al

      call PackZ_UnpackR(rhs(:,:,kstart:kend),rhst(:,:,jstart:jend))
      call PackZ_UnpackR(forclo(:,:,kstart:kend,ind),forclot(:,:,jstart:jend))


!  ********* compute the dq3* sweeping in the x3 direction
!
      ! Set common tridiagonal system (with all nonzeros)
#ifdef USE_CUDA
      !$cuf kernel do(1)
#endif
      do kc=1,n3m
        ackl_b=1.0d0/(1.0d0-ac3sk(kc)*betadx)
        ackl(kc)=1.0d0
        if (kc .eq. 1) then
          amkl(kc)=0.d0
          apkl(kc)=0.d0
        else
          amkl(kc)=-am3sk(kc)*betadx*ackl_b
          apkl(kc)=-ap3sk(kc)*betadx*ackl_b
        endif
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
              rhst(kc,ic,jc)=dqb3s(ic,jc)
            else if (kc .eq. n3m) then
              ! Subtraction of BC to avoid solving n3 equation. Recomputing apkl(n3m) term
              ! here to properly handle forclot term
              rhst(kc,ic,jc)=rhst(kc,ic,jc)*ackl_b - (-ap3sk(kc)*betadx*ackl_b*forclot(kc,ic,jc)) * dqb3n(ic,jc)
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
!     JR note: There is no need to solve all n3 equations, as last equation is just setting
!     known boundary condition (which is set explicitly in final loop in file).
!     We deal with it by modifying the n3m RHS term in previous loop.
!     TODO: Should the boundary at kc=1 also be excluded and set explicitly?
#ifdef USE_CUDA
      call trisolve_gpu(n3m,n1m*(jend-jstart+1),amkl,ackl,apkl,rhst,n3m,forclot,rwork1,rwork2)
#else
      call trisolve_cpu(n3m,n1m*(jend-jstart+1),amkl,ackl,apkl,rhst,n3m,forclot,rwork1)
#endif
!@cuf istat = cudaDeviceSynchronize() ! JDR TMP

      call PackR_UnpackZ(rhst(:,:,jstart:jend),rhs(:,:,kstart:kend))

      ! Set upper limit of next loop to hit n3 if required
      if (kend .eq. n3m) then
        ke = n3
      else
        ke = kend
      endif

#ifdef USE_CUDA
      !$cuf kernel do(3)
#endif
      do kc=kstart,ke
#ifndef USE_CUDA
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ic,jc)
#endif
        do jc=1,n2m
          do ic=1,n1m
            if (kc .eq. n3) then
              q3(ic,jc,kc) = q3(ic,jc,kc) + dqb3n(ic,jc)
            else
              q3(ic,jc,kc) = q3(ic,jc,kc) + rhs(ic,jc,kc)
            endif
          enddo
        enddo
#ifndef USE_CUDA
!$OMP END PARALLEL DO
#endif
      enddo

!@cuf istat = cudaDeviceSynchronize() ! JDR TMP

      return
      end
