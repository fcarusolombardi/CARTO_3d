
!***********************************************************************
!
      subroutine velbc
      use param
      use mpih
      use local_arrays, only: q2,q3,q1
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,ic

      real(DP) dl1q, dl2q, dlf1, dlf2

!=========================================================
!  periodic boundary condition
!
#ifdef USE_CUDA
      !$cuf kernel do (2)
#endif
      do kc=kstart-1,kend+1
#ifndef USE_CUDA
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc)
#endif
        do jc=1,n2m
          q1(n1,jc,kc) = q1(1,jc,kc)
          q2(n1,jc,kc) = q2(1,jc,kc)
          q3(n1,jc,kc) = q3(1,jc,kc)
        enddo
#ifndef USE_CUDA
!$OMP END PARALLEL DO
#endif
      end do

#ifdef USE_CUDA
      !$cuf kernel do (2)
#endif
      do kc=kstart-1,kend+1
#ifndef USE_CUDA
     !$OMP PARALLEL DO &
     !$OMP DEFAULT(SHARED) &
     !$OMP PRIVATE(ic)
#endif
        do ic=1,n1
          q1(ic,n2,kc) = q1(ic,1,kc)
          q2(ic,n2,kc) = q2(ic,1,kc)
          q3(ic,n2,kc) = q3(ic,1,kc)
        enddo
#ifndef USE_CUDA
!$OMP END PARALLEL DO
#endif
      enddo

      return
      end

