!***********************************************************************
      subroutine divgck(qmaxc,qmaxr)
      use param
      use local_arrays, only: q2,q3,q1
      use mpih
      use mpi_param, only: kstart,kend
!@cuf use cudafor
      implicit none
      real(DP),intent(out) :: qmaxc, qmaxr
      integer :: jc,kc,kp,jp,ic,ip
      real(DP)    :: dqcap,my_qmaxc
!@cuf integer :: istat

!========================================================
      my_qmaxc =-100.d0

#ifdef USE_CUDA
     !$cuf kernel do(3)
      do kc=kstart,kend
        kp=kc+1
        do jc=1,n2m
          do ic=1,n1m
            jp=jpv(jc)
            ip=ipv(ic)
#else
      do kc=kstart,kend
        kp=kc+1
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,ic,jp,ip,dqcap) &
!$OMP REDUCTION(max:my_qmaxc)
        do jc=1,n2m
          jp=jpv(jc)
          do ic=1,n1m
            ip=ipv(ic)
#endif
            dqcap= (q1(ip,jc,kc)-q1(ic,jc,kc))*dx1 &
                  +(q2(ic,jp,kc)-q2(ic,jc,kc))*dx2 &
                  +(q3(ic,jc,kp)-q3(ic,jc,kc))*udx3m(kc)
            my_qmaxc = max(abs(dqcap),my_qmaxc)
          enddo
        enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      enddo

!@cuf istat = cudaDeviceSynchronize !JDR TMP

      qmaxc=0.d0
      call MPI_ALLREDUCE(my_qmaxc,qmaxc,1,MDP,MPI_MAX,MPI_COMM_WORLD,ierr)

! qmaxr should be removed from all the cals
      qmaxr=0.d0


      return
      end
