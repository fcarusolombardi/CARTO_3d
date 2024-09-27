!   This routine calculates the maximum velocities
      subroutine vmaxv
      use param
      use local_arrays, only: q2,q3,q1
      use mpi_param, only: kstart,kend
      use mpih
!@cuf use cudafor
      implicit none
      real(DP)    :: my_vmax1,my_vmax2,my_vmax3
      integer :: jc,kc,kp,ic
!@cuf integer :: istat

      my_vmax1=real(-100.0,DP)
      my_vmax2=real(-100.0,DP)
      my_vmax3=real(-100.0,DP)

#ifdef USE_CUDA
      !$cuf kernel do(3) <<<*,*>>>
      do kc=kstart,kend
        do jc=1,n2m
          do ic=1,n1m
            my_vmax1 = max(my_vmax1,abs(q1(ic,jc,kc)))
            my_vmax2 = max(my_vmax2,abs(q2(ic,jc,kc)))
            my_vmax3 = max(my_vmax3,abs(q3(ic,jc,kc)))
          enddo
        enddo
      enddo
#else
      do kc=kstart,kend
!$OMP  PARALLEL DO &
!$OMP  DEFAULT(SHARED) &
!$OMP  PRIVATE(jc,ic) &
!$OMP  REDUCTION(max: my_vmax1,my_vmax2,my_vmax3)
        do jc=1,n2m
          do ic=1,n1m
            my_vmax1 = max(my_vmax1,abs(q1(ic,jc,kc)))
            my_vmax2 = max(my_vmax2,abs(q2(ic,jc,kc)))
            my_vmax3 = max(my_vmax3,abs(q3(ic,jc,kc)))
          enddo
        enddo
!$OMP  END PARALLEL DO
      enddo
#endif

!@cuf istat = cudaDeviceSynchronize !JDR TMP

      call MPI_ALLREDUCE(my_vmax1,vmax(1),1,MDP,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_vmax2,vmax(2),1,MDP,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_vmax3,vmax(3),1,MDP,MPI_MAX,MPI_COMM_WORLD,ierr)

!      if(myid.eq.0)then 
!        anutin= anutin*vol/kpt - (denstop-densbot)
!        anusin= anusin*volr/kps - (dsaltop-dsalbot)
!        !write(95,510) time, dabs(anutin), dabs(anusin)
!      endif
!510   format(1x,f10.4,2(1x,ES20.8))

      return
      end    
!
