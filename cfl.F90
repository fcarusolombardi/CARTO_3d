      subroutine cfl(cflm)
      use param
      use local_arrays, only: q1,q2,q3
      use mpih
      use mpi_param, only: kstart,kend
!@cuf use cudafor
      implicit none
      real(DP),intent(inout)    :: cflm
      real(DP)    :: my_cflm
      integer :: j,k,jp,kp,i,ip
      real(DP) :: qcf,udx3
!@cuf integer :: istat
      
      my_cflm=real(1.e-8,DP)


#ifdef USE_CUDA
     !$cuf kernel do(3) <<<*,*>>>
      do k=kstart,kend
        do j=1,n2m
          do i=1,n1m
            udx3=udx3m(k)
            kp=k+1
            jp=j+1
            ip=i+1
            qcf=( abs((q1(i,j,k)+q1(ip,j,k))*0.5_DP*dx1) &
                 +abs((q2(i,j,k)+q2(i,jp,k))*0.5_DP*dx2) &
                 +abs((q3(i,j,k)+q3(i,j,kp))*0.5_DP*udx3))

            my_cflm = dmax1(my_cflm,qcf)
          enddo
        enddo
      enddo
#else
      do k=kstart,kend
        udx3=udx3m(k)
        kp=k+1
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(j,i,jp,ip,qcf) &
!$OMP REDUCTION(max: my_cflm)
        do j=1,n2m
          jp=j+1
          do i=1,n1m
            ip=i+1
            qcf=( dabs((q1(i,j,k)+q1(ip,j,k))*0.5d0*dx1) &
                +dabs((q2(i,j,k)+q2(i,jp,k))*0.5d0*dx2) &
                +dabs((q3(i,j,k)+q3(i,j,kp))*0.5d0*udx3))

            my_cflm = dmax1(my_cflm,qcf)
          enddo
        enddo
!$OMP  END PARALLEL DO
      enddo
#endif

!@cuf istat = cudaDeviceSynchronize !JDR TMP
            
      call MPI_ALLREDUCE(my_cflm,cflm,1,MDP,MPI_MAX, &
             MPI_COMM_WORLD,ierr)

      return  
      end                                                               
