!***********************************************************************
!
!                                                                       
!***********************************************************************
      subroutine toteng(aengc,aengr)
      use param
      use local_arrays, only: q1,q2,q3,dens,dsal
      use mpih
      use mpi_param, only: kstart,kend,kstartr,kendr
!@cuf use cudafor
      implicit none
      real(DP),intent(out) :: aengc, aengr
      integer :: jc,kc,ic
      real(DP)    :: dvolf,dvolc,my_aengc,my_aengr
      real(DP)    :: engdens,engdsal, my_engdens, my_engdsal
!@cuf integer :: istat
        
!==================
      my_aengc = 0.d0
      my_engdens = 0.d0
#ifdef USE_CUDA
      !$cuf kernel do (3) <<<*,*>>>
      do kc=kstart,kend
       do jc=1,n2m
        do ic=1,n1m
         dvolc=(zc(kc+1)-zc(kc))/dx1/dx2
         if(kc.eq.1)then
           dvolf=(zm(kc)-zc(kc))/dx1/dx2
         else
           dvolf=(zm(kc)-zm(kc-1))/dx1/dx2
         endif
         my_aengc = my_aengc + &
          (q1(ic,jc,kc)**2 + q2(ic,jc,kc)**2)*dvolc + &
          (q3(ic,jc,kc)**2)*dvolf
         !my_engdens = my_engdens + dens(ic,jc,kc)*zm(kc)*dvolc
        enddo
       enddo
      enddo
#else
      do kc=kstart,kend
       dvolc=(zc(kc+1)-zc(kc))/dx1/dx2
       if(kc.eq.1)then
         dvolf=(zm(kc)-zc(kc))/dx1/dx2
       else
         dvolf=(zm(kc)-zm(kc-1))/dx1/dx2
       endif
!$OMP  PARALLEL DO &
!$OMP  DEFAULT(SHARED) &
!$OMP  PRIVATE(jc,ic) &
!$OMP  REDUCTION(+:my_aengc,my_engdens)
       do jc=1,n2m
        do ic=1,n1m
         my_aengc = my_aengc + &
          (q1(ic,jc,kc)**2 + q2(ic,jc,kc)**2)*dvolc + &
          (q3(ic,jc,kc)**2)*dvolf
         !my_engdens = my_engdens + dens(ic,jc,kc)*zm(kc)*dvolc
        enddo
       enddo
!$OMP  END PARALLEL DO
      enddo
#endif

!@cuf istat = cudaDeviceSynchronize !JDR TMP

!====================
      aengc = 0.d0
      call MPI_ALLREDUCE(my_aengc,aengc,1,MDP,MPI_SUM, MPI_COMM_WORLD,ierr)

      engdens = 0.d0
      !call MPI_ALLREDUCE(my_engdens,engdens,1,MDP,MPI_SUM, MPI_COMM_WORLD,ierr)

      ! Removing refined grid energy computations
      aengr=0.d0
      engdsal=0.d0

!============================
      if(myid.eq.0)then
        write(96,796) time,aengc*0.5d0+engdsal-engdens*rhop, &
                     aengc*0.5d0, engdsal-engdens*rhop, aengr*0.5d0
      endif
796   format(1x,f10.4,4(1x,ES20.8))
      
      return     
      end         





      subroutine totengQ(aengcQ1,aengcQ2,aengcQ3,aengr)
      use param
      use local_arrays, only: q1,q2,q3,dens,dsal
      use mpih
      use mpi_param, only: kstart,kend,kstartr,kendr
!@cuf use cudafor
      implicit none
      real(DP),intent(out) :: aengcQ1,aengcQ2,aengcQ3, aengr
      integer :: jc,kc,ic
      real(DP)    :: dvolf,dvolc,my_aengcQ1,my_aengcQ2,my_aengcQ3,my_aengr
      real(DP)    :: engdens,engdsal, my_engdens, my_engdsal
!@cuf integer :: istat
        
!==================
      my_aengcQ1 = 0.d0
      my_aengcQ2 = 0.d0
      my_aengcQ3 = 0.d0
      my_engdens = 0.d0
#ifdef USE_CUDA
      !$cuf kernel do (3) <<<*,*>>>
      do kc=kstart,kend
       do jc=1,n2m
        do ic=1,n1m
         dvolc=(zc(kc+1)-zc(kc))/dx1/dx2
         if(kc.eq.1)then
           dvolf=(zm(kc)-zc(kc))/dx1/dx2
         else
           dvolf=(zm(kc)-zm(kc-1))/dx1/dx2
         endif
         my_aengcQ1 = my_aengcQ1 + abs(q1(ic,jc,kc))
         my_aengcQ2 = my_aengcQ2 + abs(q2(ic,jc,kc))
         my_aengcQ3 = my_aengcQ3 + abs(q3(ic,jc,kc))
         !my_engdens = my_engdens + dens(ic,jc,kc)*zm(kc)*dvolc
        enddo
       enddo
      enddo
#else
      do kc=kstart,kend
       dvolc=(zc(kc+1)-zc(kc))/dx1/dx2
       if(kc.eq.1)then
         dvolf=(zm(kc)-zc(kc))/dx1/dx2
       else
         dvolf=(zm(kc)-zm(kc-1))/dx1/dx2
       endif
!$OMP  PARALLEL DO &
!$OMP  DEFAULT(SHARED) &
!$OMP  PRIVATE(jc,ic) &
!$OMP  REDUCTION(+:my_aengcQ1,my_aengcQ2,my_aengcQ3,my_engdens)
       do jc=1,n2m
        do ic=1,n1m
         my_aengcQ1 = my_aengcQ1 + abs(q1(ic,jc,kc))
         my_aengcQ2 = my_aengcQ2 + abs(q2(ic,jc,kc))
         my_aengcQ3 = my_aengcQ3 + abs(q3(ic,jc,kc))
        enddo
       enddo
!$OMP  END PARALLEL DO
      enddo
#endif

!@cuf istat = cudaDeviceSynchronize !JDR TMP

!====================
      aengcQ1 = 0.d0
      aengcQ2 = 0.d0
      aengcQ3 = 0.d0
      call MPI_ALLREDUCE(my_aengcQ1,aengcQ1,1,MDP,MPI_SUM, MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_aengcQ2,aengcQ2,1,MDP,MPI_SUM, MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_aengcQ3,aengcQ3,1,MDP,MPI_SUM, MPI_COMM_WORLD,ierr)

      engdens = 0.d0
      !call MPI_ALLREDUCE(my_engdens,engdens,1,MDP,MPI_SUM, MPI_COMM_WORLD,ierr)

      ! Removing refined grid energy computations
      aengr=0.d0
      engdsal=0.d0

!============================
!       if(myid.eq.0)then
!         write(96,796) time,aengc*0.5d0+engdsal-engdens*rhop, &
!                      aengc*0.5d0, engdsal-engdens*rhop, aengr*0.5d0
!       endif
! 796   format(1x,f10.4,4(1x,ES20.8))
      
      return     
      end         
