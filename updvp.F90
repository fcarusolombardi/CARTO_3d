
!***********************************************************************
!
!     this subroutine calculates the solenoidal vel field.
!       q(n+1)=qhat-grad(dph)*dt ,  pr=dph
!    third order runge-kutta is used.
!
      subroutine updvp
      use param
      use local_arrays, only: q2,q3,dph,q1
      use mpi_param, only: kstart,kend
!@cuf use cudafor
      implicit none
      integer :: jc,jm,kc,km,ic,im
      real(DP)    :: usukm,udx2,udx1,locdph
!@cuf integer :: istat

      udx1 = al*dt*dx1
      udx2 = al*dt*dx2
#ifdef USE_CUDA
      !$cuf kernel do (3)
      do kc=kstart,kend
        do jc=1,n2m
          do ic=1,n1m
          km=kmv(kc)
          usukm = al*dt*udx3c(kc)
          jm=jmv(jc)
          im=imv(ic)
#else
      do kc=kstart,kend
        km=kmv(kc)
        usukm = al*dt*udx3c(kc)
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,jm,ic,im,locdph)
        do jc=1,n2m
          jm=jmv(jc)
          do ic=1,n1m
          im=imv(ic)
#endif
          locdph=dph(ic,jc,kc)
          q1(ic,jc,kc)=q1(ic,jc,kc)- &
           (locdph-dph(im,jc,kc))*udx1
          q2(ic,jc,kc)=q2(ic,jc,kc)- &
           (locdph-dph(ic,jm,kc))*udx2
          q3(ic,jc,kc)=q3(ic,jc,kc)- &
           (locdph-dph(ic,jc,km))*usukm
        enddo 
       enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      enddo

!@cuf istat = cudaDeviceSynchronize() !JDR TMP
      return
      end

