!
!     this subroutine calculates divg(q).
!
      subroutine divg
      use param
      use local_arrays, only: q1,q2,q3,dph
      use mpi_param, only: kstart,kend
!@cuf use cudafor
      implicit none
      integer :: jc,jp,kc,kp,ic,ip
      real(DP)    :: usdtal,dqcap   
!@cuf integer :: istat

      usdtal = real(1.0,DP)/(dt*al)

#ifdef USE_CUDA
     !$cuf kernel do(3)
#else
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(kc,jc,ic,kp,jp,ip,dqcap) 
#endif
      do kc=kstart,kend
        do jc=1,n2m
          do ic=1,n1m
            kp=kc+1
            jp=jpv(jc)
            ip=ipv(ic)
            dqcap= (q1(ip,jc,kc)-q1(ic,jc,kc))*dx1 &
                 +(q2(ic,jp,kc)-q2(ic,jc,kc))*dx2  &
                 +(q3(ic,jc,kp)-q3(ic,jc,kc))*udx3m(kc)
            dph(ic,jc,kc)=dqcap*usdtal
          enddo
        enddo
      enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif

!@cuf istat = cudaDeviceSynchronize() !JDR TMP
      return
      end
