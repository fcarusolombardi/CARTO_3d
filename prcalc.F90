
!***********************************************************************
!   this subroutine performs the calculation of the pressure.
!   this depends on the fractional step
!
      subroutine prcalc
      use param
      use local_arrays, only: pr,dph
      use mpi_param, only: kstart,kend
      implicit none
      integer :: kp,km,jm,jp,jc,kc,ic,ip,im
      real(DP)    :: be,amm,acc,app
!
!    the pressure is evaluated at the center of the box.
!
!     p^{n+1} = p^{n} + phi^{n+1} - b * Nabla^2 phi^{n+1}
!
      be=al*beta

#ifdef USE_CUDA
      !$cuf kernel do (3)
      do kc=kstart,kend
        do jc=1,n2m
          do ic=1,n1m
            kp=kpv(kc)
            km=kmv(kc)
            amm=amphk(kc)
            acc=acphk(kc)
            app=apphk(kc)
            jm=jmv(jc)
            jp=jpv(jc)
            im=imv(ic)
            ip=ipv(ic)
#else
      do kc=kstart,kend
        kp=kpv(kc)
        km=kmv(kc)
        amm=amphk(kc)
        acc=acphk(kc)
        app=apphk(kc)
       !$OMP PARALLEL DO DEFAULT(SHARED) &
       !$OMP PRIVATE(jc,jm,jp,ic,im,ip)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
#endif
            pr(ic,jc,kc)=pr(ic,jc,kc)+dph(ic,jc,kc)-be*( &
              (dph(ip,jc,kc)-2.d0*dph(ic,jc,kc)+dph(im,jc,kc))*dx1q &
             +(dph(ic,jp,kc)-2.d0*dph(ic,jc,kc)+dph(ic,jm,kc))*dx2q &
             +(dph(ic,jc,kp)*app+dph(ic,jc,kc)*acc+dph(ic,jc,km)*amm))
          enddo
        enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      enddo
      return
      end
!
!
