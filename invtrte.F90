
!***********************************************************************
!                       SUBROUTINE INVTRRO
!   This subroutine performs the computation of he scalar field.
!   For details see the introduction of INVTR1
!
      subroutine invtrte
      use param
      use local_arrays, only: dens,hro,ruro,rhs
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,ip,im
      real(DP)    :: dq32,dq33,dcq3,dq31
      real(DP)    :: app,acc,amm
      real(DP)    :: alpec,udx1q,udx2q
      real(DP)    :: del1,del2, fcder

      alpec=al * kpt
      udx1q=dx1q
      udx2q=dx2q
!
!  compute the rhs of the factored equation
!  everything at i+1/2,j+1/2,k+1/2
!
      do kc=kstart,kend
      if( (kc.ge.2) .and. (kc.le.n3m-1) ) then
        km=kmv(kc)
        kp=kpv(kc)
        app=ap3ssk(kc)
        acc=ac3ssk(kc)
        amm=am3ssk(kc)
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,jm,jp,ic,ip,im) &
!$OMP PRIVATE(dq31,dq32,dq33,dcq3)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
!
!   11 second derivatives of dens
!
            dq31=(dens(ip,jc,kc)      &
                -2.d0*dens(ic,jc,kc)  &
                +dens(im,jc,kc))*udx1q
      
!
!   22 second derivatives of dens
!
            dq32=(dens(ic,jp,kc)     &
                -2.d0*dens(ic,jc,kc) &
                +dens(ic,jm,kc))*udx2q
!
!   33 second derivatives of dens
!
            dq33= dens(ic,jc,kp)*app  &
                +dens(ic,jc,kc)*acc   &
                +dens(ic,jc,km)*amm

            dcq3=dq32+dq33+dq31

!    right hand side of the density equation
!
!===========================================================
            rhs(ic,jc,kc)=(ga*hro(ic,jc,kc)+ro*ruro(ic,jc,kc) &
                   +alpec*dcq3)*dt
!===========================================================
!
!    updating of the non-linear terms
!
            ruro(ic,jc,kc)=hro(ic,jc,kc)
        enddo
        enddo
!$OMP  END PARALLEL DO
      endif
      
      if(kc.eq.1) then
        del1 = zm(1)-zc(1)
        del2 = zm(2)-zm(1)
        fcder = 2.d0/(del1*del2*(del1+del2))
        kp = kc + 1
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,jm,jp,ic,ip,im) &
!$OMP PRIVATE(dq31,dq32,dq33,dcq3)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
            dq31=(dens(ip,jc,kc)-2.d0*dens(ic,jc,kc) &
               +dens(im,jc,kc))*udx1q

            dq32=(dens(ic,jp,kc)-2.d0*dens(ic,jc,kc) &
               +dens(ic,jm,kc))*udx2q

            dq33=( dens(ic,jc,kp)*del1                     &
                 -dens(ic,jc,kc)*(del1+del2*dble(tbcbot))  &
                 +densbs(ic,jc)*del2*dble(tbcbot)          &
                )*fcder

            dcq3=dq32+dq33+dq31
!=======================================================
            rhs(ic,jc,kc)=(ga*hro(ic,jc,kc)+ro*ruro(ic,jc,kc) &
                   +alpec*dcq3)*dt
!=======================================================
            ruro(ic,jc,kc)=hro(ic,jc,kc)
          enddo
        enddo
!$OMP  END PARALLEL DO
      endif
!
!       UPPER WALL
!     
      if(kc.eq.n3m) then
        del1 = zc(n3)-zm(n3m)
        del2 = zm(n3m)-zm(n3m-1)
        fcder = 2.d0/(del1*del2*(del1+del2))
        km = kc - 1 
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,jm,jp,ic,ip,im) &
!$OMP PRIVATE(dq31,dq32,dq33,dcq3)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
            dq31=( dens(ip,jc,kc)-2.d0*dens(ic,jc,kc) &
                 +dens(im,jc,kc))*udx1q

            dq32=( dens(ic,jp,kc)-2.d0*dens(ic,jc,kc) &
                 +dens(ic,jm,kc))*udx2q

            dq33=( dens(ic,jc,km)*del1                    &
                 -dens(ic,jc,kc)*(del1+del2*dble(tbctop)) &
                 +densbn(ic,jc)*del2*dble(tbctop)         &
                )*fcder

            dcq3=dq32+dq33+dq31
!========================================================
            rhs(ic,jc,kc)=(ga*hro(ic,jc,kc)+ro*ruro(ic,jc,kc) &
                   +alpec*dcq3)*dt
!========================================================
            ruro(ic,jc,kc)=hro(ic,jc,kc)
          enddo
        enddo
!$OMP  END PARALLEL DO
      endif
      enddo

      call solxi( alpec*dt*0.5d0*dx1q )

      call solxj( alpec*dt*0.5d0*dx2q )

      call soltek

!=========================================================
!  periodic boundary condition
!  
      do kc=kstart,kend
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc)
      do jc=1,n2m
      dens(n1,jc,kc) = dens(1,jc,kc)
      enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO & 
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ic)
      do ic=1,n1
      dens(ic,n2,kc) = dens(ic,1,kc)
      enddo
!$OMP END PARALLEL DO
      enddo

      return
      end
