!                       SUBROUTINE INVTRSAL
!   This subroutine performs the computation of the salinity field.
!   For details see the introduction of INVTR1
!
      subroutine invtrsa
      use param
      use local_arrays, only: dsal,hsal,rusal,rhsr
      use mpi_param, only: kstartr,kendr
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,ip,im
      real(DP)    :: dq32,dq33,dcq3,dq31
      real(DP)    :: app,acc,amm
      real(DP)    :: alpecl,udx1qr,udx2qr
      real(DP)    :: del1, del2, fcder
!
      alpecl=al * kps
      udx1qr=dx1qr
      udx2qr=dx2qr
!
!  compute the rhs of the factored equation
!  everything at i+1/2,j+1/2,k+1/2
!
      do kc=kstartr,kendr
      if( (kc.ge.2) .and. (kc.le.n3mr-1) ) then
        km=kmvr(kc)
        kp=kpvr(kc)
        app=ap3sskr(kc)
        acc=ac3sskr(kc)
        amm=am3sskr(kc)
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,jm,jp,ic,ip,im) &
!$OMP PRIVATE(dq31,dq32,dq33,dcq3)
        do jc=1,n2mr
          jm=jmvr(jc)
          jp=jpvr(jc)
          do ic=1,n1mr
            im=imvr(ic)
            ip=ipvr(ic)
!
!   11 second derivatives of dsal
!
            dq31=(dsal(ip,jc,kc)     &
                -2.d0*dsal(ic,jc,kc) &
                +dsal(im,jc,kc))*udx1qr
      
!
!   22 second derivatives of dsal
!
            dq32=(dsal(ic,jp,kc)       &
                -2.d0*dsal(ic,jc,kc)   &
                +dsal(ic,jm,kc))*udx2qr
!
!   33 second derivatives of dsal
!
            dq33= dsal(ic,jc,kp)*app  &
                +dsal(ic,jc,kc)*acc   &
                +dsal(ic,jc,km)*amm

            dcq3=dq32+dq33+dq31
!
!    right hand side of the density equation
!
!===========================================================
            rhsr(ic,jc,kc)=(ga*hsal(ic,jc,kc)+ro*rusal(ic,jc,kc) &
                   +alpecl*dcq3)*dts
!===========================================================
!
!    updating of the non-linear terms
!
            rusal(ic,jc,kc)=hsal(ic,jc,kc)
          enddo
        enddo
!$OMP  END PARALLEL DO
      endif
      if(kc.eq.1) then
        del1 = zmr(1)-zcr(1)
        del2 = zmr(2)-zmr(1)
        fcder = 2.d0/(del1*del2*(del1+del2))
        kp = kc + 1
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,jm,jp,ic,ip,im) &
!$OMP PRIVATE(dq31,dq32,dq33,dcq3)
        do jc=1,n2mr
          jm=jmvr(jc)
          jp=jpvr(jc)
          do ic=1,n1mr
            im=imvr(ic)
            ip=ipvr(ic)
            dq31=(dsal(ip,jc,kc)-2.d0*dsal(ic,jc,kc) &
               +dsal(im,jc,kc))*udx1qr
            dq32=(dsal(ic,jp,kc)-2.d0*dsal(ic,jc,kc) &
               +dsal(ic,jm,kc))*udx2qr
            dq33=(dsal(ic,jc,kp)*del1                     &
                -dsal(ic,jc,kc)*(del1+del2*dble(sbcbot))  &
                 +dsalbs(ic,jc)*del2*dble(sbcbot))*fcder

            dcq3=dq32+dq33+dq31
!=======================================================
            rhsr(ic,jc,kc)=(ga*hsal(ic,jc,kc)+ro*rusal(ic,jc,kc) &
                   +alpecl*dcq3)*dts
!=======================================================
            rusal(ic,jc,kc)=hsal(ic,jc,kc)
          enddo
        enddo
!$OMP  END PARALLEL DO
      endif
!
!       UPPER COLD WALL
!     
      if(kc.eq.n3mr) then
        del1 = zcr(n3r)-zmr(n3mr)
        del2 = zmr(n3mr)-zmr(n3mr-1)
        fcder = 2.d0/(del1*del2*(del1+del2))
        km = kc - 1 
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,jm,jp,ic,ip,im) &
!$OMP PRIVATE(dq31,dq32,dq33,dcq3)
        do jc=1,n2mr
          jm=jmvr(jc)
          jp=jpvr(jc)
          do ic=1,n1mr
            im=imvr(ic)
            ip=ipvr(ic)
            dq31=(dsal(ip,jc,kc)-2.d0*dsal(ic,jc,kc)         &
               +dsal(im,jc,kc))*udx1qr                       
           dq32=(dsal(ic,jp,kc)-2.d0*dsal(ic,jc,kc)          &
               +dsal(ic,jm,kc))*udx2qr                       
           dq33=(dsal(ic,jc,km)*del1                         &
                -dsal(ic,jc,kc)*(del1+del2*dble(sbctop))     &
                 +dsalbn(ic,jc)*del2*dble(sbctop))*fcder

            dcq3=dq32+dq33+dq31
!========================================================
            rhsr(ic,jc,kc)=(ga*hsal(ic,jc,kc)+ro*rusal(ic,jc,kc) &
                   +alpecl*dcq3)*dts
!========================================================
            rusal(ic,jc,kc)=hsal(ic,jc,kc)
          enddo
        enddo
!$OMP  END PARALLEL DO
      endif
      enddo

      call solxri( alpecl*dts*0.5d0*dx1qr )

      call solxrj( alpecl*dts*0.5d0*dx2qr )

      call solsak

!=========================================================
! periodic boundary condition
!  
      do kc=kstartr,kendr
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc)
      do jc=1,n2mr
      dsal(n1r,jc,kc) = dsal(1,jc,kc)
      enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ic)
      do ic=1,n1r
      dsal(ic,n2r,kc) = dsal(ic,1,kc)
      enddo
!$OMP END PARALLEL DO
      enddo

      return
      end

