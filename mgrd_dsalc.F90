      subroutine mgrd_dsalc
      use param
      use local_arrays, only: dsal 
      use mgrd_arrays, only: dsalc
      use mpi_param
      use mpih
      implicit none
       
      integer ic,jc,kc,icr,jcr,kcr

      real(DP) ldzr, ldz, dmrefxy, dsalloc

!m=========================================================
!  single mesh

      IF(mref1.eq.1 .and. mref2.eq.1 .and. mref3.eq.1)then

        do kc=kstart-1,kend+1
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,ic)
         do jc=1,n2
         do ic=1,n1
           dsalc(ic,jc,kc) = dsal(ic,jc,kc)
         enddo
         enddo
!$OMP END PARALLEL DO
        enddo

!m=========================================================
!  double mesh

      ELSE

        dmrefxy = 1.d0/dble(mref1*mref2)

        do kc=kstart,kend
         ldz = 1.d0/(zc(kc+1)-zc(kc))
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,ic,kcr,jcr,icr) &
!$OMP PRIVATE(ldzr,dsalloc)
         do jc=1,n2m
         do ic=1,n1m
           dsalloc = 0.d0
           do kcr=(kc-1)*mref3+1,kc*mref3
             ldzr = zcr(kcr+1)-zcr(kcr)
             do jcr=(jc-1)*mref2+1,jc*mref2
             do icr=(ic-1)*mref1+1,ic*mref1
               dsalloc=dsalloc+dsal(icr,jcr,kcr)*ldzr
             enddo
             enddo
           enddo
           dsalc(ic,jc,kc) = dsalloc*ldz*dmrefxy
         enddo
         enddo
!$OMP END PARALLEL DO
        enddo

!m==========================================================
!   periodic B.C.

        do kc=kstart,kend
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc)
         do jc=1,n2m
           dsalc(n1,jc,kc) = dsalc(1,jc,kc)
         enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ic)
         do ic=1,n1
           dsalc(ic,n2,kc) = dsalc(ic,1,kc)
         enddo
!$OMP END PARALLEL DO
        enddo

!m==========================================================
!  boundary value 

        if(kstart.eq.1)then
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ic,jc)
         do jc=1,n2
          do ic=1,n1
            dsalc(ic,jc,0) = dsalbot 
          enddo
         enddo
!$OMP END PARALLEL DO
        endif

        if(kend.eq.n3m)then
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ic,jc)
         do jc=1,n2
          do ic=1,n1
           dsalc(ic,jc,n3) = dsaltop
          enddo
         enddo
!$OMP END PARALLEL DO
        endif

!m=====================================================

        call update_both_ghosts(n1,n2,dsalc,kstart,kend)

      ENDIF

      return
      end  subroutine mgrd_dsalc
