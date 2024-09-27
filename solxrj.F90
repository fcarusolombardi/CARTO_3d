      subroutine solxrj(betadx)
!   Solves tridiagonal system in j direction for refined grid
      use param
      use local_arrays, only : rhsr
      use mpi_param, only: kstartr,kendr
      implicit none
      integer :: jc,kc,ic,info,ipkv(m2mr)
      real(DP),intent(in) :: betadx
      real(DP) :: amjl(m2mr-1),apjl(m2mr-1),acjl(m2mr),appk(m2mr-2)
      real(DP) :: ackl_b, aoffdiag, fjl(m2mr,m1mr)
      real(DP) :: qpd(m2mr), vtq, ytq   ! for periodic modification

      ackl_b = 1.d0/(1.d0+2.d0*betadx)
      aoffdiag = -betadx*ackl_b

      amjl(1:n2mr-1)=aoffdiag
      apjl(1:n2mr-1)=aoffdiag

      acjl(1) = 2.d0
      acjl(2:n2mr-1) = 1.d0
      acjl(n2mr) = 1.d0 + aoffdiag*aoffdiag


#ifndef MKLNOTAVAIL !modifica FV controlla
      call ddttrfb(n2mr,amjl,acjl,apjl,info)
#else
      call dgttrf(n2mr,amjl,acjl,apjl,appk,ipkv,info)
#endif


      qpd(1) = -1.d0
      qpd(2:n2mr-1) = 0.d0
      qpd(n2mr) = aoffdiag


#ifndef MKLNOTAVAIL !modifica FV controlla                                                                                    
      call ddttrsb('N',n2mr,1,amjl,acjl,apjl,qpd,n2mr,info) 
#else
       call dgttrs('N',n2mr,1,amjl,acjl,apjl,appk,ipkv,qpd,n2mr,info)
#endif


      vtq = 1.d0/(1.d0+qpd(1)-aoffdiag*qpd(n2mr))
      qpd(1:n2mr) = qpd(1:n2mr)*vtq

      do kc=kstartr,kendr
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,ic)
        do ic=1,n1mr
          do jc=1,n2mr
            fjl(jc,ic)=rhsr(ic,jc,kc)*ackl_b
          enddo
        end do
!$OMP  END PARALLEL DO

#ifndef MKLNOTAVAIL !modifica FV controlla                                                                                    
       call ddttrsb('N',n2mr,n1mr,amjl,acjl,apjl,fjl,n2mr,info)
#else
       call dgttrs('N',n2mr,n1mr,amjl,acjl,apjl,appk,ipkv,fjl,n2mr,info)
#endif

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,ic,ytq)
        do ic=1,n1mr
          ytq = fjl(1,ic) - aoffdiag*fjl(n2mr,ic)
          do jc=1,n2mr
            rhsr(ic,jc,kc) = fjl(jc,ic) - ytq*qpd(jc)
          enddo
        end do
!$OMP  END PARALLEL DO
      end do 

      return
      end
