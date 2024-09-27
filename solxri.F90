      subroutine solxri(betadx)
!   Solves tridiagonal system in i direction for refined grid
      use param
      use local_arrays, only : rhsr
      use mpi_param, only: kstartr,kendr
      implicit none
      integer :: jc,kc,ic,info,ipkv(m1m)!FV
      real(DP),intent(in) :: betadx
      real(DP) :: amil(m1mr-1),acil(m1mr),apil(m1mr-1),appk(m1m-2)!FV
      real(DP) :: ackl_b,aoffdiag,fil(m1mr,m2mr)
      real(DP) :: qpd(m1mr), vtq, ytq   ! for periodic modification

      ackl_b = 1.d0/(1.d0+2.d0*betadx)
      aoffdiag = -betadx*ackl_b

      amil(1:n1mr-1)=aoffdiag
      apil(1:n1mr-1)=aoffdiag

      acil(1) = 2.d0
      acil(2:n1mr-1) = 1.d0
      acil(n1mr) = 1.d0 + aoffdiag*aoffdiag

#ifndef MKLNOTAVAIL !modifica FV controlla                                                                                    
      call ddttrfb(n1mr,amil,acil,apil,info)
#else
      call dgttrf(n1mr,amil,acil,apil,appk,ipkv,info)   
#endif



      qpd(1) = -1.d0
      qpd(2:n1mr-1) = 0.d0
      qpd(n1mr) = aoffdiag

#ifndef MKLNOTAVAIL !modifica FV controlla
      call ddttrsb('N',n1mr,1,amil,acil,apil,qpd,n1mr,info)
#else
       call dgttrs('N',n1mr,1,amil,acil,apil,appk,ipkv,qpd,n1mr,info)
#endif

      vtq = 1.d0/(1.d0+qpd(1)-aoffdiag*qpd(n1mr))
      qpd(1:n1mr) = qpd(1:n1mr)*vtq

      do kc=kstartr,kendr
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,ic)
        do jc=1,n2mr
          do ic=1,n1mr
            fil(ic,jc)=rhsr(ic,jc,kc)*ackl_b
          enddo
        end do
!$OMP  END PARALLEL DO



#ifndef MKLNOTAVAIL !modifica FV controlla
        call ddttrsb('N',n1mr,n2mr,amil,acil,apil,fil,n1mr,info)
#else
       call dgttrs('N',n1mr,n2mr,amil,acil,apil,appk,ipkv,fil,n1mr,info)   
#endif


!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,ic,ytq)
        do jc=1,n2mr
          ytq = fil(1,jc)-aoffdiag*fil(n1mr,jc)
          do ic=1,n1mr
            rhsr(ic,jc,kc) = fil(ic,jc) - ytq*qpd(ic)
          enddo
        end do
!$OMP  END PARALLEL DO
      end do 

      return
      end
