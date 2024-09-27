!***********************************************************************
!   this subroutine performs the inversion of the q3 momentum equation
!   by a factored implicit scheme, only the derivatives 11,22,33 of q3
!   are treated implicitly
!       direction x3
!
      subroutine solsak
      use param
      use local_arrays, only : dsal,rhsr
      use mpi_param
      use mpih
      implicit none
      integer :: jc,kc,ic,info,ipkv(m3mr)
      real(DP) :: betadx,ackl_b(m3mr),fkl(m3mr,m1mr)
      real(DP) :: amkT(m3mr-1),ackT(m3mr),apkT(m3mr-1),appk(m3mr-2)
      real(DP), allocatable, dimension(:,:,:) :: rhst

      allocate(rhst(1:n3mr,1:n1mr,jstartr:jendr))

      call PackZ_UnpackR_refi(rhsr,rhst)

!      dimension amkl(m3),apkl(m3),ackl(m3),fkl(m1,m3) 
!  ********* compute the dq3* sweeping in the x3 direction
!
!
      betadx=0.5d0*al*dts*kps
     
      do kc=1,n3mr
        ackl_b(kc)=1.d0/(1.d0-ac3sskr(kc)*betadx)
      enddo

      amkT(1:n3mr-1)=-betadx*am3sskr(2:n3mr)*ackl_b(2:n3mr)
      apkT(1:n3mr-1)=-betadx*ap3sskr(1:n3mr-1)*ackl_b(1:n3mr-1)
      ackT=1.d0

#ifndef MKLNOTAVAIL !modifica FV controlla
      call ddttrfb(n3mr,amkT,ackT,apkT,info)
#else
     call dgttrf(n3mr,amkT,ackT,apkT,appk,ipkv,info)
#endif

      do jc=jstartr,jendr
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ic,kc)
        do ic=1,n1mr
          do kc=1,n3mr
            fkl(kc,ic)=rhst(kc,ic,jc)*ackl_b(kc)
          end do
        enddo
!$OMP END PARALLEL DO

#ifndef MKLNOTAVAIL !modifica FV controlla
       call ddttrsb('N',n3mr,n1mr,amkT,ackT,apkT,fkl,n3mr,info)
#else
       call dgttrs('N',n3mr,n1mr,amkT,ackT,apkT,appk,ipkv,fkl,n3mr,info)
#endif
          
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ic,kc)
        do ic=1,n1mr
          do kc=1,n3mr
            rhst(kc,ic,jc)= fkl(kc,ic)
          end do
        enddo
!$OMP END PARALLEL DO
      end do

      call PackR_UnpackZ_refi(rhst,rhsr)

      do kc=kstartr,kendr
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,ic)
        do jc=1,n2mr
          do ic=1,n1mr
            dsal(ic,jc,kc) = dsal(ic,jc,kc) + rhsr(ic,jc,kc)
          enddo
        enddo
!$OMP  END PARALLEL DO
      enddo

      if(allocated(rhst)) deallocate(rhst)

      return
      end
 
