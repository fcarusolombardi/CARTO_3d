!***********************************************************************
!   this subroutine performs the inversion of the q3 momentum equation
!   by a factored implicit scheme, only the derivatives 11,22,33 of q3
!   are treated implicitly
!       direction x3
!
      subroutine soltek
      use param
      use local_arrays, only : dens,rhs
      use mpi_param
      implicit none
      integer :: jc,kc,ic,info,ipkv(m3m)
      real(DP) :: betadx,ackl_b(m3m),fkl(m3m,m1m)
      real(DP) :: amkT(m3m-1),ackT(m3m),apkT(m3m-1),appk(m3m-2)
      real(DP), allocatable, dimension(:,:,:) :: rhst

      allocate(rhst(1:n3m,1:n1m,jstart:jend))

      call PackZ_UnpackR(rhs,rhst)

!      dimension amkl(m3),apkl(m3),ackl(m3),fkl(m1,m3) 
!  ********* compute the dq3* sweeping in the x3 direction
!
!
      betadx=0.5d0*al*dt*kpt
     
      do kc=1,n3m
        ackl_b(kc)=1.d0/(1.d0-ac3ssk(kc)*betadx)
      enddo

      amkT(1:n3m-1)=-betadx*am3ssk(2:n3m)*ackl_b(2:n3m)
      apkT(1:n3m-1)=-betadx*ap3ssk(1:n3m-1)*ackl_b(1:n3m-1)
      ackT(1:n3m)=1.d0


#ifndef MKLNOTAVAIL !modifica FV controlla
      call ddttrfb(n3m,amkT,ackT,apkT,info)
#else
      call dgttrf(n3m,amkT,ackT,apkT,appk,ipkv,info)
#endif

      do jc=jstart,jend
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ic,kc)
        do ic=1,n1m
          do kc=1,n3m
            fkl(kc,ic)=rhst(kc,ic,jc)*ackl_b(kc)
          end do
        enddo
!$OMP END PARALLEL DO

#ifndef MKLNOTAVAIL !modifica FV controlla          
       call ddttrsb('N',n3m,n1m,amkT,ackT,apkT,fkl,n3m,info)
#else
       call dgttrs('N',n3m,n1m,amkT,ackT,apkT,appk,ipkv,fkl,n3m,info)
#endif
          
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ic,kc)
        do ic=1,n1m
          do kc=1,n3m
            rhst(kc,ic,jc)= fkl(kc,ic)
          end do
        enddo
!$OMP END PARALLEL DO
      end do

      call PackR_UnpackZ(rhst,rhs)

      do kc=kstart,kend
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,ic)
        do jc=1,n2m
          do ic=1,n1m
            dens(ic,jc,kc) = dens(ic,jc,kc) + rhs(ic,jc,kc)
          enddo
        enddo
!$OMP  END PARALLEL DO
      enddo

      if(allocated(rhst)) deallocate(rhst)

      return
      end
