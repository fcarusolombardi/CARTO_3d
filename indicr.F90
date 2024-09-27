      subroutine indicr
      use param
      use mpih
      implicit none
      integer :: jc,kc,ic
   
!
!   endices for hor1 direction
!
      do ic=1,n1mr
        imvr(ic)=ic-1
        ipvr(ic)=ic+1
        if(ic.eq.1) imvr(ic)=n1mr
        if(ic.eq.n1mr) ipvr(ic)=1
      enddo

!
!   indices for hor2 direction
!

!   Threaded for CPU affinity

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc)
      do jc=1,n2mr
        jmvr(jc)=jc-1
        jpvr(jc)=jc+1
        if(jc.eq.1) jmvr(jc)=n2mr
        if(jc.eq.n2mr) jpvr(jc)=1
      enddo
!$OMP  END PARALLEL DO

!
!     direction normal to non-slip walls
!
      do kc=1,n3mr
        kmvr(kc)=kc-1
        kpvr(kc)=kc+1
        if(kc.eq.1) kmvr(kc)=kc
        if(kc.eq.n3mr) kpvr(kc)=kc
      end do


#ifdef DEBUG
      if(myid.eq.0)then
        open(401,file='fact/indicir.out')
        do kc=1, n1mr
          write(401,*) ipvr(kc), imvr(kc)
        enddo
        close(401)
        open(401,file='fact/indicjr.out')
        do kc=1, n2mr
          write(401,*) jpvr(kc), jmvr(kc)
        enddo
        close(401)
        open(401,file='fact/indickr.out')
        do kc=1, n3mr
          write(401,*) kpvr(kc), kmvr(kc)
        enddo
        close(401)
      endif
#endif

      return
      end

