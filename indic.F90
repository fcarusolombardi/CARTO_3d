      subroutine indic
      use param
      use mpih
      implicit none
      integer :: jc,kc,ic,k
!
!     direction normal to non-slip walls
!
      do kc=1,n3m
        kmv(kc)=kc-1
        kpv(kc)=kc+1
        if(kc.eq.1) kmv(kc)=kc
        if(kc.eq.n3m) kpv(kc)=kc
      enddo
      do kc=1,n3m
        kpc(kc)=kpv(kc)-kc
        kmc(kc)=kc-kmv(kc)
        kup(kc)=1-kpc(kc)
        kum(kc)=1-kmc(kc)
      enddo

!
!   indices for hor1 direction
!
      do ic=1,n1m
        imv(ic)=ic-1
        ipv(ic)=ic+1
        if(ic.eq.1) imv(ic)=n1m
        if(ic.eq.n1m) ipv(ic)=1
      enddo

!
!   indices for hor2 direction
!

!   Threaded for CPU affinity

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc)
      do jc=1,n2m
        jmv(jc)=jc-1
        jpv(jc)=jc+1
        if(jc.eq.1) jmv(jc)=n2m
        if(jc.eq.n2m) jpv(jc)=1
      enddo
!$OMP  END PARALLEL DO

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc)
      do jc = 1,n2+1
       jmhv(jc) = mod(jc,n2m/2+1)
       if(jmhv(jc).eq.0) jmhv(jc) = n2m/2 + 1
      enddo
!$OMP  END PARALLEL DO


#ifdef DEBUG
      if(myid.eq.0)then
        open(401,file='fact/indici.out')
        do k=1, n1m
          write(401,*) ipv(k), imv(k)
        enddo
        close(401)
        open(401,file='fact/indicj.out')
        do k=1, n2m
          write(401,*) jpv(k), jmv(k)
        enddo
        close(401)
        open(401,file='fact/indick.out')
        do k=1, n3m
          write(401,*) kpv(k), kmv(k)
        enddo
        close(401)
      endif
#endif

      return
      end

