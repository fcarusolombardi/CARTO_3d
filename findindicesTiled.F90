!------------------------------------------------------------------
!     Routine to compute indices of all centroids on all particles
!     Only starting and ending indices of the support domain are
!     reqiured. Then loop over it.
!------------------------------------------------------------------

      subroutine findindicesTiled
      USE mpih
      USE param
      USE mls_param
      USE tile_arrays
!@cuf USE cudafor

      IMPLICIT NONE

      real(DP) post1, post2, post3
      integer i1,j1,k1,ist,jst,kst,inp,ntr
      integer ii,jj,kk
      real(DP) tstr,rcdp,x2dp
      real(DP) a1,a2,alp,alp1
!     --------INNER AND OUTER VARIABLES----------
      integer i1i,j1i,k1i,isti,jsti,ksti
      integer i1o,j1o,k1o,isto,jsto,ksto

      integer n_tile, rowtile, tile, tilex
      integer v1,v2,v3

      real(DP) tstri,rcdpi,x2dpi
      real(DP) tstro,rcdpo,x2dpo

      real(DP) a1i,a2i,alpi,alp1i
      real(DP) a1o,a2o,alpo,alp1o
!     --------------------------------------------
      real(DP) disdum,distsdum,ztmp,offmin,off
      integer ind_pal,nclipr,n3mor,k
      integer kok,kini,kfin,kmid
!@cuf integer :: istat

!FV multiple if or elseif inside cuda loop  does not work well!
      if (istr3.EQ.0) then

#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do tile = 1,nttot

        ntr = faceid_t(tile)

        n_tile = tri_tiling(ntr,1)
        rowtile = tri_tiling(ntr,2)

        tilex = tile - tstart(ntr) + 1

        post1 = albegaBar(rowtile,tilex,1)*tri_ver(1,ntr) + &
                albegaBar(rowtile,tilex,2)*tri_ver(4,ntr) + &
                albegaBar(rowtile,tilex,3)*tri_ver(7,ntr)
        post2 = albegaBar(rowtile,tilex,1)*tri_ver(2,ntr) + &
                albegaBar(rowtile,tilex,2)*tri_ver(5,ntr) + &
                albegaBar(rowtile,tilex,3)*tri_ver(8,ntr)
        post3 = albegaBar(rowtile,tilex,1)*tri_ver(3,ntr) + &
                albegaBar(rowtile,tilex,2)*tri_ver(6,ntr) + &
                albegaBar(rowtile,tilex,3)*tri_ver(9,ntr)

!       ++++++++Indices of the marker++++++++++++++++++++++++
!       X - indices
        i1=NINT((post1+xyz_tr(1))*dx1) +1
        ist=FLOOR((post1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m

!       Y - indices
        j1=NINT((post2+xyz_tr(2))*dx2) +1
        jst=FLOOR((post2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m

!       Z - indices
        k1=NINT((post3+xyz_tr(3))*dx3) +1
        kst=FLOOR((post3+xyz_tr(3))*dx3) + 1

!       OUTER and INNER INDICES
        isto = ist+1 ; isti = ist-1
        jsto = jst+1 ; jsti = jst-1
        ksto = kst+1 ; ksti = kst-1
!       -------------------------------------------------------------
        pindt(1,tile)=i1 ; pindt(2,tile)=j1 ; pindt(3,tile)=k1
        pindt(4,tile)=ist; pindt(5,tile)=jst; pindt(6,tile)=kst
!       -------------------------------------------------------------
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP


!        elseif (istr3.EQ.6) then !Chebyshev !FV not always correct

! #ifdef USE_CUDA
!       !$cuf kernel do (1)
! #endif
!       do tile = 1,nttot

!         ntr = faceid_t(tile)

!         n_tile = tri_tiling(ntr,1)
!         rowtile = tri_tiling(ntr,2)

!         tilex = tile - tstart(ntr) + 1

!         post1 = albegaBar(rowtile,tilex,1)*tri_ver(1,ntr) + &
!                 albegaBar(rowtile,tilex,2)*tri_ver(4,ntr) + &
!                 albegaBar(rowtile,tilex,3)*tri_ver(7,ntr)
!         post2 = albegaBar(rowtile,tilex,1)*tri_ver(2,ntr) + &
!                 albegaBar(rowtile,tilex,2)*tri_ver(5,ntr) + &
!                 albegaBar(rowtile,tilex,3)*tri_ver(8,ntr)
!         post3 = albegaBar(rowtile,tilex,1)*tri_ver(3,ntr) + &
!                 albegaBar(rowtile,tilex,2)*tri_ver(6,ntr) + &
!                 albegaBar(rowtile,tilex,3)*tri_ver(9,ntr)

! !       ++++++++Indices of the marker++++++++++++++++++++++++
! !       X - indices
!         i1=NINT((post1+xyz_tr(1))*dx1) +1
!         ist=FLOOR((post1+xyz_tr(1))*dx1) +1      
!         if(ist.eq.0)ist=n1m

! !       Y - indices
!         j1=NINT((post2+xyz_tr(2))*dx2) +1
!         jst=FLOOR((post2+xyz_tr(2))*dx2) + 1
!         if(jst.eq.0)jst=n2m

!  !       Z - indices
!         !a1=etazm(1+int(str3)); a2=etazm(n3+int(str3))
!         nclipr = int(str3)*mref3
!         n3mor = n3r+nclipr+nclipr
!         a1 = dcos(PI*(dble(1+int(str3))-0.5d0)/dble(n3mor))
!         a2 = dcos(PI*(dble(n3+int(str3))-0.5d0)/dble(n3mor))
        
!         alp=(post3+xyz_tr(3))/(0.5D0*alx3)
!         alp1=(1.0D0-alp)*(0.5D0*(a1-a2))
        
!         ! k1=0.5D0+((1.0D0/PI)*(float(n3)+str3+str3)*acos(alp1))
!         k1=0.5D0+((1.0D0/PI)*(real(n3)+str3+str3)*acos(alp1))          
!         k1=int(k1-str3)
!         !staggered
!         kst=k1
! !        if((post3).gt.zc(k1))kst=k1   !fv
!         if((post3).le.zc(k1))kst=k1-1 !fv
        
! !       OUTER and INNER INDICES
!         isto = ist+1 ; isti = ist-1
!         jsto = jst+1 ; jsti = jst-1
!         ksto = kst+1 ; ksti = kst-1
! !       -------------------------------------------------------------
!         pindt(1,tile)=i1 ; pindt(2,tile)=j1 ; pindt(3,tile)=k1
!         pindt(4,tile)=ist; pindt(5,tile)=jst; pindt(6,tile)=kst
! !       -------------------------------------------------------------
!       end do
! !@cuf istat = cudaDeviceSynchronize !JDR TMP
      

      else !General

#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do tile = 1,nttot

        ntr = faceid_t(tile)

        n_tile = tri_tiling(ntr,1)
        rowtile = tri_tiling(ntr,2)

        tilex = tile - tstart(ntr) + 1

        post1 = albegaBar(rowtile,tilex,1)*tri_ver(1,ntr) + &
                albegaBar(rowtile,tilex,2)*tri_ver(4,ntr) + &
                albegaBar(rowtile,tilex,3)*tri_ver(7,ntr)
        post2 = albegaBar(rowtile,tilex,1)*tri_ver(2,ntr) + &
                albegaBar(rowtile,tilex,2)*tri_ver(5,ntr) + &
                albegaBar(rowtile,tilex,3)*tri_ver(8,ntr)
        post3 = albegaBar(rowtile,tilex,1)*tri_ver(3,ntr) + &
                albegaBar(rowtile,tilex,2)*tri_ver(6,ntr) + &
                albegaBar(rowtile,tilex,3)*tri_ver(9,ntr)

!       ++++++++Indices of the marker++++++++++++++++++++++++
!       X - indices
        i1=NINT((post1+xyz_tr(1))*dx1) +1
        ist=FLOOR((post1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m

!       Y - indices
        j1=NINT((post2+xyz_tr(2))*dx2) +1
        jst=FLOOR((post2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m

!       Z - indices
        k1=0
        kini=1
        kfin=n3        
        do while (k1.eq.0)
           if ((kfin-kini).eq.1) then
              ztmp=0.5D0*(zc(kfin)+zc(kini))
              if (post3.le.ztmp) then
                 k1=kini
              else
                 k1=kfin
              endif
           endif
           kmid=(kini+kfin)/2
           if (post3.gt.zc(kini).and.post3.le.zc(kmid)) then
              kini=kini
              kfin=kmid
           else
              kini=kmid
              kfin=kfin
           endif
        enddo
        kst=k1
!        if((post3).gt.zc(k1))kst=k1   !fv
        if((post3).le.zc(k1))kst=k1-1 !fv
        
!       OUTER and INNER INDICES
        isto = ist+1 ; isti = ist-1
        jsto = jst+1 ; jsti = jst-1
        ksto = kst+1 ; ksti = kst-1
!       -------------------------------------------------------------
        pindt(1,tile)=i1 ; pindt(2,tile)=j1 ; pindt(3,tile)=k1
        pindt(4,tile)=ist; pindt(5,tile)=jst; pindt(6,tile)=kst
!       -------------------------------------------------------------
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP
      
       endif


      
      return
      end


      ! subroutine findindicesTiledOff(dOff)
      subroutine findindicesTiledOff
      USE mpih
      USE param
      USE mls_param
      USE tile_arrays
!@cuf USE cudafor

      IMPLICIT NONE

      real(DP) post1, post2, post3
      integer i1,j1,k1,ist,jst,kst,inp,ntr
      integer ii,jj,kk
      real(DP) tstr,rcdp,x2dp
      real(DP) a1,a2,alp,alp1
!     --------INNER AND OUTER VARIABLES----------
      integer i1i,j1i,k1i,isti,jsti,ksti
      integer i1o,j1o,k1o,isto,jsto,ksto

      integer n_tile, rowtile, tile, tilex
      integer v1,v2,v3

      real(DP) tstri,rcdpi,x2dpi
      real(DP) tstro,rcdpo,x2dpo

      real(DP) a1i,a2i,alpi,alp1i
      real(DP) a1o,a2o,alpo,alp1o
!     --------------------------------------------
      real(DP) disdum,distsdum,ztmp,offmin,off
      real(DP) Offx,Offy,Offz
      integer ind_pal,nclipr,n3mor,k
      integer kok,kini,kfin,kmid
!@cuf integer :: istat

!FV multiple if or elseif inside cuda loop  does not work well!
      if (istr3.EQ.0) then

#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
        do tile = ntOffstart,ntOffend       

        ntr = faceid_t(tile)

        Offx = dOff*tri_nor(1,ntr)
        Offy = dOff*tri_nor(2,ntr)
        Offz = dOff*tri_nor(3,ntr)
        
        rowtile = tri_tiling(ntr,2)

        tilex = tile - tstart(ntr) + 1

        post1 = albegaBar(rowtile,tilex,1)*tri_ver(1,ntr) + &
                albegaBar(rowtile,tilex,2)*tri_ver(4,ntr) + &
                albegaBar(rowtile,tilex,3)*tri_ver(7,ntr) + Offx
        post2 = albegaBar(rowtile,tilex,1)*tri_ver(2,ntr) + &
                albegaBar(rowtile,tilex,2)*tri_ver(5,ntr) + &
                albegaBar(rowtile,tilex,3)*tri_ver(8,ntr) + Offy
        post3 = albegaBar(rowtile,tilex,1)*tri_ver(3,ntr) + &
                albegaBar(rowtile,tilex,2)*tri_ver(6,ntr) + &
                albegaBar(rowtile,tilex,3)*tri_ver(9,ntr) + Offz

!       ++++++++Indices of the marker++++++++++++++++++++++++
!       X - indices
        i1=NINT((post1+xyz_tr(1))*dx1) +1
        ist=FLOOR((post1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m

!       Y - indices
        j1=NINT((post2+xyz_tr(2))*dx2) +1
        jst=FLOOR((post2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m

!       Z - indices
        k1=NINT((post3+xyz_tr(3))*dx3) +1
        kst=FLOOR((post3+xyz_tr(3))*dx3) + 1

!       OUTER and INNER INDICES
        isto = ist+1 ; isti = ist-1
        jsto = jst+1 ; jsti = jst-1
        ksto = kst+1 ; ksti = kst-1
!       -------------------------------------------------------------
        pindtOff(1,tile)=i1 ; pindtOff(2,tile)=j1 ; pindtOff(3,tile)=k1
        pindtOff(4,tile)=ist; pindtOff(5,tile)=jst; pindtOff(6,tile)=kst
!       -------------------------------------------------------------
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP

      else !General

#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do tile = ntOffstart,ntOffend   

         ntr = faceid_t(tile)
         
         Offx = dOff*tri_nor(1,ntr)
         Offy = dOff*tri_nor(2,ntr)
         Offz = dOff*tri_nor(3,ntr)
         
         rowtile = tri_tiling(ntr,2)

         tilex = tile - tstart(ntr) + 1

         post1 = albegaBar(rowtile,tilex,1)*tri_ver(1,ntr) + &
              albegaBar(rowtile,tilex,2)*tri_ver(4,ntr) + &
              albegaBar(rowtile,tilex,3)*tri_ver(7,ntr) + Offx
         post2 = albegaBar(rowtile,tilex,1)*tri_ver(2,ntr) + &
              albegaBar(rowtile,tilex,2)*tri_ver(5,ntr) + &
              albegaBar(rowtile,tilex,3)*tri_ver(8,ntr) + Offy
         post3 = albegaBar(rowtile,tilex,1)*tri_ver(3,ntr) + &
              albegaBar(rowtile,tilex,2)*tri_ver(6,ntr) + &
              albegaBar(rowtile,tilex,3)*tri_ver(9,ntr) + Offz

!       ++++++++Indices of the marker++++++++++++++++++++++++
!       X - indices
        i1=NINT((post1+xyz_tr(1))*dx1) +1
        ist=FLOOR((post1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m

!       Y - indices
        j1=NINT((post2+xyz_tr(2))*dx2) +1
        jst=FLOOR((post2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m

!       Z - indices
        k1=0
        kini=1
        kfin=n3        
        do while (k1.eq.0)
           if ((kfin-kini).eq.1) then
              ztmp=0.5D0*(zc(kfin)+zc(kini))
              if (post3.le.ztmp) then
                 k1=kini
              else
                 k1=kfin
              endif
           endif
           kmid=(kini+kfin)/2
           if (post3.gt.zc(kini).and.post3.le.zc(kmid)) then
              kini=kini
              kfin=kmid
           else
              kini=kmid
              kfin=kfin
           endif
        enddo
        kst=k1
!        if((post3).gt.zc(k1))kst=k1   !fv
        if((post3).le.zc(k1))kst=k1-1 !fv
        
!       OUTER and INNER INDICES
        isto = ist+1 ; isti = ist-1
        jsto = jst+1 ; jsti = jst-1
        ksto = kst+1 ; ksti = kst-1
!       -------------------------------------------------------------
        pindtOff(1,tile)=i1 ; pindtOff(2,tile)=j1 ; pindtOff(3,tile)=k1
        pindtOff(4,tile)=ist; pindtOff(5,tile)=jst; pindtOff(6,tile)=kst
!       -------------------------------------------------------------
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP
      
       endif


      
      return
      end
 

! !! findminimum without bisection !slower
! #ifdef USE_CUDA
!       !$cuf kernel do (1)
! #endif
!       do tile = 1,nttot

!         ntr = faceid_t(tile)

!         n_tile = tri_tiling(ntr,1)
!         rowtile = tri_tiling(ntr,2)

!         tilex = tile - tstart(ntr) + 1

!         post1 = albegaBar(rowtile,tilex,1)*tri_ver(1,ntr) + &
!                 albegaBar(rowtile,tilex,2)*tri_ver(4,ntr) + &
!                 albegaBar(rowtile,tilex,3)*tri_ver(7,ntr)
!         post2 = albegaBar(rowtile,tilex,1)*tri_ver(2,ntr) + &
!                 albegaBar(rowtile,tilex,2)*tri_ver(5,ntr) + &
!                 albegaBar(rowtile,tilex,3)*tri_ver(8,ntr)
!         post3 = albegaBar(rowtile,tilex,1)*tri_ver(3,ntr) + &
!                 albegaBar(rowtile,tilex,2)*tri_ver(6,ntr) + &
!                 albegaBar(rowtile,tilex,3)*tri_ver(9,ntr)

! !       ++++++++Indices of the marker++++++++++++++++++++++++
! !       X - indices
!         i1=NINT((post1+xyz_tr(1))*dx1) +1
!         ist=FLOOR((post1+xyz_tr(1))*dx1) +1      
!         if(ist.eq.0)ist=n1m

! !       Y - indices
!         j1=NINT((post2+xyz_tr(2))*dx2) +1
!         jst=FLOOR((post2+xyz_tr(2))*dx2) + 1
!         if(jst.eq.0)jst=n2m

! !       Z - indices
!         offmin=1E10
!         do k=1,n3
!            off=abs(post3-zc(k))
!            if (off.LT.offmin) then
!               offmin=off
!               k1=k
!            endif
!         enddo
!         kst=k1
! !        if((post3).gt.zc(k1))kst=k1   !fv
!         if((post3).le.zc(k1))kst=k1-1 !fv

        
! !       OUTER and INNER INDICES
!         isto = ist+1 ; isti = ist-1
!         jsto = jst+1 ; jsti = jst-1
!         ksto = kst+1 ; ksti = kst-1
! !       -------------------------------------------------------------
!         pindt(1,tile)=i1 ; pindt(2,tile)=j1 ; pindt(3,tile)=k1
!         pindt(4,tile)=ist; pindt(5,tile)=jst; pindt(6,tile)=kst
! !       -------------------------------------------------------------
!       end do
! !@cuf istat = cudaDeviceSynchronize !JDR TMP


      
