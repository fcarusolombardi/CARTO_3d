!------------------------------------------------------------------
!     Routine to compute indices of all centroids on all particles
!     Only starting and ending indices of the support domain are
!     reqiured. Then loop over it.
!------------------------------------------------------------------

      subroutine findindices
      USE mpih
      USE param
      USE mls_param
!@cuf USE cudafor
      IMPLICIT NONE

      real(DP) pos1,pos2,pos3
      integer i1,j1,k1,ist,jst,kst,inp,ntr
      integer ii,jj,kk
      real(DP) tstr,rcdp,x2dp
      real(DP) a1,a2,alp,alp1
!     --------INNER AND OUTER VARIABLES----------
      integer i1i,j1i,k1i,isti,jsti,ksti
      integer i1o,j1o,k1o,isto,jsto,ksto

      real(DP) tstri,rcdpi,x2dpi
      real(DP) tstro,rcdpo,x2dpo

      real(DP) a1i,a2i,alpi,alp1i
      real(DP) a1o,a2o,alpo,alp1o
!     --------------------------------------------
      real(DP) disdum,distsdum,ztmp
      integer ind_pal,nclipr,n3mor
      integer kok,kini,kfin,kmid,k
!@cuf integer :: istat


      !FV multiple if or elseif inside cuda loop  does not work well!
      if (istr3.EQ.0) then
      
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do ntr=1,nftot

        pos1=tri_bar(1,ntr)
        pos2=tri_bar(2,ntr)
        pos3=tri_bar(3,ntr)


!       ++++++++Indices of the marker++++++++++++++++++++++++	
!       X - indices
        i1=NINT((pos1+xyz_tr(1))*dx1) +1
        ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m

!       Y - indices
        j1=NINT((pos2+xyz_tr(2))*dx2) +1
        jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m

!       Z - indices
!
        k1=NINT((pos3+xyz_tr(3))*dx3) +1
        kst=FLOOR((pos3+xyz_tr(3))*dx3) + 1

!       OUTER and INNER INDICES
        isto = ist+1 ; isti = ist-1
        jsto = jst+1 ; jsti = jst-1
        ksto = kst+1 ; ksti = kst-1
!       dummies
        i1i = 1; i1o = 1
        j1i = 1; j1o = 1
        k1i = 1; k1o = 1
!       -------------------------------------------------------------
        pind(1,ntr)=i1 ; pind(2,ntr)=j1 ; pind(3,ntr)=k1
        pind(4,ntr)=ist; pind(5,ntr)=jst; pind(6,ntr)=kst
!       -------------------------------------------------------------
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP

      
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do ntr=1,nvtot

        pos1=xyz(1,ntr)
        pos2=xyz(2,ntr)
        pos3=xyz(3,ntr)

!       ++++++++Indices of the marker++++++++++++++++++++++++	
!       X - indices
        i1=NINT((pos1+xyz_tr(1))*dx1) +1
        ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m 

!       Y - indices
        j1=NINT((pos2+xyz_tr(2))*dx2) +1
        jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m 

!       Z - indices
!       uniform grid : wall-normal direction
        k1=NINT((pos3+xyz_tr(3))*dx3) +1
        kst=FLOOR((pos3+xyz_tr(3))*dx3) + 1

!       -------------------------------------------------------------
        pindv(1,ntr)=i1 ; pindv(2,ntr)=j1 ; pindv(3,ntr)=k1
        pindv(4,ntr)=ist; pindv(5,ntr)=jst; pindv(6,ntr)=kst
!       -------------------------------------------------------------

      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP

      else !GENERAL CASE
         
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do ntr=1,nftot
        pos1=tri_bar(1,ntr)
        pos2=tri_bar(2,ntr)
        pos3=tri_bar(3,ntr)

!       ++++++++Indices of the marker++++++++++++++++++++++++	
!       X - indices
        i1=NINT((pos1+xyz_tr(1))*dx1) +1
        ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m

!       Y - indices
        j1=NINT((pos2+xyz_tr(2))*dx2) +1
        jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m

!       Z - indices
        k1=0
        kini=1
        kfin=n3        
        do while (k1.eq.0)
           if ((kfin-kini).eq.1) then
              ztmp=0.5D0*(zc(kfin)+zc(kini))
              if (pos3.le.ztmp) then
                 k1=kini
              else
                 k1=kfin
              endif
           endif
           kmid=(kini+kfin)/2
           if (pos3.gt.zc(kini).and.pos3.le.zc(kmid)) then
              kini=kini
              kfin=kmid
           else
              kini=kmid
              kfin=kfin
           endif
        enddo
        kst=k1
!        if((pos3).gt.zc(k1))kst=k1   !fv
        if((pos3).le.zc(k1))kst=k1-1 !fv
     
!       OUTER and INNER INDICES
        isto = ist+1 ; isti = ist-1
        jsto = jst+1 ; jsti = jst-1
        ksto = kst+1 ; ksti = kst-1
!       dummies
        i1i = 1; i1o = 1
        j1i = 1; j1o = 1
        k1i = 1; k1o = 1
!       -------------------------------------------------------------
        pind(1,ntr)=i1 ; pind(2,ntr)=j1 ; pind(3,ntr)=k1
        pind(4,ntr)=ist; pind(5,ntr)=jst; pind(6,ntr)=kst
!       -------------------------------------------------------------
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP

      
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do ntr=1,nvtot

        pos1=xyz(1,ntr)
        pos2=xyz(2,ntr)
        pos3=xyz(3,ntr)

!       ++++++++Indices of the marker++++++++++++++++++++++++	
!       X - indices
        i1=NINT((pos1+xyz_tr(1))*dx1) +1
        ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m 

!       Y - indices
        j1=NINT((pos2+xyz_tr(2))*dx2) +1
        jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m 

!       Z - indices
        k1=0
        kini=1
        kfin=n3        
        do while (k1.eq.0)
           if ((kfin-kini).eq.1) then
              ztmp=0.5D0*(zc(kfin)+zc(kini))
              if (pos3.le.ztmp) then
                 k1=kini
              else
                 k1=kfin
              endif
           endif
           kmid=(kini+kfin)/2
           if (pos3.gt.zc(kini).and.pos3.le.zc(kmid)) then
              kini=kini
              kfin=kmid
           else
              kini=kmid
              kfin=kfin
           endif
        enddo
        kst=k1
!        if((pos3).gt.zc(k1))kst=k1   !fv
        if((pos3).le.zc(k1))kst=k1-1 !fv

!       -------------------------------------------------------------
        pindv(1,ntr)=i1 ; pindv(2,ntr)=j1 ; pindv(3,ntr)=k1
        pindv(4,ntr)=ist; pindv(5,ntr)=jst; pindv(6,ntr)=kst
!       -------------------------------------------------------------

      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP

      endif

      ! pos3=tri_bar(3,1)
      ! write(*,*) 'a',pind(3,1),pind(6,1),n3r
      ! write(*,*) pos3,zc(pind(3,1)),abs(pos3-zc(pind(3,1)))
      ! write(*,*) pos3,zc(pind(3,1)-1),abs(pos3-zc(pind(3,1)-1))
      ! write(*,*) pos3,zc(pind(3,1)+1),abs(pos3-zc(pind(3,1)+1))

      return
      end




!------------------------------------------------------------------
!     Routine to compute indices of all centroids on all particles
!     Only starting and ending indices of the support domain are
!     reqiured. Then loop over it.
!------------------------------------------------------------------

      subroutine findindicesFaces
      USE mpih
      USE param
      USE mls_param
!@cuf USE cudafor
      IMPLICIT NONE

      real(DP) pos1,pos2,pos3
      integer i1,j1,k1,ist,jst,kst,inp,ntr
      integer ii,jj,kk
      real(DP) tstr,rcdp,x2dp
      real(DP) a1,a2,alp,alp1
!     --------INNER AND OUTER VARIABLES----------
      integer i1i,j1i,k1i,isti,jsti,ksti
      integer i1o,j1o,k1o,isto,jsto,ksto

      real(DP) tstri,rcdpi,x2dpi
      real(DP) tstro,rcdpo,x2dpo

      real(DP) a1i,a2i,alpi,alp1i
      real(DP) a1o,a2o,alpo,alp1o
!     --------------------------------------------
      real(DP) disdum,distsdum,ztmp
      integer ind_pal,nclipr,n3mor
      integer kok,kini,kfin,kmid,k
!@cuf integer :: istat

      !FV multiple if or elseif inside cuda loop  does not work well!
      if (istr3.EQ.0) then
      
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do ntr=1,nftot

        pos1=tri_bar(1,ntr)
        pos2=tri_bar(2,ntr)
        pos3=tri_bar(3,ntr)


!       ++++++++Indices of the marker++++++++++++++++++++++++	
!       X - indices
        i1=NINT((pos1+xyz_tr(1))*dx1) +1
        ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m

!       Y - indices
        j1=NINT((pos2+xyz_tr(2))*dx2) +1
        jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m

!       Z - indices
!
        k1=NINT((pos3+xyz_tr(3))*dx3) +1
        kst=FLOOR((pos3+xyz_tr(3))*dx3) + 1

!       OUTER and INNER INDICES
        isto = ist+1 ; isti = ist-1
        jsto = jst+1 ; jsti = jst-1
        ksto = kst+1 ; ksti = kst-1
!       dummies
        i1i = 1; i1o = 1
        j1i = 1; j1o = 1
        k1i = 1; k1o = 1
!       -------------------------------------------------------------
        pind(1,ntr)=i1 ; pind(2,ntr)=j1 ; pind(3,ntr)=k1
        pind(4,ntr)=ist; pind(5,ntr)=jst; pind(6,ntr)=kst
!       -------------------------------------------------------------
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP
      
      else !GENERAL CASE
         
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do ntr=1,nftot
        pos1=tri_bar(1,ntr)
        pos2=tri_bar(2,ntr)
        pos3=tri_bar(3,ntr)

!       ++++++++Indices of the marker++++++++++++++++++++++++	
!       X - indices
        i1=NINT((pos1+xyz_tr(1))*dx1) +1
        ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m

!       Y - indices
        j1=NINT((pos2+xyz_tr(2))*dx2) +1
        jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m

!       Z - indices
        k1=0
        kini=1
        kfin=n3        
        do while (k1.eq.0)
           if ((kfin-kini).eq.1) then
              ztmp=0.5D0*(zc(kfin)+zc(kini))
              if (pos3.le.ztmp) then
                 k1=kini
              else
                 k1=kfin
              endif
           endif
           kmid=(kini+kfin)/2
           if (pos3.gt.zc(kini).and.pos3.le.zc(kmid)) then
              kini=kini
              kfin=kmid
           else
              kini=kmid
              kfin=kfin
           endif
        enddo
        kst=k1
!        if((pos3).gt.zc(k1))kst=k1   !fv
        if((pos3).le.zc(k1))kst=k1-1 !fv
     
!       OUTER and INNER INDICES
        isto = ist+1 ; isti = ist-1
        jsto = jst+1 ; jsti = jst-1
        ksto = kst+1 ; ksti = kst-1
!       dummies
        i1i = 1; i1o = 1
        j1i = 1; j1o = 1
        k1i = 1; k1o = 1
!       -------------------------------------------------------------
        pind(1,ntr)=i1 ; pind(2,ntr)=j1 ; pind(3,ntr)=k1
        pind(4,ntr)=ist; pind(5,ntr)=jst; pind(6,ntr)=kst
!       -------------------------------------------------------------
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP
   
      endif

      return
      end




!------------------------------------------------------------------
!     Routine to compute indices of all centroids on all particles
!     Only starting and ending indices of the support domain are
!     reqiured. Then loop over it.
!------------------------------------------------------------------

      subroutine findindicesNodes
      USE mpih
      USE param
      USE mls_param
!@cuf USE cudafor
      IMPLICIT NONE

      real(DP) pos1,pos2,pos3
      integer i1,j1,k1,ist,jst,kst,inp,ntr
      integer ii,jj,kk
      real(DP) tstr,rcdp,x2dp
      real(DP) a1,a2,alp,alp1
!     --------INNER AND OUTER VARIABLES----------
      integer i1i,j1i,k1i,isti,jsti,ksti
      integer i1o,j1o,k1o,isto,jsto,ksto

      real(DP) tstri,rcdpi,x2dpi
      real(DP) tstro,rcdpo,x2dpo

      real(DP) a1i,a2i,alpi,alp1i
      real(DP) a1o,a2o,alpo,alp1o
!     --------------------------------------------
      real(DP) disdum,distsdum,ztmp
      integer ind_pal,nclipr,n3mor
      integer kok,kini,kfin,kmid,k
!@cuf integer :: istat


      !FV multiple if or elseif inside cuda loop  does not work well!
      if (istr3.EQ.0) then     
      
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do ntr=1,nvtot

        pos1=xyz(1,ntr)
        pos2=xyz(2,ntr)
        pos3=xyz(3,ntr)

!       ++++++++Indices of the marker++++++++++++++++++++++++	
!       X - indices
        i1=NINT((pos1+xyz_tr(1))*dx1) +1
        ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m 

!       Y - indices
        j1=NINT((pos2+xyz_tr(2))*dx2) +1
        jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m 

!       Z - indices
!       uniform grid : wall-normal direction
        k1=NINT((pos3+xyz_tr(3))*dx3) +1
        kst=FLOOR((pos3+xyz_tr(3))*dx3) + 1

!       -------------------------------------------------------------
        pindv(1,ntr)=i1 ; pindv(2,ntr)=j1 ; pindv(3,ntr)=k1
        pindv(4,ntr)=ist; pindv(5,ntr)=jst; pindv(6,ntr)=kst
!       -------------------------------------------------------------

      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP

      else !GENERAL CASE
         
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do ntr=1,nvtot

        pos1=xyz(1,ntr)
        pos2=xyz(2,ntr)
        pos3=xyz(3,ntr)

!       ++++++++Indices of the marker++++++++++++++++++++++++	
!       X - indices
        i1=NINT((pos1+xyz_tr(1))*dx1) +1
        ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m 

!       Y - indices
        j1=NINT((pos2+xyz_tr(2))*dx2) +1
        jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m 

!       Z - indices
        k1=0
        kini=1
        kfin=n3        
        do while (k1.eq.0)
           if ((kfin-kini).eq.1) then
              ztmp=0.5D0*(zc(kfin)+zc(kini))
              if (pos3.le.ztmp) then
                 k1=kini
              else
                 k1=kfin
              endif
           endif
           kmid=(kini+kfin)/2
           if (pos3.gt.zc(kini).and.pos3.le.zc(kmid)) then
              kini=kini
              kfin=kmid
           else
              kini=kmid
              kfin=kfin
           endif
        enddo
        kst=k1
!        if((pos3).gt.zc(k1))kst=k1   !fv
        if((pos3).le.zc(k1))kst=k1-1 !fv

!       -------------------------------------------------------------
        pindv(1,ntr)=i1 ; pindv(2,ntr)=j1 ; pindv(3,ntr)=k1
        pindv(4,ntr)=ist; pindv(5,ntr)=jst; pindv(6,ntr)=kst
!       -------------------------------------------------------------

      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP
      endif

      return
      end


      


      subroutine findindicesFacesProbeP
      USE mpih
      USE param
      USE mls_param
!@cuf USE cudafor
      IMPLICIT NONE

      real(DP) pos1,pos2,pos3
      integer i1,j1,k1,ist,jst,kst,inp,ntr
      integer ii,jj,kk
      real(DP) tstr,rcdp,x2dp
      real(DP) a1,a2,alp,alp1
!     --------INNER AND OUTER VARIABLES----------
      integer i1i,j1i,k1i,isti,jsti,ksti
      integer i1o,j1o,k1o,isto,jsto,ksto

      real(DP) tstri,rcdpi,x2dpi
      real(DP) tstro,rcdpo,x2dpo

      real(DP) a1i,a2i,alpi,alp1i
      real(DP) a1o,a2o,alpo,alp1o
!     --------------------------------------------
      real(DP) disdum,distsdum,ztmp
      integer ind_pal,nclipr,n3mor,tmpa
      integer kok,kini,kfin,kmid,k
!@cuf integer :: istat

      !FV multiple if or elseif inside cuda loop  does not work well!
      if (istr3.EQ.0) then
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do ntr=1,nftot
        tmpa=ntr
        pos1=tri_bar(1,ntr)+ h_probe * tri_nor(1,ntr) * Hboxx(pind(6,ntr))
        pos2=tri_bar(2,ntr)+ h_probe * tri_nor(2,ntr) * Hboxx(pind(6,ntr))
        pos3=tri_bar(3,ntr)+ h_probe * tri_nor(3,ntr) * Hboxx(pind(6,ntr))


!       ++++++++Indices of the marker++++++++++++++++++++++++	
!       X - indices
        i1=NINT((pos1+xyz_tr(1))*dx1) +1
        ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m

!       Y - indices
        j1=NINT((pos2+xyz_tr(2))*dx2) +1
        jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m

!       Z - indices
!
        k1=NINT((pos3+xyz_tr(3))*dx3) +1
        kst=FLOOR((pos3+xyz_tr(3))*dx3) + 1

!       OUTER and INNER INDICES
        isto = ist+1 ; isti = ist-1
        jsto = jst+1 ; jsti = jst-1
        ksto = kst+1 ; ksti = kst-1
!       dummies
        i1i = 1; i1o = 1
        j1i = 1; j1o = 1
        k1i = 1; k1o = 1
!       -------------------------------------------------------------
        pind_probeP(1,ntr)=i1 ; pind_probeP(2,ntr)=j1 ; pind_probeP(3,ntr)=k1
        pind_probeP(4,ntr)=ist; pind_probeP(5,ntr)=jst; pind_probeP(6,ntr)=kst
!       -------------------------------------------------------------
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP
      
      else !GENERAL CASE
         
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
        do ntr=1,nftot
        pos1=tri_bar(1,ntr)+ h_probe * tri_nor(1,ntr) * Hboxx(pind(6,ntr))
        pos2=tri_bar(2,ntr)+ h_probe * tri_nor(2,ntr) * Hboxx(pind(6,ntr))
        pos3=tri_bar(3,ntr)+ h_probe * tri_nor(3,ntr) * Hboxx(pind(6,ntr))

!       ++++++++Indices of the marker++++++++++++++++++++++++	
!       X - indices
        i1=NINT((pos1+xyz_tr(1))*dx1) +1
        ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m

!       Y - indices
        j1=NINT((pos2+xyz_tr(2))*dx2) +1
        jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m

!       Z - indices
        k1=0
        kini=1
        kfin=n3        
        do while (k1.eq.0)
           if ((kfin-kini).eq.1) then
              ztmp=0.5D0*(zc(kfin)+zc(kini))
              if (pos3.le.ztmp) then
                 k1=kini
              else
                 k1=kfin
              endif
           endif
           kmid=(kini+kfin)/2
           if (pos3.gt.zc(kini).and.pos3.le.zc(kmid)) then
              kini=kini
              kfin=kmid
           else
              kini=kmid
              kfin=kfin
           endif
        enddo
        kst=k1
!        if((pos3).gt.zc(k1))kst=k1   !fv
        if((pos3).le.zc(k1))kst=k1-1 !fv
     
!       OUTER and INNER INDICES
        isto = ist+1 ; isti = ist-1
        jsto = jst+1 ; jsti = jst-1
        ksto = kst+1 ; ksti = kst-1
!       dummies
        i1i = 1; i1o = 1
        j1i = 1; j1o = 1
        k1i = 1; k1o = 1
!       -------------------------------------------------------------
        pind_probeP(1,ntr)=i1 ; pind_probeP(2,ntr)=j1 ; pind_probeP(3,ntr)=k1
        pind_probeP(4,ntr)=ist; pind_probeP(5,ntr)=jst; pind_probeP(6,ntr)=kst
!       -------------------------------------------------------------
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP
   
      endif

      return
      end


      subroutine findindicesFacesProbeN
      USE mpih
      USE param
      USE mls_param
!@cuf USE cudafor
      IMPLICIT NONE

      real(DP) pos1,pos2,pos3
      integer i1,j1,k1,ist,jst,kst,inp,ntr
      integer ii,jj,kk
      real(DP) tstr,rcdp,x2dp
      real(DP) a1,a2,alp,alp1
!     --------INNER AND OUTER VARIABLES----------
      integer i1i,j1i,k1i,isti,jsti,ksti
      integer i1o,j1o,k1o,isto,jsto,ksto

      real(DP) tstri,rcdpi,x2dpi
      real(DP) tstro,rcdpo,x2dpo

      real(DP) a1i,a2i,alpi,alp1i
      real(DP) a1o,a2o,alpo,alp1o
!     --------------------------------------------
      real(DP) disdum,distsdum,ztmp
      integer ind_pal,nclipr,n3mor
      integer kok,kini,kfin,kmid,k
!@cuf integer :: istat

      !FV multiple if or elseif inside cuda loop  does not work well!
      if (istr3.EQ.0) then
      
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do ntr=1,nftot

        pos1=tri_bar(1,ntr)- h_probe * tri_nor(1,ntr) * Hboxx(pind(6,ntr))
        pos2=tri_bar(2,ntr)- h_probe * tri_nor(2,ntr) * Hboxx(pind(6,ntr))
        pos3=tri_bar(3,ntr)- h_probe * tri_nor(3,ntr) * Hboxx(pind(6,ntr))


!       ++++++++Indices of the marker++++++++++++++++++++++++	
!       X - indices
        i1=NINT((pos1+xyz_tr(1))*dx1) +1
        ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m

!       Y - indices
        j1=NINT((pos2+xyz_tr(2))*dx2) +1
        jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m

!       Z - indices
!
        k1=NINT((pos3+xyz_tr(3))*dx3) +1
        kst=FLOOR((pos3+xyz_tr(3))*dx3) + 1

!       OUTER and INNER INDICES
        isto = ist+1 ; isti = ist-1
        jsto = jst+1 ; jsti = jst-1
        ksto = kst+1 ; ksti = kst-1
!       dummies
        i1i = 1; i1o = 1
        j1i = 1; j1o = 1
        k1i = 1; k1o = 1
!       -------------------------------------------------------------
        pind_probeN(1,ntr)=i1 ; pind_probeN(2,ntr)=j1 ; pind_probeN(3,ntr)=k1
        pind_probeN(4,ntr)=ist; pind_probeN(5,ntr)=jst; pind_probeN(6,ntr)=kst
!       -------------------------------------------------------------
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP
      
      else !GENERAL CASE
         
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
        do ntr=1,nftot
        pos1=tri_bar(1,ntr)- h_probe * tri_nor(1,ntr) * Hboxx(pind(6,ntr))
        pos2=tri_bar(2,ntr)- h_probe * tri_nor(2,ntr) * Hboxx(pind(6,ntr))
        pos3=tri_bar(3,ntr)- h_probe * tri_nor(3,ntr) * Hboxx(pind(6,ntr))

!       ++++++++Indices of the marker++++++++++++++++++++++++	
!       X - indices
        i1=NINT((pos1+xyz_tr(1))*dx1) +1
        ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m

!       Y - indices
        j1=NINT((pos2+xyz_tr(2))*dx2) +1
        jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m

!       Z - indices
        k1=0
        kini=1
        kfin=n3        
        do while (k1.eq.0)
           if ((kfin-kini).eq.1) then
              ztmp=0.5D0*(zc(kfin)+zc(kini))
              if (pos3.le.ztmp) then
                 k1=kini
              else
                 k1=kfin
              endif
           endif
           kmid=(kini+kfin)/2
           if (pos3.gt.zc(kini).and.pos3.le.zc(kmid)) then
              kini=kini
              kfin=kmid
           else
              kini=kmid
              kfin=kfin
           endif
        enddo
        kst=k1
!        if((pos3).gt.zc(k1))kst=k1   !fv
        if((pos3).le.zc(k1))kst=k1-1 !fv
     
!       OUTER and INNER INDICES
        isto = ist+1 ; isti = ist-1
        jsto = jst+1 ; jsti = jst-1
        ksto = kst+1 ; ksti = kst-1
!       dummies
        i1i = 1; i1o = 1
        j1i = 1; j1o = 1
        k1i = 1; k1o = 1
!       -------------------------------------------------------------
        pind_probeN(1,ntr)=i1 ; pind_probeN(2,ntr)=j1 ; pind_probeN(3,ntr)=k1
        pind_probeN(4,ntr)=ist; pind_probeN(5,ntr)=jst; pind_probeN(6,ntr)=kst
!       -------------------------------------------------------------
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP
   
      endif

      return
      end
      





      


! ! !! findminimum without bisection !slower
! #ifdef USE_CUDA
!       !$cuf kernel do (1)
! #endif
!       do ntr=1,nftot
!         pos1=tri_bar(1,ntr)
!         pos2=tri_bar(2,ntr)
!         pos3=tri_bar(3,ntr)

! !       ++++++++Indices of the marker++++++++++++++++++++++++	
! !       X - indices
!         i1=NINT((pos1+xyz_tr(1))*dx1) +1
!         ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
!         if(ist.eq.0)ist=n1m

! !       Y - indices
!         j1=NINT((pos2+xyz_tr(2))*dx2) +1
!         jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
!         if(jst.eq.0)jst=n2m

! !       Z - indices
!         offmin=1E10
!         do k=1,n3
!            off=abs(pos3-zc(k))
!            if (off.LT.offmin) then
!               offmin=off
!               k1=k
!            endif
!         enddo
!         kst=k1
!        ! if((pos3).gt.zc(k1))kst=k1   !fv
!         if((pos3).le.zc(k1))kst=k1-1 !fv
     
! !       OUTER and INNER INDICES
!         isto = ist+1 ; isti = ist-1
!         jsto = jst+1 ; jsti = jst-1
!         ksto = kst+1 ; ksti = kst-1
! !       dummies
!         i1i = 1; i1o = 1
!         j1i = 1; j1o = 1
!         k1i = 1; k1o = 1
! !       -------------------------------------------------------------
!         pind(1,ntr)=i1 ; pind(2,ntr)=j1 ; pind(3,ntr)=k1
!         pind(4,ntr)=ist; pind(5,ntr)=jst; pind(6,ntr)=kst
! !       -------------------------------------------------------------
!       end do
! !@cuf istat = cudaDeviceSynchronize !JDR TMP

      
! #ifdef USE_CUDA
!       !$cuf kernel do (1)
! #endif
!       do ntr=1,nvtot

!         pos1=xyz(1,ntr)
!         pos2=xyz(2,ntr)
!         pos3=xyz(3,ntr)

! !       ++++++++Indices of the marker++++++++++++++++++++++++	
! !       X - indices
!         i1=NINT((pos1+xyz_tr(1))*dx1) +1
!         ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
!         if(ist.eq.0)ist=n1m 

! !       Y - indices
!         j1=NINT((pos2+xyz_tr(2))*dx2) +1
!         jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
!         if(jst.eq.0)jst=n2m 

! !       Z - indices
!         offmin=1E10
!         do k=1,n3
!            off=abs(pos3-zc(k))
!            if (off.LT.offmin) then
!               offmin=off
!               k1=k
!            endif
!         enddo
!         kst=k1
! !        if((pos3).gt.zc(k1))kst=k1   !fv
!         if((pos3).le.zc(k1))kst=k1-1 !fv

! !       -------------------------------------------------------------
!         pindv(1,ntr)=i1 ; pindv(2,ntr)=j1 ; pindv(3,ntr)=k1
!         pindv(4,ntr)=ist; pindv(5,ntr)=jst; pindv(6,ntr)=kst
! !       -------------------------------------------------------------

!       end do
! !@cuf istat = cudaDeviceSynchronize !JDR TMP








      
!analytical formula for Chebyshev grid      
!       elseif (istr3.EQ.6) then !FV not always correct 
! #ifdef USE_CUDA
!       !$cuf kernel do (1)
! #endif
!       do ntr=1,nftot

!         pos1=tri_bar(1,ntr)
!         pos2=tri_bar(2,ntr)
!         pos3=tri_bar(3,ntr)


! !       ++++++++Indices of the marker++++++++++++++++++++++++	
! !       X - indices
!         i1=NINT((pos1+xyz_tr(1))*dx1) +1
!         ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
!         if(ist.eq.0)ist=n1m

! !       Y - indices
!         j1=NINT((pos2+xyz_tr(2))*dx2) +1
!         jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
!         if(jst.eq.0)jst=n2m

! !       Z - indices
!         !a1=etazm(1+int(str3)); a2=etazm(n3+int(str3))
!         nclipr = int(str3)*mref3
!         n3mor = n3r+nclipr+nclipr 
!         a1 = dcos(PI*(dble(1+int(str3))-0.5d0)/dble(n3mor))
!         a2 = dcos(PI*(dble(n3+int(str3))-0.5d0)/dble(n3mor))
        
!         alp=(pos3+xyz_tr(3))/(0.5D0*alx3)
!         alp1=(1.0D0-alp)*(0.5D0*(a1-a2))
        
!         ! k1=0.5D0+((1.0D0/PI)*(float(n3)+str3+str3)*acos(alp1))
!         k1=0.5D0+((1.0D0/PI)*(real(n3)+str3+str3)*acos(alp1))          
!         k1=int(k1-str3)
!         !staggered
!         kst=k1
! !        if((pos3).gt.zc(k1))kst=k1   !fv
!         if((pos3).le.zc(k1))kst=k1-1 !fv

        
        
! !       OUTER and INNER INDICES
!         isto = ist+1 ; isti = ist-1
!         jsto = jst+1 ; jsti = jst-1
!         ksto = kst+1 ; ksti = kst-1
! !       dummies
!         i1i = 1; i1o = 1
!         j1i = 1; j1o = 1
!         k1i = 1; k1o = 1
! !       -------------------------------------------------------------
!         pind(1,ntr)=i1 ; pind(2,ntr)=j1 ; pind(3,ntr)=k1
!         pind(4,ntr)=ist; pind(5,ntr)=jst; pind(6,ntr)=kst
! !       -------------------------------------------------------------
!       end do
! !@cuf istat = cudaDeviceSynchronize !JDR TMP

      
! #ifdef USE_CUDA
!       !$cuf kernel do (1)
! #endif
!       do ntr=1,nvtot

!         pos1=xyz(1,ntr)
!         pos2=xyz(2,ntr)
!         pos3=xyz(3,ntr)

! !       ++++++++Indices of the marker++++++++++++++++++++++++	
! !       X - indices
!         i1=NINT((pos1+xyz_tr(1))*dx1) +1
!         ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
!         if(ist.eq.0)ist=n1m 

! !       Y - indices
!         j1=NINT((pos2+xyz_tr(2))*dx2) +1
!         jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
!         if(jst.eq.0)jst=n2m 

! !       Z - indices
!         !a1=etazm(1+int(str3)); a2=etazm(n3+int(str3))
!         nclipr = int(str3)*mref3
!         n3mor = n3r+nclipr+nclipr
!         a1 = dcos(PI*(dble(1+int(str3))-0.5d0)/dble(n3mor))
!         a2 = dcos(PI*(dble(n3+int(str3))-0.5d0)/dble(n3mor))
        
!         alp=(pos3+xyz_tr(3))/(0.5D0*alx3)
!         alp1=(1.0D0-alp)*(0.5D0*(a1-a2))
        
!         ! k1=0.5D0+((1.0D0/PI)*(float(n3)+str3+str3)*acos(alp1))
!         k1=0.5D0+((1.0D0/PI)*(real(n3)+str3+str3)*acos(alp1))          
!         k1=int(k1-str3)
!         !staggered
!         kst=k1
!  !       if((pos3).gt.zc(k1))kst=k1   !fv
!         if((pos3).le.zc(k1))kst=k1-1 !fv
      
! !       -------------------------------------------------------------
!         pindv(1,ntr)=i1 ; pindv(2,ntr)=j1 ; pindv(3,ntr)=k1
!         pindv(4,ntr)=ist; pindv(5,ntr)=jst; pindv(6,ntr)=kst
! !       -------------------------------------------------------------
!       end do
! !@cuf istat = cudaDeviceSynchronize !JDR TMP

