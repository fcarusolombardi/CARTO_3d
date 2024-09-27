
!***********************************************************************
      subroutine initia
      use param
      use local_arrays
      use stat_arrays
      use ibm_param
      implicit none


      write(*,*) " HELLOINITIA"
! new dynamic allocations here
      allocate(xm(0:m1), ym(m2), zm(m3))
      allocate(xc(m1), yc(m2), zc(m3))
      allocate(udx3c(m3), udx3m(m3))
      allocate(imv(m1), ipv(m1))
      allocate(jmv(m2), jpv(m2))
      allocate(kmv(m3), kpv(m3))
      allocate(jmhv(m2+1))
      allocate(ak1(m1), ak2(m2))
      allocate(ap3sk(m3), ac3sk(m3), am3sk(m3))
      allocate(ap3ck(m3), ac3ck(m3), am3ck(m3))
      allocate(amphk(m3), acphk(m3), apphk(m3))

      allocate(indgeo(4,mpun,3), indgeoe(4,mpun,3), indgeoee(4,mpun,3))
      allocate(indgeoAS(4,mpun,3), indgeoeAS(4,mpun,3), indgeoeeAS(4,mpun,3))

      allocate(distb(4,mpun),distbAS(4,mpun))
      
      allocate(aml(max(n3, max(n1,n2))))
      allocate(acl(max(n3, max(n1,n2))))
      allocate(apl(max(n3, max(n1,n2))))

!=================================

      q1=0.d0
      q2=0.d0
      q3=0.d0
      dens=0.d0
      pr=0.d0
      dph=0.d0

      dq=0.d0
      rhs=0.d0
      rhs_t=0.d0
      forclo=0.d0
      forclo_t=0.d0
      ru1=0.d0
      ru2=0.d0
      ru3=0.d0
      ruro=0.d0
      qcap=0.d0
      hro=0.d0

      dsal=0.d0

      hsal=0.d0
      rusal=0.d0
      rhsr=0.d0

!=====================================

      q1_me=0.d0
      q2_me=0.d0
      q3_me=0.d0
      q1_rms=0.d0
      q2_rms=0.d0
      q3_rms=0.d0
      dens_me=0.d0
      dens_rms=0.d0
      dsal_me=0.d0
      dsal_rms=0.d0

      q3dens_me=0.d0
      q3dsal_me=0.d0

      dissuc=0.d0
      disste=0.d0
      disssa=0.d0
      dissur=0.d0
      
!==================================         

      xc = 0.d0
      xm = 0.d0
      yc = 0.d0
      ym = 0.d0
      ap3ck = 0.d0 
      ac3ck = 0.d0 
      am3ck = 0.d0
      ap3sk = 0.d0
      ac3sk = 0.d0
      am3sk = 0.d0
      am3ssk = 0.d0
      ac3ssk = 0.d0
      ap3ssk = 0.d0
      zc = 0.d0
      zm = 0.d0
      g3rc = 0.d0
      g3rm = 0.d0
           
!==================================         

      xcr = 0.d0
      xmr = 0.d0
      ycr = 0.d0
      ymr = 0.d0
      am3sskr = 0.d0
      ac3sskr = 0.d0
      ap3sskr = 0.d0
      zcr = 0.d0
      zmr = 0.d0
      g3rcr = 0.d0
      g3rmr = 0.d0

     
      return 
      end   
