!***********************************************************************
      subroutine continua_str
 
      use param
      use mls_param
      use ibm_param

      IMPLICIT NONE
 
      integer i,j,k,inp
      integer ndims,itime,ntr

      real(DP) :: tprfi
      real(DP) :: OH,HP,OP,A,B,C
      real(DP) :: x0,y0,z0,xH,yH,zH,xP,yP,zP
      character*70 filnam1,filnam2,filnamxdm
      character*70 strucfilename
      character*100 ipfi,ipfip,ipfihr

      ! itime=nint(1000.*time)          !FV                                                                                  
      itime=nint(1000.*time*TSTAR) !cosi salva in ms  
      itime=max(itime,0)     
      write(ipfi,199) itime
 199   format(i8.8) 
      filnam1='restart/continua_str'//trim(ipfi)//'.dat'
      open(111,file=filnam1,form='unformatted')
      rewind(111)
        write(111)xyz(1:3,1:nvtot)
        write(111)xyzv(1:3,1:nvtot)
        write(111)xyzv0(1:3,1:nvtot)
        write(111)xyza(1:3,1:nvtot)
        write(111)xyza0(1:3,1:nvtot)
        write(111)fxyz(1:3,1:nvtot)
        write(111)HR
        write(111)xyz_3d(1:3,1:nvtot_3d)
        write(111)xyzv_3d(1:3,1:nvtot_3d)
        write(111)xyzv0_3d(1:3,1:nvtot_3d)
        write(111)xyza_3d(1:3,1:nvtot_3d)
        write(111)xyza0_3d(1:3,1:nvtot_3d)
        write(111)fxyz_3d(1:3,1:nvtot_3d)
        write(111)astressEFcell_3d(1:nctot_3d)
        write(111)EFtstart_3d(1:nvtot_3d)
#ifdef ELECTRO
        write(111)XEF_1d(1:21,1:nvtotEF_1d)
        write(111)XEF_2d(1:20,1:nftotEF_2d)
#ifdef TP06
        write(111)XEF_3d(1:22,1:nctot_3d)
#endif
#ifdef MINIMAL_MODEL
        write(111)XEF_3d(1:4,1:nctot_3d)
#endif
        write(111)IstimEF_3dS1(1:nctot_3d)
        write(111)IstimEF_3dS2(1:nctot_3d)
        write(111)Purk_tstim(1:nPurkVentr)
        write(111)Bund_tstim(1:nBundAtri)
#endif
      close(111)

      !Save ipfi for bash launcher
      open(222,file='restart/TR_HR.txt')
      write(222,*)itime
      close(222)
      
      return
      end


!***********************************************************************
      subroutine continua_ecg
 
      use param
      use mls_param
      use ibm_param

      IMPLICIT NONE
 
      integer i,j,k,inp
      integer ndims,itime,ntr

      real(DP) :: tprfi
      real(DP) :: OH,HP,OP,A,B,C
      real(DP) :: x0,y0,z0,xH,yH,zH,xP,yP,zP
      character*70 filnam1,filnam2,filnamxdm
      character*70 strucfilename
      character*100 ipfi,ipfip

      ! itime=nint(1000.*time)          !FV                                                                                  
      itime=nint(1000.*time*TSTAR) !cosi salva in ms  
      itime=max(itime,0)
      write(ipfi,199) itime
 199   format(i8.8) 
      filnam1='restart/continua_ecg'//trim(ipfi)//'.dat'
      open(111,file=filnam1,form='unformatted')
      rewind(111)
      write(111)time
      write(111)xyz(1:3,1:nvtot)
      write(111)xyz_3d(1:3,1:nvtot_3d)
      write(111)potEFnode_1d(1:nvtotEF_1d)
      write(111)potEFface_2d(1:nftotEF_2d)
      write(111)potEFcell_3d(1:nctot_3d)
      close(111)


      return
      end


