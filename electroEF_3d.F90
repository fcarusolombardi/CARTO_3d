!===================================================
      subroutine ElectroRunEF_3d(timex,nBeat)
      use constants
      use mls_param
      USE ieee_arithmetic
!@cuf   use cudafor, cuf_minloc => minloc, cuf_maxloc => maxloc 
      implicit none

      !-------------------------------------------
      ! Fibrillation Variables Declaration

      character(len=4):: dummy
      real(DP):: xS1,yS1,zS1
      real(DP):: xS2,yS2,zS2
      real(DP):: amplitudeS1,amplitudeS2
      real(DP):: duration_signal_S1,duration_signal_S2,delayS2
      real(DP):: pot_apd,pot_old,treshold_apd
      integer:: nBeat,step
      !-------------------------------------------

      real(DP):: timex,timeEFv,alpv1,timeMS,dtMS
      real(DP):: timeVentr,timeAtrio,errorEFbido
      real(DP):: rhspot,num,den,numext,dvvjj,potsur,potextsur,volcell
      real(DP):: distCF,g1interp,corrsten,bordo,potv,tempov
      real(DP):: xC1,yC1,zC1,xV1,yV1,zV1
      real(DP):: xgrad,ygrad,zgrad
      real(DP):: xgradext,ygradext,zgradext
      real(DP):: xvers,yvers,zvers
      real(DP):: xnormale,ynormale,znormale
      real(DP):: mint11,mint12,mint13,mint21,mint22,mint23,mint31,mint32,mint33
      real(DP):: mext11,mext12,mext13,mext21,mext22,mext23,mext31,mext32,mext33
      real(DP):: msum11,msum12,msum13,msum21,msum22,msum23,msum31,msum32,msum33

      real(DP):: Volt,Cai,CaSR,CaSS
      real(DP):: Nai,Ki,INa_m,INa_h 
      real(DP):: INa_j,IKr_xr1,IKr_xr2,IKs_xs 
      real(DP):: Ito_r,Ito_s,ICaL_d,ICaL_f
      real(DP):: ICaL_f2,ICaL_fCaSS,RR

      real(DP):: R,F,T,RTONF
      real(DP):: Ko,Cao,Nao,Vc,Vsr,Vss
      real(DP):: Bufc,Kbufc,Bufsr,Kbufsr,Bufss,Kbufss
      real(DP):: Vmaxup,KupEF,Vrel,k1_,k2_,k3,k4,EC,maxsr,minsr,Vleak,Vxfer
      real(DP):: CAPACITANCE,Gkr,pKNa,Gks,GK1,Gto,GNa,GbNA,KmK,KmNa,knak
      real(DP):: GCaL,GbCa,knaca,KmNai,KmCa,ksat,n,GpCa,KpCa,GpK,Istim

      real(DP):: IKr_xr1_inf,IKr_xr1_a,IKr_xr1_b,IKr_xr1_t,IKr_xr1P1
      real(DP):: IKr_xr2_inf,IKr_xr2_a,IKr_xr2_b,IKr_xr2_t,IKr_xr2P1
      real(DP):: Ek,sxr1,sxr2,IKr

      real(DP):: IKs_xs_inf,IKs_xs_a,IKs_xs_b,IKs_xs_t,IKs_xsP1,Eks,IKs
      real(DP):: IK1_a,IK1_b,rec_iK1,IK1
      real(DP):: Ito_r_inf,Ito_s_inf,Ito_r_t,Ito_s_t,Ito_rP1,Ito_sP1,Ito

      real(DP):: Ena,INa_m_a,INa_m_b,INa_m_t,INa_m_inf,INa_mP1
      real(DP)::     INa_h_a,INa_h_b,INa_h_t,INa_h_inf,INa_hP1
      real(DP)::     INa_j_a,INa_j_b,INa_j_t,INa_j_inf,INa_jP1
      real(DP):: INa,IbNa,rec_iNaK,INaK

      real(DP):: ICaL_d_inf,ICaL_d_a,ICaL_d_b,ICaL_d_c,ICaL_d_t,CaL_dP1
      real(DP):: ICaL_f_inf,ICaL_f_a,ICaL_f_b,ICaL_f_c,ICaL_f_t,CaL_fP1
      real(DP):: ICaL_f2_inf,ICaL_f2_a,ICaL_f2_b,ICaL_f2_c,ICaL_f2_t,CaL_f2P1
      real(DP):: ICaL_fCaSS_inf,ICaL_fCaSS_t,CaL_fCaSSP1,ICaL

      real(DP):: ECa,IbCa,INaCa,IpCa,rec_ipK,IpK,dv,inverseVcF2,inverseVcF,inversevssF2
      real(DP):: kCaSR,k1,k2,alphaRR,RRP1,OO,Irel,Ileak,Iup,Ixfer
      real(DP):: CaCSQN,dCaSR,bjsr,cjsr,CaSSBuf,dCaSS,bcss,ccss
      real(DP):: CaBuf,dCai,bc,cc,dNai,dKi
!      real(DP),dimension(3,1:nctot_3d) :: gradcell_3d
      real(DP):: sta1,sta2,sta3,sta4,sta5,sta6,sta7,sta8,sta9
      real(DP):: sta10,sta11,sta12,sta13,sta14,sta15,sta16,sta17,sta18,sta19,sta20,sta21
      real(DP):: sta1np1,sta2np1,sta3np1,sta4np1,sta5np1,sta6np1,sta7np1,sta8np1,sta9np1
      real(DP):: sta10np1,sta11np1,sta12np1,sta13np1,sta14np1,sta15np1,sta16np1,sta17np1,sta18np1,sta19np1,sta20np1,sta21np1
      real(DP):: rat1,rat2,rat3,rat4,rat5,rat6,rat7,rat8,rat9
      real(DP):: rat10,rat11,rat12,rat13,rat14,rat15,rat16,rat17,rat18,rat19,rat20,rat21
      real(DP):: const1c,const2c,const3c,const4c
      real(DP):: const10c,const11c,const12c,const13c,const14c,const15c,const16c,const17c,const18c,const19c
      real(DP):: const20c,const21c,const22c,const23c,const24c,const25c,const26c,const27c,const28c,const29c
      real(DP):: const30c,const31c,const32c,const33c,const34c,const35c,const36c,const37c,const38c,const39c
      real(DP):: const40c,const41c,const42c,const43c,const44c,const45c,const46c,const47c,const48c,const49c
      real(DP):: alg1,alg2,alg3,alg4,alg5,alg6,alg7,alg8,alg9
      real(DP):: alg10,alg11,alg12,alg13,alg14,alg15,alg16,alg17,alg18,alg19
      real(DP):: alg20,alg21,alg22,alg23,alg24,alg25,alg26,alg27,alg28,alg29
      real(DP):: alg30,alg31,alg32,alg33,alg34,alg35,alg36,alg37,alg38,alg39
      real(DP):: alg40,alg41,alg42,alg43,alg44,alg45,alg46,alg47,alg48,alg49
      real(DP):: alg50,alg51,alg52,alg53,alg54,alg55,alg56,alg57,alg58,alg59
      real(DP):: alg60,alg61,alg62,alg63,alg64,alg65,alg66,alg67,alg68,alg69
      real(DP):: alg70,alg71,alg72,alg73,alg74,alg75,alg76,alg77,alg78,alg79
      real(DP):: xEc,yEc,zEc,xBc,yBc,zBc,Rm1,ecg_gradrm1_x,ecg_gradrm1_y,ecg_gradrm1_z
      real(DP):: xLV,yLV,zLV,xRV,yRV,zRV,xRA,yRA,zRA
      real(DP):: aoo, duration_signal,testECG,torto,timev,timeLead_LV,timeLead_RV,timeLead_RA,distveG2
      integer:: i,j,inp,chamb,chamb1,chamb2,c1,c2,v1,v2,v3,v4,f1,cslave,chambcheck,chambcheckbc
      integer:: vsi,vei,fsi,fei,esi,eei,csi,cei
      character*150 :: stri

      ! Conduction velocity
      !------------------------------------------------------------------------------
      integer:: indmin,indmax,v_max,v_min
      real(DP)::LATx,LATy,LATz,cvmod,num1,num2,num3,latf,cvloc,cvlocx,cvlocy,cvlocz
      real(DP)::cvlocmax,cvlocmin
      real(DP)::xgradlat,ygradlat,zgradlat
      real(DP)::cvf1,cvf2,cvf3
      real(DP)::grad11,grad12,grad13
      real(DP)::grad21,grad22,grad23
      real(DP)::grad31,grad32,grad33
      real(DP)::xV_max,yV_max,zV_max
      real(DP)::xV_min,yV_min,zV_min
      real(DP)::dxlat,dylat,dzlat,modgradlat
      real(DP),dimension(4)::lat_cell_vert
      ! minimal model variables
      !------------------------------------------------------------------------------
      integer::CONST1,CONST2,CONST3,epi_endo_mur
      real(DP):: const4,const5,const6,const7,const8,const9
      real(DP):: const10,const11,const12,const13,const14,const15,const16,const17,const18,const19
      real(DP):: const20,const21,const22,const23,const24,const25,const26,const27,const28,const29
      real(DP):: const30,const31,const32,app1,app2,app3,app4,app5
      real(DP):: remod

      !------------------------------------------------------------------------------
#ifdef MITCHELL_SCHAEFFER
      ! Mitchell-Schaeffer
      !------------------------------------------------------------------------------
      real(DP)::invtau_in,invtau_out,invtau_open,invtau_close,u_gate
      !------------------------------------------------------------------------------
#endif
!@cuf   integer :: istat
!#ifdef USE_CUDA
!        attributes(managed) :: gradcell_3d
!#endif 
#ifdef USE_CUDA
      attributes(managed) :: lat_cell_vert
#endif 

!NOTAFV per considerare l'effetto delle piccole variazione geometriche durante la
!depolarizzazione dovremmo aggiornare (ag ogni time step o multiplo) anche le informazioni
!geometriche del dominio chiamando le routines che sono alla fine di read_geoSingleBody_3d
!tutto è predisposto, manca solo un modo per aggiornare/seguire l'orientamento delle fibre 


!#ifdef S1S2-3D
      open(unit=15,file='S1S2_3dStim.in',status='old')
      read(15,*) dummy
      read(15,*) dummy
      read(15,*) dummy
      read(15,*) dummy
      read(15,*) dummy
      read(15,*) dummy
      read(15,*) dummy
      read(15,*) dummy
      read(15,*) dummy
      read(15,*) duration_signal_S1,duration_signal_S2,delayS2
      close(15)
!#endif
      treshold_apd = -79.0D0
!      delayS2=60.0D0/delayS2*1000.0D0
#ifdef TP06      
!Parameters, should i move it to param?
    !%Constants
    R = 8314.472;
    F = 96485.3415;
    T = 310.0;
    RTONF = (R*T) / F;
    !%External concentrations
    Ko = 5.4;
    Cao = 2.0;
    Nao = 140.0;
    !%Intracellular volumes
    Vc = 0.016404;
    Vsr = 0.001094;
    Vss = 0.00005468;
    !%Calcium buffering dynamics
    Bufc = 0.2;
    Kbufc = 0.001;
    Bufsr = 10.;
    Kbufsr = 0.3;
    Bufss = 0.4;
    Kbufss = 0.00025;
    !%Intracellular calcium flux dynamics
    Vmaxup = 0.006375;
    KupEF = 0.00025;
    Vrel = 0.102;
    k1_ = 0.15;
    k2_ = 0.045;
    k3 = 0.060;
    k4 = 0.005;!%0.000015;
    EC = 1.5;
    maxsr = 2.5;
    minsr = 1.;
    Vleak = 0.00036;
    Vxfer = 0.0038;
    !%Cellular capacitance         
    CAPACITANCE = 0.185;
    !%Parameters for currents
    !%Parameters for IKr
    Gkr = 0.153;
    !%Parameters for Iks
    pKNa = 0.03;
    Gks = 0.392;  !Epi
    ! Gks = Gks*3.D0 !ricambia
    ! molt = 1.5 --> -27ms
    ! molt = 2.0 --> -46ms
    ! molt = 2.5 --> -60ms
    ! molt = 3.0 --> -72ms
    ! molt = 3.5 --> -80ms
    ! molt = 4.0 --> -88ms
    ! molt = 4.5 --> -95ms
    ! molt = 5.0 --> -101ms
!    Gks = 0.392;  !ENDO
!    Gks = 0.098;  !else, myocardium?
    !% //Parameters for Ik1
    GK1 = 5.405;
    !% //Parameters for Ito
    Gto = 0.294;
    ! Gto = Gto*0.5D0 !ricambia
!    Gto = 0.073; !ENDO
!    Gto = 0.294; !else,myocardium?
    !% //Parameters for INa
    GNa = 14.838;
    !% //Parameters for IbNa
    GbNa = 0.00029;
    !% //Parameters for INaK
    KmK = 1.0;
    KmNa = 40.0;
    knak = 2.724;
    !% //Parameters for ICaL
    GCaL = 0.00003980;
    !% //Parameters for IbCa
    GbCa = 0.000592;
    !% //Parameters for INaCa
    knaca = 1000;
    KmNai = 87.5;
    KmCa = 1.38;
    ! KmCa = KmCa*3.D0 !ricambia
    ksat = 0.1;
    n = 0.35;
    !% //Parameters for IpCa
    GpCa = 0.1238;
    KpCa = 0.0005;
    !% //Parameters for IpK;
    GpK = 0.0146;
#endif
#ifdef COURTEMANCHE    
!Courtemanche
    CONST1c = 8.3143;
    CONST2c = 310;
    CONST3c = 96.4867;
    CONST4c = 100;

    CONST10c = 7.8;
    CONST11c = 140;
    CONST12c = 0.09;
    CONST13c = 5.4;
    CONST14c = 3;
    CONST15c = 0.1652;
    CONST16c = 0.029411765;
    CONST17c = 0.12941176;
    CONST18c = 0.12375;
    CONST19c = 10;
    CONST20c = 1.5;
    CONST21c = 0.59933874;
    CONST22c = 0.0006744375;
    CONST23c = 0.001131;
    CONST24c = 0;
    CONST25c = 1.8;
    CONST26c = 1600;
    CONST27c = 87.5;
    CONST28c = 1.38;
    CONST29c = 0.1;
    CONST30c = 0.35;
    CONST31c = 0.275;
    CONST32c = 30;
    CONST33c = 180;      
    CONST34c = 0.005;
    CONST35c = 0.00092;
    CONST36c = 15;
    CONST37c = 0.05;
    CONST38c = 0.07;
    CONST39c = 10;
    CONST40c = 0.00238;
    CONST41c = 0.0005;
    CONST42c = 0.8;
    CONST43c = 20100;
    CONST44c =  CONST43c*0.680000;
    CONST45c = 2.00000;
    CONST46c =  (1.00000/7.00000)*(exp(CONST11c/67.3000) - 1.00000);
    CONST47c = 8.00000;
    CONST48c =  0.00480000*CONST43c;
    CONST49c =  0.0552000*CONST43c;
#endif
#ifdef MITCHELL_SCHAEFFER
    invtau_open  = 1.0D0/100.0D0
    invtau_close = 1.0D0/140.0D0
    invtau_in    = 1.0D0/0.13D0
    invtau_out   = 1.0D0/2.60D0
    u_gate    = 0.10
#endif

    timeMS=timex*TSTAR*1000.d0
    dtMS=dt*TSTAR*1000.d0*real(stepEF_3d)
    ! dtMS=dt*real(stepEF_3d)

! ! !!STIMOLO                
!     duration_signal = 5 !ms ricambia 
!     duration_signal = 5/(1000.0d0*TSTAR) 
!     torto = 5 !ms                                                              
!     torto = torto/(1000.0d0*TSTAR)
! !Stimolo Purkinje
! !coupling Purkinje -> Ventr
!        IstimEF_3d(:)=0.0D0
!        do i=1,nPurkVentr  
!           timev=Purk_tstim(i)
!           if (  (timex.GT.(timev+torto)).AND.(timex.LT.(timev+torto+duration_signal)) ) then
! !             fmaster=fVentrmaster(i)
!              cslave=cVentrslave(i)
!              IstimEF_3d(cslave)=-60.D0
!           endif
!        enddo              

!        do i=1,nBundAtri  
!           timev=Bund_tstim(i)
!           if (  (timex.GT.(timev)).AND.(timex.LT.(timev+duration_signal)) ) then
! !             vmaster=vAtrimaster(i)
!              cslave=cAtrislave(i)
!              IstimEF_3d(cslave)=-60.D0
! !             IstimEF_3d(cslave)=-120.D0
!           endif
!        enddo                      

! !!LEADS
! if (leads.EQ.1) then
! !2-ELETTRODI  !NOTA QUI I TEMI SONO IN ms                                                                                                                                                     !tempi elettrodi                                                                                                                                                                              !   duration_signal = 5.0D0 !ms ricambia                                                                                                                                                      
! duration_signal = 0.5D0 !ms ricambia                                                                                                                                                          
! timeLead_LV=200.0D0!/(TSTAR*1000) ms  !5                                                                                                                                                      
! timeLead_LV=timeLead_LV-30.0D0
! timeLead_RV=timeLead_LV
! timeLead_RA=0!/(TSTAR*1000)  ms !0.5 
! xLV = cell_bar(1,lead_LV)
! yLV = cell_bar(2,lead_LV)
! zLV = cell_bar(3,lead_LV)

! xRV = cell_bar(1,lead_RV)
! yRV = cell_bar(2,lead_RV)
! zRV = cell_bar(3,lead_RV)

! xRA = cell_bar(1,lead_RA)
! yRA = cell_bar(2,lead_RA)
! zRA = cell_bar(3,lead_RA)

! if ((timeMS.GT.timeLead_LV).AND.(timeMS.LT.(timeLead_LV+duration_signal))) then
! #ifdef USE_CUDA
!         !$cuf kernel do (1)                                                                                                                                                                    
! #endif
!       do i=cstart_3d(1),cend_3d(1)
!          chamb = cell_to_chamb_3d(i)
!          xBc = cell_bar(1,i)
!          yBc = cell_bar(2,i)
!          zBc = cell_bar(3,i)

!          if (chamb.EQ.1) then
!             distveG2=sqrt((xBc-xLV)**2+(yBc-yLV)**2+(zBc-zLV)**2)
!             if (distveG2.LT.(5.D0/(1000.D0*LSTAR))) then
!             ! if (distveG2.LT.5.D0) then
!                IstimEF_3d(i)=-60.D0
!             endif
!          endif
!       enddo
! !@cuf istat = cudaDeviceSynchronize !JDR TMP                                                                                                                                                   
! endif


! if ((timeMS.GT.timeLead_RV).AND.(timeMS.LT.(timeLead_RV+duration_signal))) then
! #ifdef USE_CUDA
!         !$cuf kernel do (1)                                                                                                                                                                    
! #endif
!       do i=cstart_3d(1),cend_3d(1)
!          chamb = cell_to_chamb_3d(i)
!          xBc = cell_bar(1,i)
!          yBc = cell_bar(2,i)
!          zBc = cell_bar(3,i)

!          if (chamb.EQ.3) then
!             distveG2=sqrt((xBc-xRV)**2+(yBc-yRV)**2+(zBc-zRV)**2)
!             if (distveG2.LT.(5.D0/(1000.D0*LSTAR))) then
!             ! if (distveG2.LT.5.D0) then
!                IstimEF_3d(i)=-60.D0
!             endif
!          endif
!       enddo
! !@cuf istat = cudaDeviceSynchronize !JDR TMP                                                                                                                                                   
! endif


! if ((timeMS.GT.timeLead_RA).AND.(timeMS.LT.(timeLead_RA+duration_signal))) then
! #ifdef USE_CUDA
!    !$cuf kernel do (1)                                                     
! #endif
!    do i=cstart_3d(1),cend_3d(1)
!       chamb = cell_to_chamb_3d(i)
!       xBc = cell_bar(1,i)
!       yBc = cell_bar(2,i)
!       zBc = cell_bar(3,i)

!       if (chamb.EQ.4) then
!          distveG2=sqrt((xBc-xRA)**2+(yBc-yRA)**2+(zBc-zRA)**2)
!          if (distveG2.LT.(5.D0/(1000.D0*LSTAR))) then
!             ! if (distveG2.LT.5.D0) then
!             IstimEF_3d(i)=-60.D0
!          endif
!       endif
!    enddo
!    !@cuf istat = cudaDeviceSynchronize !JDR TMP                                 
! endif
! endif !if leads

    ! !!!!!!! fine ISTIM!!!!!!

    
! #ifdef USE_CUDA
!         !$cuf kernel do (1)
! #endif
!     do i=cstart_3d(1),cend_3d(1)
!        !if (cell_bar(1,i).lt.0.5D0) then
!        !if (LabelSmooth_cell(i).lt.0.7D0) then 
!           cell_to_chamb_3d(i)=1
!        !endif
!     enddo
!     !@cuf istat = cudaDeviceSynchronize !JDR TMP

!!COMPUTE GRADIENT
!find nodal values
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
!do i=vstart_3d(1),vend_3d(1)
       do i=vstart_3d(1),vend_3d(1)
        num=0
#ifdef BIDOMAIN
        numext=0
#endif
        den=0
#ifdef ELEGEO0
       xV1 = xyz0_3d(1,i)
       yV1 = xyz0_3d(2,i)
       zV1 = xyz0_3d(3,i)
#else
       xV1 = xyz_3d(1,i)
       yV1 = xyz_3d(2,i)
       zV1 = xyz_3d(3,i)
#endif
       
       do j=1,n_cell_of_vert_3d(i)
            c1=cell_of_vert_3d(j,i)
            if (cell_to_chamb_3d(c1).lt.5) then!ne.9) then !MOD_BC
               xC1 = cell_bar(1,c1)
               yC1 = cell_bar(2,c1)
               zC1 = cell_bar(3,c1)

               dvvjj = sqrt( (xV1-xC1)**2+(yV1-yC1)**2+(zV1-zC1)**2)
               num = num + potEFcell_3d(c1)/dvvjj
               den = den + 1.D0/dvvjj
#ifdef BIDOMAIN
               numext = numext + potextEFcell_3d(c1)/dvvjj
#endif
            endif !MOD_BC
        enddo
        if (den.NE.0.0D0) then
           potEFnode_3d(i)=num/den
        else
           potEFnode_3d(i)=-83.0D0
        endif
#ifdef BIDOMAIN
        if (den.NE.0.0D0) then
           potextEFnode_3d(i)=numext/den
        else
           potextEFnode_3d(i)=-83.0D0
        endif
#endif
    enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!find face values
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
    do i=fstart_3d(1),fend_3d(1)
       v1=vert_of_face_3d(1,i)
       v2=vert_of_face_3d(2,i)
       v3=vert_of_face_3d(3,i)
       potEFface_3d(i)=(potEFnode_3d(v1)+potEFnode_3d(v2)+potEFnode_3d(v3))/3.D0;
#ifdef BIDOMAIN
       potextEFface_3d(i)=(potextEFnode_3d(v1)+potextEFnode_3d(v2)+potextEFnode_3d(v3))/3.D0;
#endif
    enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

! !find gradient on cells (gauss-green theorem)
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=cstart_3d(1),cend_3d(1)
        xgrad=0.0D0
        ygrad=0.0D0
        zgrad=0.0D0
#ifdef BIDOMAIN
        xgradext=0.0D0
        ygradext=0.0D0
        zgradext=0.0D0
#endif
        do j=1,4
           f1=face_of_cell_3d(j,i)
           xnormale=normalfaceofcells_3d(1,j,i)
           ynormale=normalfaceofcells_3d(2,j,i)
           znormale=normalfaceofcells_3d(3,j,i)

#ifdef ELEGEO0
          potsur = potEFface_3d(f1)*sur0_3d(f1)
#else
          potsur = potEFface_3d(f1)*sur_3d(f1)
#endif
           xgrad = xgrad + potsur*xnormale
           ygrad = ygrad + potsur*ynormale
           zgrad = zgrad + potsur*znormale
#ifdef BIDOMAIN
#ifdef ELEGEO0           
           potextsur = potextEFface_3d(f1)*sur0_3d(f1)
#else
           potextsur = potextEFface_3d(f1)*sur_3d(f1)
#endif           
           xgradext = xgradext + potextsur*xnormale
           ygradext = ygradext + potextsur*ynormale
           zgradext = zgradext + potextsur*znormale
#endif
        enddo
#ifdef ELEGEO0
        volcell = vol0_3d(i)
#else
        volcell = vol_3d(i)
#endif
        gradcell_3d(1,i) =  xgrad/volcell
        gradcell_3d(2,i) =  ygrad/volcell
        gradcell_3d(3,i) =  zgrad/volcell
#ifdef BIDOMAIN
        gradextcell_3d(1,i) =  xgradext/volcell
        gradextcell_3d(2,i) =  ygradext/volcell
        gradextcell_3d(3,i) =  zgradext/volcell
#endif

     enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!find gradient on faces
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
     do i=fstart_3d(1),fend_3d(1)
        c1 = cell_of_face_3d(1,i)
        c2 = cell_of_face_3d(2,i)
        
        g1interp = g1interpface_3d(i)
        xvers  = versCFface_3d(1,i)
        yvers  = versCFface_3d(2,i)
        zvers  = versCFface_3d(3,i)
        distCF = distCFface_3d(i)
        if ((c1.NE.0).AND.(c2.NE.0)) then
           chamb1 = cell_to_chamb_3d(c1)
           chamb2 = cell_to_chamb_3d(c2)
           if ((chamb1.lt.5).AND.(chamb2.lt.5)) then
              !grandiente medio
              xgrad = g1interp*gradcell_3d(1,c1) + (1-g1interp)*gradcell_3d(1,c2)
              ygrad = g1interp*gradcell_3d(2,c1) + (1-g1interp)*gradcell_3d(2,c2)
              zgrad = g1interp*gradcell_3d(3,c1) + (1-g1interp)*gradcell_3d(3,c2)
              !gradiente corretto
              !gradientecorretto=gradientemedio+(-dot(gradientemedio,versCF)+(pot(c2)-pot(c1))/distCF)*versCF;
              corrsten= -(xgrad*xvers+ygrad*yvers+zgrad*zvers) + (potEFcell_3d(c2)-potEFcell_3d(c1))/distCF
              xgrad = xgrad + corrsten*xvers
              ygrad = ygrad + corrsten*yvers
              zgrad = zgrad + corrsten*zvers           
              
#ifdef BIDOMAIN
              xgradext = g1interp*gradextcell_3d(1,c1) + (1-g1interp)*gradextcell_3d(1,c2)
              ygradext = g1interp*gradextcell_3d(2,c1) + (1-g1interp)*gradextcell_3d(2,c2)
              zgradext = g1interp*gradextcell_3d(3,c1) + (1-g1interp)*gradextcell_3d(3,c2)
              !gradiente corretto
              corrsten= -(xgradext*xvers+ygradext*yvers+zgradext*zvers) + (potextEFcell_3d(c2)-potextEFcell_3d(c1))/distCF
              xgradext = xgradext + corrsten*xvers
              ygradext = ygradext + corrsten*yvers
              zgradext = zgradext + corrsten*zvers           
#endif
           else
              xgrad = 0.0D0 !for Neumann bcs
              ygrad = 0.0D0 
              zgrad = 0.0D0 
#ifdef BIDOMAIN
              xgradext = 0.0D0 !for Neumann bcs
              ygradext = 0.0D0 
              zgradext = 0.0D0 
#endif
           endif
        elseif ((c1.NE.0).AND.(c2.EQ.0)) then
           ! xgrad = gradcell_3d(1,c1)
           ! ygrad = gradcell_3d(2,c1)
           ! zgrad = gradcell_3d(3,c1)
           xgrad = 0.0D0 !for Neumann bcs
           ygrad = 0.0D0 
           zgrad = 0.0D0 
#ifdef BIDOMAIN
           xgradext = 0.0D0 !for Neumann bcs
           ygradext = 0.0D0 
           zgradext = 0.0D0 
#endif
        elseif ((c1.EQ.0).AND.(c2.NE.0)) then
           ! xgrad = gradcell_3d(1,c2)
           ! ygrad = gradcell_3d(2,c2)
           ! zgrad = gradcell_3d(3,c2)
           xgrad = 0.0D0 !for Neumann bcs
           ygrad = 0.0D0 
           zgrad = 0.0D0 
#ifdef BIDOMAIN
           xgradext = 0.0D0 !for Neumann bcs
           ygradext = 0.0D0 
           zgradext = 0.0D0 
#endif
        endif
        gradface_3d(1,i) = xgrad
        gradface_3d(2,i) = ygrad
        gradface_3d(3,i) = zgrad
#ifdef BIDOMAIN
        gradextface_3d(1,i) = xgradext
        gradextface_3d(2,i) = ygradext
        gradextface_3d(3,i) = zgradext
#endif
     enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!CELL MODEL
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=cstart_3d(1),cend_3d(1)
         !         write(*,*) "LS", LabelStenosi(i)
         chamb = cell_to_chamb_3d(i) !ricambiaFV
         if (LabelStenosi(i).EQ.0) then
         !Stimulus (Bundles+Purkinje+S1S2_Stim)
         ! if ((timeMS.GT.0).AND.(timeMS.LT.duration_signal)) then
         !    Istim = IstimLV_3d(i)
         ! else
         !    Istim = 0.0D0
         ! endif
         !!Istim = 0.0D0 !ricambia
         !Istim=IstimEF_3d(i)
         if (nBeat.eq.0)then
           !  if ((timeMS.GT.0.0D0).AND.(timeMS.LT.duration_signal_S1)) then
         !       !Istim = IstimEF_3d(i)+IstimEF_3dS1(i)
         !       Istim = IstimEF_3dS1(i)
         !    ! elseif ((timeMS.GT.delayS2).AND.(timeMS.LT.delayS2+duration_signal_S2)) then
         !    !     !Istim = IstimEF_3d(i)+IstimEF_3dS2(i)
         !    !     !Istim = IstimEF_3dS2(i)
         !    !    Istim = IstimEF_3dS2(i)
         !    else
         !       !Istim = IstimEF_3d(i)
         !       Istim = 0.0D0
         !    endif
         ! elseif (nBeat.eq.1)then
            !if ((timeMS.GT.0.0D0).AND.(timeMS.LT.0.0D0+duration_signal_S1)) then
               !Istim = IstimEF_3datria(i)!+IstimEF_3dS1(i)
               !Istim = IstimEF_3dS1(i)
            if ((timeMS.GT.0.0D0).AND.(timeMS.LT.0.0D0+duration_signal_S1)) then
               !Istim = IstimEF_3d(i)!+IstimEF_3dS1(i)
               Istim = IstimEF_3dS1(i)
            elseif ((timeMS.GT.delayS2+0.0D0).AND.(timeMS.LT.delayS2+duration_signal_S2+0.0D0)) then
                !Istim = IstimEF_3d(i)+IstimEF_3dS2(i)
                !Istim = IstimEF_3dS2(i)
               Istim = IstimEF_3dS2(i)
            else
               !Istim = IstimEF_3d(i)
               Istim = 0.0D0
            endif
         else
            !Istim = IstimEF_3d(i)
            !Istim = 0.0D0
            if ((timeMS.GT.0.0D0).AND.(timeMS.LT.0.0D0+duration_signal_S1)) then
               !Istim = IstimEF_3d(i)!+IstimEF_3dS1(i)
               Istim = IstimEF_3dS1(i)
            else
               !Istim = IstimEF_3d(i)
               Istim = 0.0D0
            endif
         endif
         
!spatial term
        bordo = 0.0D0
        do j=1,4
           f1=face_of_cell_3d(j,i)

           xnormale=normalfaceofcells_3d(1,j,i)
           ynormale=normalfaceofcells_3d(2,j,i)
           znormale=normalfaceofcells_3d(3,j,i)
   
           xgrad = gradface_3d(1,f1) 
           ygrad = gradface_3d(2,f1) 
           zgrad = gradface_3d(3,f1) 
#ifdef BIDOMAIN
           xgradext = gradextface_3d(1,f1) 
           ygradext = gradextface_3d(2,f1) 
           zgradext = gradextface_3d(3,f1) 
#endif

           mint11 = Mintfaces_3d(1,1,f1)
           mint12 = Mintfaces_3d(1,2,f1)
           mint13 = Mintfaces_3d(1,3,f1)
           mint21 = Mintfaces_3d(2,1,f1)
           mint22 = Mintfaces_3d(2,2,f1)
           mint23 = Mintfaces_3d(2,3,f1)
           mint31 = Mintfaces_3d(3,1,f1)
           mint32 = Mintfaces_3d(3,2,f1)
           mint33 = Mintfaces_3d(3,3,f1)
#ifdef CARTO
!              bordo = bordo + &
!                   CARTO_Dcell3d(i)*(  xnormale*(mint11*xgrad + mint12*ygrad + mint13*zgrad) &
!                    + ynormale*(mint21*xgrad + mint22*ygrad + mint23*zgrad) &
! #ifdef ELEGEO0                
!                    + znormale*(mint31*xgrad + mint32*ygrad + mint33*zgrad) )*sur0_3d(f1)
! #else
!               + znormale*(mint31*xgrad + mint32*ygrad + mint33*zgrad) )*sur_3d(f1)
! #endif
             bordo = bordo + &
                  CARTO_Dface3d(f1)*(  xnormale*xgrad + ynormale*ygrad &
#ifdef ELEGEO0                
                  + znormale*zgrad )*sur0_3d(f1)
#else
                  + znormale*zgrad )*sur_3d(f1)
#endif
              
#else
           if ((cell_to_chamb_3d(c1).lt.5).and.(cell_to_chamb_3d(c2).lt.5)) then !MOD_BC ((cell_to_chamb_3d(c1).ne.9).and.(cell_to_chamb_3d(c2).ne.9)) 
              bordo = bordo + &
                   (  xnormale*(mint11*xgrad + mint12*ygrad + mint13*zgrad) &
                   + ynormale*(mint21*xgrad + mint22*ygrad + mint23*zgrad) &
#ifdef ELEGEO0                
                   + znormale*(mint31*xgrad + mint32*ygrad + mint33*zgrad) )*sur0_3d(f1)
#else
              + znormale*(mint31*xgrad + mint32*ygrad + mint33*zgrad) )*sur_3d(f1)
#endif           
           !           bordo=bordo+dot(Mintfiber*gradiente,normalface)*sur(f1);      
           
#ifdef BIDOMAIN
              bordo = bordo + &
                   (  xnormale*(mint11*xgradext + mint12*ygradext + mint13*zgradext) &
                   + ynormale*(mint21*xgradext + mint22*ygradext + mint23*zgradext) &
#ifdef ELEGEO0
                   + znormale*(mint31*xgradext + mint32*ygradext + mint33*zgradext) )*sur0_3d(f1)
#else
              + znormale*(mint31*xgradext + mint32*ygradext + mint33*zgradext) )*sur_3d(f1)
#endif
#endif
           endif!MOD_BC
#endif
        enddo
        
#ifdef ELEGEO0
         rhspot = bordo/vol0_3d(i)
#else
         rhspot = bordo/vol_3d(i)
#endif
        rhspot = rhspot/((1000.D0*LSTAR)**2) !devo dimens lo spazio in mm  

        if ((chamb.EQ.1).OR.(chamb.EQ.3)) then !Ventricle, multiple cell models available
           ! if (XEF_3d(22,i).eq.2) then
           !    cycle
           ! end if
#ifdef TP06
! !all the rest        
        Volt = XEF_3d(1,i);
        Cai  = XEF_3d(2,i);
        CaSR = XEF_3d(3,i);
        CaSS = XEF_3d(4,i);
        Nai  = XEF_3d(5,i);
        Ki   = XEF_3d(6,i);
        INa_m = XEF_3d(7,i);
        INa_h = XEF_3d(8,i);
        INa_j = XEF_3d(9,i);
        IKr_xr1 = XEF_3d(10,i);
        IKr_xr2 = XEF_3d(11,i);
        IKs_xs  = XEF_3d(12,i);
        Ito_r   = XEF_3d(13,i);
        Ito_s   = XEF_3d(14,i);
        ICaL_d  = XEF_3d(15,i);
        ICaL_f  = XEF_3d(16,i);
        ICaL_f2 = XEF_3d(17,i);
        ICaL_fCaSS = XEF_3d(18,i);
        RR         = XEF_3d(19,i);

    !%IKr
        IKr_xr1_inf = 1. / (1. + exp((-26. - Volt) / 7.));
        IKr_xr1_a = 450. / (1. + exp((-45. - Volt) / 10.));
        IKr_xr1_b = 6. / (1. + exp((Volt - (-30.)) / 11.5));
        IKr_xr1_t = IKr_xr1_a * IKr_xr1_b;
        IKr_xr2_inf = 1. / (1. + exp((Volt - (-88.)) / 24.));
        IKr_xr2_a = 3. / (1. + exp((-60. - Volt) / 20.));
        IKr_xr2_b = 1.12 / (1. + exp((Volt - 60.) / 20.));
        IKr_xr2_t = IKr_xr2_a*IKr_xr2_b;
        !%     dIKr_xr1=(IKr_xr1_inf-IKr_xr1)/IKr_xr1_t; %gate1
        IKr_xr1P1 = IKr_xr1_inf + (IKr_xr1 - IKr_xr1_inf)*exp(-dtMS/IKr_xr1_t);
        !%     dIKr_xr2=(IKr_xr2_inf-IKr_xr2)/IKr_xr2_t; %gate2
        IKr_xr2P1 = IKr_xr2_inf + (IKr_xr2 - IKr_xr2_inf)*exp(-dtMS/IKr_xr2_t);
        Ek = RTONF*(log((Ko / Ki)));
        sxr1 = IKr_xr1;
        sxr2 = IKr_xr2;
        IKr = Gkr*sqrt(Ko / 5.4)*sxr1*sxr2*(Volt - Ek);

    !%IKs
        IKs_xs_inf = 1. / (1. + exp((-5. - Volt) / 14.));
        IKs_xs_a = (1400. / (sqrt(1. + exp((5. - Volt) / 6))));
        IKs_xs_b = (1. / (1. + exp((Volt - 35.) / 15.)));
        IKs_xs_t = IKs_xs_a*IKs_xs_b + 80;
        !%     dIKs_xs = (IKs_xs_inf-IKs_xs)/IKs_xs_t; %gate3
        IKs_xsP1 = IKs_xs_inf + (IKs_xs - IKs_xs_inf)*exp(-dtMS/IKs_xs_t);
        Eks = RTONF*(log((Ko + pKNa*Nao) / (Ki + pKNa*Nai)));
        IKs = Gks*IKs_xs*IKs_xs*(Volt - Eks);
    
    !%IK1
        IK1_a = 0.1 / (1. + exp(0.06*(Volt - Ek - 200)));
        IK1_b = (3.*exp(0.0002*(Volt - Ek + 100))+exp(0.1*(Volt - Ek - 10))) / (1. + exp(-0.5*(Volt - Ek)));
        rec_iK1 = IK1_a / (IK1_a + IK1_b);
        IK1 = GK1*rec_iK1*(Volt - Ek);

    !%Ito
        Ito_r_inf = 1. / (1. + exp((20 - Volt) / 6.)); !Epi
        Ito_s_inf = 1. / (1. + exp((Volt + 20) / 5.));
        Ito_r_t = 9.5*exp(-(Volt + 40.)*(Volt + 40.) / 1800.) + 0.8;
        Ito_s_t = 85.*exp(-(Volt + 45.)*(Volt + 45.) / 320.) + 5. / (1. + exp((Volt - 20.) / 5.)) + 3.;
!%ENDO         Ito_r_inf = 1. / (1. + exp((20 - Volt) / 6.));
!%         Ito_s_inf = 1. / (1. + exp((Volt + 28) / 5.));
!%         Ito_r_t = 9.5*exp(-(Volt + 40.)*(Volt + 40.) / 1800.) + 0.8;
!%         Ito_s_t = 1000.*exp(-(Volt + 67)*(Volt + 67) / 1000.) + 8.;
!%else     Ito_r_inf = 1. / (1. + exp((20 - Volt) / 6.));
!%         Ito_s_inf = 1. / (1. + exp((Volt + 20) / 5.));
!%         Ito_r_t = 9.5*exp(-(Volt + 40.)*(Volt + 40.) / 1800.) + 0.8;
!%         Ito_s_t = 85.*exp(-(Volt + 45.)*(Volt + 45.) / 320.) + 5. / (1. + exp((Volt - 20.) / 5.)) + 3.;
!%     end
!%     dIto_r = (Ito_r_inf-Ito_r)/Ito_r_t; %gate7
        Ito_rP1 = Ito_r_inf + (Ito_r-Ito_r_inf)*exp(-dtMS/Ito_r_t);
!%     dIto_s = (Ito_s_inf-Ito_s)/Ito_s_t; %gate8
        Ito_sP1 = Ito_s_inf + (Ito_s-Ito_s_inf)*exp(-dtMS/Ito_s_t);
        Ito = Gto*Ito_r*Ito_s*(Volt - Ek);

    !%INa
        Ena = RTONF*(log((Nao / Nai)));
        INa_m_a = 1. / (1. + exp((-60. - Volt) / 5.));
        INa_m_b = 0.1 / (1. + exp((Volt + 35.) / 5.)) + 0.10 / (1. + exp((Volt - 50.) / 200.));
        INa_m_t = INa_m_a*INa_m_b;
        INa_m_inf = 1. / ((1. + exp((-56.86 - Volt) / 9.03))*(1. + exp((-56.86 - Volt) / 9.03)));
        if (Volt.GE.-40.d0) then
           INa_h_a = 0.;
           INa_h_b = (0.77 / (0.13*(1. + exp(-(Volt + 10.66) / 11.1))));
           INa_h_t = 1.0 / (INa_h_a + INa_h_b);
        else
           INa_h_a = (0.057*exp(-(Volt + 80.) / 6.8));
           INa_h_b = (2.7*exp(0.079*Volt) + (3.1e5)*exp(0.3485*Volt));
           INa_h_t = 1.0 / (INa_h_a + INa_h_b);
        endif
        INa_h_inf = 1. / ((1. + exp((Volt + 71.55) / 7.43))*(1. + exp((Volt + 71.55) / 7.43)));
        if (Volt.GE.-40.d0) then
           INa_j_a = 0.;
           INa_j_b = (0.6*exp((0.057)*Volt) / (1. + exp(-0.1*(Volt + 32.))));
           INa_j_t = 1.0 / (INa_j_a + INa_j_b);
        else
           INa_j_a = (((-2.5428e4)*exp(0.2444*Volt)-(6.948e-6)*exp(-0.04391*Volt))*(Volt + 37.78)/(1. + exp(0.311*(Volt + 79.23))));
           INa_j_b = (0.02424*exp(-0.01052*Volt) / (1. + exp(-0.1378*(Volt + 40.14))));
           INa_j_t = 1.0 / (INa_j_a + INa_j_b);
        endif
        INa_j_inf = INa_h_inf;
!%     dINa_m = (INa_m_inf-INa_m)/INa_m_t; %gate4
        INa_mP1 = INa_m_inf + (INa_m-INa_m_inf)*exp(-dtMS/INa_m_t);  
!%    dINa_h = (INa_h_inf-INa_h)/INa_h_t; %gate5
        INa_hP1 = INa_h_inf + (INa_h-INa_h_inf)*exp(-dtMS/INa_h_t);
!%     dINa_j = (INa_j_inf-INa_j)/INa_j_t; %gate6
        INa_jP1 = INa_j_inf + (INa_j-INa_j_inf)*exp(-dtMS/INa_j_t);    
        INa = GNa*(INa_m**3)*INa_h*INa_j*(Volt - Ena);
    
    !%IbNa
        IbNa = GbNa*(Volt - Ena);
    
    !%INaK
        rec_iNaK = (1. / (1. + 0.1245*exp(-0.1*Volt*F / (R*T)) + 0.0353*exp(-Volt*F / (R*T))));
        INaK = knak*(Ko / (Ko + KmK))*(Nai / (Nai + KmNa))*rec_iNaK;

    !%ICaL
        ICaL_d_inf = 1. / (1. + exp((-8 - Volt) / 7.5));
        ICaL_d_a = 1.4 / (1. + exp((-35 - Volt) / 13)) + 0.25;
        ICaL_d_b = 1.4 / (1. + exp((Volt + 5) / 5));
        ICaL_d_c = 1. / (1. + exp((50 - Volt) / 20));
        ICaL_d_t = ICaL_d_a*ICaL_d_b + ICaL_d_c;
    
        ICaL_f_inf = 1. / (1. + exp((Volt + 20) / 7));
        ICaL_f_a = 1102.5*exp(-(Volt + 27)*(Volt + 27) / 225);
        ICaL_f_b = 200. / (1 + exp((13 - Volt) / 10.));
        ICaL_f_c = (180. / (1 + exp((Volt + 30) / 10))) + 20;
        ICaL_f_t = ICaL_f_a + ICaL_f_b + ICaL_f_c;
        
        ICaL_f2_inf = 0.67 / (1. + exp((Volt + 35) / 7)) + 0.33;
        ICaL_f2_a = 600 * exp(-(Volt + 25)*(Volt + 25) / 170);
        ICaL_f2_b = 31 / (1. + exp((25 - Volt) / 10));
        ICaL_f2_c = 16 / (1. + exp((Volt + 30) / 10));
        ICaL_f2_t = ICaL_f2_a + ICaL_f2_b + ICaL_f2_c;
        
        ICaL_fCaSS_inf = 0.6 / (1 + (CaSS / 0.05)*(CaSS / 0.05)) + 0.4;
        ICaL_fCaSS_t = 80. / (1 + (CaSS / 0.05)*(CaSS / 0.05)) + 2.;
        
        !%     dCaL_d = (ICaL_d_inf-ICaL_d)/ICaL_d_t; %gate9
        CaL_dP1 = ICaL_d_inf + (ICaL_d-ICaL_d_inf)*exp(-dtMS/ICaL_d_t);    
        !%     dCaL_f = (ICaL_f_inf-ICaL_f)/ICaL_f_t; %gate10
        CaL_fP1 = ICaL_f_inf + (ICaL_f-ICaL_f_inf)*exp(-dtMS/ICaL_f_t);    
        !%     dCaL_f2 = (ICaL_f2_inf-ICaL_f2)/ICaL_f2_t; %gate11
        CaL_f2P1 = ICaL_f2_inf + (ICaL_f2-ICaL_f2_inf)*exp(-dtMS/ICaL_f2_t);        
        !%    dCaL_fCaSS = (ICaL_fCaSS_inf-ICaL_fCaSS)/ICaL_fCaSS_t; %gate12
        CaL_fCaSSP1 = ICaL_fCaSS_inf + (ICaL_fCaSS-ICaL_fCaSS_inf)*exp(-dtMS/ICaL_fCaSS_t);            
        ICaL = GCaL*ICaL_d*ICaL_f*ICaL_f2*ICaL_fCaSS * 4 * (Volt - 15)*(F**2 / (R*T))*(0.25*exp(2 * (Volt - 15)*F / (R*T))*CaSS - Cao) / (exp(2 * (Volt - 15)*F / (R*T)) - 1.);

        !%IbCa
        ECa = 0.5*RTONF*(log((Cao / Cai)));
        IbCa = GbCa*(Volt - ECa);
    
        !%INaCa
        INaCa = knaca*(1. / (KmNai**3 + Nao**3))*(1. / (KmCa + Cao))*(1. / (1 + ksat*exp((n - 1)*Volt*F / (R*T))))*(exp(n*Volt*F / (R*T))*Nai**3*Cao -exp((n - 1)*Volt*F / (R*T))*Nao**3*Cai*2.5);
    
        !%IpCa
        IpCa = GpCa*Cai / (KpCa + Cai);
    
        !%IpK
        rec_ipK = 1. / (1. + exp((25 - Volt) / 5.98));
        IpK = GpK*rec_ipK*(Volt - Ek);
    
        !%update the membrane voltage
        dv = -(IKr+IKs+IK1+Ito+INa+IbNa+INaK+ICaL+IbCa+INaCa+IpCa+IpK+Istim)
    
        !%Ca transient and intracellular concentrations
        inverseVcF2 = 1 / (2 * Vc*F);
        inverseVcF = 1. / (Vc*F);
        inversevssF2 = 1 / (2 * Vss*F);
    
        kCaSR = maxsr - ((maxsr - minsr) / (1 + (EC / CaSR)*(EC / CaSR)));
        k1 = k1_ / kCaSR;
        k2 = k2_ * kCaSR;
        !% dRR = k4 * (1 - RR) - k2*CaSS*RR;
        alphaRR = k4+k2*CaSS;
        RRP1 = (RR-k4/alphaRR)*exp(-alphaRR*dtMS)+k4/alphaRR;    
        OO = k1*CaSS*CaSS*RR / (k3 + k1*CaSS*CaSS);


        Irel = Vrel*OO*(CaSR - CaSS);
        Ileak = Vleak*(CaSR - Cai);
        Iup = Vmaxup / (1. + ((KupEF**2) / (Cai**2)));
        Ixfer = Vxfer*(CaSS - Cai);


        CaCSQN = Bufsr*CaSR / (CaSR + Kbufsr);
        dCaSR = (Iup - Irel - Ileak);
        bjsr = Bufsr - CaCSQN - dCaSR - CaSR + Kbufsr;
        cjsr = Kbufsr*(CaCSQN + dCaSR + CaSR);
        dCaSR = (sqrt(bjsr*bjsr + 4 * cjsr) - bjsr) / 2-CaSR;
        
        
        CaSSBuf = Bufss*CaSS / (CaSS + Kbufss);
        dCaSS = (-Ixfer*(Vc / Vss) + Irel*(Vsr / Vss) + (-ICaL*inversevssF2*CAPACITANCE));
        bcss = Bufss - CaSSBuf - dCaSS - CaSS + Kbufss;
        ccss = Kbufss*(CaSSBuf + dCaSS + CaSS);
        dCaSS = (sqrt(bcss*bcss + 4 * ccss) - bcss) / 2-CaSS;
        
        CaBuf = Bufc*Cai / (Cai + Kbufc);
        dCai = ((-(IbCa + IpCa - 2 * INaCa)*inverseVcF2*CAPACITANCE) - (Iup - Ileak)*(Vsr / Vc) + Ixfer);
        bc = Bufc - CaBuf - dCai - Cai + Kbufc;
        cc = Kbufc*(CaBuf + dCai + Cai);
        dCai = (sqrt(bc*bc + 4 * cc) - bc) / 2-Cai;
        
        dNai = -(INa + IbNa + 3 * INaK + 3 * INaCa)*inverseVcF*CAPACITANCE;
        
        dKi = -(Istim + IK1 + Ito + IKr + IKs - 2 * INaK + IpK)*inverseVcF*CAPACITANCE; 

!AGGIORNA        
        XEF_3d(1,i) = XEF_3d(1,i) + dtMS*(dv+rhspot);
        XEF_3d(2,i) = XEF_3d(2,i) + dtMS*dCai;
        XEF_3d(3,i) = XEF_3d(3,i) + dtMS*dCaSR;
        XEF_3d(4,i) = XEF_3d(4,i) + dtMS*dCaSS;
        XEF_3d(5,i) = XEF_3d(5,i) + dtMS*dNai;
        XEF_3d(6,i) = XEF_3d(6,i) + dtMS*dKi;
        XEF_3d(7,i) = INa_mP1; !%gate4
        XEF_3d(8,i) = INa_hP1; !%gate5
        XEF_3d(9,i) = INa_jP1; !%gate6
        XEF_3d(10,i) = IKr_xr1P1; !%gate1
        XEF_3d(11,i) = IKr_xr2P1; !%gate2
        XEF_3d(12,i) = IKs_xsP1;  !%gate3
        XEF_3d(13,i) = Ito_rP1;   !%gate7
        XEF_3d(14,i) = Ito_sP1;   !%gate8
        XEF_3d(15,i) = CaL_dP1;   !%gate9
        XEF_3d(16,i) = CaL_fP1;   !%gate10
        XEF_3d(17,i) = CaL_f2P1;  !%gate11
        XEF_3d(18,i) = CaL_fCaSSP1; !%gate12
        XEF_3d(19,i) = RRP1; !%gate13

        potEFcell_3d(i) = XEF_3d(1,i)
#endif !TP06
#ifdef MINIMAL_MODEL_AAA
           !minimal model
!           Dij = 0.1171D0 !mm^2/ms
!           Dij = 0.168D0 !mm^2/ms !Parameter Fig. 3c
           !CONSTANTSmv1 = 0; !epi
           !CONSTANTSmv2 = 1; !endo    
           !CONSTANTSmv3 = 0; !m 
           ! if (myotag(i).EQ.1) then
           !    !CONSTANTSmv2 = 1; !endo
           !    CONSTANTSmv1 = 1; !endo (epi in bueno-orovio)
           ! elseif (myotag(i).EQ.2) then
           !    CONSTANTSmv3 = 1; !m
           ! else
           !    !CONSTANTSmv1 = 1; !epi
           !    CONSTANTSmv2 = 1; !endo
        ! endif
! #ifdef CARTO
!         epi_endo_mur = 1
! #else      
!         epi_endo_mur = myotag(i)
        ! #endif
        epi_endo_mur = 1 
        CONST4 = -83;
        CONST5 = 2.7;
        ! if (epi_endo_mur.EQ.2) then !Mur
      !    CONST6 = 0.3;
      !    CONST7 = 0.13;
      !    CONST8 = 1.45;
      !    CONST9 = 2.7342;
      !    CONST10 = 2.0994;
      !    CONST11 = 0.9087;
           
      !    CONST12 =4.00;
      !    CONST13 =1.45;
      !    CONST14 =0.1;
      !    CONST15 =410.0;
      !    CONST16 =0.005;
      !    CONST17 =0.078;
      !    CONST18 =1.61;
      !    CONST19 =91.00;
      !    CONST20 =0.80;
      !    CONST21 =2.1;
      !    CONST22 =0.6;
      !    CONST23 =3.3849;
      !    CONST24 =0.01;
      !    CONST25 =0.5;
      !    CONST26 =70.00;
      !    CONST27 =8.00;
      !    CONST28 =200.00;
      !    CONST29 =0.016;
      !    CONST30 =280.;
      !    CONST31 =80.00;
      !    CONST32 =7.00;
      !   elseif (epi_endo_mur.EQ.1) then !endo
      !      ! CONST12 = 2.0D0;
      !      ! CONST13 = 10.0D0;
      !      ! CONST14 = 0.024D0!0.2D0!0.02D0;
      !      ! CONST15 = 470.0D0;
      !      ! CONST16 = 0.006D0;
      !      ! CONST17 = 0.1040D0;
      !      ! CONST18 = 1.56D0;
      !      ! CONST19 = 40.0D0;
      !      ! CONST20 = 1.2D0;
      !      ! CONST21 = 2.0D0;
      !      ! CONST22 = 0.65D0;
      !      ! CONST23 = 2.9013D0;
      !      ! CONST24 = 0.0273D0;
      !      ! CONST25 = 0.78D0;
      !      ! CONST26 = 6.0D0;
      !      ! CONST27 = 140.0D0;
      !      ! CONST28 = 200.0D0;
      !      ! CONST29 = 0.016D0;
      !      ! CONST30 = 280.0D0;
      !      ! CONST31 = 75.0D0;
      !      ! CONST32 = 6.0D00;

           
      !    ! From supplmenetary material, Set 4V-3 (Fenton Science Adv.)
      !    CONST6 = 0.218655D0;
      !    CONST7 = 0.15D0;
      !    CONST8 = 7.188259D0;
      !    CONST9 = 2.031406D0;
      !    CONST10 = 2.545600D0;
      !    CONST11 = 0.856429D0;
         
      !    CONST12 = 13.090045D0;
      !    CONST13 = 573.569884D0;
      !    CONST14 = 0.044570D0
      !    CONST15 = 284.169705D0;
      !    CONST16 = 0.007950D0;
      !    CONST17 = 0.168277D0;
      !    CONST18 = 1.451815D0;
      !    CONST19 = 28.726329D0;
      !    CONST20 = 0.281850D0;
      !    CONST21 = 1.937860D0;
      !    CONST22 = 0.682362D0;
      !    CONST23 = 2.347223D0;
      !    CONST24 = 0.049502D0;
      !    CONST25 = 0.651285D0;
      !    CONST26 = 19.751628D0;
      !    CONST27 = 101.167275D0;
      !    CONST28 = 57.162563D0;
      !    CONST29 = 0.028506D0;
      !    CONST30 = 213.296779D0;
      !    CONST31 = 19.147136D0;
      !    CONST32 = 9.999739D0;
         
      !  ! ! From supplmenetary material, Set 4V-4 (Fenton Science Adv.)
      !    ! CONST6 = 0.33020D0;
      !      ! CONST7 = 0.15D0;
      !      ! CONST8 = 14.737856D0;
      !      ! CONST9 = 2.016938D0;
      !      ! CONST10 = 1.791970D0;
      !      ! CONST11 = 0.857440D0;
           
      !      ! CONST12 = 7.005449D0;
      !      ! CONST13 = 107.183777D0;
      !      ! CONST14 = 0.153932D0
      !      ! CONST15 = 365.098651D0;
      !      ! CONST16 = 0.004032D0;
      !      ! CONST17 = 0.129816D0;
      !      ! CONST18 = 1.45D0;
      !      ! CONST19 = 46.850196D0;
      !      ! CONST20 = 0.404662D0;
      !      ! CONST21 = 1.80D0;
      !      ! CONST22 = 0.70D0;
      !      ! CONST23 = 2.119200D0;
      !      ! CONST24 = 0.011585D0;
      !      ! CONST25 = 0.40D0;
      !      ! CONST26 = 87.085277D0;
      !      ! CONST27 = 105.663454D0;
      !      ! CONST28 = 106.727040D0;
      !      ! CONST29 = 0.037092D0;
      !      ! CONST30 = 105.848052D0;
      !      ! CONST31 = 18.754498D0;
      !      ! CONST32 = 8.195328D0;
        
      ! else !epi
      !      CONST6 = 0.3;
      !      CONST7 = 0.13;
      !      CONST8 = 1.45;
      !      CONST9 = 2.7342;
      !      CONST10 = 2.0994;
      !      CONST11 = 0.9087;
           
      !      CONST12 =16.00;
      !      CONST13 =1150.00;
      !      CONST14 =0.006;
      !      CONST15 =400.00;
      !      CONST16 =0.006;
      !      CONST17 =0.11;
      !      CONST18 =1.55;
      !      CONST19 =30.02;
      !      CONST20 =0.996;
      !      CONST21 =2.046;
      !      CONST22 =0.65;
      !      CONST23 =1.8875;
      !      CONST24 =0.07;
      !      CONST25 =0.94;
      !      CONST26 =60.0;
      !      CONST27 =15.00;
      !      CONST28 =65.00;
      !      CONST29 =0.03;
      !      CONST30 =200.00;
      !      CONST31 =60.00;
      !      CONST32 =6.00;
      ! endif
      
      CONST6 = 0.3d0;
      CONST7 = 0.13d0;
      CONST8 = 1.4506d0;
      CONST9 = 2.7342d0;
      CONST10 = 2.0994d0;
      CONST11 = 0.9087d0;
      
      CONST12 =2.0d0;
      CONST13 =40.0d0;
      CONST14 =0.20d0;
      CONST15 =470.0d0;
      CONST16 =0.006d0;
      CONST17 =0.10d0;
      CONST18 =1.56d0;
      CONST19 =40.0d0;
      CONST20 =1.2d0;
      CONST21 =2.0d0;
      CONST22 =0.65d0;
      CONST23 =2.9013d0;
      CONST24 =0.0273d0;
      CONST25 =0.78d0;
      CONST26 =40.0d0;
      CONST27 =115.0d0;
      CONST28 =20.0d0;
      CONST29 =0.00615d0;
      CONST30 =8.0d0;
      CONST31 =55.0d0;
      
      !Added for tau_w+
      CONST32 =6.0d0;
      CONST33 =0.0005d0;
      CONST34 =175.0d0;
      CONST35 =230.0d0;
      
      STA1 = XEF_3d(1,i); !nondimensional here
      STA2 = XEF_3d(2,i);
      STA3 = XEF_3d(3,i);
      STA4 = XEF_3d(4,i);
      
      if (STA1.LT.CONST7) then
         ALG5=0.0;
      else
         ALG5=1.0;
      end if
      
   
      ALG8 =  (1.00000 - ALG5)*CONST9+ ALG5*CONST12;
      STA4np1 =(1.00000+ tanh( CONST10*(STA1 - CONST11)))/2.00000 + (STA4-(1.00000+ tanh( CONST10*(STA1 - CONST11)))/2.00000)*exp(-dtMS/ALG8);
      if (STA1.LT.CONST6) then
         ALG3=0.0;
      else
         ALG3=1.0;
      end if
    
      if (STA1.LT.CONST14) then
         ALG7=1.0;
      else
         ALG7=0.0;
      end if
      if (STA1.LT.CONST14) then
         ALG4=0.0;
      else
         ALG4=1.0;
      end if
      ALG10 =  ALG4*CONST13+ (1.00000 - ALG4)*CONST31;
      RAT2 = ( (1.00000 - ALG3)*(ALG7 - STA2))/ALG10 - ( ALG3*STA2)/CONST8;
      app1=CONST8*ALG7*(1.0-ALG3);
      app2=(CONST8-CONST8*ALG3+ALG10*ALG3);
      app3=ALG10*CONST8;
      app4=app1/app2;
      app5=app2/app3;
      STA2np1=app4 + (STA2 -app4) *exp(-app5*dtMS);
      
      if (STA1.LT.CONST16)  then
         ALG6=0.00;
      else
         ALG6=1.00;
      endif
      ALG11 =  (1.00000 - ALG6)*(1.00000 - ( STA1*1.00000)/CONST24)+ ALG6*CONST25;
      ALG13 = CONST26+( (CONST27 - CONST26)*(1.00000+ tanh( CONST28*(STA1 - CONST29))))/2.00000;!tau_w-
      ALG17 = CONST33+( (CONST34 - CONST33)*(1.00000+ tanh( CONST30*(STA1 - CONST33))))/2.00000;!tau_w+
      app1=ALG17*ALG11*(1.0-ALG6);
      app2=(ALG17-ALG17*ALG6+ALG13*ALG6);
      app3=ALG13*ALG17;
      app4=app1/app2;
      app5=app2/app3;
      STA3np1=app4 + (STA3 -app4) *exp(-app5*dtMS);
      ALG9 = (  - ALG3*STA2*(STA1 - CONST6)*(CONST18 - STA1))/CONST17;
      ALG12 =  (1.00000 - ALG6)*CONST15+ ALG6*CONST32;
      ALG14 = CONST19+( (CONST20 - CONST19)*(1.00000+ tanh( CONST21*(STA1 - CONST22))))/2.00000;
      ALG15 = ( STA1*(1.00000 - ALG5))/ALG12+ALG5/ALG14;
      ALG16 = (  - ALG5*STA3*STA4)/CONST23;
      
           
      RAT1 =  - (ALG9+ALG15+ALG16) -Istim; 
      
      CONST4 = -83;
      CONST5 = 2.7;

      remod=1.0D0!5.0D0*CARTO_Dcell3d(i)
      
      XEF_3d(1,i) = STA1 + dtMS*(remod*RAT1+rhspot/(CONST5-CONST4)); 
      XEF_3d(2,i) = STA2np1
      XEF_3d(3,i) = STA3np1
      XEF_3d(4,i) = STA4np1
      
      potEFcell_3d(i) = CONST4+ (CONST5-CONST4)*XEF_3d(1,i)

#endif
#ifdef MITCHELL_SCHAEFFER
      remod=1.0D0/0.09D0*CARTO_Dcell3d(i)!10.0D0/1.65D0*CARTO_Dcell3d(i)
      STA1 = XEF_3d(1,i); !nondimensional here
      STA2 = XEF_3d(2,i);
      if (STA1.lt.u_gate)then
         STA2np1=STA2+dtMS*(1.0D0-STA2)*invtau_open
      else
         STA2np1=STA2-dtMS*STA2*invtau_close
      endif

      CONST4 = -83.0D0;
      CONST5 = 28.58D0;!2.7;

      XEF_3d(1,i) =STA1 +dtMS*(remod*(STA2np1*STA1*STA1*(1-STA1)*invtau_in-STA1*invtau_out-Istim)+rhspot/(CONST5-CONST4))
      XEF_3d(2,i) =STA2np1

      potEFcell_3d(i) = CONST4+ (CONST5-CONST4)*XEF_3d(1,i)
#endif
   elseif ((chamb.EQ.2).OR.(chamb.EQ.4)) then !Courtemanche
#ifdef COURTEMANCHE
           STA1 = XEF_3d(1,i); !Volt
           STA2 = XEF_3d(2,i);
           STA3 = XEF_3d(3,i);
           STA4 = XEF_3d(4,i);
           STA5 = XEF_3d(5,i);
           STA6 = XEF_3d(6,i);
           STA7 = XEF_3d(7,i);
           STA8 = XEF_3d(8,i);
           STA9 = XEF_3d(9,i);
           STA10 = XEF_3d(10,i);
           STA11 = XEF_3d(11,i);
           STA12 = XEF_3d(12,i);
           STA13 = XEF_3d(13,i);
           STA14 = XEF_3d(14,i);
           STA15 = XEF_3d(15,i);
           STA16 = XEF_3d(16,i);
           STA17 = XEF_3d(17,i);
           STA18 = XEF_3d(18,i);
           STA19 = XEF_3d(19,i);
           STA20 = XEF_3d(20,i);
           STA21 = XEF_3d(21,i);

           RAT1=0.0D0;RAT2=0.0D0;RAT6=0.0D0; !mi sembra non necessario
           RAT13=0.0D0;RAT17=0.0D0;RAT21=0.0D0;


           ALG13 = 1/(1.00000+STA13/0.000350000);
           STA16np1 = ALG13 + (STA16 - ALG13)*exp(-dtMS/CONST45c);
           ALG11 = 1/(1.00000+exp((STA1+10.0000)/ - 8.00000));
           if (abs(STA1+10.0000).LT.1.00000e-10) then
              ALG28 = 4.57900/(1.00000+exp((STA1+10.0000)/ - 6.24000));
           else
              ALG28 = (1.00000 - exp((STA1+10.0000)/ - 6.24000))/( 0.0350000*(STA1+10.0000)*(1.00000+exp((STA1+10.0000)/ - 6.24000)));
           end if
           STA14np1 = ALG11 + (STA14 - ALG11)*exp(-dtMS/ALG28);
           ALG12 = exp( - (STA1+28.0000)/6.90000)/(1.00000+exp( - (STA1+28.0000)/6.90000));
           ALG29 =  9.00000/( 0.0197000*exp(  (- 0.0337000**2.00000)*( (STA1+10.0000)**2.00000))+0.0200000);
           STA15np1 = ALG12 + (STA15 - ALG12)*exp(-dtMS/ALG29);
    
           if (abs(STA1 - 7.90000).LT.1.00000e-10) then
              ALG14 = ( 6.00000*0.200000)/1.30000;
           else
              ALG14 = ( 6.00000*(1.00000 - exp( - (STA1 - 7.90000)/5.00000)))/( (1.00000+ 0.300000*exp( - (STA1 - 7.90000)/5.00000))*1.00000*(STA1 - 7.90000));
           end if
           ALG30 = 1.00000 - 1/(1.00000+exp( - (STA1 - 40.0000)/17.0000));
           STA20np1 = ALG30 + (STA20 - ALG30)*exp(-dtMS/ALG14);
           
           if (abs(STA1 - 47.13).LT.1.00000e-10) then
              ALG2 = 3.20000;
           else
              ALG2 = ( 0.320000*(STA1+47.1300))/(1.00000 - exp(  - 0.100000*(STA1+47.1300)));
           endif
!    %QUESTO NON MI TORNA; LO MODIFICO
!    %ALG2 = piecewise({STA1== - 47.1300, 3.20000 }, ( 0.320000*(STA1+47.1300))/(1.00000 - exp(  - 0.100000*(STA1+47.1300))));
           ALG19 =  0.0800000*exp( - STA1/11.0000);
           ALG32 = ALG2/(ALG2+ALG19);
           ALG42 = 1.00000/(ALG2+ALG19);
           STA3np1 = ALG32 + (STA3 - ALG32)*exp(-dtMS/ALG42);
           if (STA1.LT.-40.0000) then
              ALG3 = 0.135000*exp((STA1+80.0000)/ - 6.80000);
              ALG20 = 3.56000*exp( 0.0790000*STA1)+ 310000.*exp( 0.350000*STA1);
              ALG4 = ( (  - 127140.*exp( 0.244400*STA1) -  3.47400e-05*exp(  - 0.0439100*STA1))*(STA1+37.7800))/(1.00000+exp( 0.311000*(STA1+79.2300))) ;
              ALG21 = ( 0.121200*exp(  - 0.0105200*STA1))/(1.00000+exp(  - 0.137800*(STA1+40.1400)));
           else
              ALG3 = 0.0;
              ALG20 = 1.00000/( 0.130000*(1.00000+exp((STA1+10.6600)/ - 11.1000)));
              ALG4 =  0.0;
              ALG21 = ( 0.300000*exp(  - 2.53500e-07*STA1))/(1.00000+exp(  - 0.100000*(STA1+32.0000)));
           endif

           ALG33 = ALG3/(ALG3+ALG20);
           ALG43 = 1.00000/(ALG3+ALG20);
           STA4np1 = ALG33 + (STA4 - ALG33)*exp(-dtMS/ALG43);
           ALG34 = ALG4/(ALG4+ALG21);
           ALG44 = 1.00000/(ALG4+ALG21);
           STA5np1 = ALG34 + (STA5 - ALG34)*exp(-dtMS/ALG44);
           ALG5 =  0.650000*1/(exp((STA1 -  - 10.0000)/ - 8.50000)+exp(((STA1 -  - 10.0000) - 40.0000)/ - 59.0000));
           ALG22 =  0.650000/(2.50000+exp(((STA1 -  - 10.0000)+72.0000)/17.0000));
           ALG35 = 1/(CONST14c*(ALG5+ALG22));
           ALG45 = 1/(1.00000+exp(((STA1 -  - 10.0000)+10.4700)/ - 17.5400));
           STA7np1 = ALG45 + (STA7 - ALG45)*exp(-dtMS/ALG35);
           ALG6 = 1/(18.5300+ 1.00000*exp(((STA1 -  - 10.0000)+103.700)/10.9500));
           ALG23 = 1/(35.5600+ 1.00000*exp(((STA1 -  - 10.0000) - 8.74000)/ - 7.44000));
           ALG36 = 1/(CONST14c*(ALG6+ALG23));
           ALG46 = 1/(1.00000+exp(((STA1 -  - 10.0000)+33.1000)/5.30000));
           STA8np1 = ALG46 + (STA8 - ALG46)*exp(-dtMS/ALG36);
           ALG7 =  0.650000/(exp((STA1 -  - 10.0000)/ - 8.50000)+exp(((STA1 -  - 10.0000) - 40.0000)/ - 59.0000));
           ALG24 =  0.650000/(2.50000+exp(((STA1 -  - 10.0000)+72.0000)/17.0000));
           ALG37 = 1/(CONST14c*(ALG7+ALG24));
           ALG47 = 1/(1.00000+exp(((STA1 -  - 10.0000)+20.3000)/ - 9.60000));
           STA9np1 = ALG47 + (STA9 - ALG47)*exp(-dtMS/ALG37);
           ALG8 = 1/(21.0000+ 1.00000*exp(((STA1 -  - 10.0000) - 195.000)/ - 28.0000));
           ALG25 = 1.00000/exp(((STA1 -  - 10.0000) - 168.000)/ - 16.0000);
           ALG38 = 1/(CONST14c*(ALG8+ALG25));
           ALG48 = 1/(1.00000+exp(((STA1 -  - 10.0000) - 109.450)/27.4800));
           STA10np1 = ALG48 + (STA10 - ALG48)*exp(-dtMS/ALG38);
           if (abs(STA1+14.1000).LT.1.00000e-10) then
              ALG9 = 0.00150000;
           else
              ALG9 = ( 0.000300000*(STA1+14.1000))/(1.00000 - exp((STA1+14.1000)/ - 5.00000));
           end if
           if (abs(STA1 - 3.33280).LT.1.00000e-10) then
              ALG26 = 0.000378361;
           else
              ALG26 = ( 7.38980e-05*(STA1 - 3.33280))/(exp((STA1 - 3.33280)/5.12370) - 1.00000);
           end if
           ALG39 = 1/(ALG9+ALG26);
           ALG49 = 1/(1.00000+exp((STA1+14.1000)/ - 6.50000));
           STA11np1 = ALG49 + (STA11 - ALG49)*exp(-dtMS/ALG39);

           if (abs(STA1 - 19.9000).LT.1.00000e-10) then
              ALG10 = 0.000680000;
              ALG27 = 0.000315000;
           else
              ALG10 = ( 4.00000e-05*(STA1 - 19.9000))/(1.00000 - exp((STA1 - 19.9000)/ - 17.0000));
              ALG27 = ( 3.50000e-05*(STA1 - 19.9000))/(exp((STA1 - 19.9000)/9.00000) - 1.00000);
           endif
           ALG40 =  0.500000/(ALG10+ALG27);
           ALG50 = 1/sqrt(1.00000+exp((STA1 - 19.9000)/ - 12.7000));
           STA12np1 = ALG50 + (STA12 - ALG50)*exp(-dtMS/ALG40);

           ALG41 =  (( CONST1c*CONST2c)/CONST3c)*log(CONST13c/STA6);
           ALG51 = ( CONST4c*CONST12c*(STA1 - ALG41))/(1.00000+exp( 0.0700000*(STA1+80.0000)));
           ALG52 =  CONST4c*CONST15c*(STA7**3)*STA8*(STA1 - ALG41);
           ALG53 = 0.00500000+0.0500000/(1.00000+exp((STA1 - 15.0000)/ - 13.0000));
           ALG54 =  CONST4c*ALG53*(STA9**3)*STA10*(STA1 - ALG41);
           ALG55 = ( CONST4c*CONST16c*STA11*(STA1 - ALG41))/(1.00000+exp((STA1+15.0000)/22.4000));
           ALG56 =  CONST4c*CONST17c*(STA12**2)*(STA1 - ALG41);
           ALG58 = 1/(1.00000+ 0.124500*exp((  - 0.100000*CONST3c*STA1)/( CONST1c*CONST2c))+ 0.0365000*CONST46c*exp((  - CONST3c*STA1)/( CONST1c*CONST2c)));
           ALG59 = ( (( CONST4c*CONST21c*ALG58*1.00000)/(1.00000+(CONST19c/STA2)**(1.5)))*CONST13c)/(CONST13c+CONST20c);
           ALG61 =  CONST4c*CONST24c*(STA1 - ALG41);
           RAT6 = ( 2.00000*ALG59 - (ALG51+ALG52+ALG54+ALG55+ALG56+ALG61))/( CONST44c*CONST3c);
           ALG18 =  (( CONST1c*CONST2c)/CONST3c)*log(CONST11c/STA2);
           ALG31 =  CONST4c*CONST10c*(STA3**3)*STA4*STA5*(STA1 - ALG18);
           ALG64 = ( CONST4c*CONST26c*( exp(( CONST30c*CONST3c*STA1)/( CONST1c*CONST2c))*(STA2**3)*CONST25c -  exp(( (CONST30c - 1.00000)*CONST3c*STA1)/( CONST1c*CONST2c))*(CONST11c**3)*STA13))/( ((CONST27c**3.00000)+(CONST11c**3.00000))*(CONST28c+CONST25c)*(1.00000+ CONST29c*exp(( (CONST30c - 1.00000)*STA1*CONST3c)/( CONST1c*CONST2c))));
           ALG62 =  CONST4c*CONST22c*(STA1 - ALG18);
           RAT2 = (  - 3.00000*ALG59 - ( 3.00000*ALG64+ALG62+ALG31))/( CONST44c*CONST3c);

           ALG57 =  CONST4c*CONST18c*STA14*STA15*STA16*(STA1 - 65.0000);
           ALG65 = ( CONST4c*CONST31c*STA13)/(0.000500000+STA13);
           ALG60 =  (( CONST1c*CONST2c)/( 2.00000*CONST3c))*log(CONST25c/STA13);
           ALG63 =  CONST4c*CONST23c*(STA1 - ALG60);
           
!           RAT1 =  - (ALG31+ALG51+ALG52+ALG54+ALG55+ALG56+ALG62+ALG63+ALG59+ALG65+ALG64+ALG57)/CONST4c +Istim; 
!           RAT1 =  - (ALG31+ALG51+ALG52+ALG54+ALG55+ALG56+ALG62+ALG63+ALG59+ALG65+ALG64+ALG57)/CONST4c ;
           RAT1 =  - (ALG31+ALG51+ALG52+ALG54+ALG55+ALG56+ALG62+ALG63+ALG59+ALG65+ALG64+ALG57)/CONST4c -Istim;

           ALG66 =  CONST32c*(STA18**2.00000)*STA19*STA20*(STA17 - STA13);
           ALG68 = (STA21 - STA17)/CONST33c;
           RAT17 =  (ALG68 - ALG66)/(1.00000+( CONST39c*CONST42c)/(STA17+(CONST42c**2.00000)));
           ALG67 =  1000.00*( 1.00000e-15*CONST48c*ALG66 -  (1.00000e-15/( 2.00000*CONST3c))*( 0.500000*ALG57 -  0.200000*ALG64));
           ALG69 = 1/(1.00000+exp( - (ALG67 - 3.41750e-13)/1.36700e-15));
           STA18np1 = ALG69 + (STA18 - ALG69)*exp(-dtMS/CONST47c);
           ALG70 = 1.91000+ 2.09000/(1.00000+exp( - (ALG67 - 3.41750e-13)/1.36700e-15));
           ALG72 = 1.00000 - 1/(1.00000+exp( - (ALG67 - 6.83500e-14)/1.36700e-15));
           STA19np1 = ALG72 + (STA19 - ALG72)*exp(-dtMS/ALG70);
           ALG71 = CONST34c/(1.00000+CONST35c/STA13);
           ALG73 = ( CONST34c*STA21)/CONST36c;
           RAT21 = ALG71 - (ALG73+( ALG68*CONST48c)/CONST49c);
           ALG74 = ( 2.00000*ALG64 - (ALG65+ALG57+ALG63))/( 2.00000*CONST44c*CONST3c)+( CONST49c*(ALG73 - ALG71)+ ALG66*CONST48c)/CONST44c;
           ALG75 = 1.00000+( CONST38c*CONST41c)/((STA13+CONST41c)**2.00000)+( CONST37c*CONST40c)/((STA13+CONST40c)**2.00000);
           RAT13 = ALG74/ALG75;


!           XEF_3d(1,i) = STA1 + dtMS*(RAT1+rhspot+Istim);
           XEF_3d(1,i) = STA1 + dtMS*(RAT1+rhspot); !l'ho spostata in RAT1, come nel modello ventricolare     
           XEF_3d(2,i) = STA2 + dtMS*RAT2
           XEF_3d(3,i) = STA3np1
           XEF_3d(4,i) = STA4np1
           XEF_3d(5,i) = STA5np1
           XEF_3d(6,i) = STA6 + dtMS*RAT6
           XEF_3d(7,i) = STA7np1
           XEF_3d(8,i) = STA8np1
           XEF_3d(9,i) = STA9np1
           XEF_3d(10,i) = STA10np1
           XEF_3d(11,i) = STA11np1
           XEF_3d(12,i) = STA12np1
           XEF_3d(13,i) = STA13 + dtMS*RAT13
           XEF_3d(14,i) = STA14np1
           XEF_3d(15,i) = STA15np1
           XEF_3d(16,i) = STA16np1
           XEF_3d(17,i) = STA17 + dtMS*RAT17
           XEF_3d(18,i) = STA18np1
           XEF_3d(19,i) = STA19np1
           XEF_3d(20,i) = STA20np1
           XEF_3d(21,i) = STA21 + dtMS*RAT21
           
           potEFcell_3d(i) = XEF_3d(1,i)
           !          potEFcell_3d(i) = 0.0D0
#endif!Courtemanche
#ifdef MINIMAL_MODEL !ATRIA

           !minimal model atria ( Parameters from Anatomical and spiral 
           !                      wave reentry in a simplified model for atrial electrophysiology) Richter et. al 2017
           CONST6 = 0.3D0;
           CONST7 = 0.18171D0;
           CONST8 = 1.7026D0;
           CONST9 = 9.876D0;
           CONST10 = 2.2268D0;
           CONST11 = 0.81568D0;
           
           CONST12 =4.2036D0;
           CONST13 =1150.0D0;
           CONST14 =0.1007D0;
           CONST15 =250.03D0;
           CONST16 =0.015473D0;
           CONST17 =0.083536D0;
           CONST18 =1.0089D0;
           CONST19 =73.675D0;
           CONST20 =6.5537D0;
           CONST21 =2.9748D0;
           CONST22 =0.592093D0;
           CONST23 =10.699D0;
           CONST24 =0.2233D0;
           CONST25 =0.902D0;
           CONST26 =79.963D0;
           CONST27 =28.136D0;
           CONST28 =60.219D0;
           CONST29 =0.009991D0;
           CONST30 =213.55D0;
           CONST31 =16.3D0;
           CONST32 =16.632D0;
      
      

           STA1 = XEF_3d(1,i); !nondimensional here
           STA2 = XEF_3d(2,i);
           STA3 = XEF_3d(3,i);
           STA4 = XEF_3d(4,i);
           
           if (STA1.LT.CONST7) then
              ALG5=0.0;
           else
              ALG5=1.0;
           end if
           
           
           ALG8 =  (1.00000 - ALG5)*CONST9+ ALG5*CONST12;
           STA4np1 =(1.00000+ tanh( CONST10*(STA1 - CONST11)))/2.00000 + (STA4-(1.00000+ tanh( CONST10*(STA1 - CONST11)))/2.00000)*exp(-dtMS/ALG8);
           if (STA1.LT.CONST6) then
              ALG3=0.0;
           else
              ALG3=1.0;
           end if
           
           if (STA1.LT.CONST14) then
              ALG7=1.0;
           else
              ALG7=0.0;
           end if
           if (STA1.LT.CONST14) then
              ALG4=0.0;
           else
              ALG4=1.0;
           end if
           ALG10 =  ALG4*CONST13+ (1.00000 - ALG4)*CONST31;
           RAT2 = ( (1.00000 - ALG3)*(ALG7 - STA2))/ALG10 - ( ALG3*STA2)/CONST8;
           app1=CONST8*ALG7*(1.0-ALG3);
           app2=(CONST8-CONST8*ALG3+ALG10*ALG3);
           app3=ALG10*CONST8;
           app4=app1/app2;
           app5=app2/app3;
           STA2np1=app4 + (STA2 -app4) *exp(-app5*dtMS);
           
           if (STA1.LT.CONST16)  then
              ALG6=0.00;
           else
              ALG6=1.00;
           endif
           ALG11 =  (1.00000 - ALG6)*(1.00000 - ( STA1*1.00000)/CONST24)+ ALG6*CONST25;
           ALG13 = CONST26+( (CONST27 - CONST26)*(1.00000+ tanh( CONST28*(STA1 - CONST29))))/2.00000;
           app1=CONST30*ALG11*(1.0-ALG6);
           app2=(CONST30-CONST30*ALG6+ALG13*ALG6);
           app3=ALG13*CONST30;
           app4=app1/app2;
           app5=app2/app3;
           STA3np1=app4 + (STA3 -app4) *exp(-app5*dtMS);
           ALG9 = (  - ALG3*STA2*(STA1 - CONST6)*(CONST18 - STA1))/CONST17;
           ALG12 =  (1.00000 - ALG6)*CONST15+ ALG6*CONST32;
           ALG14 = CONST19+( (CONST20 - CONST19)*(1.00000+ tanh( CONST21*(STA1 - CONST22))))/2.00000;
           ALG15 = ( STA1*(1.00000 - ALG5))/ALG12+ALG5/ALG14;
           ALG16 = (  - ALG5*STA3*STA4)/CONST23;
           
           
           RAT1 =  - (ALG9+ALG15+ALG16) -Istim; 
           
           CONST4 = -83;
           CONST5 = 2.7;
           
           XEF_3d(1,i) = STA1 + dtMS*(RAT1+rhspot/(CONST5-CONST4)); 
           XEF_3d(2,i) = STA2np1
           XEF_3d(3,i) = STA3np1
           XEF_3d(4,i) = STA4np1
           
           potEFcell_3d(i) = CONST4+ (CONST5-CONST4)*XEF_3d(1,i)
           
#endif
        endif

     endif

     !Three-dimensional APD data
     pot_apd = potEFcell_3d(i)
     pot_old = old_1(i)
     if ((pot_apd.gt.treshold_apd.and.pot_old.lt.treshold_apd).and.&
          f_apd(i).eq.0) then
        t_apd_3d(1,i) = timeMS! (dtms/(pot_apd-pot_old))*(treshold_apd-pot_apd+((pot_apd-pot_old)/dtms)*timeMS)
        f_apd(i) = 1
     elseif ((pot_apd.lt.treshold_apd.and.pot_old.gt.treshold_apd).and.&
          f_apd(i).eq.1) then
        t_apd_3d(2,i) = timeMS!(dtms/(pot_apd-pot_old))*(treshold_apd-pot_apd+((pot_apd-pot_old)/dtms)*timeMS)
        f_apd(i) = 2
     elseif ((pot_apd.gt.treshold_apd.and.pot_old.lt.treshold_apd).and.&
          f_apd(i).eq.2) then
        t_apd_3d(3,i) = timeMS!(dtms/(pot_apd-pot_old))*(treshold_apd-pot_apd+((pot_apd-pot_old)/dtms)*timeMS)
        f_apd(i) = 3
     elseif ((pot_apd.lt.treshold_apd.and.pot_old.gt.treshold_apd).and.&
          f_apd(i).eq.3) then
        t_apd_3d(4,i) = timeMS!(dtms/(pot_apd-pot_old))*(treshold_apd-pot_apd+((pot_apd-pot_old)/dtms)*timeMS)
        f_apd(i) = 0
     endif
     old_1(i) = potEFcell_3d(i)

   enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP


   
#ifdef USE_CUDA 
       !$cuf kernel do (1)
#endif
       do i=1,nctot_3d
          if (LabelStenosi(i).EQ.1) then
             v1=vert_of_cell_3d(1,i)
             v2=vert_of_cell_3d(2,i)
             v3=vert_of_cell_3d(3,i)
             v4=vert_of_cell_3d(4,i)
             
             potEFcell_3d(i)=0.25D0*(potEFnode_3d(v1)+potEFnode_3d(v2)+potEFnode_3d(v3)+potEFnode_3d(v4))

          endif
       enddo
       !@cuf istat = cudaDeviceSynchronize !JDR TMP
     !find nodal values CV
#ifdef USE_CUDA
     !$cuf kernel do (1)
#endif
     do i=vstart_3d(1),vend_3d(1)
        num      = 0.0d0
        den      = 0.0d0
#ifdef ELEGEO0
        xV1 = xyz0_3d(1,i)
        yV1 = xyz0_3d(2,i)
        zV1 = xyz0_3d(3,i)
#else
        xV1 = xyz_3d(1,i)
        yV1 = xyz_3d(2,i)
        zV1 = xyz_3d(3,i)
#endif
        do j=1,n_cell_of_vert_3d(i)
           c1=cell_of_vert_3d(j,i)
           if (cell_to_chamb_3d(c1).lt.5) then!ne.9) then !MOD_BC
              xC1 = cell_bar(1,c1)
              yC1 = cell_bar(2,c1)
              zC1 = cell_bar(3,c1)
              dvvjj = sqrt( (xV1-xC1)**2+(yV1-yC1)**2+(zV1-zC1)**2)
              num = num + t_apd_3d(1,c1)/dvvjj
              den = den + 1.D0/dvvjj
           endif !MOD_BC
        enddo
        if (den.NE.0.0D0) then
           LATnode_3d(i)=num/den
        else
           LATnode_3d(i)=0.0d0
        endif
     enddo
     !@cuf istat = cudaDeviceSynchronize !JDR TMP

       
       !Conduction velocity calculation
       !find gradient on faces
#ifdef USE_CUDA
       !$cuf kernel do (1)
#endif
       do i=fstart_3d(1),fend_3d(1)
          c1 = cell_of_face_3d(1,i)
          c2 = cell_of_face_3d(2,i)
          
          g1interp = g1interpface_3d(i)
          xvers  = versCFface_3d(1,i)
          yvers  = versCFface_3d(2,i)
          zvers  = versCFface_3d(3,i)
          distCF = distCFface_3d(i)
          if ((c1.NE.0).AND.(c2.NE.0)) then
             chamb1 = cell_to_chamb_3d(c1)
             chamb2 = cell_to_chamb_3d(c2)
             if ((chamb1.lt.5).AND.(chamb2.lt.5)) then
                !grandiente medio
                latf = g1interp*t_apd_3d(1,c1) + (1-g1interp)*t_apd_3d(1,c2)
             else
                latf = 0.0D0 !for Neumann bcs
             endif
          elseif ((c1.NE.0).AND.(c2.EQ.0)) then
             latf = t_apd_3d(1,c1)
          elseif ((c1.EQ.0).AND.(c2.NE.0)) then
             latf = t_apd_3d(1,c2)
          endif
          LATface_3d(i) = latf
       enddo
       !@cuf istat = cudaDeviceSynchronize !JDR TMP

       !Find grad(LAT) on cells
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=cstart_3d(1),cend_3d(1)
        xgradlat=0.0D0
        ygradlat=0.0D0
        zgradlat=0.0D0
        do j=1,4
           f1=face_of_cell_3d(j,i)
           xnormale=normalfaceofcells_3d(1,j,i)
           ynormale=normalfaceofcells_3d(2,j,i)
           znormale=normalfaceofcells_3d(3,j,i)

#ifdef ELEGEO0
          latf = LATface_3d(f1)*sur0_3d(f1)
#else
          latf = LATface_3d(f1)*sur_3d(f1)
#endif
          xgradlat = xgradlat + latf*xnormale
          ygradlat = ygradlat + latf*ynormale
          zgradlat = zgradlat + latf*znormale
        enddo
#ifdef ELEGEO0
        volcell = vol0_3d(i)
#else
        volcell = vol_3d(i)
#endif  
        LATx  =  xgradlat/volcell/(1000.D0*LSTAR)
        LATy  =  ygradlat/volcell/(1000.D0*LSTAR)
        LATz  =  zgradlat/volcell/(1000.D0*LSTAR)
        if ((LATx*LATx+LATy*LATy+LATz*LATz).gt.0.01d0) then

           v1 = vert_of_cell_3d(1,i)
           v2 =	vert_of_cell_3d(2,i)
           v3 =	vert_of_cell_3d(3,i)
           v4 =	vert_of_cell_3d(4,i)
           lat_cell_vert(1) = LATnode_3d(v1)
           lat_cell_vert(2) = LATnode_3d(v2)
           lat_cell_vert(3) = LATnode_3d(v3)
           lat_cell_vert(4) = LATnode_3d(v4)
           indmax = maxloc(lat_cell_vert(:),1)
           indmin = minloc(lat_cell_vert(:),1)
           v_max = vert_of_cell_3d(indmax,i)
           v_min = vert_of_cell_3d(indmin,i)

#ifdef ELEGEO0
        xV_max = xyz0_3d(1,v_max)
        yV_max = xyz0_3d(2,v_max)
        zV_max = xyz0_3d(3,v_max)
        xV_min = xyz0_3d(1,v_min)
        yV_min = xyz0_3d(2,v_min)
        zV_min = xyz0_3d(3,v_min)
#else   
        xV_max = xyz_3d(1,v_max)
        yV_max = xyz_3d(2,v_max)
        zV_max = xyz_3d(3,v_max)  
        xV_min = xyz_3d(1,v_min)
        yV_min = xyz_3d(2,v_min)
        zV_min = xyz_3d(3,v_min)
#endif 
        dxlat= (xV_max-xV_min)*LATx/modgradlat
        dylat= (yV_max-yV_min)*LATy/modgradlat
        dzlat= (zV_max-zV_min)*LATz/modgradlat
        cvmod =  sqrt(dxlat*dxlat+dylat*dylat+dzlat*dzlat)/&
                 (LATnode_3d(v_max)-LATnode_3d(v_min))/modgradlat
        else
           cvmod = 0.0d0
        endif

        CVcell(1,i)=LATx*cvmod
        CVcell(2,i)=LATy*cvmod
        CVcell(3,i)=LATz*cvmod

        
     enddo
     !@cuf istat = cudaDeviceSynchronize !JDR TMP
       !find CV on faces
#ifdef USE_CUDA
       !$cuf kernel do (1)
#endif
       do i=fstart_3d(1),fend_3d(1)
          c1 = cell_of_face_3d(1,i)
          c2 = cell_of_face_3d(2,i)
          
          g1interp = g1interpface_3d(i)
          xvers  = versCFface_3d(1,i)
          yvers  = versCFface_3d(2,i)
          zvers  = versCFface_3d(3,i)
          distCF = distCFface_3d(i)
          if ((c1.NE.0).AND.(c2.NE.0)) then
             chamb1 = cell_to_chamb_3d(c1)
             chamb2 = cell_to_chamb_3d(c2)
             if ((chamb1.lt.5).AND.(chamb2.lt.5)) then
                !grandiente medio
                cvf1 = g1interp*CVcell(1,c1) + (1-g1interp)*CVcell(1,c2)
                cvf2 = g1interp*CVcell(2,c1) + (1-g1interp)*CVcell(2,c2)
                cvf3 = g1interp*CVcell(3,c1) + (1-g1interp)*CVcell(3,c2)
             else
                cvf1 = 0.0d0
                cvf2 = 0.0d0
                cvf3 = 0.0d0
             endif
          elseif ((c1.NE.0).AND.(c2.EQ.0)) then
             cvf1 = CVcell(1,c1)
             cvf2 = CVcell(2,c1)
             cvf3 = CVcell(3,c1)
          elseif ((c1.EQ.0).AND.(c2.NE.0)) then
             cvf1 = CVcell(1,c2)
             cvf2 = CVcell(2,c2)
             cvf3 = CVcell(3,c2)
          endif
          CVface_3d(1,i) = cvf1
          CVface_3d(2,i) = cvf2
          CVface_3d(3,i) = cvf3
       enddo
       !@cuf istat = cudaDeviceSynchronize !JDR TMP

!Find grad(CV) on cells
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=cstart_3d(1),cend_3d(1)
        grad11=0.0D0
        grad12=0.0D0
        grad13=0.0D0
        grad21=0.0D0
        grad22=0.0D0
        grad23=0.0D0
        grad31=0.0D0
        grad32=0.0D0
        grad33=0.0D0
        do j=1,4
           f1=face_of_cell_3d(j,i)
           xnormale=normalfaceofcells_3d(1,j,i)
           ynormale=normalfaceofcells_3d(2,j,i)
           znormale=normalfaceofcells_3d(3,j,i)

#ifdef ELEGEO0
          cvf1 = CVface_3d(1,f1)*sur0_3d(f1)
          cvf2 = CVface_3d(2,f1)*sur0_3d(f1)
          cvf3 = CVface_3d(3,f1)*sur0_3d(f1)
#else
          cvf1 = CVface_3d(1,f1)*sur_3d(f1)
          cvf2 = CVface_3d(2,f1)*sur_3d(f1)
          cvf3 = CVface_3d(3,f1)*sur_3d(f1)
#endif
          grad11 = grad11 + cvf1*xnormale
          grad12 = grad12 + cvf2*xnormale
          grad13 = grad13 + cvf3*xnormale
          
          grad21 = grad21 + cvf1*ynormale
          grad22 = grad22 + cvf2*ynormale
          grad23 = grad23 + cvf3*ynormale
          
          grad31 = grad31 + cvf1*znormale
          grad32 = grad32 + cvf2*znormale
          grad33 = grad33 + cvf3*znormale
          
        enddo
#ifdef ELEGEO0
        volcell = vol0_3d(i)
#else
        volcell = vol_3d(i)
#endif

        CVgrad_cell(1,1,i)=grad11
        CVgrad_cell(1,2,i)=grad12
        CVgrad_cell(1,3,i)=grad13

        CVgrad_cell(2,1,i)=grad21
        CVgrad_cell(2,2,i)=grad22
        CVgrad_cell(2,3,i)=grad23
        
        CVgrad_cell(3,1,i)=grad31
        CVgrad_cell(3,2,i)=grad32
        CVgrad_cell(3,3,i)=grad33

        
     enddo
     !@cuf istat = cudaDeviceSynchronize !JDR TMP

     
     !find nodal values CV
#ifdef USE_CUDA
     !$cuf kernel do (1)
#endif
     do i=vstart_3d(1),vend_3d(1)
        num1     = 0.0d0
        num2     = 0.0d0
        num3     = 0.0d0
        den      = 0.0d0
        cvlocmax = 0.0d0
        cvlocmin = 0.0D0
#ifdef ELEGEO0
        xV1 = xyz0_3d(1,i)
        yV1 = xyz0_3d(2,i)
        zV1 = xyz0_3d(3,i)
#else
        xV1 = xyz_3d(1,i)
        yV1 = xyz_3d(2,i)
        zV1 = xyz_3d(3,i)
#endif
        do j=1,n_cell_of_vert_3d(i)
           c1=cell_of_vert_3d(j,i)
           if (cell_to_chamb_3d(c1).lt.5) then!ne.9) then !MOD_BC
              xC1 = cell_bar(1,c1)
              yC1 = cell_bar(2,c1)
              zC1 = cell_bar(3,c1)
              dvvjj = sqrt( (xV1-xC1)**2+(yV1-yC1)**2+(zV1-zC1)**2)
              num1 = num1 + CVcell(1,c1)/dvvjj
              num2 = num2 + CVcell(2,c1)/dvvjj
              num3 = num3 + CVcell(3,c1)/dvvjj
              den = den + 1.D0/dvvjj
           endif !MOD_BC
        enddo
        cvloc =1.0d0
        if (den.NE.0.0D0) then
           CVnode(1,i)=num1/den*cvloc
           CVnode(2,i)=num2/den*cvloc
           CVnode(3,i)=num3/den*cvloc
        else
           CVnode(1,i)=0.0D0
           CVnode(2,i)=0.0D0
           CVnode(3,i)=0.0D0
        endif
     enddo
     !@cuf istat = cudaDeviceSynchronize !JDR TMP

     !find nodal values grad(CV)
#ifdef USE_CUDA
     !$cuf kernel do (1)
#endif
     do i=vstart_3d(1),vend_3d(1)

        grad11=0.0D0
        grad12=0.0D0
        grad13=0.0D0

        grad21=0.0D0
        grad22=0.0D0
        grad23=0.0D0

        grad31=0.0D0
        grad32=0.0D0
        grad33=0.0D0
        
        den   =0.0d0
#ifdef ELEGEO0
        xV1 = xyz0_3d(1,i)
        yV1 = xyz0_3d(2,i)
        zV1 = xyz0_3d(3,i)
#else
        xV1 = xyz_3d(1,i)
        yV1 = xyz_3d(2,i)
        zV1 = xyz_3d(3,i)
#endif
        do j=1,n_cell_of_vert_3d(i)
           c1=cell_of_vert_3d(j,i)
           if (cell_to_chamb_3d(c1).lt.5) then!ne.9) then !MOD_BC
              xC1 = cell_bar(1,c1)
              yC1 = cell_bar(2,c1)
              zC1 = cell_bar(3,c1)
              dvvjj = sqrt( (xV1-xC1)**2+(yV1-yC1)**2+(zV1-zC1)**2)

              grad11 = grad11 + CVgrad_cell(1,1,c1)/dvvjj
              grad12 = grad12 +	CVgrad_cell(1,2,c1)/dvvjj
              grad13 = grad13 +	CVgrad_cell(1,3,c1)/dvvjj
              
              grad21 = grad21 + CVgrad_cell(2,1,c1)/dvvjj
              grad22 = grad22 +	CVgrad_cell(2,2,c1)/dvvjj
              grad23 = grad23 +	CVgrad_cell(2,3,c1)/dvvjj

              grad31 = grad31 + CVgrad_cell(3,1,c1)/dvvjj
              grad32 = grad32 +	CVgrad_cell(3,2,c1)/dvvjj
              grad33 = grad33 +	CVgrad_cell(3,3,c1)/dvvjj

              den = den + 1.D0/dvvjj
           endif !MOD_BC
        enddo

        CVgrad_node(1,1,i)=grad11
        CVgrad_node(1,2,i)=grad12
        CVgrad_node(1,3,i)=grad13

        CVgrad_node(2,1,i)=grad21
        CVgrad_node(2,2,i)=grad22
        CVgrad_node(2,3,i)=grad23
        
        CVgrad_node(3,1,i)=grad31
        CVgrad_node(3,2,i)=grad32
        CVgrad_node(3,3,i)=grad33

        CVdiv(i)   = grad11 + grad22 + grad33

        CVrot(1,i) = grad32 - grad23
        CVrot(2,i) = grad13 - grad31
        CVrot(3,i) = grad21 - grad12
        
     enddo
     !@cuf istat = cudaDeviceSynchronize !JDR TMP


!find nodal values of dV/dt
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
!do i=vstart_3d(1),vend_3d(1)
       do i=vstart_3d(1),vend_3d(1)
        num=0
#ifdef BIDOMAIN
        numext=0
#endif
        den=0
#ifdef ELEGEO0
       xV1 = xyz0_3d(1,i)
       yV1 = xyz0_3d(2,i)
       zV1 = xyz0_3d(3,i)
#else
       xV1 = xyz_3d(1,i)
       yV1 = xyz_3d(2,i)
       zV1 = xyz_3d(3,i)
#endif
       
       do j=1,n_cell_of_vert_3d(i)
            c1=cell_of_vert_3d(j,i)
            if (cell_to_chamb_3d(c1).lt.5) then!ne.9) then !MOD_BC
               xC1 = cell_bar(1,c1)
               yC1 = cell_bar(2,c1)
               zC1 = cell_bar(3,c1)

               dvvjj = sqrt( (xV1-xC1)**2+(yV1-yC1)**2+(zV1-zC1)**2)
               num = num + potEFcell_3d(c1)/dvvjj
               den = den + 1.D0/dvvjj
#ifdef BIDOMAIN
               numext = numext + potextEFcell_3d(c1)/dvvjj
#endif
            endif !MOD_BC
        enddo
        if (den.NE.0.0D0) then
           XEF_3ddt(i)=(num/den-potEFnode_3d(i))/dtMS
        else
           XEF_3ddt(i)=-100000
        endif
#ifdef BIDOMAIN
        if (den.NE.0.0D0) then
           potextEFnode_3d(i)=numext/den
        else
           potextEFnode_3d(i)=-83.0D0
        endif
#endif
    enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP
       
       
#ifdef BIDOMAIN

      call boxEFbido(checkEFbido(1:nctot_3d,1),potEFcell_3d(1:nctot_3d),0) !check vero sistemalin
      call boxEFbido(checkEFbido(1:nctot_3d,2),potextEFcell_3d(1:nctot_3d),1)  !check vero sistemalin

      countjgmres=0
      do j=1,maxrestEFbido
          call SolveVextGmres(errorEFbido)
          countjgmres=countjgmres+1
!          write(*,*) 'GMRES', ntime, j, errorEFbido 
          if (errorEFbido.LT.tolEFbido) go to 667
       enddo
667 continue 

      call boxEFbido(checkEFbido(1:nctot_3d,3),potextEFcell_3d(1:nctot_3d),1) !check vero sistemalin

       scarto1gmres=0.0!check vero sistemalin
       scarto2gmres=0.0!check vero sistemalin
#ifdef USE_CUDA 
       !$cuf kernel do (1)
#endif
       do i=1,nctot_3d
          scarto1gmres = scarto1gmres + ( checkEFbido(i,1) - checkEFbido(i,2) )**2
          scarto2gmres = scarto2gmres + ( checkEFbido(i,1) - checkEFbido(i,3) )**2
       enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

#endif

       
! !TEMPI
! #ifdef USE_CUDA 
!        !$cuf kernel do (1)
! #endif
!        do i=1,nvtot_3d
!           chamb=vert_to_chamb_3d(i)
!           potv = potEFnode_3d(i)
!           if ((chamb.LE.4).AND.(potv.GT.-30.0)) then
!              tempov=EFtstart_3d(i)
!              if (nBeat.GT.floor(tempov/period)) then     
!                 EFtstart_3d(i)=time
!              endif
!           endif
!        enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP
  
      return
      end

!===================================================
#ifdef BIDOMAIN
      subroutine SolveVextGmres(errorEFbido)
      use constants
      use mls_param
!@cuf   use cudafor

      implicit none
      real(DP):: bEFbido_norm,rEFbido_norm,errorEFbido
      real(DP):: bEFbido_norm2,x1,happo,qppEFbido_norm,tempEFbido
      real(DP),dimension(maxitEFbido):: snEFbido,csEFbido
      real(DP),dimension(maxitEFbido+1):: eEFbido,betaEFbido
      integer:: i,iEF,kEF,inp,kEFcount,INFOLS
      integer:: vsi,vei,esi,eei,fsi,fei,csi,cei

      real(DP),dimension(nctot_3d) :: bEFbido,yEFbido,rEFbido,qpEFbido,qppEFbido
!@cuf   integer :: istat
#ifdef USE_CUDA
      attributes(managed) :: bEFbido,yEFbido,rEFbido,qpEFbido,qppEFbido,bEFbido_norm,rEFbido_norm,qppEFbido_norm,betaEFbido
!      attributes(device) :: pippo,bEFbido_norm
#endif


      ! write(*,*) "tol=",tolEFbido
      ! write(*,*) "mit=",maxitEFbido

      ! call boxEFbido(bEFbido(1:nctot_3d))
      ! write(*,*) 'b1',ntime,bEFbido(1345)
      call boxEFbido(bEFbido(1:nctot_3d),potEFcell_3d,0)
      call boxEFbido(yEFbido(1:nctot_3d),potextEFcell_3d,1)

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=cstart_3d(1),cend_3d(1)
         rEFbido(i)=bEFbido(i) - yEFbido(i)
      enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!      call calculate_norm(bEFbido_norm,1,nctot_3d,bEFbido)       
      inp=1;
      csi = cstart_3d(inp) ; cei = cend_3d(inp)
      call calculate_norm(bEFbido_norm,csi,cei,bEFbido(csi:cei))
      call calculate_norm(rEFbido_norm,csi,cei,rEFbido(csi:cei))

      errorEFbido = rEFbido_norm / bEFbido_norm
      ! bEFbido_norm2=0.0
      ! do i=1,nctot_3d
      !    bEFbido_norm2=bEFbido_norm2+bEFbido(i)**2
      ! enddo
      ! bEFbido_norm2 = sqrt(bEFbido_norm2)
      ! write(*,*) 'b1',ntime,bEFbido_norm
      ! write(*,*) 'b2',ntime,bEFbido_norm2

      snEFbido=0.0
      csEFbido=0.0
      eEFbido=0.0
      eEFbido(1)=errorEFbido
      betaEFbido=0.0
      betaEFbido(1)=rEFbido_norm


      x1 = rEFbido_norm !strano
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=csi,cei
         QEFbido(i,1)= rEFbido(i)/x1
      enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

      !LOOP
      kEFcount=0
      do kEF=1,maxitEFbido
      
!####### START ARNOLDI FUNCTION ################  
         call boxEFbido(qpEFbido(1:nctot_3d),QEFbido(1:nctot_3d,kEF),1)

         do iEF=1,kEF
            !h(i) = q' * Q(:, i);
            happo=0.0
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
            do i=csi,cei
               happo = happo + qpEFbido(i)*QEFbido(i,iEF)
            enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP
 
            HEFbido(iEF,kEF)=happo

            !q = q - h(i) * Q(:, i);
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
            do i=cstart_3d(1),cend_3d(1)
               qppEFbido(i) = qpEFbido(i) - happo*QEFbido(i,iEF)
            enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP
         enddo !-> iEF=1,kEF

         ! h(k + 1) = norm(q)
            call calculate_norm(qppEFbido_norm,csi,cei,qppEFbido(csi:cei))
            HEFbido(kEF+1,kEF)=qppEFbido_norm

         !q = q / h(k + 1);
            x1 = qppEFbido_norm !strano
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
            do i=cstart_3d(1),cend_3d(1)
               QEFbido(i,kEF+1) = qppEFbido(i)/x1
            enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!####### END ARNOLDI FUNCTION ################


!####### START APPLY ROTATION ################
            do iEF=1,kEF-1
               tempEFbido = csEFbido(iEF) * HEFbido(iEF,kEF) + snEFbido(iEF) * HEFbido(iEF + 1,kEF)
               HEFbido(iEF + 1,kEF) = -snEFbido(iEF) * HEFbido(iEF,kEF) + csEFbido(iEF) * HEFbido(iEF + 1,kEF)
               HEFbido(iEF,kEF) = tempEFbido 
            enddo

            !terza funzione nested nella seconda
            tempEFbido = sqrt( HEFbido(kEF,kEF)**2 + HEFbido(kEF+1,kEF)**2 ) 
            csEFbido(kEF) = HEFbido(kEF,kEF)/tempEFbido
            snEFbido(kEF) = HEFbido(kEF+1,kEF)/tempEFbido
            
            !eliminate H(i + 1, i)
            HEFbido(kEF,kEF) = csEFbido(kEF) * HEFbido(kEF,kEF) + snEFbido(kEF) * HEFbido(kEF + 1,kEF);
            HEFbido(kEF + 1,kEF) = 0.0;
!####### END APPLY ROTATION ################
            
            ! update the residual vector
            betaEFbido(kEF + 1) = -snEFbido(kEF) * betaEFbido(kEF);
            betaEFbido(kEF)     =  csEFbido(kEF) * betaEFbido(kEF);
            errorEFbido = abs(betaEFbido(kEF + 1)) / bEFbido_norm

            ! save the error
            eEFbido(kEF+1)=errorEFbido

            kEFcount=kEFcount+1
            ! if(errorEFbido.lt.tolEFbido) then
            !    write(*,*) "gmres converged",kEF,errorEFbido
            !    go to 666
            ! endif
            if(errorEFbido.lt.tolEFbido) go to 666


 !           write(*,*) ntime,kEF,errorEFbido

      enddo  !kEF=1,maxitEFbido
      
666 continue
      !solve hessenberg system
      CALL DGESV( kEFcount, 1, HEFbido(1:kEFcount,1:kEFcount), kEFcount, IPIVEFbido(1:kEFcount), betaEFbido(1:kEFcount), kEFcount, INFOLS )

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=csi,cei
         do kEF=1,KEFcount
            potextEFcell_3d(i)= potextEFcell_3d(i) + QEFbido(i,kEF)*betaEFbido(kEF)
         enddo
      enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP



      return 
      end

!===================================================
      subroutine boxEFbido(vettoreout_cell,vettorein_cell,Mflag)
      use constants
      use mls_param
!@cuf   use cudafor

      implicit none

      real(DP):: num,den,xv1,yv1,zv1,xc1,yc1,zc1,dvvjj
      real(DP):: xnormale,ynormale,znormale,potsur,volcell
      real(DP):: xvers,yvers,zvers,Istim,xgrad,ygrad,zgrad
      real(DP):: distCF,g1interp,corrsten,bordo
      real(DP):: mint11,mint12,mint13,mint21,mint22,mint23,mint31,mint32,mint33
      integer:: i,j,c1,v1,v2,v3,f1,c2,Mflag
      character*150 :: stri
      real(DP),dimension(nctot_3d),intent(out) :: vettoreout_cell
      real(DP),dimension(nctot_3d),intent(in) :: vettorein_cell
      real(DP),dimension(nvtot_3d) :: vettorein_node
      real(DP),dimension(nftot_3d) :: vettorein_face
!@cuf   integer :: istat
#ifdef USE_CUDA
      attributes(managed) :: vettoreout_cell,vettorein_cell,vettorein_node,vettorein_face
#endif


!!COMPUTE GRADIENT
!find nodal values
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
    do i=vstart_3d(1),vend_3d(1)
        num=0
        den=0
#ifdef ELEGEO0        
       xV1 = xyz0_3d(1,i)
       yV1 = xyz0_3d(2,i)
       zV1 = xyz0_3d(3,i)
#else
       xV1 = xyz_3d(1,i)
       yV1 = xyz_3d(2,i)
       zV1 = xyz_3d(3,i)
#endif
        do j=1,n_cell_of_vert_3d(i)
            c1=cell_of_vert_3d(j,i)
            xC1 = cell_bar(1,c1)
            yC1 = cell_bar(2,c1)
            zC1 = cell_bar(3,c1)

            dvvjj = sqrt( (xV1-xC1)**2+(yV1-yC1)**2+(zV1-zC1)**2)
            num = num + vettorein_cell(c1)/dvvjj
            den = den + 1.D0/dvvjj
        enddo
        vettorein_node(i)=num/den
    enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!find face values
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
    do i=fstart_3d(1),fend_3d(1)
       v1=vert_of_face_3d(1,i)
       v2=vert_of_face_3d(2,i)
       v3=vert_of_face_3d(3,i)
       vettorein_face(i)=(vettorein_node(v1)+vettorein_node(v2)+vettorein_node(v3))/3.D0;
    enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!find gradient on cells (gauss-green theorem)
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=cstart_3d(1),cend_3d(1)
        xgrad=0.0D0
        ygrad=0.0D0
        zgrad=0.0D0
        do j=1,4
           f1=face_of_cell_3d(j,i)
           xnormale=normalfaceofcells_3d(1,j,i)
           ynormale=normalfaceofcells_3d(2,j,i)
           znormale=normalfaceofcells_3d(3,j,i)

#ifdef ELEGEO0           
           potsur = vettorein_face(f1)*sur0_3d(f1)
#else
           potsur = vettorein_face(f1)*sur_3d(f1)
#endif
           xgrad = xgrad + potsur*xnormale
           ygrad = ygrad + potsur*ynormale
           zgrad = zgrad + potsur*znormale
        enddo
#ifdef ELEGEO0
        volcell = vol0_3d(i)
#else
        volcell = vol_3d(i)
#endif
        gradcell_3d(1,i) =  xgrad/volcell
        gradcell_3d(2,i) =  ygrad/volcell
        gradcell_3d(3,i) =  zgrad/volcell
     enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!find gradient on faces
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
     do i=fstart_3d(1),fend_3d(1)
        c1 = cell_of_face_3d(1,i)
        c2 = cell_of_face_3d(2,i)

        chambcheckbc = cell_to_chamb_3d(c1)-cell_to_chamb_3d(c2)
        g1interp = g1interpface_3d(i)
        xvers  = versCFface_3d(1,i)
        yvers  = versCFface_3d(2,i)
        zvers  = versCFface_3d(3,i)
        distCF = distCFface_3d(i)
        if ((c1.NE.0).AND.(c2.NE.0)) then
           !grandiente medio
           xgrad = g1interp*gradcell_3d(1,c1) + (1-g1interp)*gradcell_3d(1,c2)
           ygrad = g1interp*gradcell_3d(2,c1) + (1-g1interp)*gradcell_3d(2,c2)
           zgrad = g1interp*gradcell_3d(3,c1) + (1-g1interp)*gradcell_3d(3,c2)
           !gradiente corretto
            !gradientecorretto=gradientemedio+(-dot(gradientemedio,versCF)+(pot(c2)-pot(c1))/distCF)*versCF;
           corrsten= -(xgrad*xvers+ygrad*yvers+zgrad*zvers) + (vettorein_cell(c2)-vettorein_cell(c1))/distCF
           xgrad = xgrad + corrsten*xvers
           ygrad = ygrad + corrsten*yvers
           zgrad = zgrad + corrsten*zvers           

        elseif ((c1.NE.0).AND.(c2.EQ.0)) then
           ! xgrad = gradcell_3d(1,c1)
           ! ygrad = gradcell_3d(2,c1)
           ! zgrad = gradcell_3d(3,c1)
           xgrad = 0.0D0 !for Neumann bcs
           ygrad = 0.0D0 
           zgrad = 0.0D0 
        elseif ((c1.EQ.0).AND.(c2.NE.0)) then
           ! xgrad = gradcell_3d(1,c2)
           ! ygrad = gradcell_3d(2,c2)
           ! zgrad = gradcell_3d(3,c2)
           xgrad = 0.0D0 !for Neumann bcs
           ygrad = 0.0D0 
           zgrad = 0.0D0
        elseif (chambcheckbc.NE.0)
           xgrad = 0.0D0 !for Neumann bcs
           ygrad = 0.0D0 
           zgrad = 0.0D0
           
        endif

        !tolto switch bidomain
        gradface_3d(1,i) = xgrad
        gradface_3d(2,i) = ygrad
        gradface_3d(3,i) = zgrad

     enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!CELL MODEL
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=cstart_3d(1),cend_3d(1)
!spatial term
        bordo = 0
        do j=1,4
           f1=face_of_cell_3d(j,i)
           xnormale=normalfaceofcells_3d(1,j,i)
           ynormale=normalfaceofcells_3d(2,j,i)
           znormale=normalfaceofcells_3d(3,j,i)

           xgrad = gradface_3d(1,f1) 
           ygrad = gradface_3d(2,f1) 
           zgrad = gradface_3d(3,f1) 

           if (Mflag.EQ.0) then
              mint11 = - Mintfaces_3d(1,1,f1)
              mint12 = - Mintfaces_3d(1,2,f1)
              mint13 = - Mintfaces_3d(1,3,f1)
              mint21 = - Mintfaces_3d(2,1,f1)
              mint22 = - Mintfaces_3d(2,2,f1)
              mint23 = - Mintfaces_3d(2,3,f1)
              mint31 = - Mintfaces_3d(3,1,f1)
              mint32 = - Mintfaces_3d(3,2,f1)
              mint33 = - Mintfaces_3d(3,3,f1)
           else 
              mint11 = Mintfaces_3d(1,1,f1)
              mint12 = Mintfaces_3d(1,2,f1)
              mint13 = Mintfaces_3d(1,3,f1)
              mint21 = Mintfaces_3d(2,1,f1)
              mint22 = Mintfaces_3d(2,2,f1)
              mint23 = Mintfaces_3d(2,3,f1)
              mint31 = Mintfaces_3d(3,1,f1)
              mint32 = Mintfaces_3d(3,2,f1)
              mint33 = Mintfaces_3d(3,3,f1)

              mint11 = mint11 + Mextfaces_3d(1,1,f1)
              mint12 = mint12 + Mextfaces_3d(1,2,f1)
              mint13 = mint13 + Mextfaces_3d(1,3,f1)
              mint21 = mint21 + Mextfaces_3d(2,1,f1)
              mint22 = mint22 + Mextfaces_3d(2,2,f1)
              mint23 = mint23 + Mextfaces_3d(2,3,f1)
              mint31 = mint31 + Mextfaces_3d(3,1,f1)
              mint32 = mint32 + Mextfaces_3d(3,2,f1)
              mint33 = mint33 + Mextfaces_3d(3,3,f1)
           endif

           bordo = bordo + &
                (  xnormale*(mint11*xgrad + mint12*ygrad + mint13*zgrad) &
                + ynormale*(mint21*xgrad + mint22*ygrad + mint23*zgrad) &
#ifdef ELEGEO0
                + znormale*(mint31*xgrad + mint32*ygrad + mint33*zgrad) )*sur_3d(f1)
#else
           + znormale*(mint31*xgrad + mint32*ygrad + mint33*zgrad) )*sur0_3d(f1)
#endif
        enddo
           vettoreout_cell(i)=bordo
      enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP


      return 
      end
#endif
!===================================================

!===================================================
      subroutine ElectroStartEF_3d
      use constants
      use mls_param
      use scar_tagging
!@cuf   use cudafor

      implicit none

      character(len=4):: dummy
       
      real(DP):: timex,timeEFv,alpv1
      real(DP):: Volt,Cai,CaSR,CaSS
      real(DP):: Nai,Ki,INa_m,INa_h 
      real(DP):: INa_j,IKr_xr1,IKr_xr2,IKs_xs 
      real(DP):: Ito_r,Ito_s,ICaL_d,ICaL_f
      real(DP):: ICaL_f2,ICaL_fCaSS,RR


      !real(DP):: xBc,yBc,zBc
      real(DP):: off,off0,offminLV,offminRV,offminRA
      real(DP):: xLV,yLV,zLV,xLA,yLA,zLA 
      real(DP):: xRV,yRV,zRV,xRA,yRA,zRA 
      real(DP):: distveG2,px,py,pz
      real(DP):: Chi,Cm,L_signal,A_signal,amplitude_signal
      integer:: i,indv1,indv2,act,inp,chamb,c1,c2
      integer, dimension(nctot_3d) :: scar_id
      ! integer:: lead_LV,lead_RV,lead_RA
      character*150 :: stri

      !Label Smooth extension
      real(DP):: xC1,yC1,zC1
      real(DP):: xV1,yV1,zV1
      integer:: step,j
      real(DP):: num,den,dvvjj
      
      
      !-------------------------------------------------------
      ! S1S2-Stimulation protocol variables
      real(DP):: tmp_2dS1,alpha,amplitudeS1,amplitudeS2
      real(DP):: S11,S12,S13,S14
      real(DP):: S21,S22,S23,S24
      real(DP):: xS1,yS1,zS1
      real(DP):: xS2,yS2,zS2
      real(DP):: pxS1,pyS1,pzS1
      real(DP):: pxS2,pyS2,pzS2
      real(DP):: distS1,distS2
      real(DP):: st_dev3D,MFstim
      integer::v1,v2,v3,v4
      !-------------------------------------------------------

      !Stenosi FV
      integer:: countSTE,countVent,countSTElv,f1,f2,f3,f4
      real(DP)::s1,s2,s3,s4,supmin,supmax,ratios
      real(DP)::x1,y1,z1
      real(DP)::x2,y2,z2
      real(DP)::x3,y3,z3
      real(DP)::x4,y4,z4
      real(DP)::e1,e2,e3,e4
      real(DP)::rhs1,rhs2,rhs3,rhs4
      real(DP)::rin,rout,a,b,c,d,r1,r2,r3,r4
      real(DP)::a0,b0,c0
      !real(DP)::a1,b1,c1
      real(DP)::xf,yf,zf
      real(DP)::xBc,yBc,zBc
      real(DP):: xnormale,ynormale,znormale
      real(DP)::sum_vol,avg_vol,avg_skew

      real(DP)::tmpDcart,scaleD,g1interp
      
      real(DP)::npt_abl
      real(DP),allocatable,dimension(:,:)::data_abl

#ifdef USE_CUDA
      attributes(managed)::data_abl
#endif
!@cuf   integer :: istat
    
    !  write(*,*)'Entering scaray'
    !  call scaray(0,0,1) !Tag cells in the Gray_zone/Scar_core                                    
    !  call write_scar_vtk

      ! !%Initial values of state variables
#ifdef TP06
      Volt = -85.23; !%membrane potential
      Cai = 0.000126; !%intracellular calcium
      CaSR = 3.64; !%sarcoplasmic reticulum calcium
      CaSS = 0.00036; !%subspace calcium
      Nai = 8.604; !%intracellular sodium
      Ki = 136.89; !%intracellular potassium
      INa_m = 0.00172; !%fast sodium current m gate
      INa_h = 0.7444; !%fast sodium current h gate
      INa_j = 0.7045; !%fast sodium current j gate
      IKr_xr1 = 0.00621; !%rapid time-dependent potassium current Xr1 gate
      IKr_xr2 = 0.4712;  !%rapid time-dependent potassium current Xr2 gate
      IKs_xs = 0.0095; !%rapid time-dependent potassium current Xs gate
      Ito_r = 2.42e-08; !%transient outward current r gate
      Ito_s = 0.999998; !%transient outward current s gate
      ICaL_d = 3.373e-05; !%L-type Ca current d gate
      ICaL_f = 0.7888; !%L-type Ca current f gate
      ICaL_f2 = 0.9755; !%L-type Ca current f2 gate
      ICaL_fCaSS = 0.9953; !%L-type Ca current fCass gate
      RR = 0.9073; !%ryanodine receptor R_prime
#endif
      !scar_id(:) = XEF_3d(22,:)
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=cstart_3d(1),cend_3d(1)

         chamb = cell_to_chamb_3d(i)
         
       if ((chamb.EQ.1).OR.(chamb.EQ.3)) then
#ifdef TP06
            XEF_3d(1,i) = Volt 
            XEF_3d(2,i) = Cai 
            XEF_3d(3,i) = CaSR 
            XEF_3d(4,i) = CaSS 
            XEF_3d(5,i) = Nai 
            XEF_3d(6,i) = Ki 
            XEF_3d(7,i) = INa_m 
            XEF_3d(8,i) = INa_h 
            XEF_3d(9,i) = INa_j
            XEF_3d(10,i) = IKr_xr1 
            XEF_3d(11,i) = IKr_xr2 
            XEF_3d(12,i) = IKs_xs 
            XEF_3d(13,i) = Ito_r 
            XEF_3d(14,i) = Ito_s 
            XEF_3d(15,i) = ICaL_d 
            XEF_3d(16,i) = ICaL_f 
            XEF_3d(17,i) = ICaL_f2 
            XEF_3d(18,i) = ICaL_fCaSS 
            XEF_3d(19,i) = RR

            potEFcell_3d(i) = XEF_3d(1,i)
#endif
#ifdef MINIMAL_MODEL_AAA
            XEF_3d(1,i) = 0.0D0;
            XEF_3d(2,i) = 1.0D0;
            XEF_3d(3,i) = 1.0D0;
            XEF_3d(4,i) = 0.0D0;

            potEFcell_3d(i) = -83.0D0!XEF_3d(1,i)!+XEF_3d(1,i)*(2.7D0 - 83.0D0)!-83.0d0!XEF_3d(1,i)

            minpotLV=-83.0D0
            minpotRV=-83.0D0
#endif
#ifdef MITCHELL_SCHAEFFER
            XEF_3d(1,i) = 0.0D0;
            XEF_3d(2,i) = 1.0D0;
            XEF_3d(3,i) = 0.0D0;
            XEF_3d(4,i) = 0.0D0;

            potEFcell_3d(i) = -83.0D0!XEF_3d(1,i)!+XEF_3d(1,i)*(2.7D0 - 83.0D0)!-83.0d0!XEF_3d(1,i)

            minpotLV=-83.0D0
            minpotRV=-83.0D0
#endif


            !potEFcell_3d(i) = Volt
         elseif ((chamb.EQ.2).OR.(chamb.EQ.4)) then
         ! else !devo inserire un valore 'decente' per evitare gradienti
#ifdef COURTEMANCHE
            XEF_3d(1,i) =  -81.18;
            XEF_3d(2,i) = 1.117e+01;
            XEF_3d(3,i) = 2.908e-3;
            XEF_3d(4,i) = 9.649e-1;
            XEF_3d(5,i) = 9.775e-1;
            XEF_3d(6,i) = 1.39e+02;
            XEF_3d(7,i) = 3.043e-2;
            XEF_3d(8,i) = 9.992e-1;
            XEF_3d(9,i) = 4.966e-3;
            XEF_3d(10,i) = 9.986e-1;
            XEF_3d(11,i) = 3.296e-5;
            XEF_3d(12,i) = 1.869e-2;
            XEF_3d(13,i) = 1.013e-4;
            XEF_3d(14,i) = 1.367e-4;
            XEF_3d(15,i) = 9.996e-1;
            XEF_3d(16,i) = 7.755e-1;
            XEF_3d(17,i) = 1.488;
            XEF_3d(18,i) = 2.35e-112;
            XEF_3d(19,i) = 1;
            XEF_3d(20,i) = 0.9992;
            XEF_3d(21,i) = 1.488;

            potEFcell_3d(i) = XEF_3d(1,i)

            minpotLA=-81.18D0
            minpotRA=-81.18D0
#endif
#ifdef MINIMAL_MODEL
            XEF_3d(1,i) = 0.0D0;
            XEF_3d(2,i) = 1.0D0;
            XEF_3d(3,i) = 1.0D0;
            XEF_3d(4,i) = 0.0D0;

            potEFcell_3d(i) = -83.0D0!XEF_3d(1,i)!+XEF_3d(1,i)*(2.7D0 - 83.0D0)!-83.0d0!XEF_3d(1,i)

            minpotLA=-83.0D0
            minpotRA=-83.0D0
#endif
         else

            XEF_3d(1,i) = 0.0D0;
            XEF_3d(2,i) = 1.0D0;
            XEF_3d(3,i) = 1.0D0;
            XEF_3d(4,i) = 0.0D0;

            potEFcell_3d(i) = -81.18D0

         endif
         !if (cell_bar(1,i).lt.0.5D0) cell_to_chamb_3d(i)=9 !MOD_BC
      enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

      !call write_scar_vtk      
!#ifdef SCAR_IZ
!       write(*,*)'Entering scaray'
!       call scaray(0,0,1) !Tag cells in the Gray_zone/Scar_core
!       call write_scar_vtk
      
! ! Put fiber field isotropic in the scar
! #ifdef USE_CUDA
!         !$cuf kernel do (1)
! #endif
!         do i=cstart_3d(1),cend_3d(1)
! !            if (scar_cell(i).NE.0) then     
!               AmatrFibers_cell_3d(1,1,i) = 1.0D0
!               AmatrFibers_cell_3d(2,1,i) = 0.0D0
!               AmatrFibers_cell_3d(3,1,i) = 0.0D0
              
!               AmatrFibers_cell_3d(1,2,i) = 0.0D0
!               AmatrFibers_cell_3d(2,2,i) = 1.0D0
!               AmatrFibers_cell_3d(3,2,i) = 0.0D0
           
!               AmatrFibers_cell_3d(1,3,i) = 0.0D0
!               AmatrFibers_cell_3d(2,3,i) = 0.0D0
!               AmatrFibers_cell_3d(3,3,i) = 1.0D0
! !           endif
              
              
!         enddo
!  !@cuf   istat = cudaDeviceSynchronize !JDR TMP

      
!#else
!      XEF_3d(22,:)=0. !All cells are in the healthy tissue; no scar
!#endif
      
      ! minpotLV=Volt
      ! minpotLA=-81.18
      ! minpotRV=Volt
      ! minpotRA=-81.18

      meshquality_3d(:) = 0.0D0
      astressEFcell_3d(:)=0.D0
      astressEFnode_3d(:)=0.D0
      astressEFedge_3d(:)=0.D0
!#ifdef S1S2-3D
      open(unit=15,file='S1S2_3dStim.in',status='old')
      read(15,*) dummy
      read(15,*) xS1,yS1,zS1
      read(15,*) dummy
      read(15,*) xS2,yS2,zS2
      read(15,*) dummy
      read(15,*) amplitudeS1,amplitudeS2
      read(15,*) dummy
      read(15,*) st_dev3D,MFstim
      close(15)
      
      write(*,*) '3D S1S2 Stimulus'
      write(*,*) "xS1 yS1 zS1", xS1, yS1, zS1
      write(*,*) "xS2 yS2 zS2", xS2, yS2, zS2

      !LABEL SKEWED TRET

      !AVG SUR
      countVent=0
      sum_vol=0.0D0
      ! do i=cstart_3d(1),cend_3d(1)
      !    chamb=cell_to_chamb_3d(i)
      !    if ((chamb.EQ.1).or.(chamb.EQ.3)) then
      !       sum_vol=sum_vol+vol0_3d(i)
      !       countVent=countVent+1
      !    endif
      ! enddo
      ! avg_vol = sum_vol/countVent

!       avg_vol=0.d0
!       do i=cstart_3d(1),cend_3d(1)
!          avg_vol=avg_vol+vol_3d(i)

! !!!!!!!!!!!!
!          f1=face_of_cell_3d(1,i)
!          f2=face_of_cell_3d(2,i)
!          f3=face_of_cell_3d(3,i)
!          f4=face_of_cell_3d(4,i)
!          xBc=cell_bar(1,i)
!          yBc=cell_bar(2,i)
!          zBc=cell_bar(3,i)

!          rin = 3.0D0*vol0_3d(i)/(sur_3d(f1)+sur_3d(f2)+sur_3d(f3)+sur_3d(f4))
         
!          v1 = vert_of_cell_3d(1,i)
!          v2 = vert_of_cell_3d(2,i)
!          v3 = vert_of_cell_3d(3,i)
!          v4 = vert_of_cell_3d(4,i)

!          a0 = sqrt((xyz_3d(1,v2)-xyz_3d(1,v1))**2+&
!                    (xyz_3d(2,v2)-xyz_3d(2,v1))**2+&
!                    (xyz_3d(3,v2)-xyz_3d(3,v1))**2)
!          a1 = sqrt((xyz_3d(1,v3)-xyz_3d(1,v4))**2+&
!                    (xyz_3d(2,v3)-xyz_3d(2,v4))**2+&
!                    (xyz_3d(3,v3)-xyz_3d(3,v4))**2)
         
!          b0 = sqrt((xyz_3d(1,v1)-xyz_3d(1,v4))**2+&
!                    (xyz_3d(2,v1)-xyz_3d(2,v4))**2+&
!                    (xyz_3d(3,v1)-xyz_3d(3,v4))**2)
!          b1 = sqrt((xyz_3d(1,v2)-xyz_3d(1,v3))**2+&
!                    (xyz_3d(2,v2)-xyz_3d(2,v3))**2+&
!                    (xyz_3d(3,v2)-xyz_3d(3,v3))**2)
         
!          c0 = sqrt((xyz_3d(1,v1)-xyz_3d(1,v3))**2+&
!                    (xyz_3d(2,v1)-xyz_3d(2,v3))**2+&
!                    (xyz_3d(3,v1)-xyz_3d(3,v3))**2)
!          c1 = sqrt((xyz_3d(1,v2)-xyz_3d(1,v4))**2+&
!                    (xyz_3d(2,v2)-xyz_3d(2,v4))**2+&
!                    (xyz_3d(3,v2)-xyz_3d(3,v4))**2)

!          rout = sqrt(( a0*a1 + b0*b1 + c0*c1 )*&
!                   (-a0*a1 + b0*b1 + c0*c1 )*&
!                   ( a0*a1 - b0*b1 + c0*c1 )*&
!                   ( a0*a1 + b0*b1 - c0*c1 )/(576.0D0*vol0_3d(i)**2))
         
!          ! s1=sur_3d(f1)
!          ! s2=sur_3d(f2)
!          ! s3=sur_3d(f3)
!          ! s4=sur_3d(f4)
         
!          ! supmin=min(s1,s2,s3,s4)
!          ! supmax=max(s1,s2,s3,s4)
         
!          !ratios=supmin/supmax
!          ratios=(0.3D0-(rin/rout))/0.3D0 !Normalized with the optimal ratios for a regular tetrahedron
         
!          avg_skew = avg_skew + ratios
!       enddo
!       avg_vol  = avg_vol/nctot_3d
!       avg_skew = avg_skew/nctot_3d

      
      countSTE=0
      countSTElv=0
      !IstimEF_3datria(:) = 0.0D0
      IstimEF_3dS1(:)    = 0.0D0
      IstimEF_3dS2(:)    = 0.0D0
      CVcell(:,:)        = 0.0D0
      CVnode(:,:)        = 0.0D0
      CVrot(:,:)         = 0.0D0
      CVdiv(:)           = 0.0D0
      CVgrad_cell(:,:,:) = 0.0D0
      CVgrad_node(:,:,:) = 0.0D0
      LATface_3d(:)      = 0.0D0
      LATnode_3d(:)      = 0.0D0
      !write(*,*) "countSTE = ", countSTE, countSTElv, countVent, nctot_3d
      write(*,*) "countSTE = ", countSTE, nctot_3d, real(countSTE)/real((cend_3d(1)-cstart_3d(1)))*100.d0
#ifdef CARTO
      scaleD=0.9D0/5.0D0
      open(unit=15,file='meshes/CARTO_diffusivity_3d.txt',status='old')
      do i=vstart_3d(1),vend_3d(1)
         read(15,*) tmpDcart
         if (vert_to_chamb_3d(i).eq.2) then
            CARTO_Dnode3d(i)=0.18D0
         elseif (vert_to_chamb_3d(i).ge.5) then
            CARTO_Dnode3d(i)=0.000001D0
         else
            CARTO_Dnode3d(i)=scaleD*tmpDcart 
         endif
      enddo
      close(15)
      write(*,*) 'CARTO diffusivity 3D loaded'
      
      open(unit=15,file='meshes/ABLATION_CARTO.txt',status='old')
      read(15,*) npt_abl
      allocate(data_abl(4,npt_abl))
      do i=1,npt_abl
         read(15,*) data_abl(1,i),data_abl(2,i),data_abl(3,i),data_abl(4,i)
      enddo
      close(15)
      write(*,*) 'CARTO ABLATION PATH LOADED'

      !$cuf kernel do (1)
      do i=cstart_3d(1),cend_3d(1)
         v1 = vert_of_cell_3d(1,i)
         v2 = vert_of_cell_3d(2,i)
         v3 = vert_of_cell_3d(3,i)
         v4 = vert_of_cell_3d(4,i)
         CARTO_Dcell3d(i)=0.25D0*(CARTO_Dnode3d(v1)+CARTO_Dnode3d(v2)+CARTO_Dnode3d(v3)+CARTO_Dnode3d(v4))
      enddo
      !@cuf istat = cudaDeviceSynchronize !JDR TMP

      !INTERPOLAZIONI VARIE TENSORE CONDUTTIVITA
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
    do i=fstart_3d(1),fend_3d(1)
       c1=cell_of_face_3d(1,i)
       c2=cell_of_face_3d(2,i)        
       g1interp=g1interpface_3d(i)
       if ((c1.NE.0).AND.(c2.NE.0)) then
          CARTO_Dface3d(i)=g1interp*CARTO_Dcell3d(c1)+(1-g1interp)*CARTO_Dcell3d(c2)
       elseif ((c1.NE.0).AND.(c2.EQ.0)) then
          CARTO_Dface3d(i)=CARTO_Dcell3d(c1)
       elseif ((c1.EQ.0).AND.(c2.NE.0)) then
          CARTO_Dface3d(i)=CARTO_Dcell3d(c2)
       endif    
    enddo
!@cuf   istat = cudaDeviceSynchronize !JDR TMP

#else
      
      open(unit=15,file='meshes/myotag.txt',status='old')
      do i=cstart_3d(1),cend_3d(1)
         read(15,*) myotag(i)
      enddo
      close(15)
      write(*,*) 'Myocardium tagging loaded'
#endif
      
#ifdef USE_CUDA
      !$cuf kernel do (1) 
#endif

      do i=cstart_3d(1),cend_3d(1)
         
         chamb = cell_to_chamb_3d(i)
         if (chamb.eq.1.OR.chamb.eq.3) then
            countVent = countVent+1
         endif
         if (vol0_3d(i).LT.0.15D0*avg_vol) then!.and.(ratios.LT.0.1D0*avg_skew)) then
            LabelStenosi(i)=1
            countSTE=countSTE+1
            ! if ((chamb.EQ.1).or.(chamb.EQ.3)) then
            !    countSTElv=countSTElv+1
            ! endif
         endif
      
         if ((chamb.eq.1).or.(chamb.eq.3)) then
            xBc = cell_bar(1,i)
            yBc = cell_bar(2,i)
            zBc = cell_bar(3,i)
            distveG2=sqrt((xBc-xS1)**2+(yBc-yS1)**2+(zBc-zS1)**2)
            if ((distveG2.LT.(5.D0/(1000.D0*LSTAR))).and.(scar_cell(i).ne.2)) then
!               IstimEF_3d(i)=-60.D0
               IstimEF_3dS1(i)=-0.30D0 !minimal stimolo
            endif
            distveG2=sqrt((xBc-xS2)**2+(yBc-yS2)**2+(zBc-zS2)**2)
            if ((distveG2.LT.(15.0D0/(1000.D0*LSTAR))).and.(scar_cell(i).ne.2)) then 
!               IstimEF_3d(i)=-60.D0
               IstimEF_3dS2(i)=-0.30D0 !minimal stimolo
            endif            
         endif
         
      enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP
! #endif

#ifdef CARTO      
#ifdef USE_CUDA
      !$cuf kernel do (2) 
#endif
      do j=1,npt_abl
         do i=vstart_3d(1),vend_3d(1)
            xS1 = data_abl(1,j)
            yS1 = data_abl(2,j)
            zS1 = data_abl(3,j)
            
            xBc = xyz_3d(1,i)
            yBc = xyz_3d(2,i)
            zBc = xyz_3d(3,i)
            distveG2=sqrt((xBc-xS1)**2+(yBc-yS1)**2+(zBc-zS1)**2)
            if (distveG2.LT.(1.050D0*data_abl(4,j)/(1000.0D0*LSTAR))) then
               !CARTO_Dnode3d(i)=0.00000001D0
            endif
     
         enddo
      enddo
      !@cuf istat = cudaDeviceSynchronize !JDR TMP

      !$cuf kernel do (1)
      do i=cstart_3d(1),cend_3d(1)
         v1 = vert_of_cell_3d(1,i)
         v2 = vert_of_cell_3d(2,i)
         v3 = vert_of_cell_3d(3,i)
         v4 = vert_of_cell_3d(4,i)
         CARTO_Dcell3d(i)=0.25D0*(CARTO_Dnode3d(v1)+CARTO_Dnode3d(v2)+CARTO_Dnode3d(v3)+CARTO_Dnode3d(v4))
      enddo
      !@cuf istat = cudaDeviceSynchronize !JDR TMP
#endif
      
      !write(*,*) 'Istim calculated'
      ! Write vtk with scar tags and S1 S2 stimuli profiles 
     ! call write_scar_vtk
     ! write(*,*) 'Istim calculated'
    
!posizionio elettrodi
!posiziono elettrodi
!LMV1-1     
      xLV = 2.37
      yLV =-0.59
      zLV = 2.54

      xRV = 0.61
      yRV =-0.30
      zRV = 0.067
      
      xRA = -1.03
      yRA = -0.065
      zRA =  3.75
    
      off0 = 100000.
      offminLV = off0
      offminRV = off0
      offminRA = off0
      do i=cstart_3d(1),cend_3d(1)
         chamb = cell_to_chamb_3d(i)
         xBc = cell_bar(1,i)         
         yBc = cell_bar(2,i)
         zBc = cell_bar(3,i)

         if (chamb.EQ.1) then
            off = sqrt( (xLV-xBc)**2 + (yLV-yBc)**2 + (zLV-zBc)**2 )
            if (off.LT.offminLV) then
               lead_LV=i
               offminLV = off
            endif
         elseif (chamb.EQ.3) then
            off = sqrt( (xRV-xBc)**2 + (yRV-yBc)**2 + (zRV-zBc)**2 )
            if (off.LT.offminRV) then
               lead_RV=i
               offminRV = off
            endif
         elseif (chamb.EQ.4) then
            off = sqrt( (xRA-xBc)**2 + (yRA-yBc)**2 + (zRA-zBc)**2 )
            if (off.LT.offminRA) then
               lead_RA=i
               offminRA = off
            endif
         endif

      enddo
#ifdef CARTO
      lead_RA  = lead_RV
      offminRA = offminRV
      lead_LV  = lead_RV
      offminLV = offminRV

      
#endif

      ! write(*,*) "lead_RA",xRA,yRA,zRA
      ! write(*,*) "lead_RA",cell_bar(1:3,lead_RA)
      ! write(*,*) "lead_RV",xRV,yRV,zRV
      ! write(*,*) "lead_RV",cell_bar(1:3,lead_RV)
      ! write(*,*) "lead_LV",xLV,yLV,zLV
      ! write(*,*) "lead_LV",cell_bar(1:3,lead_LV)

      call calculate_MintMext_cells


#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i = 1,nvtot_3d      
         f_apd(i) = 0
         old_1(i) = 0.0D0
      enddo
      !@cuf istat = cudaDeviceSynchronize !JDR TMP
      
!       !-------------------------------------------
!       !Extend label_smooth between atria and veins
! #ifdef USE_CUDA
!         !$cuf kernel do (1)
! #endif
!       do i = 1,nvtot_3d      
!          if (vert_to_chamb_3d(i).lt.5) then
!             LabelSmooth(i)=1.0D0
!          else
!             LabelSmooth(i)=0.0D0
!          endif

!       enddo
! !@cuf istat = cudaDeviceSynchronize !JDR TMP

!       do step=1,2
! #ifdef USE_CUDA
!         !$cuf kernel do (1)
! #endif
!          do i = 1,nvtot_3d
            
!             num=0
!             den=0
! #ifdef ELEGEO0
!             xV1 = xyz0_3d(1,i)
!             yV1 = xyz0_3d(2,i)
!             zV1 = xyz0_3d(3,i)
! #else
!             xV1 = xyz_3d(1,i)
!             yV1 = xyz_3d(2,i)
!             zV1 = xyz_3d(3,i)
! #endif
!             do j=1,n_edge_of_vert_3d(i)
!                v1=vert_of_vert_3d(j,i)
! #ifdef ELEGEO0
!                xC1 = xyz0_3d(1,v1)
!                yC1 = xyz0_3d(2,v1)
!                zC1 = xyz0_3d(3,v1)
! #else
!                xC1 = xyz_3d(1,v1)
!                yC1 = xyz_3d(2,v1)
!                zC1 = xyz_3d(3,v1)
! #endif               
!                dvvjj = sqrt( (xV1-xC1)**2+(yV1-yC1)**2+(zV1-zC1)**2)
!                num = num + LabelSmooth(v1)/dvvjj
!                den = den + 1.D0/dvvjj
!             enddo
!             LabelSmooth(i)=num/den
!          enddo
!          !@cuf istat = cudaDeviceSynchronize !JDR TMP
!       enddo
!       !$cuf kernel do (1) 
!       do i=cstart_3d(1),cend_3d(1)
!          v1 = vert_of_cell_3d(1,i)
!          v2 = vert_of_cell_3d(2,i)
!          v3 = vert_of_cell_3d(3,i)
!          v4 = vert_of_cell_3d(4,i)
         
!          LabelSmooth_cell(i)=0.25*(LabelSmooth(v1)+LabelSmooth(v2)+LabelSmooth(v3)+LabelSmooth(v4))
!       enddo
      !@cuf istat = cudaDeviceSynchronize !JDR TMP


      deallocate(data_abl)
      
      return
      end

!------------------------------------------------------
        subroutine calculate_MintMext_cells
        use constants
        use mls_param
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snc,enc,v1,v2,v3,v4,i,c1,c2,chamb,c_inf_loc
        real(DP) :: Chi,Cm,g1interp,parco
        real(DP) :: sigma_if,sigma_is,sigma_in
        real(DP) :: sigma_ef,sigma_es,sigma_en
        real(DP) :: sigma_f,sigma_s,sigma_n
        real(DP) :: xfv,yfv,zfv,xsv,ysv,zsv,xnv,ynv,znv
        real(DP) :: scalaVentr,scalaAtri,scalaElse,valM
        real(DP) :: x_inf_loc,y_inf_loc,z_inf_loc,size_inf
        real(DP) :: xcell,ycell,zcell,campo_inf,off,offmin

        !------------------------------------------------
        ! Variables for infarcted zone
        character(len=4):: dummy
        real(DP) :: dr,alpha,Px,Py,Pz
        real(DP) :: distPv1,distPv2,distPv3,distPv4
        !-----------------------------------------------
!@cuf   integer :: istat
! #ifdef USE_CUDA
!         attributes(managed) :: vert_of_cell, xyz, vol
! #endif

        
!        scalaVentr=1.0D0
        scalaVentr=2.0D0

        scalaAtri=2.0D0!1.0D0
!        scalaElse=0.5D0 !non 0 altrimenti nan
        scalaElse=0.0D0!0.00001D0!0.0D0 !non 0 altrimenti nan MESSO 0


        !MONO/BIDOMAIN MODEL
        ! if (XEF_3d(22,i).eq.1.0)then
        !    Chi=140; !%mm^-1
        !    Cm=0.01; !%mF mm^?2 
        !    sigma_if = 0.17/(Cm*Chi);  !% mS / mm
        !    sigma_is = 0.019/(Cm*Chi); !% mS / mm
        !    sigma_in = 0.019/(Cm*Chi); !% mS / mm
        !    sigma_ef = 0.62/(Cm*Chi);  !% mS / mm
        !    sigma_es = 0.62/(Cm*Chi);  !% mS / mm
        !    sigma_en = 0.62/(Cm*Chi);  !% mS / mm      
        !    sigma_f = sigma_if*sigma_ef/(sigma_if+sigma_ef);
        !    sigma_s = sigma_is*sigma_es/(sigma_is+sigma_es);
        !    sigma_n = sigma_in*sigma_en/(sigma_in+sigma_en);
        ! else
           Chi=140; !%mm^-1                     
           Cm=0.01; !%mF mm^?2                  
           sigma_if = 0.17/(Cm*Chi);  !% mS / mm
           sigma_is = 0.019/(Cm*Chi); !% mS / mm
           sigma_in = 0.019/(Cm*Chi); !% mS / mm
           sigma_ef = 0.62/(Cm*Chi);  !% mS / mm
           sigma_es = 0.24/(Cm*Chi);  !% mS / mm
           sigma_en = 0.24/(Cm*Chi);  !% mS / mm
           sigma_f = sigma_if*sigma_ef/(sigma_if+sigma_ef);
           sigma_s = sigma_is*sigma_es/(sigma_is+sigma_es);
           sigma_n = sigma_in*sigma_en/(sigma_in+sigma_en);
!         endif
        !       queste vanno sovrascritte piu in basso
      
#ifdef USE_CUDA  
        !$cuf kernel do (1)
#endif
      do i=cstart_3d(1),cend_3d(1)
         chamb = cell_to_chamb_3d(i)
         if ((chamb.EQ.1).OR.(chamb.EQ.3)) then !ventricoli
            parco = scalaVentr
         elseif ((chamb.EQ.2).OR.(chamb.EQ.4)) then !atri
            parco = scalaAtri
         else
            parco =scalaElse 
         endif

         xfv =AmatrFibers_cell_3d(1,1,i) 
         yfv =AmatrFibers_cell_3d(2,1,i) 
         zfv =AmatrFibers_cell_3d(3,1,i) 
         xsv =AmatrFibers_cell_3d(1,2,i) 
         ysv =AmatrFibers_cell_3d(2,2,i) 
         zsv =AmatrFibers_cell_3d(3,2,i) 
         xnv =AmatrFibers_cell_3d(1,3,i) 
         ynv =AmatrFibers_cell_3d(2,3,i) 
         znv =AmatrFibers_cell_3d(3,3,i)

         Mintcells_3d(1,1,i) = parco*(sigma_if*xfv**2  + sigma_in*xnv**2  + sigma_is*xsv**2)
         Mintcells_3d(1,2,i) = parco*(sigma_if*xfv*yfv + sigma_in*xnv*ynv + sigma_is*xsv*ysv)
         Mintcells_3d(1,3,i) = parco*(sigma_if*xfv*zfv + sigma_in*xnv*znv + sigma_is*xsv*zsv)
         Mintcells_3d(2,1,i) = parco*(sigma_if*xfv*yfv + sigma_in*xnv*ynv + sigma_is*xsv*ysv)
         Mintcells_3d(2,2,i) = parco*(sigma_if*yfv**2  + sigma_in*ynv**2  + sigma_is*ysv**2)
         Mintcells_3d(2,3,i) = parco*(sigma_if*yfv*zfv + sigma_in*ynv*znv + sigma_is*ysv*zsv)
         Mintcells_3d(3,1,i) = parco*(sigma_if*xfv*zfv + sigma_in*xnv*znv + sigma_is*xsv*zsv)
         Mintcells_3d(3,2,i) = parco*(sigma_if*yfv*zfv + sigma_in*ynv*znv + sigma_is*ysv*zsv)
         Mintcells_3d(3,3,i) = parco*(sigma_if*zfv**2  + sigma_in*znv**2  + sigma_is*zsv**2)

         Mextcells_3d(1,1,i) = parco*(sigma_ef*xfv**2  + sigma_en*xnv**2  + sigma_es*xsv**2)
         Mextcells_3d(1,2,i) = parco*(sigma_ef*xfv*yfv + sigma_en*xnv*ynv + sigma_es*xsv*ysv)
         Mextcells_3d(1,3,i) = parco*(sigma_ef*xfv*zfv + sigma_en*xnv*znv + sigma_es*xsv*zsv)
         Mextcells_3d(2,1,i) = parco*(sigma_ef*xfv*yfv + sigma_en*xnv*ynv + sigma_es*xsv*ysv)
         Mextcells_3d(2,2,i) = parco*(sigma_ef*yfv**2  + sigma_en*ynv**2  + sigma_es*ysv**2)
         Mextcells_3d(2,3,i) = parco*(sigma_ef*yfv*zfv + sigma_en*ynv*znv + sigma_es*ysv*zsv)
         Mextcells_3d(3,1,i) = parco*(sigma_ef*xfv*zfv + sigma_en*xnv*znv + sigma_es*xsv*zsv)
         Mextcells_3d(3,2,i) = parco*(sigma_ef*yfv*zfv + sigma_en*ynv*znv + sigma_es*ysv*zsv)
         Mextcells_3d(3,3,i) = parco*(sigma_ef*zfv**2  + sigma_en*znv**2  + sigma_es*zsv**2)
      end do
!@cuf   istat = cudaDeviceSynchronize !JDR TMP


!      open(unit=15,file='S1S2_3D_Stim.in',status='old')
!      read(15,*) dummy
!      read(15,*) dummy
!      read(15,*) dummy
!      read(15,*) dummy
!      read(15,*) dummy
!      read(15,*) dummy
!      read(15,*) dummy
!      read(15,*) dummy
!      read(15,*) dr
!      close(15)
      
#ifndef BIDOMAIN
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=cstart_3d(1),cend_3d(1)
         
         chamb = cell_to_chamb_3d(i)
         
         if (((chamb.EQ.1).OR.(chamb.EQ.3)).and.(scar_cell(i).EQ.0)) then !ventricoli 
            parco = scalaVentr
            sigma_f= 1.263609/(Cm*Chi);
         elseif (((chamb.EQ.1).OR.(chamb.EQ.3)).and.(scar_cell(i).EQ.1)) then !Scar_BZ vent 
            parco = 0.5D0
            sigma_f= 1.263609/(Cm*Chi);
         elseif (((chamb.EQ.1).OR.(chamb.EQ.3)).and.(scar_cell(i).EQ.2)) then !Scar_core vent
            parco = 0.25D0
            sigma_f= 1.263609/(Cm*Chi);
         elseif ((chamb.EQ.2)) then !atrio sinistro
            parco = scalaAtri
            sigma_f= 1.459190/(Cm*Chi);
         elseif ((chamb.EQ.4)) then !atrio destro
            parco = scalaAtri
            sigma_f= 1.603532/(Cm*Chi);
         else
            parco =scalaElse 
            sigma_f= 1.263609/(Cm*Chi); !un valore a caso, viene molt per zero
         endif
            sigma_n=sigma_f/7.5779D0
            sigma_s=sigma_f/7.5779D0
#ifdef CARTO
            sigma_f=1.0D0
            sigma_s=1.0D0
            sigma_n=1.0D0
#endif
            if (minf.GE.1) then
               v1=vert_of_cell_3d(1,i)
               v2=vert_of_cell_3d(2,i)
               v3=vert_of_cell_3d(3,i)
               v4=vert_of_cell_3d(4,i)   
               xcell=0.25D0*(xyz0_3d(1,v1)+xyz0_3d(1,v2)+xyz0_3d(1,v3)+xyz0_3d(1,v4))
               ycell=0.25D0*(xyz0_3d(2,v1)+xyz0_3d(2,v2)+xyz0_3d(2,v3)+xyz0_3d(2,v4))
               zcell=0.25D0*(xyz0_3d(3,v1)+xyz0_3d(3,v2)+xyz0_3d(3,v3)+xyz0_3d(3,v4))
               campo_inf=1.0D0 - exp( -( ( (xcell-x_inf_loc)**2+(ycell-y_inf_loc)**2+(zcell-z_inf_loc)**2)**4/size_inf**8 ) )
               parco=parco*campo_inf
            endif
         !-----------------------------------------------------------------------------------------
         ! Insert defect in conductivity (Infarcted Zone)
         !v1=vert_of_cell_3d(1,i)
         !v2=vert_of_cell_3d(2,i)
         !v3=vert_of_cell_3d(3,i)
         !v4=vert_of_cell_3d(4,i)

         
         !distPv1 = sqrt( (xyz0_3d(1,v1)-Px)**2 + (xyz0_3d(2,v1)-Py)**2 + (xyz0_3d(3,v1)-Pz)**2) + dr
         !distPv2 = sqrt( (xyz0_3d(1,v2)-Px)**2 + (xyz0_3d(2,v2)-Py)**2 + (xyz0_3d(3,v2)-Pz)**2) + dr
         !distPv3 = sqrt( (xyz0_3d(1,v3)-Px)**2 + (xyz0_3d(2,v3)-Py)**2 + (xyz0_3d(3,v3)-Pz)**2) + dr
         !distPv4 = sqrt( (xyz0_3d(1,v4)-Px)**2 + (xyz0_3d(2,v4)-Py)**2 + (xyz0_3d(3,v4)-Pz)**2) + dr
         
         !Mdefectv1 = sigma_f*0.1D0*((exp(alpha*(distPv1)) -1)/(exp(alpha*(distPv1))+1)+1)
         !Mdefectv2 = sigma_f*0.1D0*((exp(alpha*(distPv2)) -1)/(exp(alpha*(distPv2))+1)+1)
         !Mdefectv3 = sigma_f*0.1D0*((exp(alpha*(distPv3)) -1)/(exp(alpha*(distPv3))+1)+1)
         !Mdefectv4 = sigma_f*0.1D0*((exp(alpha*(distPv4)) -1)/(exp(alpha*(distPv4))+1)+1)
         !Mdefectv4
         !--------------------------------------------------------------------------------------------
         xfv =AmatrFibers_cell_3d(1,1,i) 
         yfv =AmatrFibers_cell_3d(2,1,i) 
         zfv =AmatrFibers_cell_3d(3,1,i) 
         xsv =AmatrFibers_cell_3d(1,2,i) 
         ysv =AmatrFibers_cell_3d(2,2,i) 
         zsv =AmatrFibers_cell_3d(3,2,i) 
         xnv =AmatrFibers_cell_3d(1,3,i) 
         ynv =AmatrFibers_cell_3d(2,3,i) 
         znv =AmatrFibers_cell_3d(3,3,i) 

         Mintcells_3d(1,1,i) = parco*(sigma_f*xfv**2  + sigma_n*xnv**2  + sigma_s*xsv**2)
         Mintcells_3d(1,2,i) = parco*(sigma_f*xfv*yfv + sigma_n*xnv*ynv + sigma_s*xsv*ysv)
         Mintcells_3d(1,3,i) = parco*(sigma_f*xfv*zfv + sigma_n*xnv*znv + sigma_s*xsv*zsv)
         Mintcells_3d(2,1,i) = parco*(sigma_f*xfv*yfv + sigma_n*xnv*ynv + sigma_s*xsv*ysv)
         Mintcells_3d(2,2,i) = parco*(sigma_f*yfv**2  + sigma_n*ynv**2  + sigma_s*ysv**2)
         Mintcells_3d(2,3,i) = parco*(sigma_f*yfv*zfv + sigma_n*ynv*znv + sigma_s*ysv*zsv)
         Mintcells_3d(3,1,i) = parco*(sigma_f*xfv*zfv + sigma_n*xnv*znv + sigma_s*xsv*zsv)
         Mintcells_3d(3,2,i) = parco*(sigma_f*yfv*zfv + sigma_n*ynv*znv + sigma_s*ysv*zsv)
         Mintcells_3d(3,3,i) = parco*(sigma_f*zfv**2  + sigma_n*znv**2  + sigma_s*zsv**2)

         ! Mintcells_3d(1,1,i) = sigma_f*xfv**2  + sigma_n*xnv**2  + sigma_s*xsv**2
         ! Mintcells_3d(1,2,i) = sigma_f*xfv*yfv + sigma_n*xnv*ynv + sigma_s*xsv*ysv
         ! Mintcells_3d(1,3,i) = sigma_f*xfv*zfv + sigma_n*xnv*znv + sigma_s*xsv*zsv
         ! Mintcells_3d(2,1,i) = sigma_f*xfv*yfv + sigma_n*xnv*ynv + sigma_s*xsv*ysv
         ! Mintcells_3d(2,2,i) = sigma_f*yfv**2  + sigma_n*ynv**2  + sigma_s*ysv**2
         ! Mintcells_3d(2,3,i) = sigma_f*yfv*zfv + sigma_n*ynv*znv + sigma_s*ysv*zsv
         ! Mintcells_3d(3,1,i) = sigma_f*xfv*zfv + sigma_n*xnv*znv + sigma_s*xsv*zsv
         ! Mintcells_3d(3,2,i) = sigma_f*yfv*zfv + sigma_n*ynv*znv + sigma_s*ysv*zsv
         ! Mintcells_3d(3,3,i) = sigma_f*zfv**2  + sigma_n*znv**2  + sigma_s*zsv**2
      enddo
#endif 
!@cuf   istat = cudaDeviceSynchronize !JDR TMP

      
!INTERPOLAZIONI VARIE TENSORE CONDUTTIVITA
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
    do i=fstart_3d(1),fend_3d(1)
       c1=cell_of_face_3d(1,i)
       c2=cell_of_face_3d(2,i)        
       g1interp=g1interpface_3d(i)
       if ((c1.NE.0).AND.(c2.NE.0)) then
          Mintfaces_3d(1,1,i)=g1interp*Mintcells_3d(1,1,c1)+(1-g1interp)*Mintcells_3d(1,1,c2)
          Mintfaces_3d(1,2,i)=g1interp*Mintcells_3d(1,2,c1)+(1-g1interp)*Mintcells_3d(1,2,c2)
          Mintfaces_3d(1,3,i)=g1interp*Mintcells_3d(1,3,c1)+(1-g1interp)*Mintcells_3d(1,3,c2)
          Mintfaces_3d(2,1,i)=g1interp*Mintcells_3d(2,1,c1)+(1-g1interp)*Mintcells_3d(2,1,c2)
          Mintfaces_3d(2,2,i)=g1interp*Mintcells_3d(2,2,c1)+(1-g1interp)*Mintcells_3d(2,2,c2)
          Mintfaces_3d(2,3,i)=g1interp*Mintcells_3d(2,3,c1)+(1-g1interp)*Mintcells_3d(2,3,c2)
          Mintfaces_3d(3,1,i)=g1interp*Mintcells_3d(3,1,c1)+(1-g1interp)*Mintcells_3d(3,1,c2)
          Mintfaces_3d(3,2,i)=g1interp*Mintcells_3d(3,2,c1)+(1-g1interp)*Mintcells_3d(3,2,c2)
          Mintfaces_3d(3,3,i)=g1interp*Mintcells_3d(3,3,c1)+(1-g1interp)*Mintcells_3d(3,3,c2)

          Mextfaces_3d(1,1,i)=g1interp*Mextcells_3d(1,1,c1)+(1-g1interp)*Mextcells_3d(1,1,c2)
          Mextfaces_3d(1,2,i)=g1interp*Mextcells_3d(1,2,c1)+(1-g1interp)*Mextcells_3d(1,2,c2)
          Mextfaces_3d(1,3,i)=g1interp*Mextcells_3d(1,3,c1)+(1-g1interp)*Mextcells_3d(1,3,c2)
          Mextfaces_3d(2,1,i)=g1interp*Mextcells_3d(2,1,c1)+(1-g1interp)*Mextcells_3d(2,1,c2)
          Mextfaces_3d(2,2,i)=g1interp*Mextcells_3d(2,2,c1)+(1-g1interp)*Mextcells_3d(2,2,c2)
          Mextfaces_3d(2,3,i)=g1interp*Mextcells_3d(2,3,c1)+(1-g1interp)*Mextcells_3d(2,3,c2)
          Mextfaces_3d(3,1,i)=g1interp*Mextcells_3d(3,1,c1)+(1-g1interp)*Mextcells_3d(3,1,c2)
          Mextfaces_3d(3,2,i)=g1interp*Mextcells_3d(3,2,c1)+(1-g1interp)*Mextcells_3d(3,2,c2)
          Mextfaces_3d(3,3,i)=g1interp*Mextcells_3d(3,3,c1)+(1-g1interp)*Mextcells_3d(3,3,c2)
       elseif ((c1.NE.0).AND.(c2.EQ.0)) then
          Mintfaces_3d(1,1,i)=Mintcells_3d(1,1,c1)
          Mintfaces_3d(1,2,i)=Mintcells_3d(1,2,c1)
          Mintfaces_3d(1,3,i)=Mintcells_3d(1,3,c1)
          Mintfaces_3d(2,1,i)=Mintcells_3d(2,1,c1)
          Mintfaces_3d(2,2,i)=Mintcells_3d(2,2,c1)
          Mintfaces_3d(2,3,i)=Mintcells_3d(2,3,c1)
          Mintfaces_3d(3,1,i)=Mintcells_3d(3,1,c1)
          Mintfaces_3d(3,2,i)=Mintcells_3d(3,2,c1)
          Mintfaces_3d(3,3,i)=Mintcells_3d(3,3,c1)

          Mextfaces_3d(1,1,i)=Mextcells_3d(1,1,c1)
          Mextfaces_3d(1,2,i)=Mextcells_3d(1,2,c1)
          Mextfaces_3d(1,3,i)=Mextcells_3d(1,3,c1)
          Mextfaces_3d(2,1,i)=Mextcells_3d(2,1,c1)
          Mextfaces_3d(2,2,i)=Mextcells_3d(2,2,c1)
          Mextfaces_3d(2,3,i)=Mextcells_3d(2,3,c1)
          Mextfaces_3d(3,1,i)=Mextcells_3d(3,1,c1)
          Mextfaces_3d(3,2,i)=Mextcells_3d(3,2,c1)
          Mextfaces_3d(3,3,i)=Mextcells_3d(3,3,c1)
       elseif ((c1.EQ.0).AND.(c2.NE.0)) then
          Mintfaces_3d(1,1,i)=Mintcells_3d(1,1,c2)
          Mintfaces_3d(1,2,i)=Mintcells_3d(1,2,c2)
          Mintfaces_3d(1,3,i)=Mintcells_3d(1,3,c2)
          Mintfaces_3d(2,1,i)=Mintcells_3d(2,1,c2)
          Mintfaces_3d(2,2,i)=Mintcells_3d(2,2,c2)
          Mintfaces_3d(2,3,i)=Mintcells_3d(2,3,c2)
          Mintfaces_3d(3,1,i)=Mintcells_3d(3,1,c2)
          Mintfaces_3d(3,2,i)=Mintcells_3d(3,2,c2)
          Mintfaces_3d(3,3,i)=Mintcells_3d(3,3,c2)

          Mextfaces_3d(1,1,i)=Mextcells_3d(1,1,c2)
          Mextfaces_3d(1,2,i)=Mextcells_3d(1,2,c2)
          Mextfaces_3d(1,3,i)=Mextcells_3d(1,3,c2)
          Mextfaces_3d(2,1,i)=Mextcells_3d(2,1,c2)
          Mextfaces_3d(2,2,i)=Mextcells_3d(2,2,c2)
          Mextfaces_3d(2,3,i)=Mextcells_3d(2,3,c2)
          Mextfaces_3d(3,1,i)=Mextcells_3d(3,1,c2)
          Mextfaces_3d(3,2,i)=Mextcells_3d(3,2,c2)
          Mextfaces_3d(3,3,i)=Mextcells_3d(3,3,c2)
       endif    
    enddo
!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_MintMext_cells
!------------------------------------------------------
!===================================================
      subroutine GradientEF_3d
      use constants
      use mls_param
      USE ieee_arithmetic
!@cuf   use cudafor

      implicit none

      real(DP):: rhspot,num,den,numext,dvvjj,potsur,potextsur,volcell
      real(DP):: distCF,g1interp,corrsten,bordo,potv,tempov
      real(DP):: xC1,yC1,zC1,xV1,yV1,zV1
      real(DP):: xgrad,ygrad,zgrad
      real(DP):: xgradext,ygradext,zgradext
      real(DP):: xvers,yvers,zvers
      real(DP):: xnormale,ynormale,znormale
      real(DP):: mint11,mint12,mint13,mint21,mint22,mint23,mint31,mint32,mint33
      real(DP):: mext11,mext12,mext13,mext21,mext22,mext23,mext31,mext32,mext33
      real(DP):: msum11,msum12,msum13,msum21,msum22,msum23,msum31,msum32,msum33
      real(DP), dimension(nvtot_3d) :: potEFnode_loc,potextEFnode_loc
      integer:: i,j,inp,chamb,c1,c2,v1,v2,v3,v4,f1,chambcheck
      integer:: vsi,vei,fsi,fei,esi,eei,csi,cei,controllo,ff,cc
      character*150 :: stri
!@cuf   integer :: istat
#ifdef USE_CUDA
       attributes(managed) :: potEFnode_loc,potextEFnode_loc
#endif 

!NOTAFV per considerare l'effetto delle piccole variazione geometriche durante la
!depolarizzazione dovremmo aggiornare (ag ogni time step o multiplo) anche le informazioni
!geometriche del dominio chiamando le routines che sono alla fine di read_geoSingleBody_3d
!tutto è predisposto, manca solo un modo per aggiornare/seguire l'orientamento delle fibre 


!!COMPUTE GRADIENT
!find nodal values
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
    do i=vstart_3d(1),vend_3d(1)
       chamb=vert_to_chamb_3d(i)
        num=0
#ifdef BIDOMAIN
        numext=0
#endif
        den=0
#ifdef ELEGEO0
        xV1 = xyz0_3d(1,i)
        yV1 = xyz0_3d(2,i)
        zV1 = xyz0_3d(3,i)
#else
        xV1 = xyz_3d(1,i)
        yV1 = xyz_3d(2,i)
        zV1 = xyz_3d(3,i)
#endif
        do j=1,n_cell_of_vert_3d(i)
            c1=cell_of_vert_3d(j,i)
        chambcheck=cell_to_chamb_3d(c1)
        if ( (chamb.EQ.chambcheck) &
       .OR. ((chamb.EQ.1).AND.(chambcheck.EQ.3)) .OR. ((chamb.EQ.3).AND.(chambcheck.EQ.1)) &
       .OR. ((chamb.EQ.2).AND.(chambcheck.EQ.4)) .OR. ((chamb.EQ.4).AND.(chambcheck.EQ.2))   )  then !MODIFICA  FV

            xC1 = cell_bar(1,c1)
            yC1 = cell_bar(2,c1)
            zC1 = cell_bar(3,c1)

            dvvjj = sqrt( (xV1-xC1)**2+(yV1-yC1)**2+(zV1-zC1)**2)
            num = num + (potEFcell_3d(c1) + 83.0D0)/(85.7D0*dvvjj)
            den = den + 1.D0/dvvjj
#ifdef BIDOMAIN
            numext = numext + potextEFcell_3d(c1)/dvvjj
#endif
        endif !chambcheck
        enddo
        if (den.NE.0.0D0) potEFnode_loc(i)=num/den
#ifdef BIDOMAIN
        if (den.NE.0.0D0) potextEFnode_loc(i)=numext/den
#endif
    enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!find face values
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
    do i=fstart_3d(1),fend_3d(1)
       v1=vert_of_face_3d(1,i)
       v2=vert_of_face_3d(2,i)
       v3=vert_of_face_3d(3,i)
       potEFface_3d(i)=(potEFnode_loc(v1)+potEFnode_loc(v2)+potEFnode_loc(v3))/3.D0;
#ifdef BIDOMAIN
       potextEFface_3d(i)=(potextEFnode_loc(v1)+potextEFnode_loc(v2)+potextEFnode_loc(v3))/3.D0;
#endif
    enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

! !find gradient on cells (gauss-green theorem)
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=cstart_3d(1),cend_3d(1)
        xgrad=0.0D0
        ygrad=0.0D0
        zgrad=0.0D0
#ifdef BIDOMAIN
        xgradext=0.0D0
        ygradext=0.0D0
        zgradext=0.0D0
#endif
        do j=1,4
           f1=face_of_cell_3d(j,i)
           xnormale=normalfaceofcells_3d(1,j,i)
           ynormale=normalfaceofcells_3d(2,j,i)
           znormale=normalfaceofcells_3d(3,j,i)

#ifdef ELEGEO0
           potsur = potEFface_3d(f1)*sur0_3d(f1)
#else
           potsur = potEFface_3d(f1)*sur_3d(f1)
#endif
           xgrad = xgrad + potsur*xnormale
           ygrad = ygrad + potsur*ynormale
           zgrad = zgrad + potsur*znormale
#ifdef BIDOMAIN
#ifdef ELEGEO0
           potextsur = potextEFface_3d(f1)*sur0_3d(f1)
#else
           potextsur = potextEFface_3d(f1)*sur_3d(f1)
#endif
           xgradext = xgradext + potextsur*xnormale
           ygradext = ygradext + potextsur*ynormale
           zgradext = zgradext + potextsur*znormale
#endif
        enddo
#ifdef ELEGEO0
               volcell = vol0_3d(i)
#else
               volcell = vol_3d(i)
#endif
        gradcell_3d(1,i) =  xgrad/volcell
        gradcell_3d(2,i) =  ygrad/volcell
        gradcell_3d(3,i) =  zgrad/volcell
#ifdef BIDOMAIN
        gradextcell_3d(1,i) =  xgradext/volcell
        gradextcell_3d(2,i) =  ygradext/volcell
        gradextcell_3d(3,i) =  zgradext/volcell
#endif

     enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP       


! !clear gradients on boundaries 
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=cstart_3d(1),cend_3d(1)
         chamb=cell_to_chamb_3d(i)
         
         controllo=0
         do j=1,4
            ff=face_of_cell_3d(j,i)
            c1=cell_of_face_3d(1,ff)
            c2=cell_of_face_3d(2,ff)
            if (c1.EQ.i) then
               cc=c2
            elseif (c2.EQ.i) then
               cc=c1
            else
               write(*,*) "error gradient ECG"
               stop
            endif

            if (cc.NE.0) then !check cc not zero !!
               chambcheck=cell_to_chamb_3d(cc)
            else
               chambcheck=-1
            endif

                    if ( (chamb.EQ.chambcheck) &
                         .OR. ((chamb.EQ.1).AND.(chambcheck.EQ.3)) .OR. ((chamb.EQ.3).AND.(chambcheck.EQ.1)) &
                         .OR. ((chamb.EQ.2).AND.(chambcheck.EQ.4)) .OR. ((chamb.EQ.4).AND.(chambcheck.EQ.2))   )  then !MODIFICA  FV
                       
                       controllo=1
                       
                    endif
                 
         enddo

         if (controllo.EQ.0) then
            gradcell_3d(1,i) =  0.0d0
            gradcell_3d(2,i) =  0.0d0
            gradcell_3d(3,i) =  0.0d0
#ifdef BIDOMAIN
            gradextcell_3d(1,i) =  0.0d0
            gradextcell_3d(2,i) =  0.0d0
            gradextcell_3d(3,i) =  0.0d0
#endif            
         endif


     enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP       

      
      return
      end subroutine GradientEF_3d
!===================================================
!===================================================
      subroutine ElectroTensionEF_3d(timex)
      use constants
      use param, only: time
      use mls_param
!@cuf   use cudafor

      implicit none

      real(DP):: timex,timeEFv,alpv1,tempov,tempovBeat
      real(DP):: ksig_LV,ksig_LA,ksig_RV,ksig_RA,ksig
      real(DP):: eps0_LV,eps0_LA,eps0_RV,eps0_RA,eps0
      real(DP):: timeMS,dtMS
      real(DP):: c1minpot,num,den,xV1,yV1,zV1,xC1,yC1,zC1
      real(DP):: dvvjj,astressold,astressnew
      real(DP) :: x_inf_loc,y_inf_loc,z_inf_loc,size_inf
      real(DP) :: campo_inf
      integer:: i,indv1,indv2,act,inp,chamb,v1,v2,j,c1
      character*150 :: stri
!@cuf   integer :: istat


      timeMS=timex*TSTAR*1000.d0
      dtMS=dt*TSTAR*1000.d0*real(stepEF_3d)
          
      ksig_LV=0.001  !stress(Pa) adimensionalizzati
      ksig_LA=0.001
      
      ksig_LV=ksig_LV*(1/0.08)*20.d0
      ksig_LA=ksig_LA*(5/0.04)*20.d0

      ksig_LV=ksig_LV*10.d0
      ! ksig_LA=ksig_LA*20.d0
      ksig_LA=ksig_LA*2.d0      
      
      ksig_RV=ksig_LV*1.5d0
      ksig_RA=ksig_LA
      
      eps0_LV =0.01; !%mV^−1
      eps0_LA =0.02; !%mV^−1
      eps0_RV =eps0_LV; !%mV^−1
      eps0_RA =eps0_LA; !%mV^−1


#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do i=1,nctot_3d
         chamb=cell_to_chamb_3d(i)
         if (chamb.EQ.1) then
            c1minpot=minpotLV
            ksig=ksig_LV*CARTO_Dcell3d(i)/0.18D0
            eps0=eps0_LV
         elseif (chamb.EQ.2) then
            c1minpot=minpotLA
            ksig=ksig_LA
            eps0=eps0_LA
         elseif (chamb.EQ.3) then
            c1minpot=minpotRV
            ksig=ksig_RV
            eps0=eps0_RV
         elseif (chamb.EQ.4) then
            c1minpot=minpotRA
            ksig=ksig_RA
            eps0=eps0_RA
         else
            ksig=0 !this should be zero to keep the astress null
            eps0=0.01 !let set a value
            c1minpot=minpotLV !let set a value
         endif
         astressold=astressEFcell_3d(i)
         astressnew=1.D0/(1.D0+dtMS*eps0)*(astressold+dtMS*eps0*ksig*(potEFcell_3d(i)-c1minpot));
         astressEFcell_3d(i)=astressnew
      enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

      if (minf.GE.1) then
         size_inf=xyz_inf(4)
         x_inf_loc=xyz_inf(1)
         y_inf_loc=xyz_inf(2)
         z_inf_loc=xyz_inf(3)
      endif
!find nodal values
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=vstart_3d(1),vend_3d(1)         
       chamb=vert_to_chamb_3d(i)
       num=0
       den=0
#ifdef ELEGEO0
       xV1 = xyz0_3d(1,i)
       yV1 = xyz0_3d(2,i)
       zV1 = xyz0_3d(3,i)
#else
       xV1 = xyz_3d(1,i)
       yV1 = xyz_3d(2,i)
       zV1 = xyz_3d(3,i)
#endif
        do j=1,n_cell_of_vert_3d(i)
           c1=cell_of_vert_3d(j,i)
           xC1 = cell_bar(1,c1)
           yC1 = cell_bar(2,c1)
           zC1 = cell_bar(3,c1)

           dvvjj = sqrt( (xV1-xC1)**2+(yV1-yC1)**2+(zV1-zC1)**2)
           num = num + astressEFcell_3d(c1)/dvvjj
           den = den + 1.D0/dvvjj
        enddo
        if (den.NE.0.0D0) astressEFnode_3d(i)=num/den

        if (minf.GE.1) then   !for infarction, TODO ADD IZ from scaray
           campo_inf=1.0D0 - exp( -( ( (xyz0_3d(1,i)-x_inf_loc)**2+(xyz0_3d(2,i)-y_inf_loc)**2+(xyz0_3d(3,i)-z_inf_loc)**2)**4/size_inf**8 ) )                 
           astressEFnode_3d(i)=astressEFnode_3d(i)*campo_inf
        endif
     enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!find face values
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
    do i=1,netot_3d
       v1=vert_of_edge_3d(1,i)
       v2=vert_of_edge_3d(2,i)
       astressEFedge_3d(i)=0.5D0*(astressEFnode_3d(v1)+astressEFnode_3d(v2))
    enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP
    

      return
      end subroutine ElectroTensionEF_3d
