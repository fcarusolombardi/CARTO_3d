!===================================================
      subroutine ElectroRunEF_1d(timex)
      use constants
      use mls_param
!@cuf   use cudafor

      implicit none

      real(DP):: timex,timeEFv,alpv1,timeMS,dtMS
      real(DP):: timeVentr,timeAtrio,celld,poti,dedge
      real(DP):: rhspot,num,den,dvvjj,potsur,volcell
      real(DP):: distCF,g1interp,corrsten,bordo,potv,tempiv
      real(DP):: xB1,yB1,zB1,xV1,yV1,zV1
      real(DP):: xgrad,ygrad,zgrad
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
      real(DP):: const1s,const2s,const3s,const4s,const5s,const6s,const7s,const8s,const9s
      real(DP):: const10s,const11s,const12s,const13s,const14s,const15s,const16s,const17s,const18s,const19s
      real(DP):: const20s,const21s,const22s,const23s,const24s,const25s,const26s,const27s,const28s,const29s
      real(DP):: const30s,const31s,const32s,const33s,const34s,const35s,const36s,const37s,const38s,const39s
      real(DP):: const40s,const41s,const42s,const43s,const44s,const45s,const46s,const47s,const48s,const49s
      real(DP):: const50s,const51s,const52s

      !Minimal Model
      ! -----------------------------------------------------
      real(DP):: const4,const5,const6,const7,const8,const9
      real(DP):: const10,const11,const12,const13,const14,const15,const16,const17,const18,const19
      real(DP):: const20,const21,const22,const23,const24,const25,const26,const27,const28,const29
      real(DP):: const30,const31,const32,app1,app2,app3,app4,app5

      ! ----------------------------------------------------- 

      real(DP):: aoo, duration_signal,sur_face,ftempo,tempov
      integer:: i,j,inp,chamb,f1,f2,v1,v2,v3,e1,vgc
      integer:: vmaster,vslave,fmaster,fslave,cslave
      character*150 :: stri
!@cuf   integer :: istat
!#ifdef USE_CUDA
!        attributes(managed) :: gradcell_3d
!#endif 

!NOTAFV per considerare l'effetto delle piccole variazione geometriche durante la
!depolarizzazione dovremmo aggiornare (ag ogni time step o multiplo) anche le informazioni
!geometriche del dominio chiamando le routines che sono alla fine di read_geoSingleBody_3d
!tutto è predisposto, manca solo un modo per aggiornare/seguire l'orientamento delle fibre 

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
!    Gks = Gks*3.D0 !ricambia
!    Gks = 0.392;  !ENDO
!    Gks = 0.098;  !else, myocardium?
    !% //Parameters for Ik1
    GK1 = 5.405;
    !% //Parameters for Ito
    Gto = 0.294;
!    Gto = Gto*0.5D0 !ricambia
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
!    KmCa = KmCa*3.D0 !ricambia
    ksat = 0.1;
    n = 0.35;
    !% //Parameters for IpCa
    GpCa = 0.1238;
    KpCa = 0.0005;
    !% //Parameters for IpK;
    GpK = 0.0146;


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


!Stewart 2009
    CONST1s = 8314.472; !%R in component membrane (joule_per_mole_kelvin)
    CONST2s = 310; !%T in component membrane (kelvin)
    CONST3s = 96485.3415; !%F in component membrane (coulomb_per_millimole)
    CONST4s = 0.185; !%Cm in component membrane (microF)
    CONST5s = 0.016404; !%V_c in component membrane (micrometre3)
    CONST6s = 0.03; !%P_kna in component reversal_potentials (dimensionless)
    CONST7s = 5.4; !%K_o in component potassium_dynamics (millimolar)
    CONST8s = 140; !%Na_o in component sodium_dynamics (millimolar)
    CONST9s = 2; !%Ca_o in component calcium_dynamics (millimolar)
    CONST10s = 0.0145654; !%g_f_Na in component hyperpolarization_activated_current (nanoS_per_picoF)
    CONST11s = 0.0234346; !%g_f_K in component hyperpolarization_activated_current (nanoS_per_picoF)
    CONST12s = 0.065; !%g_K1 in component inward_rectifier_potassium_current (nanoS_per_picoF)
    CONST13s = 0.0918; !%g_Kr in component rapid_time_dependent_potassium_current (nanoS_per_picoF)
    CONST14s = 0.2352; !%g_Kr in component rapid_time_dependent_potassium_current (nanoS_per_picoF)
    CONST15s = 130.5744; !%g_Na in component fast_sodium_current (nanoS_per_picoF)
    CONST16s = 0.00029; !%g_bna in component sodium_background_current (nanoS_per_picoF)
    CONST17s = 3.98e-5; !%g_CaL in component L_type_Ca_current (litre_per_farad_second) 
    CONST18s = 0.000592;!%g_bca in component calcium_background_current (nanoS_per_picoF)
    CONST19s = 0.08184; !%g_to in component transient_outward_current (nanoS_per_picoF)
    CONST20s = 0.0227; !%g_sus in component sustained_outward_current (nanoS_per_picoF)
    CONST21s = 2.724; !%P_NaK in component sodium_potassium_pump_current (picoA_per_picoF)
    CONST22s = 1; !%K_mk in component sodium_potassium_pump_current (millimolar)
    CONST23s = 40; !%K_mNa in component sodium_potassium_pump_current (millimolar)
    CONST24s = 1000; !%K_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF)
    CONST25s = 0.1; !%K_sat in component sodium_calcium_exchanger_current (dimensionless)
    CONST26s = 2.5; !%alpha in component sodium_calcium_exchanger_current (dimensionless)
    CONST27s = 0.35; !%gamma in component sodium_calcium_exchanger_current (dimensionless)
    CONST28s = 1.38; !%Km_Ca in component sodium_calcium_exchanger_current (millimolar)
    CONST29s = 87.5; !%Km_Nai in component sodium_calcium_exchanger_current (millimolar)
    CONST30s = 0.1238;! %g_pCa in component calcium_pump_current (picoA_per_picoF)
    CONST31s = 0.0005; !%K_pCa in component calcium_pump_current (millimolar)
    CONST32s = 0.0146; !%g_pK in component potassium_pump_current (nanoS_per_picoF)
    CONST33s = 0.15; !%k1_prime in component calcium_dynamics (per_millimolar2_per_millisecond)
    CONST34s = 0.045; !%k2_prime in component calcium_dynamics (per_millimolar_per_millisecond)
    CONST35s = 0.06; !%k3 in component calcium_dynamics (per_millisecond)
    CONST36s = 0.005; !%k4 in component calcium_dynamics (per_millisecond)
    CONST37s = 1.5; !%EC in component calcium_dynamics (millimolar)
    CONST38s = 2.5; !%max_sr in component calcium_dynamics (dimensionless)
    CONST39s = 1; !%min_sr in component calcium_dynamics (dimensionless)
    CONST40s = 0.102; !%V_rel in component calcium_dynamics (per_millisecond)
    CONST41s = 0.0038; !%V_xfer in component calcium_dynamics (per_millisecond)
    CONST42s = 0.00025; !%K_up in component calcium_dynamics (millimolar)
    CONST43s = 0.00036; !%V_leak in component calcium_dynamics (per_millisecond)
    CONST44s = 0.006375; !%Vmax_up in component calcium_dynamics (millimolar_per_millisecond)
    CONST45s = 0.2; !%Buf_c in component calcium_dynamics (millimolar)
    CONST46s = 0.001; !%K_buf_c in component calcium_dynamics (millimolar)
    CONST47s = 10; !%Buf_sr in component calcium_dynamics (millimolar)
    CONST48s = 0.3; !%K_buf_sr in component calcium_dynamics (millimolar)
    CONST49s = 0.4; !%Buf_ss in component calcium_dynamics (millimolar)
    CONST50s = 0.00025; !%K_buf_ss in component calcium_dynamics (millimolar
    CONST51s = 0.001094; !%V_sr in component calcium_dynamics (micrometre3)
    CONST52s = 5.468e-5; !%V_ss in component calcium_dynamics (micrometre3)




    timeMS=timex*TSTAR*1000.d0
    dtMS=dt*TSTAR*1000.d0/real(stepEF_1d)       !CONTROLLA

!!----->PARTE SPAZIALE
! !!STIMOLO
    duration_signal = 5 !ms ricambia
!    duration_signal = 10 !ms ricambia


!CELL MODEL
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=1,nvtotEF_1d 
        inp = vert_to_partEF_1d(i)

!stimulus
        if ((timeMS.GT.0).AND.(timeMS.LT.duration_signal)) then
           Istim = IstimEF_1d(i)
        else
           Istim = 0.0
        endif
!spatial term
        bordo = 0.D0
        celld = 0.D0
        poti = potEFnode_1d(i)

        do j=1,n_edge_of_vertEF_1d(i)  
           e1=edge_of_vertEF_1d(j,i)
#ifdef ELEGEO0
           dedge=dist0EF_1d(e1)
#else           
           dedge=distEF_1d(e1) !FV ricambia
#endif
           celld=celld+dedge/2.D0;
           v1=vert_of_edgeEF_1d(1,e1)
           v2=vert_of_edgeEF_1d(2,e1) 
           if ((v1.NE.i).AND.(v2.EQ.i)) then
               vgc=v1
            elseif ((v2.NE.i).AND.(v1.EQ.i)) then
               vgc=v2
            else
               write(*,*) "error EF1d"
               stop
            endif

           bordo = bordo + MintedgesEF_1d(e1)*(potEFnode_1d(vgc)-poti)/dedge
        enddo
        rhspot = bordo/celld
        rhspot = rhspot/((1000.D0*LSTAR)**2) !CONTROLLA devo dimens lo spazio in mm                   
        ! rhspot = rhspot/((LSTAR)**2) !devo dimens lo spazio in mm                            
        if (inp.EQ.1) then
#ifdef COURTEMANCHE
           STA1 = XEF_1d(1,i); !Volt
           STA2 = XEF_1d(2,i);
           STA3 = XEF_1d(3,i);
           STA4 = XEF_1d(4,i);
           STA5 = XEF_1d(5,i);
           STA6 = XEF_1d(6,i);
           STA7 = XEF_1d(7,i);
           STA8 = XEF_1d(8,i);
           STA9 = XEF_1d(9,i);
           STA10 = XEF_1d(10,i);
           STA11 = XEF_1d(11,i);
           STA12 = XEF_1d(12,i);
           STA13 = XEF_1d(13,i);
           STA14 = XEF_1d(14,i);
           STA15 = XEF_1d(15,i);
           STA16 = XEF_1d(16,i);
           STA17 = XEF_1d(17,i);
           STA18 = XEF_1d(18,i);
           STA19 = XEF_1d(19,i);
           STA20 = XEF_1d(20,i);
           STA21 = XEF_1d(21,i);

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
           RAT1 =  - (ALG31+ALG51+ALG52+ALG54+ALG55+ALG56+ALG62+ALG63+ALG59+ALG65+ALG64+ALG57)/CONST4c ;
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


           XEF_1d(1,i) = STA1 + dtMS*(RAT1+rhspot+Istim);
           XEF_1d(2,i) = STA2 + dtMS*RAT2
           XEF_1d(3,i) = STA3np1
           XEF_1d(4,i) = STA4np1
           XEF_1d(5,i) = STA5np1
           XEF_1d(6,i) = STA6 + dtMS*RAT6
           XEF_1d(7,i) = STA7np1
           XEF_1d(8,i) = STA8np1
           XEF_1d(9,i) = STA9np1
           XEF_1d(10,i) = STA10np1
           XEF_1d(11,i) = STA11np1
           XEF_1d(12,i) = STA12np1
           XEF_1d(13,i) = STA13 + dtMS*RAT13
           XEF_1d(14,i) = STA14np1
           XEF_1d(15,i) = STA15np1
           XEF_1d(16,i) = STA16np1
           XEF_1d(17,i) = STA17 + dtMS*RAT17
           XEF_1d(18,i) = STA18np1
           XEF_1d(19,i) = STA19np1
           XEF_1d(20,i) = STA20np1
           XEF_1d(21,i) = STA21 + dtMS*RAT21
#endif !Courtemanche
            
#ifdef MINIMAL_MODEL !ATRIA !TODO :-> Parametrizations for different bundles?
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
      

           STA1 = XEF_1d(1,i); !nondimensional here
           STA2 = XEF_1d(2,i);
           STA3 = XEF_1d(3,i);
           STA4 = XEF_1d(4,i);
           
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
           
           XEF_1d(1,i) = STA1 + dtMS*(RAT1+rhspot/(CONST5-CONST4)); 
           XEF_1d(2,i) = STA2np1
           XEF_1d(3,i) = STA3np1
           XEF_1d(4,i) = STA4np1
           
           potEFnode_1d(i) = CONST4+ (CONST5-CONST4)*XEF_1d(1,i)
#endif !minimal model atria
        elseif (inp.EQ.2) then
#ifdef COURTEMANCHE
              Volt = XEF_1d(1,i);
              Cai  = XEF_1d(2,i);
              CaSR = XEF_1d(3,i);
              CaSS = XEF_1d(4,i);
              Nai  = XEF_1d(5,i);
              Ki   = XEF_1d(6,i);
              INa_m = XEF_1d(7,i);
              INa_h = XEF_1d(8,i);
              INa_j = XEF_1d(9,i);
              IKr_xr1 = XEF_1d(10,i);
              IKr_xr2 = XEF_1d(11,i);
              IKs_xs  = XEF_1d(12,i);
              Ito_r   = XEF_1d(13,i);
              Ito_s   = XEF_1d(14,i);
              ICaL_d  = XEF_1d(15,i);
              ICaL_f  = XEF_1d(16,i);
              ICaL_f2 = XEF_1d(17,i);
              ICaL_fCaSS = XEF_1d(18,i);
              RR         = XEF_1d(19,i);

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
        dv = -(IKr+IKs+IK1+Ito+INa+IbNa+INaK+ICaL+IbCa+INaCa+IpCa+IpK)
    
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
        
        dKi = -(Istim + IK1 + Ito + IKr + IKs - 2 * INaK + IpK)*inverseVcF*CAPACITANCE; !ricambia
        
!AGGIORNA        
        XEF_1d(1,i) = XEF_1d(1,i) + dtMS*(dv+rhspot+Istim);
        XEF_1d(2,i) = XEF_1d(2,i) + dtMS*dCai;
        XEF_1d(3,i) = XEF_1d(3,i) + dtMS*dCaSR;
        XEF_1d(4,i) = XEF_1d(4,i) + dtMS*dCaSS;
        XEF_1d(5,i) = XEF_1d(5,i) + dtMS*dNai;
        XEF_1d(6,i) = XEF_1d(6,i) + dtMS*dKi;
        XEF_1d(7,i) = INa_mP1; !%gate4
        XEF_1d(8,i) = INa_hP1; !%gate5
        XEF_1d(9,i) = INa_jP1; !%gate6
        XEF_1d(10,i) = IKr_xr1P1; !%gate1
        XEF_1d(11,i) = IKr_xr2P1; !%gate2
        XEF_1d(12,i) = IKs_xsP1;  !%gate3
        XEF_1d(13,i) = Ito_rP1;   !%gate7
        XEF_1d(14,i) = Ito_sP1;   !%gate8
        XEF_1d(15,i) = CaL_dP1;   !%gate9
        XEF_1d(16,i) = CaL_fP1;   !%gate10
        XEF_1d(17,i) = CaL_f2P1;  !%gate11
        XEF_1d(18,i) = CaL_fCaSSP1; !%gate12
        XEF_1d(19,i) = RRP1; !%gate13
#endif!Courtemanche
    
#ifdef MINIMAL_MODEL !Atria, here ENDO parametrization

        !minimal model endo, Set 4V-3 (Fenton Science Adv.) 
        CONST6 = 0.218655D0;
        CONST7 = 0.15D0;
        CONST8 = 7.188259D0;
        CONST9 = 2.031406D0;
        CONST10 = 2.545600D0;
        CONST11 = 0.856429D0;
        
        CONST12 = 13.090045D0;
        CONST13 = 573.569884D0;
        CONST14 = 0.044570D0
        CONST15 = 284.169705D0;
        CONST16 = 0.007950D0;
        CONST17 = 0.168277D0;
        CONST18 = 1.451815D0;
        CONST19 = 28.726329D0;
        CONST20 = 0.281850D0;
        CONST21 = 1.937860D0;
        CONST22 = 0.682362D0;
        CONST23 = 2.347223D0;
        CONST24 = 0.049502D0;
        CONST25 = 0.651285D0;
        CONST26 = 19.751628D0;
        CONST27 = 101.167275D0;
        CONST28 = 57.162563D0;
        CONST29 = 0.028506D0;
        CONST30 = 213.296779D0;
        CONST31 = 19.147136D0;
        CONST32 = 9.999739D0;
         
      

           STA1 = XEF_1d(1,i); !nondimensional here
           STA2 = XEF_1d(2,i);
           STA3 = XEF_1d(3,i);
           STA4 = XEF_1d(4,i);
           
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
           
           
           
           XEF_1d(1,i) = STA1 + dtMS*(RAT1+rhspot/(CONST5-CONST4)); 
           XEF_1d(2,i) = STA2np1
           XEF_1d(3,i) = STA3np1
           XEF_1d(4,i) = STA4np1
           
           potEFnode_1d(i) = CONST4+ (CONST5-CONST4)*XEF_1d(1,i)
        
#endif !minimal model atria
     elseif (inp.EQ.3) then !Here Purkinjie parametrization 
#ifdef COURTEMANCHE
           STA1 = XEF_1d(1,i); !Volt                                      
           STA2 = XEF_1d(2,i);
           STA3 = XEF_1d(3,i);
           STA4 = XEF_1d(4,i);
           STA5 = XEF_1d(5,i);
           STA6 = XEF_1d(6,i);
           STA7 = XEF_1d(7,i);
           STA8 = XEF_1d(8,i);
           STA9 = XEF_1d(9,i);
           STA10 = XEF_1d(10,i);
           STA11 = XEF_1d(11,i);
           STA12 = XEF_1d(12,i);
           STA13 = XEF_1d(13,i);
           STA14 = XEF_1d(14,i);
           STA15 = XEF_1d(15,i);
           STA16 = XEF_1d(16,i);
           STA17 = XEF_1d(17,i);
           STA18 = XEF_1d(18,i);
           STA19 = XEF_1d(19,i);
           STA20 = XEF_1d(20,i);

           
    ALG9 = 1.00000/(1.00000+exp((STA1+20.0000)/7.00000));
    ALG23 =  1102.50*exp( - ((STA1+27.0000)**2.00000)/225.000)+200.000/(1.00000+exp((13.0000 - STA1)/10.0000))+180.000/(1.00000+exp((STA1+30.0000)/10.0000))+20.0000;
    STA14np1 = ALG9 + (STA14-ALG9)*exp(-dtMS/ALG23);  
    ALG10 = 0.670000/(1.00000+exp((STA1+35.0000)/7.00000))+0.330000;
    ALG24 =  562.000*exp( - ((STA1+27.0000)**2.00000)/240.000)+31.0000/(1.00000+exp((25.0000 - STA1)/10.0000))+80.0000/(1.00000+exp((STA1+30.0000)/10.0000));
    STA15np1 = ALG10 + (STA15-ALG10)*exp(-dtMS/ALG24);  
    ALG11 = 0.600000/(1.00000+((STA12/0.0500000)**2.00000))+0.400000;
    ALG25 = 80.0000/(1.00000+((STA12/0.0500000)**2.00000))+2.00000;
    STA16np1 = ALG11 + (STA16-ALG11)*exp(-dtMS/ALG25);
    ALG12 = 1.00000/(1.00000+exp((STA1+27.0000)/13.0000));
    ALG26 =  85.0000*exp( - ((STA1+25.0000)**2.00000)/320.000)+5.00000/(1.00000+exp((STA1 - 40.0000)/5.00000))+42.0000;
    STA17np1 = ALG12 + (STA17-ALG12)*exp(-dtMS/ALG26);
    ALG13 = 1.00000/(1.00000+exp((20.0000 - STA1)/13.0000));
    ALG27 =  10.4500*exp( - ((STA1+40.0000)**2.00000)/1800.00)+7.30000;
    STA18np1 = ALG13 + (STA18-ALG13)*exp(-dtMS/ALG27);
    ALG1 = 1.00000/(1.00000+exp((STA1+80.6000)/6.80000));
    ALG15 =  1.00000*exp( - 2.90000 -  0.0400000*STA1);
    ALG29 =  1.00000*exp(3.60000+ 0.110000*STA1);
    ALG38 = 4000.00/(ALG15+ALG29);
    STA5np1 = ALG1 + (STA5-ALG1)*exp(-dtMS/ALG38);
    ALG2 = 1.00000/(1.00000+exp(( - 26.0000 - STA1)/7.00000));
    ALG16 = 450.000/(1.00000+exp(( - 45.0000 - STA1)/10.0000));
    ALG30 = 6.00000/(1.00000+exp((STA1+30.0000)/11.5000));
    ALG39 =  1.00000*ALG16*ALG30;
    STA6np1 = ALG2 + (STA6-ALG2)*exp(-dtMS/ALG39);
    ALG3 = 1.00000/(1.00000+exp((STA1+88.0000)/24.0000));
    ALG17 = 3.00000/(1.00000+exp(( - 60.0000 - STA1)/20.0000));
    ALG31 = 1.12000/(1.00000+exp((STA1 - 60.0000)/20.0000));
    ALG40 =  1.00000*ALG17*ALG31;
    STA7np1 = ALG3 + (STA7-ALG3)*exp(-dtMS/ALG40);
    ALG4 = 1.00000/(1.00000+exp(( - 5.00000 - STA1)/14.0000));
    ALG18 = 1400.00/sqrt((1.00000+exp((5.00000 - STA1)/6.00000)));
    ALG32 = 1.00000/(1.00000+exp((STA1 - 35.0000)/15.0000));
    ALG41 =  1.00000*ALG18*ALG32+80.0000;
    STA8np1 = ALG4 + (STA8-ALG4)*exp(-dtMS/ALG41);
    ALG5 = 1.00000/((1.00000+exp(( - 56.8600 - STA1)/9.03000))**2.00000);
    ALG19 = 1.00000/(1.00000+exp(( - 60.0000 - STA1)/5.00000));
    ALG33 = 0.100000/(1.00000+exp((STA1+35.0000)/5.00000))+0.100000/(1.00000+exp((STA1 - 50.0000)/200.000));
    ALG42 =  1.00000*ALG19*ALG33;
    STA9np1 = ALG5 + (STA9-ALG5)*exp(-dtMS/ALG42);
    ALG6 = 1.00000/( (1.00000+exp((STA1+71.5500)/7.43000))**2.00000);
    if (STA1.LT.- 40.000) then
        ALG20 = 0.0570000*exp( - (STA1+80.0000)/6.80000);
        ALG34 = 2.70000*exp( 0.0790000*STA1)+ 310000*exp( 0.348500*STA1);
        
        ALG21 =(( (  - 25428.0*exp( 0.244400*STA1) -  6.94800e-06*exp(  - 0.0439100*STA1))*(STA1+37.7800))/1.00000)/(1.00000+exp( 0.311000*(STA1+79.2300)));
        ALG35 = ( 0.0242400*exp(  - 0.0105200*STA1))/(1.00000+exp(  - 0.137800*(STA1+40.1400)));
    else
        ALG20 = 0.0;
        ALG34 =0.770000/( 0.130000*(1.00000+exp((STA1+10.6600)/ - 11.1000)));
        
        ALG21 =0.0;
        ALG35 = ( 0.600000*exp( 0.0570000*STA1))/(1.00000+exp(  - 0.100000*(STA1+32.0000)));
    endif     
    ALG43 = 1.00000/(ALG20+ALG34);
    STA10np1 = ALG6 + (STA10-ALG6)*exp(-dtMS/ALG43);
    ALG7 = 1.00000/((1.00000+exp((STA1+71.5500)/7.43000))**2.00000);
    ALG44 = 1.00000/(ALG21+ALG35);
    STA11np1 = ALG7 + (STA11-ALG7)*exp(-dtMS/ALG44);
    ALG8 = 1.00000/(1.00000+exp(( - 8.00000 - STA1)/7.50000));
    ALG22 = 1.40000/(1.00000+exp(( - 35.0000 - STA1)/13.0000))+0.250000;
    ALG36 = 1.40000/(1.00000+exp((STA1+5.00000)/5.00000));
    ALG45 = 1.00000/(1.00000+exp((50.0000 - STA1)/20.0000));
    ALG47 =  1.00000*ALG22*ALG36+ALG45;
    STA13np1 = ALG8 + (STA13-ALG8)*exp(-dtMS/ALG47);
    ALG62 = (( (( CONST21s*CONST7s)/(CONST7s+CONST22s))*STA3)/(STA3+CONST23s))/(1.00000+ 0.124500*exp((  - 0.100000*STA1*CONST3s)/( CONST1s*CONST2s))+ 0.0353000*exp((  - STA1*CONST3s)/( CONST1s*CONST2s)));
    ALG14 =  (( CONST1s*CONST2s)/CONST3s)*log(CONST8s/STA3);
    ALG55 =  CONST15s*(STA9**3.00000)*STA10*STA11*(STA1 - ALG14);
    ALG56 =  CONST16s*(STA1 - ALG14);
    ALG63 = ( CONST24s*( exp(( CONST27s*STA1*CONST3s)/( CONST1s*CONST2s))*(STA3**3.00000)*CONST9s -  exp(( (CONST27s - 1.00000)*STA1*CONST3s)/( CONST1s*CONST2s))*(CONST8s**3.00000)*STA4*CONST26s))/( ((CONST29s**3.00000)+(CONST8s**3.00000))*(CONST28s+CONST9s)*(1.00000+ CONST25s*exp(( (CONST27s - 1.00000)*STA1*CONST3s)/( CONST1s*CONST2s))));
    ALG48 =  STA5*CONST10s*(STA1 - ALG14);
    RAT3 =  ((  - 1.00000*(ALG55+ALG56+ALG48+ 3.00000*ALG62+ 3.00000*ALG63))/( 1.00000*CONST5s*CONST3s))*CONST4s;
    ALG28 =  (( CONST1s*CONST2s)/CONST3s)*log(CONST7s/STA2);
    ALG51 = 1.00000/(1.00000+exp( 0.100000*(STA1+75.4400)));
    ALG52 =  CONST12s*ALG51*((STA1 - 8.00000) - ALG28);
    ALG59 =  CONST19s*STA18*STA17*(STA1 - ALG28);
    ALG60 = 1.00000/(1.00000+exp((5.00000 - STA1)/17.0000));
    ALG61 =  CONST20s*ALG60*(STA1 - ALG28);
    ALG53 =  CONST13s* sqrt(CONST7s/5.40000)*STA6*STA7*(STA1 - ALG28);
    ALG37 =  (( CONST1s*CONST2s)/CONST3s)*log((CONST7s+ CONST6s*CONST8s)/(STA2+ CONST6s*STA3));
    ALG54 =  CONST14s*(STA8**2.00000)*(STA1 - ALG37);
    ALG57 = ( (( CONST17s*STA13*STA14*STA15*STA16*4.00000*(STA1 - 15.0000)*(CONST3s**2.00000))/( CONST1s*CONST2s))*( 0.250000*STA12*exp(( 2.00000*(STA1 - 15.0000)*CONST3s)/( CONST1s*CONST2s)) - CONST9s))/(exp(( 2.00000*(STA1 - 15.0000)*CONST3s)/( CONST1s*CONST2s)) - 1.00000);
    ALG46 =  (( 0.500000*CONST1s*CONST2s)/CONST3s)*log(CONST9s/STA4);
    ALG58 =  CONST18s*(STA1 - ALG46);
    ALG65 = ( CONST32s*(STA1 - ALG28))/(1.00000+exp((25.0000 - STA1)/5.98000));
    ALG64 = ( CONST30s*STA4)/(STA4+CONST31s);
    ALG49 =  STA5*CONST11s*(STA1 - ALG28);
!    ALG50 = ALG48+ALG49;
    ALG50 = 0.0D0*ALG48+ALG49; !Modifica FV inibire self-pacing
    RAT1 =  -(ALG52+ALG59+ALG61+ALG53+ALG54+ALG57+ALG62+ALG55+ALG56+ALG63+ALG58+ALG65+ALG64+ALG50) ;
    RAT2 =  ((  - 1.00000*((ALG52+ALG59+ALG49+ALG61+ALG53+ALG54+ALG65) -  2.00000*ALG62))/( 1.00000*CONST5s*CONST3s))*CONST4s;
    ALG66 = CONST44s/(1.00000+(CONST42s**2.00000)/(STA4**2.00000));
    ALG67 =  CONST43s*(STA19 - STA4);
    ALG68 =  CONST41s*(STA12 - STA4);
    ALG70 = 1.00000/(1.00000+( CONST45s*CONST46s)/((STA4+CONST46s)**2.00000));
    RAT4 =  ALG70*((( (ALG67 - ALG66)*CONST51s)/CONST5s+ALG68) - ( 1.00000*((ALG58+ALG64) -  2.00000*ALG63)*CONST4s)/( 2.00000*1.00000*CONST5s*CONST3s));
    ALG69 = CONST38s - (CONST38s - CONST39s)/(1.00000+ (CONST37s/STA19)**2.00000 );
    ALG72 =  CONST34s*ALG69;
    RAT20 =   - ALG72*STA12*STA20+ CONST36s*(1.00000 - STA20);
    ALG71 = CONST33s/ALG69;
    ALG73 = ( ALG71*(STA12**2.00000)*STA20)/(CONST35s+ ALG71*(STA12**2.00000));
    ALG74 =  CONST40s*ALG73*(STA19 - STA12);
    ALG75 = 1.00000/(1.00000+( CONST47s*CONST48s)/((STA19+CONST48s)**2.00000));
    RAT19 =  ALG75*(ALG66 - (ALG74+ALG67));
    ALG76 = 1.00000/(1.00000+( CONST49s*CONST50s)/( (STA12+CONST50s)**2.00000));
    RAT12 =  ALG76*(((  - 1.00000*ALG57*CONST4s)/( 2.00000*1.00000*CONST52s*CONST3s)+( ALG74*CONST51s)/CONST52s) - ( ALG68*CONST5s)/CONST52s);


    XEF_1d(1,i) = STA1 + dtMS*(RAT1+rhspot+Istim); 
    XEF_1d(2,i) = STA2 + dtMS*RAT2
    XEF_1d(3,i) = STA3 + dtMS*RAT3
    XEF_1d(4,i) = STA4 + dtMS*RAT4
    XEF_1d(5,i) = STA5np1
    XEF_1d(6,i) = STA6np1
    XEF_1d(7,i) = STA7np1
    XEF_1d(8,i) = STA8np1
    XEF_1d(9,i) = STA9np1
    XEF_1d(10,i) = STA10np1
    XEF_1d(11,i) = STA11np1
    XEF_1d(12,i) = STA12 + dtMS*RAT12
    XEF_1d(13,i) = STA13np1
    XEF_1d(14,i) = STA14np1
    XEF_1d(15,i) = STA15np1
    XEF_1d(16,i) = STA16np1
    XEF_1d(17,i) = STA17np1
    XEF_1d(18,i) = STA18np1
    XEF_1d(19,i) = STA19 + dtMS*RAT19
    XEF_1d(20,i) = STA20 + dtMS*RAT20
#endif !Courtemanche
    
#ifdef MINIMAL_MODEL !ATRIA !TODO :-> Here put parametrization Minimal Model Purkinjie -----> Find in literature
           !minimal model, endo provisional TODO :-> change parametrization 
         CONST6 = 0.218655D0;
         CONST7 = 0.15D0;
         CONST8 = 7.188259D0;
         CONST9 = 2.031406D0;
         CONST10 = 2.545600D0;
         CONST11 = 0.856429D0;
         
         CONST12 = 13.090045D0;
         CONST13 = 573.569884D0;
         CONST14 = 0.044570D0
         CONST15 = 284.169705D0;
         CONST16 = 0.007950D0;
         CONST17 = 0.168277D0;
         CONST18 = 1.451815D0;
         CONST19 = 28.726329D0;
         CONST20 = 0.281850D0;
         CONST21 = 1.937860D0;
         CONST22 = 0.682362D0;
         CONST23 = 2.347223D0;
         CONST24 = 0.049502D0;
         CONST25 = 0.651285D0;
         CONST26 = 19.751628D0;
         CONST27 = 101.167275D0;
         CONST28 = 57.162563D0;
         CONST29 = 0.028506D0;
         CONST30 = 213.296779D0;
         CONST31 = 19.147136D0;
         CONST32 = 9.999739D0;
         
      

           STA1 = XEF_1d(1,i); !nondimensional here
           STA2 = XEF_1d(2,i);
           STA3 = XEF_1d(3,i);
           STA4 = XEF_1d(4,i);
           
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
           
           XEF_1d(1,i) = STA1 + dtMS*(RAT1+rhspot/(CONST5-CONST4)); 
           XEF_1d(2,i) = STA2np1
           XEF_1d(3,i) = STA3np1
           XEF_1d(4,i) = STA4np1
           
           potEFnode_1d(i) = CONST4+ (CONST5-CONST4)*XEF_1d(1,i)
           
#endif !minimal model atria
        endif !inp

       enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP


#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=1,nvtotEF_1d 
         potEFnode_1d(i) = XEF_1d(1,i)
      enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP
     
         
       !coupling atrial bundles -> AV node
       do i=1,nBundAV
          vmaster=vAVmaster(i)
          vslave=vAVslave(i)
          ! write(*,*) "fsl", vmaster,fslave
          potEFnode_1d(vslave) = potEFnode_1d(vmaster) 
          XEF_1d(1,vslave) = XEF_1d(1,vmaster)          
       enddo

       !coupling atrial bundles -> Purkinje !BSN
       if (bsn.EQ.1) then
          do i=1,nBundPurk
             vmaster=vPurkmaster(i)
             fslave=fPurkslave(i)

             if (i.EQ.1) then !i eq 1 destro
                potEFface_2d(fslave) = potEFnode_1d(vmaster)
                XEF_2d(1,fslave) = XEF_1d(1,vmaster)
             endif
          enddo
       else
          do i=1,nBundPurk
             vmaster=vPurkmaster(i)
             fslave=fPurkslave(i)

             potEFface_2d(fslave) = potEFnode_1d(vmaster)
             XEF_2d(1,fslave) = XEF_1d(1,vmaster)
          enddo
       endif

       !coupling atrial bundles -> Atri
       ! do i=1,nBundAtri
       !    vmaster=vAtrimaster(i)
       !    cslave=cAtrislave(i)
       !    ! write(*,*) "fsl", vmaster,fslave
       !    potEFcell_3d(cslave) = potEFnode_1d(vmaster) 
       !    XEF_3d(1,cslave) = XEF_1d(1,vmaster)          
       ! enddo

       do i=1,nBundAtri
          vmaster=vAtrimaster(i)
          potv=potEFnode_1d(vmaster)
          tempov=Bund_tstim(i)
          if ((potv.GT.-30.0).AND.(tempov.LT.0.d0)) then
                Bund_tstim(i)=timex
          endif          
       enddo

       ! !coupling Purkinje -> Ventr
       ! do i=1,nPurkVentr
       !    fmaster=fVentrmaster(i)
       !    cslave=cVentrslave(i)
       !    potEFcell_3d(cslave) = potEFface_2d(fmaster) 
       !    XEF_3d(1,cslave) = XEF_2d(1,fmaster)          
       ! enddo

       do i=1,nPurkVentr
          fmaster=fVentrmaster(i) 
          potv = potEFface_2d(fmaster)  
          tempov=Purk_tstim(i)
          if ((potv.GT.-30.0).AND.(tempov.LT.0.d0)) then
                Purk_tstim(i)=timex
          endif
       enddo



!TEMPI
#ifdef USE_CUDA 
       !$cuf kernel do (1)
#endif
       do i=1,nvtotEF_1d
          potv = potEFnode_1d(i)
          if ((potv.GT.-30.0)) then
             tempov=EFtstart_1d(i)
             if (nBeat.GT.floor(tempov/period)) then     
                EFtstart_1d(i)=time
             endif
          endif
       enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP



      return
      end

!===================================================
      subroutine ElectroStartEF_1d
      use constants
      use mls_param
!@cuf   use cudafor

      implicit none

      real(DP):: timex,timeEFv,alpv1
      real(DP):: Volt,Cai,CaSR,CaSS
      real(DP):: Nai,Ki,INa_m,INa_h 
      real(DP):: INa_j,IKr_xr1,IKr_xr2,IKs_xs 
      real(DP):: Ito_r,Ito_s,ICaL_d,ICaL_f
      real(DP):: ICaL_f2,ICaL_fCaSS,RR

      real(DP):: xBf,yBf,zBf
      real(DP):: xSA,ySA,zSA,distSAmin,distSA
      real(DP):: xLV,yLV,zLV,xLA,yLA,zLA 
      real(DP):: xRV,yRV,zRV,xRA,yRA,zRA 
      real(DP):: distveG2,px,py,pz,zmaxSA
      real(DP):: Chi,Cm,L_signal,A_signal,amplitude_signal
      integer:: i,indv1,indv2,act,inp,chamb,vSA,vsi,vei
      character*150 :: stri
!@cuf   integer :: istat


!%Initial values of state variables
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


#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=1,nvtotEF_1d 
        inp = vert_to_partEF_1d(i)

        if (inp.EQ.1) then !Courtemanche
#ifdef COURTEMANCHE
              XEF_1d(1,i) =  -81.18;
              XEF_1d(2,i) = 1.117e+01;
              XEF_1d(3,i) = 2.908e-3;
              XEF_1d(4,i) = 9.649e-1;
              XEF_1d(5,i) = 9.775e-1;
              XEF_1d(6,i) = 1.39e+02;
              XEF_1d(7,i) = 3.043e-2;
              XEF_1d(8,i) = 9.992e-1;
              XEF_1d(9,i) = 4.966e-3;
              XEF_1d(10,i) = 9.986e-1;
              XEF_1d(11,i) = 3.296e-5;
              XEF_1d(12,i) = 1.869e-2;
              XEF_1d(13,i) = 1.013e-4;
              XEF_1d(14,i) = 1.367e-4;
              XEF_1d(15,i) = 9.996e-1;
              XEF_1d(16,i) = 7.755e-1;
              XEF_1d(17,i) = 1.488;
              XEF_1d(18,i) = 2.35e-112;
              XEF_1d(19,i) = 1;
              XEF_1d(20,i) = 0.9992;
              XEF_1d(21,i) = 1.488;

              potEFnode_1d(i) = XEF_1d(1,i)
#endif
#ifdef MINIMAL_MODEL
              XEF_1d(1,i) = 0.0D0;
              XEF_1d(2,i) = 1.0D0;
              XEF_1d(3,i) = 1.0D0;
              XEF_1d(4,i) = 0.0D0;
              
              potEFnode_1d(i) = -73.0D0
#endif
           elseif (inp.EQ.2) then  !tTP2006
#ifdef TP06
              XEF_1d(1,i) = Volt 
              XEF_1d(2,i) = Cai 
              XEF_1d(3,i) = CaSR 
              XEF_1d(4,i) = CaSS 
              XEF_1d(5,i) = Nai 
              XEF_1d(6,i) = Ki 
              XEF_1d(7,i) = INa_m 
              XEF_1d(8,i) = INa_h 
              XEF_1d(9,i) = INa_j
              XEF_1d(10,i) = IKr_xr1 
              XEF_1d(11,i) = IKr_xr2 
              XEF_1d(12,i) = IKs_xs 
              XEF_1d(13,i) = Ito_r 
              XEF_1d(14,i) = Ito_s 
              XEF_1d(15,i) = ICaL_d 
              XEF_1d(16,i) = ICaL_f 
              XEF_1d(17,i) = ICaL_f2 
              XEF_1d(18,i) = ICaL_fCaSS 
              XEF_1d(19,i) = RR

              potEFnode_1d(i) = XEF_1d(1,i)
#endif
#ifdef MINIMAL_MODEL
              XEF_1d(1,i) = 0.0D0;
              XEF_1d(2,i) = 1.0D0;
              XEF_1d(3,i) = 1.0D0;
              XEF_1d(4,i) = 0.0D0;
              
              potEFnode_1d(i) = -73.0D0
#endif
           elseif (inp.EQ.3) then !Stewart2009
#ifdef STEWART09
              XEF_1d(1,i) = -69.1370441635924; !%V in component membrane (millivolt)
              XEF_1d(2,i) = 136.781894160227; !%K_i in component potassium_dynamics (millimolar)
              XEF_1d(3,i) = 8.80420286531673; !%Na_i in component sodium_dynamics (millimolar)
              XEF_1d(4,i) = 0.000101878186157052;!%Ca_i in component calcium_dynamics (millimolar)
              XEF_1d(5,i) = 0.0457562667986602; !%y in component hyperpolarization_activated_current_y_gate (dimensionless)
              XEF_1d(6,i) = 0.00550281999719088; !%Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless)
              XEF_1d(7,i) = 0.313213286437995;  !%Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless)
              XEF_1d(8,i) = 0.00953708522974789;!%Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless)   
              XEF_1d(9,i) = 0.0417391656294997; !%m in component fast_sodium_current_m_gate (dimensionless)INa_j
              XEF_1d(10,i) = 0.190678733735145; !%h in component fast_sodium_current_h_gate (dimensionless)
              XEF_1d(11,i) = 0.238219836154029; !%j in component fast_sodium_current_j_gate (dimensionless)
              XEF_1d(12,i) = 0.000446818714055411; !%Ca_ss in component calcium_dynamics (millimolar)
              XEF_1d(13,i) = 0.000287906256206415; !%d in component L_type_Ca_current_d_gate (dimensionless)
              XEF_1d(14,i) = 0.989328560287987; !%f in component L_type_Ca_current_f_gate (dimensionless)
              XEF_1d(15,i) = 0.995474890442185; !%f2 in component L_type_Ca_current_f2_gate (dimensionless)
              XEF_1d(16,i) = 0.999955429598213; !%fCass in component L_type_Ca_current_fCass_gate (dimensionless)
              XEF_1d(17,i) = 0.96386101799501; !%s in component transient_outward_current_s_gate (dimensionless)
              XEF_1d(18,i) = 0.00103618091196912; !%r in component transient_outward_current_r_gate (dimensionless)ICaL_fCaSS 
              XEF_1d(19,i) = 3.10836886659417; !%Ca_SR in component calcium_dynamics (millimolar)
              XEF_1d(20,i) = 0.991580051907845; !%R_prime in component calcium_dynamics (dimensionless)

              potEFnode_1d(i) = XEF_1d(1,i)
#endif
#ifdef MINIMAL_MODEL
              XEF_1d(1,i) = 0.0D0;
              XEF_1d(2,i) = 1.0D0;
              XEF_1d(3,i) = 1.0D0;
              XEF_1d(4,i) = 0.0D0;
              
              potEFnode_1d(i) = -73.0D0 
#endif
      endif !inp

      enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

      !TROVA SA NODE
      xSA=-2.0d0
      ySA= 0.42d0
      zSA=-1.48d0
      
      vsi = vstartEF_1d(1) ; vei = vendEF_1d(1)
      distSAmin=1000.d0
      do i=vsi,vei
         ! if (n_edge_of_vertEF_1d(i).EQ.3) then
         !    if (xyzEF_1d(3,i).GT.zmaxSA) then
         !       vSA = i
         !       zmaxSA=xyzEF_1d(3,i)
         !    endif
         ! endif
         distSA= (xyzEF_1d(1,i)-xSA)**2 + (xyzEF_1d(2,i)-ySA)**2 + (xyzEF_1d(3,i)-zSA)**2
         if (distSA.LT.distSAmin) then
            vSA = i
            distSAmin=distSA
         endif
      enddo
      
      ! px = -1.612
      ! py = -0.618
      ! pz =  4.236
      px = xyzEF_1d(1,vSA)
      py = xyzEF_1d(2,vSA)
      pz = xyzEF_1d(3,vSA)
      write(*,*) "vSA", vSA, px,py,pz

      IstimEF_1d(:)=0.0D0

!ISTIM
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
!      do i=vstartEF_1d(1),vendEF_1d(1) !1,nvtotEF_1d 
      do i=1,nvtotEF_1d 
         xBf = xyz0EF_1d(1,i)
         yBf = xyz0EF_1d(2,i)
         zBf = xyz0EF_1d(3,i)

         distveG2=((xBf-px)**2+(yBf-py)**2+(zBf-pz)**2)*((LSTAR*1000.d0)**2)
         distveG2=distveG2*10.0D0
#ifdef COURTEMANCHE
         IstimEF_1d(i)=20.D0*exp( -(0.02D0*distveG2) ) !CAMBIA STD?
#endif
#ifdef MINIMAL_MODEL
         IstimEF_1d(i)=1.D0*exp( -(0.02D0*distveG2) ) !CAMBIA STD?
#endif
       enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

      call calculate_MintMext_nodesEF_1d

      return
      end


!------------------------------------------------------
        subroutine calculate_MintMext_nodesEF_1d
        use constants
        use mls_param
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snf,enf,v1,v2,v3,v4,i,f1,f2,inp
        real(DP) :: Chi,Cm,g1interp,parco
        real(DP) :: sigma_if,sigma_is,sigma_in
        real(DP) :: sigma_ef,sigma_es,sigma_en
        real(DP) :: sigma_f,sigma_s,sigma_n
        real(DP) :: xfv,yfv,zfv,xsv,ysv,zsv,xnv,ynv,znv

!@cuf   integer :: istat
! #ifdef USE_CUDA
!         attributes(managed) :: vert_of_cell, xyz, vol
! #endif


        inp = 1
!MONO/BIDOMAIN MODEL
        Chi=140; !%mm^-1
        Cm=0.01; !%mF mm^?2 
        ! sigma_if = 0.17/(Cm*Chi);  !% mS / mm
        ! sigma_ef = 0.62/(Cm*Chi);  !% mS / mm
        ! ! sigma_f = sigma_if*sigma_ef/(sigma_if+sigma_ef);
        ! sigma_f = 1.29/(Cm*Chi);  !% mS / mm
        ! sigma_f = sigma_f*0.5*2.0;
        sigma_f=1.214048/(Cm*Chi); 

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=vstartEF_1d(inp),vendEF_1d(inp)
         MintnodesEF_1d(i) = sigma_f
      enddo
!@cuf   istat = cudaDeviceSynchronize !JDR TMP


#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=estartEF_1d(inp),eendEF_1d(inp)
         v1=vert_of_edgeEF_1d(1,i)
         v2=vert_of_edgeEF_1d(2,i)
         MintedgesEF_1d(i) = (MintnodesEF_1d(v1)+MintnodesEF_1d(v2))/2.0D0
      enddo
!@cuf   istat = cudaDeviceSynchronize !JDR TMP



        inp = 2
!MONO/BIDOMAIN MODEL
        Chi=140; !%mm^-1
        Cm=0.01; !%mF mm^?2 
        sigma_f = 0.0129/(Cm*Chi); 
        sigma_f = (0.0129/4.0D0);
        ! sigma_f = (0.0129/5.0D0); 

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=vstartEF_1d(inp),vendEF_1d(inp)
         MintnodesEF_1d(i) = sigma_f
      enddo
!@cuf   istat = cudaDeviceSynchronize !JDR TMP


#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=estartEF_1d(inp),eendEF_1d(inp)
         v1=vert_of_edgeEF_1d(1,i)
         v2=vert_of_edgeEF_1d(2,i)
         MintedgesEF_1d(i) = (MintnodesEF_1d(v1)+MintnodesEF_1d(v2))/2.0D0
      enddo
!@cuf   istat = cudaDeviceSynchronize !JDR TMP


        inp = 3
!MONO/BIDOMAIN MODEL
        Chi=140; !%mm^-1
        Cm=0.01; !%mF mm^?2 
!         sigma_if = 0.17/(Cm*Chi);  !% mS / mm
!         sigma_ef = 0.62/(Cm*Chi);  !% mS / mm
! !        sigma_f = sigma_if*sigma_ef/(sigma_if+sigma_ef);
! !        sigma_f = 3.95/(Cm*Chi); !FV ricambia
! !        sigma_f = 2.00/(Cm*Chi);
!         sigma_f = 1.00/(Cm*Chi);


        sigma_f=0.570839/(Cm*Chi);
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=vstartEF_1d(inp),vendEF_1d(inp)
         MintnodesEF_1d(i) = sigma_f
      enddo
!@cuf   istat = cudaDeviceSynchronize !JDR TMP


#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=estartEF_1d(inp),eendEF_1d(inp)
         v1=vert_of_edgeEF_1d(1,i)
         v2=vert_of_edgeEF_1d(2,i)
         MintedgesEF_1d(i) = (MintnodesEF_1d(v1)+MintnodesEF_1d(v2))/2.0D0
      enddo
!@cuf   istat = cudaDeviceSynchronize !JDR TMP


        return
        end subroutine calculate_MintMext_nodesEF_1d
