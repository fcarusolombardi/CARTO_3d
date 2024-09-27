!===================================================
      subroutine ElectroRunEF_2d(timex)
      use constants
      use mls_param
      USE ieee_arithmetic
!@cuf   use cudafor

      implicit none
  
      real(DP):: timex,timeEFv,alpv1,timeMS,dtMS
      real(DP):: timeVentr,timeAtrio
      real(DP):: rhspot,num,den,dvvjj,potsur,volcell
      real(DP):: distCF,g1interp,corrsten,bordo,potv,tempov
      real(DP):: xB1,yB1,zB1,xV1,yV1,zV1
      real(DP):: xgrad,ygrad,zgrad
      real(DP):: xvers,yvers,zvers
      real(DP):: xnormale,ynormale,znormale
      real(DP):: mint11,mint12,mint13,mint21,mint22,mint23,mint31,mint32,mint33
      real(DP):: mext11,mext12,mext13,mext21,mext22,mext23,mext31,mext32,mext33
      real(DP):: msum11,msum12,msum13,msum21,msum22,msum23,msum31,msum32,msum33
      real(DP):: Istim
      real(DP):: sta1,sta2,sta3,sta4,sta5,sta6,sta7,sta8,sta9
      real(DP):: sta10,sta11,sta12,sta13,sta14,sta15,sta16,sta17,sta18,sta19,sta20,sta21
      real(DP):: sta1np1,sta2np1,sta3np1,sta4np1,sta5np1,sta6np1,sta7np1,sta8np1,sta9np1
      real(DP):: sta10np1,sta11np1,sta12np1,sta13np1,sta14np1,sta15np1,sta16np1,sta17np1,sta18np1,sta19np1,sta20np1,sta21np1
      real(DP):: rat1,rat2,rat3,rat4,rat5,rat6,rat7,rat8,rat9
      real(DP):: rat10,rat11,rat12,rat13,rat14,rat15,rat16,rat17,rat18,rat19,rat20,rat21

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

      ! minimal model variables
      !------------------------------------------------------------------------------
      real(DP):: const4,const5,const6,const7,const8,const9
      real(DP):: const10,const11,const12,const13,const14,const15,const16,const17,const18,const19
      real(DP):: const20,const21,const22,const23,const24,const25,const26,const27,const28,const29
      real(DP):: const30,const31,const32,app1,app2,app3,app4,app5
      real(DP):: Dij,GINF,TAU,alpha,remod
      !------------------------------------------------------------------------------

      !-------------------------------------------
      ! Fibrillation Variables Declaration

      character(len=4):: dummy
      real(DP):: xS1,yS1,zS1
      real(DP):: xS2,yS2,zS2
      real(DP):: amplitudeS1,amplitudeS2
      real(DP):: duration_signal_S1,duration_signal_S2,delayS2
      real(DP):: pot_apd,pot_old,treshold_apd
      !CV
      real(DP):: tpq,tpr,tmin,tmmax,sign_face
      real(DP):: x1_1,x1_2,x1_3
      real(DP):: x2_1,x2_2,x2_3
      real(DP):: x3_1,x3_2,x3_3
      real(DP):: xv1_1,xv1_2,xv1_3
      real(DP):: xv2_1,xv2_2,xv2_3
      real(DP):: xv3_1,xv3_2,xv3_3
      real(DP):: dx12,dx13,dx23
      real(DP):: t12,t13,t23
      real(DP):: cvf1,cvf2,cvf3
      real(DP):: cvp1,cvp2,cvp3
      real(DP):: normx12,normx13,norm23
      real(DP):: xpq_norm,xpr_norm,xpr_norm2,xpq_norm2,xqr_norm2
      real(DP):: xpq1,xpq2,xpq3
      real(DP):: nor1,nor2,nor3
      real(DP):: cos_theta,sin_theta,tan_alpha,tmpcv
      integer :: ip,iq,ir,imin,imax,f1,f2,f3
      !integer(DP):: nBeat
      !-------------------------------------------


      real(DP):: aoo, duration_signal,sur_face,ftempo
      integer:: i,j,k,inp,chamb,v1,v2,v3,e1
      integer:: fslave,vmaster
      character*150 :: stri
!@cuf   integer :: istat
!#ifdef USE_CUDA
!        attributes(managed) :: gradcell_3d
!#endif 

      open(unit=15,file='S1S2_2dStim.in',status='old')
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

     
      treshold_apd = -65.0D0
      
!NOTAFV per considerare l'effetto delle piccole variazione geometriche durante la
!depolarizzazione dovremmo aggiornare (ag ogni time step o multiplo) anche le informazioni
!geometriche del dominio chiamando le routines che sono alla fine di read_geoSingleBody_3d
!tutto è predisposto, manca solo un modo per aggiornare/seguire l'orientamento delle fibre 
#ifdef STEWART09
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
#endif


#ifdef MV_PURK !here ENDO parametrization

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
         
      
#endif

    ! timeMS=timex!*TSTAR*1000.d0
    ! ! dtMS=dt!*TSTAR*1000.d0
    ! dtMS=dt/real(stepEF_2d)
    timeMS=timex*TSTAR*1000.d0
    dtMS=dt*TSTAR*1000.d0/real(stepEF_2d) 


!find nodal values
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
    do i=1,nvtotEF_2d 
        num=0.0D0
        den=0.0D0
#ifdef ELEGEO0
        xV1 = xyz0EF_2d(1,i)
        yV1 = xyz0EF_2d(2,i)
        zV1 = xyz0EF_2d(3,i)        
#else
        xV1 = xyzEF_2d(1,i)
        yV1 = xyzEF_2d(2,i)
        zV1 = xyzEF_2d(3,i)
#endif        
        do j=1,n_face_of_vertEF_2d(i)
            f1  = face_of_vertEF_2d(j,i)
            xB1 = tri_barEF_2d(1,f1)
            yB1 = tri_barEF_2d(2,f1)
            zB1 = tri_barEF_2d(3,f1)

            dvvjj = sqrt( (xV1-xB1)**2+(yV1-yB1)**2+(zV1-zB1)**2)
            num = num + potEFface_2d(f1)/dvvjj
            den = den + 1.D0/dvvjj            
        enddo
        potEFnode_2d(i)=num/den
    enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP
   
!!COMPUTE GRADIENT
!find pot on edges
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
    do i=1,netotEF_2d 
       v1=vert_of_edgeEF_2d(1,i)
       v2=vert_of_edgeEF_2d(2,i)
       potEFedge_2d(i)=(potEFnode_2d(v1)+potEFnode_2d(v2))*0.5D0;
    enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP
!find gradient on cells (gauss-green theorem)
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=1,nftotEF_2d 
        xgrad=0.0D0
        ygrad=0.0D0
        zgrad=0.0D0
        do j=1,3
           e1=edge_of_faceEF_2d(j,i)
           xnormale=normaledgeoffacesEF_2d(1,j,i)
           ynormale=normaledgeoffacesEF_2d(2,j,i)
           znormale=normaledgeoffacesEF_2d(3,j,i)

#ifdef ELEGEO0
           potsur = potEFedge_2d(e1)*dist0EF_2d(e1)
#else           
          potsur = potEFedge_2d(e1)*distEF_2d(e1)
#endif
           xgrad = xgrad + potsur*xnormale
           ygrad = ygrad + potsur*ynormale
           zgrad = zgrad + potsur*znormale
        enddo

#ifdef ELEGEO0        
        sur_face = sur0EF_2d(i)
#else        
        sur_face = surEF_2d(i)
#endif
        
        gradfaceEF_2d(1,i) =  xgrad/sur_face
        gradfaceEF_2d(2,i) =  ygrad/sur_face
        gradfaceEF_2d(3,i) =  zgrad/sur_face
     enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!find gradient on edges
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
     do i=1,netotEF_2d 
        f1 = face_of_edgeEF_2d(1,i)
        f2 = face_of_edgeEF_2d(2,i)
        g1interp = g1interpedgeEF_2d(i)
        xvers  = versCFedgeEF_2d(1,i)
        yvers  = versCFedgeEF_2d(2,i)
        zvers  = versCFedgeEF_2d(3,i)
        distCF = distCFedgeEF_2d(i)
        if ((f1.NE.0).AND.(f2.NE.0)) then
           !grandiente medio
           xgrad = g1interp*gradfaceEF_2d(1,f1) + (1-g1interp)*gradfaceEF_2d(1,f2)
           ygrad = g1interp*gradfaceEF_2d(2,f1) + (1-g1interp)*gradfaceEF_2d(2,f2)
           zgrad = g1interp*gradfaceEF_2d(3,f1) + (1-g1interp)*gradfaceEF_2d(3,f2)
           !gradiente corretto
            !gradientecorretto=gradientemedio+(-dot(gradientemedio,versCF)+(pot(c2)-pot(c1))/distCF)*versCF;
           corrsten= -(xgrad*xvers+ygrad*yvers+zgrad*zvers) + (potEFface_2d(f2)-potEFface_2d(f1))/distCF
           xgrad = xgrad + corrsten*xvers
           ygrad = ygrad + corrsten*yvers
           zgrad = zgrad + corrsten*zvers           
        elseif ((f1.NE.0).AND.(f2.EQ.0)) then
           ! xgrad = gradcell_3d(1,c1)
           ! ygrad = gradcell_3d(2,c1)
           ! zgrad = gradcell_3d(3,c1)
           xgrad = 0.0D0 !for Neumann bcs
           ygrad = 0.0D0 
           zgrad = 0.0D0 
        elseif ((f1.EQ.0).AND.(f2.NE.0)) then
           ! xgrad = gradcell_3d(1,c2)
           ! ygrad = gradcell_3d(2,c2)
           ! zgrad = gradcell_3d(3,c2)
           xgrad = 0.0D0 !for Neumann bcs
           ygrad = 0.0D0 
           zgrad = 0.0D0 
        endif

        gradedgeEF_2d(1,i) = xgrad
        gradedgeEF_2d(2,i) = ygrad
        gradedgeEF_2d(3,i) = zgrad
     enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!CELL MODEL
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=1,nftotEF_2d 

         !stimulus
         !        Istim = 0.0!IstimEF_2d(i)
         if (nBeat.eq.0) then
            if ((timeMS.GT.0.0D0).AND.(timeMS.LT.duration_signal_S1)) then
               Istim = IstimEF_2dS1(i)
            elseif ((timeMS.GT.delayS2).AND.(timeMS.LT.duration_signal_S2+delayS2)) then
               Istim = IstimEF_2dS2(i)
            else
               Istim = 0.0
            endif
         else
            Istim = 0.0D0
         endif
         
         !spatial term
         bordo = 0.D0
         do j=1,3
            e1=edge_of_faceEF_2d(j,i)
            xnormale=normaledgeoffacesEF_2d(1,j,i)
            ynormale=normaledgeoffacesEF_2d(2,j,i)
            znormale=normaledgeoffacesEF_2d(3,j,i)
            
            xgrad = gradedgeEF_2d(1,e1) 
            ygrad = gradedgeEF_2d(2,e1) 
            zgrad = gradedgeEF_2d(3,e1) 
            
            mint11 = MintedgesEF_2d(1,1,e1)
            mint12 = MintedgesEF_2d(1,2,e1)
            mint13 = MintedgesEF_2d(1,3,e1)
            mint21 = MintedgesEF_2d(2,1,e1)
            mint22 = MintedgesEF_2d(2,2,e1)
            mint23 = MintedgesEF_2d(2,3,e1)
            mint31 = MintedgesEF_2d(3,1,e1)
            mint32 = MintedgesEF_2d(3,2,e1)
            mint33 = MintedgesEF_2d(3,3,e1)  
            
!             bordo = bordo + &
!                  (  xnormale*(mint11*xgrad + mint12*ygrad + mint13*zgrad) &
!                  + ynormale*(mint21*xgrad + mint22*ygrad + mint23*zgrad) &
! #ifdef ELEGEO0
!                  + znormale*(mint31*xgrad + mint32*ygrad + mint33*zgrad) )*dist0EF_2d(e1)
! #else           
!             + znormale*(mint31*xgrad + mint32*ygrad + mint33*zgrad) )*distEF_2d(e1)
! #endif           
             bordo = bordo + &
                  CARTO_Dface(i)*(  xnormale*xgrad + ynormale*ygrad &
#ifdef ELEGEO0                
                  + znormale*zgrad )*dist0EF_2d(e1)
#else
                  + znormale*zgrad )*distEF_2d(e1)
#endif
              

         enddo
#ifdef ELEGEO0
         rhspot = bordo/sur0EF_2d(i)
#else       
         rhspot = bordo/surEF_2d(i)
#endif
         rhspot = rhspot/((1000.D0*LSTAR)**2) !CONTROLLA devo dimens lo spazio in mm             

!       if (abs(rhspot - 69.138).LT.0.01) then
       ! if ((i.EQ.3425).OR.(i.EQ.3426).OR.(i.EQ.3427)) then
       !    write(*,*) "sss", ntime, i, rhspot
       ! if ((ntime.EQ.8)) then
       !    stop
       ! endif
       ! endif
        ! rhspot = rhspot/((LSTAR)**2) !devo dimens lo spazio in mm                                                             
#ifdef STEWART09 
         STA1 = XEF_2d(1,i); !Volt                                      
         STA2 = XEF_2d(2,i);
         STA3 = XEF_2d(3,i);
         STA4 = XEF_2d(4,i);
         STA5 = XEF_2d(5,i);
         STA6 = XEF_2d(6,i);
         STA7 = XEF_2d(7,i);
         STA8 = XEF_2d(8,i);
         STA9 = XEF_2d(9,i);
         STA10 = XEF_2d(10,i);
         STA11 = XEF_2d(11,i);
         STA12 = XEF_2d(12,i);
         STA13 = XEF_2d(13,i);
         STA14 = XEF_2d(14,i);
         STA15 = XEF_2d(15,i);
         STA16 = XEF_2d(16,i);
         STA17 = XEF_2d(17,i);
         STA18 = XEF_2d(18,i);
         STA19 = XEF_2d(19,i);
         STA20 = XEF_2d(20,i);

           
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


!    XEF_2d(1,i) = STA1 + dtMS*(RAT1+rhspot+Istim);
         XEF_2d(1,i) = STA1 + dtMS*(RAT1+rhspot);  !ricambia
         XEF_2d(2,i) = STA2 + dtMS*RAT2
         XEF_2d(3,i) = STA3 + dtMS*RAT3
         XEF_2d(4,i) = STA4 + dtMS*RAT4
         XEF_2d(5,i) = STA5np1
         XEF_2d(6,i) = STA6np1
         XEF_2d(7,i) = STA7np1
         XEF_2d(8,i) = STA8np1
         XEF_2d(9,i) = STA9np1
         XEF_2d(10,i) = STA10np1
         XEF_2d(11,i) = STA11np1
         XEF_2d(12,i) = STA12 + dtMS*RAT12
         XEF_2d(13,i) = STA13np1
         XEF_2d(14,i) = STA14np1
         XEF_2d(15,i) = STA15np1
         XEF_2d(16,i) = STA16np1
         XEF_2d(17,i) = STA17np1
         XEF_2d(18,i) = STA18np1
         XEF_2d(19,i) = STA19 + dtMS*RAT19
         XEF_2d(20,i) = STA20 + dtMS*RAT20
    
         potEFface_2d(i) = XEF_2d(1,i)
#endif
    
#ifdef MV_PURK !here ENDO parametrization  

           STA1 = XEF_2d(1,i); !nondimensional here
           STA2 = XEF_2d(2,i);
           STA3 = XEF_2d(3,i);
           STA4 = XEF_2d(4,i);
           
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
           
           
           if (CARTO_Dface(i).gt.0.49D0) then
              remod=1.0D0
           else
              remod=CARTO_Dface(i)
           endif
           XEF_2d(1,i) = STA1 + dtMS*(remod*RAT1+rhspot/(CONST5-CONST4)); 
           XEF_2d(2,i) = STA2np1
           XEF_2d(3,i) = STA3np1
           XEF_2d(4,i) = STA4np1
           
           potEFface_2d(i) = CONST4+ (CONST5-CONST4)*XEF_2d(1,i)
        
#endif !minimal model atria
       enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP


!TEMPI
#ifdef USE_CUDA 
       !$cuf kernel do (1)
#endif
       do i=1,nvtotEF_2d
          potv = potEFnode_2d(i)
          if (potv.GT.-30.0) then
             tempov=EFtstart_2d(i)
             if (nBeat.GT.floor(tempov/period)) then     
                EFtstart_2d(i)=time
             endif
          endif

          !Two-dimensional APD data
          pot_apd = CONST4+potEFnode_2d(i)*(CONST5-CONST4)
          pot_old = CONST4+old_1_2d(i)*(CONST5-CONST4)
          if ((pot_apd.gt.treshold_apd.and.pot_old.lt.treshold_apd).and.&
               f_apd_2d(i).eq.0) then
             t_apd_2d(1,i) = (dtms/(pot_apd-pot_old))*(treshold_apd-pot_apd+((pot_apd-pot_old)/dtms)*timeMS)
             f_apd_2d(i) = 1
          elseif ((pot_apd.lt.treshold_apd.and.pot_old.gt.treshold_apd).and.&
               f_apd_2d(i).eq.1) then
             t_apd_2d(2,i) =  (dtms/(pot_apd-pot_old))*(treshold_apd-pot_apd+((pot_apd-pot_old)/dtms)*timeMS)
             f_apd_2d(i) = 2
          elseif ((pot_apd.gt.treshold_apd.and.pot_old.lt.treshold_apd).and.&
               f_apd_2d(i).eq.2) then
             t_apd_2d(3,i) = (dtms/(pot_apd-pot_old))*(treshold_apd-pot_apd+((pot_apd-pot_old)/dtms)*timeMS)
             f_apd_2d(i) = 3
          elseif ((pot_apd.lt.treshold_apd.and.pot_old.gt.treshold_apd).and.&
             f_apd_2d(i).eq.3) then
             t_apd_2d(4,i) = (dtms/(pot_apd-pot_old))*(treshold_apd-pot_apd+((pot_apd-pot_old)/dtms)*timeMS)
             f_apd_2d(i) = 4
          endif
          old_1_2d(i) = potEFnode_2d(i)

       enddo
       !@cuf istat = cudaDeviceSynchronize !JDR TMP


! !CV manifold 2d
! !!COMPUTE GRADIENT
       ! !find pot on edges
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
    do i=1,netotEF_2d 
       v1=vert_of_edgeEF_2d(1,i)
       v2=vert_of_edgeEF_2d(2,i)
       potEFedge_2d(i)=(t_apd_2d(1,v1)+t_apd_2d(1,v2))*0.5D0;
    enddo
    !@cuf istat = cudaDeviceSynchronize !JDR TMP
    
!find gradient on cells (gauss-green theorem)
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=1,nftotEF_2d 
        xgrad=0.0D0
        ygrad=0.0D0
        zgrad=0.0D0
        do j=1,3
           e1=edge_of_faceEF_2d(j,i)
           xnormale=normaledgeoffacesEF_2d(1,j,i)
           ynormale=normaledgeoffacesEF_2d(2,j,i)
           znormale=normaledgeoffacesEF_2d(3,j,i)

#ifdef ELEGEO0
           potsur = potEFedge_2d(e1)*dist0EF_2d(e1)*1000.D0*LSTAR
#else           
          potsur = potEFedge_2d(e1)*distEF_2d(e1)*1000.D0*LSTAR
#endif
           xgrad = xgrad + potsur*xnormale
           ygrad = ygrad + potsur*ynormale
           zgrad = zgrad + potsur*znormale
        enddo

#ifdef ELEGEO0        
        sur_face = sur0EF_2d(i)
#else        
        sur_face = surEF_2d(i)
#endif
        
        ! gradfaceEF_2d(1,i) =  xgrad/sur_face
        ! gradfaceEF_2d(2,i) =  ygrad/sur_face
        ! gradfaceEF_2d(3,i) =  zgrad/sur_face
        cvp1 =sqrt(xgrad**2+ygrad**2+zgrad**2)/(sur_face*(1000.D0*LSTAR)**2)
        cvp2 =1.0D0/cvp1
        ! if (cvp1.eq.0.0D0) then
        !    CVface_EF_2d(i) = 0.0D0
        ! elseif (cvp2.gt.0.90D0) then
        !    CVface_EF_2d(i) = 100000.0D0
        !else
        CVface_EF_2d(i) = cvp1
        !endif
       
     enddo
     !@cuf istat = cudaDeviceSynchronize !JDR TMP       

! #ifdef USE_CUDA
!         !$cuf kernel do (1)
! #endif
!      do i=1,nftotEF_2d

!         v1 = vert_of_faceEF_2d(1,i)
!         v2 = vert_of_faceEF_2d(2,i)
!         v3 = vert_of_faceEF_2d(3,i)
!         if (CVface_EF_2d(i).gt.1000.0D0) then

!            CVface_EF_2d(i)=0.0D0
!            tmpcv = 0.0D0
!            do k=1,3
!               v1 = vert_of_faceEF_2d(k,i)
!               num=0.0D0
!               den=0.0D0
! #ifdef ELEGEO0
!               xV1 = xyz0EF_2d(1,v1)
!               yV1 = xyz0EF_2d(2,v1)
!               zV1 = xyz0EF_2d(3,v1)        
! #else
!               xV1 = xyzEF_2d(1,v1)
!               yV1 = xyzEF_2d(2,v1)
!               zV1 = xyzEF_2d(3,v1)
! #endif        
!               do j=1,n_face_of_vertEF_2d(v1)
!                  f1  = face_of_vertEF_2d(j,v1)
!                  xB1 = tri_barEF_2d(1,f1)
!                  yB1 = tri_barEF_2d(2,f1)
!                  zB1 = tri_barEF_2d(3,f1)

!                  dvvjj = sqrt( (xV1-xB1)**2+(yV1-yB1)**2+(zV1-zB1)**2)
!                  num = num + CVface_EF_2d(f1)/dvvjj
!                  den = den + 1.D0/dvvjj            
!               enddo
!               tmpcv=tmpcv+num/den
!            enddo
!            CVface_EF_2d(i)=tmpcv/3.0D0
           
           
!         else
!            cycle
!         endif
!      enddo
!      !@cuf istat = cudaDeviceSynchronize !JDR TMP       

      return
      end

!===================================================
      subroutine ElectroStartEF_2d
      use constants
      use mls_param
!@cuf   use cudafor

      implicit none

      real(DP):: timex,timeEFv,alpv1
      ! real(DP):: Volt,Cai,CaSR,CaSS
      ! real(DP):: Nai,Ki,INa_m,INa_h 
      ! real(DP):: INa_j,IKr_xr1,IKr_xr2,IKs_xs 
      ! real(DP):: Ito_r,Ito_s,ICaL_d,ICaL_f
      ! real(DP):: ICaL_f2,ICaL_fCaSS,RR

      real(DP):: xBf,yBf,zBf
      real(DP):: xLV,yLV,zLV,xLA,yLA,zLA 
      real(DP):: xRV,yRV,zRV,xRA,yRA,zRA 
      real(DP):: distveG2,px,py,pz
      real(DP):: Chi,Cm,L_signal,A_signal,amplitude_signal
      integer:: i,indv1,indv2,act,inp,chamb
      character*150 :: stri

      !-------------------------------------------------------
      ! S1S2-Stimulation prodocol cariables
      real(DP):: tmp_2dS1,alpha,amplitudeS1,amplitudeS2
      real(DP):: S11,S12,S13,S14
      real(DP):: S21,S22,S23,S24
      real(DP):: xS1,yS1,zS1
      real(DP):: xS2,yS2,zS2
      real(DP):: pxS1,pyS1,pzS1
      real(DP):: pxS2,pyS2,pzS2
      real(DP):: distS1,distS2
      real(DP):: st_dev3D,MFstim
      real(DP):: xBc,yBc,zBc
      integer::v1,v2,v3,v4
      
      character(len=4):: dummy
      !-------------------------------------------------------

      
!@cuf   integer :: istat

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=1,nftotEF_2d 
#ifdef STEWART09
              XEF_2d(1,i) = -69.1370441635924; !%V in component membrane (millivolt)
              XEF_2d(2,i) = 136.781894160227; !%K_i in component potassium_dynamics (millimolar)
              XEF_2d(3,i) = 8.80420286531673; !%Na_i in component sodium_dynamics (millimolar)
              XEF_2d(4,i) = 0.000101878186157052;!%Ca_i in component calcium_dynamics (millimolar)
              XEF_2d(5,i) = 0.0457562667986602; !%y in component hyperpolarization_activated_current_y_gate (dimensionless)
              XEF_2d(6,i) = 0.00550281999719088; !%Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless)
              XEF_2d(7,i) = 0.313213286437995;  !%Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless)
              XEF_2d(8,i) = 0.00953708522974789;!%Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless)   
              XEF_2d(9,i) = 0.0417391656294997; !%m in component fast_sodium_current_m_gate (dimensionless)INa_j
              XEF_2d(10,i) = 0.190678733735145; !%h in component fast_sodium_current_h_gate (dimensionless)
              XEF_2d(11,i) = 0.238219836154029; !%j in component fast_sodium_current_j_gate (dimensionless)
              XEF_2d(12,i) = 0.000446818714055411; !%Ca_ss in component calcium_dynamics (millimolar)
              XEF_2d(13,i) = 0.000287906256206415; !%d in component L_type_Ca_current_d_gate (dimensionless)
              XEF_2d(14,i) = 0.989328560287987; !%f in component L_type_Ca_current_f_gate (dimensionless)
              XEF_2d(15,i) = 0.995474890442185; !%f2 in component L_type_Ca_current_f2_gate (dimensionless)
              XEF_2d(16,i) = 0.999955429598213; !%fCass in component L_type_Ca_current_fCass_gate (dimensionless)
              XEF_2d(17,i) = 0.96386101799501; !%s in component transient_outward_current_s_gate (dimensionless)
              XEF_2d(18,i) = 0.00103618091196912; !%r in component transient_outward_current_r_gate (dimensionless)ICaL_fCaSS 
              XEF_2d(19,i) = 3.10836886659417; !%Ca_SR in component calcium_dynamics (millimolar)
              XEF_2d(20,i) = 0.991580051907845; !%R_prime in component calcium_dynamics (dimensionless)

              potEFface_2d(i) = XEF_2d(1,i)
#endif

#ifdef MV_PURK
              XEF_2d(1,i) = 0.0D0;
              XEF_2d(2,i) = 1.0D0;
              XEF_2d(3,i) = 1.0D0;
              XEF_2d(4,i) = 0.0D0;

              potEFface_2d(i) = -83.0D0
#endif
         enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP


      ! Istim0EF_2d(:)=0.0D0
         IstimEF_2dS1(:)=0.0D0
         IstimEF_2dS2(:)=0.0D0

         open(unit=15,file='S1S2_2dStim.in',status='old')
         read(15,*) dummy
         read(15,*) xS1,yS1,zS1
         read(15,*) dummy
         read(15,*) xS2,yS2,zS2
         read(15,*) dummy
         read(15,*) amplitudeS1,amplitudeS2
         read(15,*) dummy
         read(15,*) st_dev3D,MFstim
         close(15)
         
         write(*,*) '2D S1S2 Stimulus'
         write(*,*) "xS1 yS1 zS1", xS1, yS1, zS1
         write(*,*) "xS2 yS2 zS2", xS2, yS2, zS2

         
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do i=1,nftotEF_2d         

         ! v1 = vert_of_faceEF_2d(1,i)
         ! v2 = vert_of_faceEF_2d(2,i)
         ! v3 = vert_of_faceEF_2d(3,i)
         
         ! ! ----------------------------------------
         ! ! V1
         ! pxS1 = xyz0EF_2d(1,v1)
         ! pyS1 = xyz0EF_2d(2,v1)
         ! pzS1 = xyz0EF_2d(3,v1)
         
         ! distS1 = sqrt( (pxS1-xS1)**2 + (pyS1-yS1)**2 + (pzS1-zS1)**2 )
         ! distS2	= sqrt( (pxS1-xS2)**2 + (pyS1-yS2)**2 + (pzS1-zS2)**2 )
         ! S11 = amplitudeS1*exp(-(st_dev3D*distS1*MFstim))
         ! S21 = amplitudeS2*exp(-(st_dev3D*distS2*MFstim)) 
         ! ! ----------------------------------------
         ! ! V2
         ! pxS1 = xyz0EF_2d(1,v2)
         ! pyS1 = xyz0EF_2d(2,v2)
         ! pzS1 = xyz0EF_2d(3,v2)

         ! distS1 = sqrt( (pxS1-xS1)**2 + (pyS1-yS1)**2 + (pzS1-zS1)**2 )
         ! distS2 = sqrt( (pxS1-xS2)**2 + (pyS1-yS2)**2 + (pzS1-zS2)**2 )
         ! S12 = amplitudeS1*exp(-(st_dev3D*distS1*MFstim)) 
         ! S22 = amplitudeS2*exp(-(st_dev3D*distS2*MFstim))
         ! ! ----------------------------------------
         ! ! V3
         ! pxS1 = xyz0EF_2d(1,v3)
         ! pyS1 = xyz0EF_2d(2,v3)
         ! pzS1 = xyz0EF_2d(3,v3)

         ! distS1 = sqrt( (pxS1-xS1)**2 + (pyS1-yS1)**2 + (pzS1-zS1)**2 )
         ! distS2 = sqrt( (pxS1-xS2)**2 + (pyS1-yS2)**2 + (pzS1-zS2)**2 )
         ! S13 = amplitudeS1*exp(-(st_dev3D*distS1*MFstim)) 
         ! S23 = amplitudeS2*exp(-(st_dev3D*distS2*MFstim)) 
         
         ! IstimEF_2dS1(i) = -(S11+S12+S13)/3.0D0
         ! IstimEF_2dS2(i) = -(S21+S22+S23)/3.0D0
         xBc=tri_barEF_2d(1,i)
         yBc=tri_barEF_2d(2,i)
         zBc=tri_barEF_2d(3,i)
         
         distveG2=sqrt((xBc-xS1)**2+(yBc-yS1)**2+(zBc-zS1)**2)
         if ((distveG2.LT.(5.D0/(1000.D0*LSTAR))).and.(scar_cell(i).ne.2)) then
            !               IstimEF_3d(i)=-60.D0
            IstimEF_2dS1(i)=-1.D0 !minimal stimolo
         endif
         distveG2=sqrt((xBc-xS2)**2+(yBc-yS2)**2+(zBc-zS2)**2)
         if ((distveG2.LT.(8.0D0/(1000.D0*LSTAR))).and.(scar_cell(i).ne.2)) then 
            !               IstimEF_3d(i)=-60.D0
            IstimEF_2dS2(i)=-1.D0 !minimal stimolo
         endif
      enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP
         
      call calculate_MintMext_facesEF_2d

      return
      end


!------------------------------------------------------
        subroutine calculate_MintMext_facesEF_2d
        use constants
        use mls_param
        use param
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snf,enf,v1,v2,v3,v4,i,j,f1,f2,latoPurk
        real(DP) :: Chi,Cm,g1interp,parco,val_mono
        real(DP) :: sigma_if,sigma_is,sigma_in
        real(DP) :: sigma_ef,sigma_es,sigma_en
        real(DP) :: sigma_f,sigma_s,sigma_n
        real(DP) :: xfv,yfv,zfv,xsv,ysv,zsv,xnv,ynv,znv,scalePurk
        real(DP):: xB1,yB1,zB1,xV1,yV1,zV1,num,den,dvvjj
        real(DP),dimension(:),allocatable :: tmp_D

!@cuf   integer :: istat
#ifdef USE_CUDA
!         attributes(managed) :: vert_of_cell, xyz, vol
        attributes(managed) :: tmp_D
#endif
        allocate(tmp_D(nvtotEF_2d))
!MONO/BIDOMAIN MODEL
!       Chi=140; !%mm^-1
!       Cm=0.01; !%mF mm^?2 

!       ! sigma_if = 2.0*3.95/(Cm*Chi);  !% mS / mm                                                                                                                                               
!       ! sigma_is = 2.0*3.95/(Cm*Chi); !% mS / mm                                                                                                                                                
!       ! sigma_in = 2.0*3.95/(Cm*Chi); !% mS / mm                                                                                                                                                
!       ! sigma_ef = 2.0*3.95/(Cm*Chi);  !% mS / mm                                                                                                                                               
!       ! sigma_es = 2.0*3.95/(Cm*Chi);  !% mS / mm                                                                                                                                               
!       ! sigma_en = 2.0*3.95/(Cm*Chi);  !% mS / mm 
      
! !      val_mono=0.75 -> 202ms
!       val_mono=2.0*1.0 !173
!       val_mono=0.001 !ricambia
!       sigma_if = 2.0D0*val_mono/(Cm*Chi);  !% mS / mm                                          
!       sigma_is = 2.0D0*val_mono/(Cm*Chi); !% mS / mm                                           
!       sigma_in = 2.0D0*val_mono/(Cm*Chi); !% mS / mm                                           
!       sigma_ef = 2.0D0*val_mono/(Cm*Chi);  !% mS / mm                                          
!       sigma_es = 2.0D0*val_mono/(Cm*Chi);  !% mS / mm                                          
!       sigma_en = 2.0D0*val_mono/(Cm*Chi);  !% mS / mm 

!       sigma_f = sigma_if*sigma_ef/(sigma_if+sigma_ef);
!       sigma_s = sigma_is*sigma_es/(sigma_is+sigma_es);
!       sigma_n = sigma_in*sigma_en/(sigma_in+sigma_en);
!       ! Mint = [sigma_if 0 0;0 sigma_is 0;0 0 sigma_in]/(Cm*Chi);
!       ! Mext = [sigma_ef 0 0;0 sigma_es 0;0 0 sigma_en]/(Cm*Chi);
!       ! Mmon = [sigma_f 0 0;0 sigma_s 0;0 0 sigma_n]/(Cm*Chi);

!       sigma_f = 3.159691/(Cm*Chi);
!       sigma_s = 3.159691/(Cm*Chi);
!       sigma_n = 3.159691/(Cm*Chi);

        CARTO_Dnode(:)=0.0D0
        tmp_D(:)=0.0D0
        CARTO_Dface(:)=0.0D0
        scar_tag_2d(:)=0
      !open(unit=15,file='meshes/CARTO_diffusivity_2d.txt',status='old')
      do i=vstartEF_2d(1),vendEF_2d(1)
         !read(15,*) CARTO_Dnode(i)!tmp_D(i)
         tmp_D(i)=0.5D0
      enddo
      close(15)
      !open(unit=15,file='meshes/scar_tag_2d.txt',status='old')
      !do i=vstartEF_2d(1),vendEF_2d(1)
         !read(15,*) scar_tag_2d(i)=0
      !   scar_tag_2d(i)=0
      !enddo
      !close(15)
      write(*,*) 'CARTO diffusivity field loaded'

! #ifdef USE_CUDA
!         !$cuf kernel do (1)
! #endif
!       do i=vstartEF_2d(1),vendEF_2d(1)

!          if (scar_tag_2d(i).eq.1) then ! Border Zone
!             CARTO_Dnode(i)=0.1*sqrt(tmp_D(i))
!          elseif (scar_tag_2d(i).eq.2) then !Scar
!             CARTO_Dnode(i)=0.005*sqrt(tmp_D(i))
!          elseif (scar_tag_2d(i).eq.3) then !Septum
!             CARTO_Dnode(i)=1.0e-6*sqrt(tmp_D(i))
!          else
!             CARTO_Dnode(i)=sqrt(tmp_D(i))
!          endif
         
!       enddo
! !@cuf istat = cudaDeviceSynchronize !JDR TMP
      deallocate(tmp_D)
      
      !find face values
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
    do i=fstartEF_2d(1),fendEF_2d(1)
        num=0.0D0
        den=0.0D0

        xB1 = tri_barEF_2d(1,i)
        yB1 = tri_barEF_2d(2,i)
        zB1 = tri_barEF_2d(3,i)
        do j=1,3
           v1 = vert_of_faceEF_2d(j,i) 
           
#ifdef ELEGEO0
           xV1 = xyz0EF_2d(1,v1)
           yV1 = xyz0EF_2d(2,v1)
           zV1 = xyz0EF_2d(3,v1)        
#else
           xV1 = xyzEF_2d(1,v1)
           yV1 = xyzEF_2d(2,v1)
           zV1 = xyzEF_2d(3,v1)
#endif
            dvvjj = sqrt( (xV1-xB1)**2+(yV1-yB1)**2+(zV1-zB1)**2)
            num = num + CARTO_Dnode(v1)/dvvjj
            den = den + 1.D0/dvvjj            
        enddo
        CARTO_Dface(i)=num/den
    enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP


      
!      k=sigma_f

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=1,nftotEF_2d 
         ! xfv =AmatrFibersEF_2d(1,1,i) 
         ! yfv =AmatrFibersEF_2d(2,1,i) 
         ! zfv =AmatrFibersEF_2d(3,1,i) 
         ! xsv =AmatrFibersEF_2d(1,2,i) 
         ! ysv =AmatrFibersEF_2d(2,2,i) 
         ! zsv =AmatrFibersEF_2d(3,2,i) 
         ! xnv =AmatrFibersEF_2d(1,3,i) 
         ! ynv =AmatrFibersEF_2d(2,3,i) 
         ! znv =AmatrFibersEF_2d(3,3,i) 

         xfv =1.0D0
         yfv =0.0D0 
         zfv =0.0D0
         xsv =0.0D0
         ysv =1.0D0
         zsv =0.0D0
         xnv =0.0D0
         ynv =0.0D0
         znv =1.0D0
         
         sigma_if=CARTO_Dface(i)
         sigma_in=sigma_if
         sigma_is=sigma_if
         
         sigma_ef=CARTO_Dface(i)
         sigma_en=sigma_ef
         sigma_es=sigma_ef

         ! latoPurk=LabelPurkRightLeft(i)
         ! if (latoPurk.EQ.1) then !right            
         !    scalePurk=scaleRightPurk
         ! elseif (latoPurk.EQ.2) then !left
         !    scalePurk=scaleLeftPurk
         ! endif
         scalePurk=0.1D0 !differenziazione left/right purk non attiva in questa versione (guarda ultima ECG solo elettrico)

         MintfacesEF_2d(1,1,i) = scalePurk*(sigma_if*xfv**2  + sigma_in*xnv**2  + sigma_is*xsv**2)
         MintfacesEF_2d(1,2,i) = scalePurk*(sigma_if*xfv*yfv + sigma_in*xnv*ynv + sigma_is*xsv*ysv)
         MintfacesEF_2d(1,3,i) = scalePurk*(sigma_if*xfv*zfv + sigma_in*xnv*znv + sigma_is*xsv*zsv)
         MintfacesEF_2d(2,1,i) = scalePurk*(sigma_if*xfv*yfv + sigma_in*xnv*ynv + sigma_is*xsv*ysv)
         MintfacesEF_2d(2,2,i) = scalePurk*(sigma_if*yfv**2  + sigma_in*ynv**2  + sigma_is*ysv**2)
         MintfacesEF_2d(2,3,i) = scalePurk*(sigma_if*yfv*zfv + sigma_in*ynv*znv + sigma_is*ysv*zsv)
         MintfacesEF_2d(3,1,i) = scalePurk*(sigma_if*xfv*zfv + sigma_in*xnv*znv + sigma_is*xsv*zsv)
         MintfacesEF_2d(3,2,i) = scalePurk*(sigma_if*yfv*zfv + sigma_in*ynv*znv + sigma_is*ysv*zsv)
         MintfacesEF_2d(3,3,i) = scalePurk*(sigma_if*zfv**2  + sigma_in*znv**2  + sigma_is*zsv**2)

         MextfacesEF_2d(1,1,i) = scalePurk*(sigma_ef*xfv**2  + sigma_en*xnv**2  + sigma_es*xsv**2)
         MextfacesEF_2d(1,2,i) = scalePurk*(sigma_ef*xfv*yfv + sigma_en*xnv*ynv + sigma_es*xsv*ysv)
         MextfacesEF_2d(1,3,i) = scalePurk*(sigma_ef*xfv*zfv + sigma_en*xnv*znv + sigma_es*xsv*zsv)
         MextfacesEF_2d(2,1,i) = scalePurk*(sigma_ef*xfv*yfv + sigma_en*xnv*ynv + sigma_es*xsv*ysv)
         MextfacesEF_2d(2,2,i) = scalePurk*(sigma_ef*yfv**2  + sigma_en*ynv**2  + sigma_es*ysv**2)
         MextfacesEF_2d(2,3,i) = scalePurk*(sigma_ef*yfv*zfv + sigma_en*ynv*znv + sigma_es*ysv*zsv)
         MextfacesEF_2d(3,1,i) = scalePurk*(sigma_ef*xfv*zfv + sigma_en*xnv*znv + sigma_es*xsv*zsv)
         MextfacesEF_2d(3,2,i) = scalePurk*(sigma_ef*yfv*zfv + sigma_en*ynv*znv + sigma_es*ysv*zsv)
         MextfacesEF_2d(3,3,i) = scalePurk*(sigma_ef*zfv**2  + sigma_en*znv**2  + sigma_es*zsv**2)



      end do
!@cuf   istat = cudaDeviceSynchronize !JDR TMP

! #ifndef BIDOMAIN !FV always bidomain for purkinje
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=1,nftotEF_2d 

         ! xfv =AmatrFibersEF_2d(1,1,i) 
         ! yfv =AmatrFibersEF_2d(2,1,i) 
         ! zfv =AmatrFibersEF_2d(3,1,i) 
         ! xsv =AmatrFibersEF_2d(1,2,i) 
         ! ysv =AmatrFibersEF_2d(2,2,i) 
         ! zsv =AmatrFibersEF_2d(3,2,i) 
         ! xnv =AmatrFibersEF_2d(1,3,i) 
         ! ynv =AmatrFibersEF_2d(2,3,i) 
         ! znv =AmatrFibersEF_2d(3,3,i) 

         xfv =1.0D0
         yfv =0.0D0 
         zfv =0.0D0
         xsv =0.0D0
         ysv =1.0D0
         zsv =0.0D0
         xnv =0.0D0
         ynv =0.0D0
         znv =1.0D0
         
         sigma_f=CARTO_Dface(i)
         sigma_n=sigma_f
         sigma_s=sigma_f
         ! latoPurk=LabelPurkRightLeft(i)
         ! if (latoPurk.EQ.1) then !right            
         !    scalePurk=scaleRightPurk
         ! elseif (latoPurk.EQ.2) then !left
         !    scalePurk=scaleLeftPurk
         ! endif
         scalePurk=0.1D0 !differenziazione left/right purk non attiva in questa versione (guarda ultima ECG solo elettrico)
         scalePurk=1.0D0 !minimal model carto 2d
         MintfacesEF_2d(1,1,i) = scalePurk*(sigma_f*xfv**2  + sigma_n*xnv**2  + sigma_s*xsv**2)
         MintfacesEF_2d(1,2,i) = scalePurk*(sigma_f*xfv*yfv + sigma_n*xnv*ynv + sigma_s*xsv*ysv)
         MintfacesEF_2d(1,3,i) = scalePurk*(sigma_f*xfv*zfv + sigma_n*xnv*znv + sigma_s*xsv*zsv)
         MintfacesEF_2d(2,1,i) = scalePurk*(sigma_f*xfv*yfv + sigma_n*xnv*ynv + sigma_s*xsv*ysv)
         MintfacesEF_2d(2,2,i) = scalePurk*(sigma_f*yfv**2  + sigma_n*ynv**2  + sigma_s*ysv**2)
         MintfacesEF_2d(2,3,i) = scalePurk*(sigma_f*yfv*zfv + sigma_n*ynv*znv + sigma_s*ysv*zsv)
         MintfacesEF_2d(3,1,i) = scalePurk*(sigma_f*xfv*zfv + sigma_n*xnv*znv + sigma_s*xsv*zsv)
         MintfacesEF_2d(3,2,i) = scalePurk*(sigma_f*yfv*zfv + sigma_n*ynv*znv + sigma_s*ysv*zsv)
         MintfacesEF_2d(3,3,i) = scalePurk*(sigma_f*zfv**2  + sigma_n*znv**2  + sigma_s*zsv**2)
      enddo
! #endif 
!@cuf   istat = cudaDeviceSynchronize !JDR TMP

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
    do i=1,netotEF_2d 
       f1=face_of_edgeEF_2d(1,i)
       f2=face_of_edgeEF_2d(2,i)        
       g1interp=g1interpedgeEF_2d(i)
       if ((f1.NE.0).AND.(f2.NE.0)) then
          MintedgesEF_2d(1,1,i)=g1interp*MintfacesEF_2d(1,1,f1)+(1-g1interp)*MintfacesEF_2d(1,1,f2)
          MintedgesEF_2d(1,2,i)=g1interp*MintfacesEF_2d(1,2,f1)+(1-g1interp)*MintfacesEF_2d(1,2,f2)
          MintedgesEF_2d(1,3,i)=g1interp*MintfacesEF_2d(1,3,f1)+(1-g1interp)*MintfacesEF_2d(1,3,f2)
          MintedgesEF_2d(2,1,i)=g1interp*MintfacesEF_2d(2,1,f1)+(1-g1interp)*MintfacesEF_2d(2,1,f2)
          MintedgesEF_2d(2,2,i)=g1interp*MintfacesEF_2d(2,2,f1)+(1-g1interp)*MintfacesEF_2d(2,2,f2)
          MintedgesEF_2d(2,3,i)=g1interp*MintfacesEF_2d(2,3,f1)+(1-g1interp)*MintfacesEF_2d(2,3,f2)
          MintedgesEF_2d(3,1,i)=g1interp*MintfacesEF_2d(3,1,f1)+(1-g1interp)*MintfacesEF_2d(3,1,f2)
          MintedgesEF_2d(3,2,i)=g1interp*MintfacesEF_2d(3,2,f1)+(1-g1interp)*MintfacesEF_2d(3,2,f2)
          MintedgesEF_2d(3,3,i)=g1interp*MintfacesEF_2d(3,3,f1)+(1-g1interp)*MintfacesEF_2d(3,3,f2)

          MextedgesEF_2d(1,1,i)=g1interp*MextfacesEF_2d(1,1,f1)+(1-g1interp)*MextfacesEF_2d(1,1,f2)
          MextedgesEF_2d(1,2,i)=g1interp*MextfacesEF_2d(1,2,f1)+(1-g1interp)*MextfacesEF_2d(1,2,f2)
          MextedgesEF_2d(1,3,i)=g1interp*MextfacesEF_2d(1,3,f1)+(1-g1interp)*MextfacesEF_2d(1,3,f2)
          MextedgesEF_2d(2,1,i)=g1interp*MextfacesEF_2d(2,1,f1)+(1-g1interp)*MextfacesEF_2d(2,1,f2)
          MextedgesEF_2d(2,2,i)=g1interp*MextfacesEF_2d(2,2,f1)+(1-g1interp)*MextfacesEF_2d(2,2,f2)
          MextedgesEF_2d(2,3,i)=g1interp*MextfacesEF_2d(2,3,f1)+(1-g1interp)*MextfacesEF_2d(2,3,f2)
          MextedgesEF_2d(3,1,i)=g1interp*MextfacesEF_2d(3,1,f1)+(1-g1interp)*MextfacesEF_2d(3,1,f2)
          MextedgesEF_2d(3,2,i)=g1interp*MextfacesEF_2d(3,2,f1)+(1-g1interp)*MextfacesEF_2d(3,2,f2)
          MextedgesEF_2d(3,3,i)=g1interp*MextfacesEF_2d(3,3,f1)+(1-g1interp)*MextfacesEF_2d(3,3,f2)
       elseif ((f1.NE.0).AND.(f2.EQ.0)) then
          MintedgesEF_2d(1,1,i)=MintfacesEF_2d(1,1,f1)
          MintedgesEF_2d(1,2,i)=MintfacesEF_2d(1,2,f1)
          MintedgesEF_2d(1,3,i)=MintfacesEF_2d(1,3,f1)
          MintedgesEF_2d(2,1,i)=MintfacesEF_2d(2,1,f1)
          MintedgesEF_2d(2,2,i)=MintfacesEF_2d(2,2,f1)
          MintedgesEF_2d(2,3,i)=MintfacesEF_2d(2,3,f1)
          MintedgesEF_2d(3,1,i)=MintfacesEF_2d(3,1,f1)
          MintedgesEF_2d(3,2,i)=MintfacesEF_2d(3,2,f1)
          MintedgesEF_2d(3,3,i)=MintfacesEF_2d(3,3,f1)

          MextedgesEF_2d(1,1,i)=MextfacesEF_2d(1,1,f1)
          MextedgesEF_2d(1,2,i)=MextfacesEF_2d(1,2,f1)
          MextedgesEF_2d(1,3,i)=MextfacesEF_2d(1,3,f1)
          MextedgesEF_2d(2,1,i)=MextfacesEF_2d(2,1,f1)
          MextedgesEF_2d(2,2,i)=MextfacesEF_2d(2,2,f1)
          MextedgesEF_2d(2,3,i)=MextfacesEF_2d(2,3,f1)
          MextedgesEF_2d(3,1,i)=MextfacesEF_2d(3,1,f1)
          MextedgesEF_2d(3,2,i)=MextfacesEF_2d(3,2,f1)
          MextedgesEF_2d(3,3,i)=MextfacesEF_2d(3,3,f1)
       elseif ((f1.EQ.0).AND.(f2.NE.0)) then
          MintedgesEF_2d(1,1,i)=MintfacesEF_2d(1,1,f2)
          MintedgesEF_2d(1,2,i)=MintfacesEF_2d(1,2,f2)
          MintedgesEF_2d(1,3,i)=MintfacesEF_2d(1,3,f2)
          MintedgesEF_2d(2,1,i)=MintfacesEF_2d(2,1,f2)
          MintedgesEF_2d(2,2,i)=MintfacesEF_2d(2,2,f2)
          MintedgesEF_2d(2,3,i)=MintfacesEF_2d(2,3,f2)
          MintedgesEF_2d(3,1,i)=MintfacesEF_2d(3,1,f2)
          MintedgesEF_2d(3,2,i)=MintfacesEF_2d(3,2,f2)
          MintedgesEF_2d(3,3,i)=MintfacesEF_2d(3,3,f2)

          MextedgesEF_2d(1,1,i)=MextfacesEF_2d(1,1,f2)
          MextedgesEF_2d(1,2,i)=MextfacesEF_2d(1,2,f2)
          MextedgesEF_2d(1,3,i)=MextfacesEF_2d(1,3,f2)
          MextedgesEF_2d(2,1,i)=MextfacesEF_2d(2,1,f2)
          MextedgesEF_2d(2,2,i)=MextfacesEF_2d(2,2,f2)
          MextedgesEF_2d(2,3,i)=MextfacesEF_2d(2,3,f2)
          MextedgesEF_2d(3,1,i)=MextfacesEF_2d(3,1,f2)
          MextedgesEF_2d(3,2,i)=MextfacesEF_2d(3,2,f2)
          MextedgesEF_2d(3,3,i)=MextfacesEF_2d(3,3,f2)
       endif    
    enddo
!@cuf   istat = cudaDeviceSynchronize !JDR TMP

104 format(9(2x,f12.8))
105 format(1(2x,f12.8))
106 format(1(2x,i8))
107 format(3(2x,f12.8))

        return
        end subroutine calculate_MintMext_facesEF_2d



!===================================================
      subroutine GradientEF_2d(timex)
      use constants
      use mls_param
!@cuf   use cudafor

      implicit none

      real(DP):: timex,timeEFv,alpv1,timeMS,dtMS
      real(DP):: timeVentr,timeAtrio
      real(DP):: rhspot,num,den,dvvjj,potsur,volcell
      real(DP):: distCF,g1interp,corrsten,bordo
      real(DP):: xB1,yB1,zB1,xV1,yV1,zV1
      real(DP):: xgrad,ygrad,zgrad
      real(DP):: xvers,yvers,zvers
      real(DP):: xnormale,ynormale,znormale
      real(DP):: mint11,mint12,mint13,mint21,mint22,mint23,mint31,mint32,mint33
      real(DP):: mext11,mext12,mext13,mext21,mext22,mext23,mext31,mext32,mext33
      real(DP):: msum11,msum12,msum13,msum21,msum22,msum23,msum31,msum32,msum33
      real(DP):: Istim

      real(DP):: aoo, duration_signal,sur_face,ftempo
      integer:: i,j,inp,chamb,f1,f2,v1,v2,v3,e1
      integer:: fslave,vmaster
      character*150 :: stri
!@cuf   integer :: istat
!#ifdef USE_CUDA
!        attributes(managed) :: gradcell_3d
!#endif 

!NOTAFV per considerare l'effetto delle piccole variazione geometriche durante la
!depolarizzazione dovremmo aggiornare (ag ogni time step o multiplo) anche le informazioni
!geometriche del dominio chiamando le routines che sono alla fine di read_geoSingleBody_3d
!tutto è predisposto, manca solo un modo per aggiornare/seguire l'orientamento delle fibre 

!find nodal values
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
    do i=1,nvtotEF_2d 
        num=0.0D0
        den=0.0D0
        ! xV1 = xyzEF_2d(1,i)
        ! yV1 = xyzEF_2d(2,i)
        ! zV1 = xyzEF_2d(3,i)
        xV1 = xyz0EF_2d(1,i)
        yV1 = xyz0EF_2d(2,i)
        zV1 = xyz0EF_2d(3,i)
        do j=1,n_face_of_vertEF_2d(i)
            f1  = face_of_vertEF_2d(j,i)
            xB1 = tri_barEF_2d(1,f1)
            yB1 = tri_barEF_2d(2,f1)
            zB1 = tri_barEF_2d(3,f1)

            dvvjj = sqrt( (xV1-xB1)**2+(yV1-yB1)**2+(zV1-zB1)**2)
            num = num + potEFface_2d(f1)/dvvjj
            den = den + 1.D0/dvvjj            
        enddo
        potEFnode_2d(i)=num/den
    enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!!COMPUTE GRADIENT
!find pot on edges
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
    do i=1,netotEF_2d 
       v1=vert_of_edgeEF_2d(1,i)
       v2=vert_of_edgeEF_2d(2,i)
       potEFedge_2d(i)=(potEFnode_2d(v1)+potEFnode_2d(v2))/2.D0;
    enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!find gradient on cells (gauss-green theorem)
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=1,nftotEF_2d 
        xgrad=0.0D0
        ygrad=0.0D0
        zgrad=0.0D0
        do j=1,3
           e1=edge_of_faceEF_2d(j,i)
           xnormale=normaledgeoffacesEF_2d(1,j,i)
           ynormale=normaledgeoffacesEF_2d(2,j,i)
           znormale=normaledgeoffacesEF_2d(3,j,i)
           
!           potsur = potEFedge_2d(e1)*distEF_2d(e1)
           potsur = potEFedge_2d(e1)*dist0EF_2d(e1)
           xgrad = xgrad + potsur*xnormale
           ygrad = ygrad + potsur*ynormale
           zgrad = zgrad + potsur*znormale
        enddo
!        sur_face = surEF_2d(i)
        sur_face = sur0EF_2d(i)
        gradfaceEF_2d(1,i) =  xgrad/sur_face
        gradfaceEF_2d(2,i) =  ygrad/sur_face
        gradfaceEF_2d(3,i) =  zgrad/sur_face

     enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

      
      return
      end
