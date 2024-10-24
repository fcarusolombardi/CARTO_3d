subroutine preambolotsch
        !@cuf use cudafor
        !@cuf use cudadevice 
      use param
      use local_arrays
      use mgrd_arrays
      use mpih
      use mpi_param, only: kstart,kend,kstartr,kendr
      use mls_param     
      use mls_local, only: for_xc,for_yc,for_zc, for_sc, coll
      use ibm_param, only: ftAS, ftAS2, ftAS3
      use probes
      use nvtx
      implicit none
!@cuf integer :: istat
      integer :: ns,inp,ntr,nsub
      integer :: j,k,i,v1,ic,jc,kc,aa,iii
      integer :: vsi, vei, esi, eei, fsi, fei
      real(DP) :: potHand, potHandND, ftAS3max
      real(DP) :: pressAO, pressLV,pressLA
      real(DP) :: pressAP, pressRV,pressRA
      real(DP) :: dinvm,alpv1,timeExt,shiftSMZ
      real(DP) :: integrale1,integrale2,integrale3,dtEFshift
      real(DP) :: Kpid_p,Kpid_i,Kpid_d,saturo,pressT
      real(DP) :: alp_LA,VolumeLV_T,VolumeLV,pesoAV,timeMS

      !EGM
#ifdef EGM
      real(DP)::EGMadd,EGMsum,invdist,invdist3,volcell
      real(DP)::Delta_egm1,Delta_egm2,Delta_egm3
#ifdef USE_CUDA
      attributes(managed)::EGMadd,invdist,invdist3,volcell
      attributes(managed)::Delta_egm1,Delta_egm2,Delta_egm3
#endif   
#endif

      
      integer :: jjj,indv1,indv2, jjEFfast,csi,cei
      character*50 filename,strucfilename
#ifdef DEBUG
      integer :: kstartp
      real(DP) :: cksum2,cksum3,mck2,mck3,cksum1,mck1
#endif

      timeMS=time*(TSTAR*1000.d0)
      
      !Find pressures for contact model
      pressOUT = 0
      !LEFT HEART
      if (myid .eq. rankAO) then
#ifdef USE_CUDA
        istat = cudaMemcpy(pressOUT(1), pr(masonAO(1), masonAO(2), masonAO(3)), 1, cudaMemcpyDeviceToHost)
#else
        pressOUT(1) = pr(masonAO(1), masonAO(2), masonAO(3))
#endif
      end if

      if (myid .eq. rankLV) then
#ifdef USE_CUDA
        istat = cudaMemcpy(pressOUT(2), pr(masonLV(1), masonLV(2), masonLV(3)), 1, cudaMemcpyDeviceToHost)
#else
        pressOUT(2) = pr(masonLV(1), masonLV(2), masonLV(3))
#endif
      end if
 
      if (myid .eq. rankLA) then
#ifdef USE_CUDA
        istat = cudaMemcpy(pressOUT(3), pr(masonLA(1), masonLA(2), masonLA(3)), 1, cudaMemcpyDeviceToHost)
#else
        pressOUT(3) = pr(masonLA(1), masonLA(2), masonLA(3))
#endif
      end if

!       !RIGHT HEART
!       if (myid .eq. rankAP) then
! #ifdef USE_CUDA
!         istat = cudaMemcpy(pressOUT(4), pr(masonAP(1), masonAP(2), masonAP(3)), 1, cudaMemcpyDeviceToHost)
! #else
!         pressOUT(4) = pr(masonAP(1), masonAP(2), masonAP(3))
! #endif
!       end if

!       if (myid .eq. rankRV) then
! #ifdef USE_CUDA
!         istat = cudaMemcpy(pressOUT(5), pr(masonRV(1), masonRV(2), masonRV(3)), 1, cudaMemcpyDeviceToHost)
! #else
!         pressOUT(5) = pr(masonRV(1), masonRV(2), masonRV(3))
! #endif
!       end if
 
!       if (myid .eq. rankRA) then
! #ifdef USE_CUDA
!         istat = cudaMemcpy(pressOUT(6), pr(masonRA(1), masonRA(2), masonRA(3)), 1, cudaMemcpyDeviceToHost)
! #else
!         pressOUT(6) = pr(masonRA(1), masonRA(2), masonRA(3))
! #endif
!       end if
      
      !call MPI_Bcast(pressAO,1,MPI_DOUBLE,rankAO,MPI_COMM_WORLD,ierr)
      !call MPI_Bcast(pressVE,1,MPI_DOUBLE,rankVE,MPI_COMM_WORLD,ierr)
      !call MPI_Bcast(pressAT,1,MPI_DOUBLE,rankAT,MPI_COMM_WORLD,ierr)
!      call MPI_Allreduce(MPI_IN_PLACE,pressOUT,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(MPI_IN_PLACE,pressOUT,6,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

      pressAO = pressOUT(1)
      pressLV = pressOUT(2)
      pressLA = pressOUT(3)
      pressAP = pressOUT(4)
      pressRV = pressOUT(5)
      pressRA = pressOUT(6)

! Set Contact flags
      !Aortic valve
      if(pressAO.gt.pressLV)then 
            ContactA = 1 
      else
            ContactA = 0
      end if
      !Mitral valve  
      if(pressLV.gt.pressLA)then
            ContactM = 1
      else
            ContactM = 0
      end if
      ! !Pulmonary valve  
      ! if(pressAP.gt.pressRV)then
      !       ContactP = 1
      ! else
      !       ContactP = 0
      ! end if
      ! !Tricuspid valve  
      ! if(pressRV.gt.pressRA)then
      !       ContactT = 1
      ! else
      !       ContactT = 0
      ! end if

!->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
      !Aorta Cap and Pulmonary Veins (posso spostarlo in gcurv)
!      ftAS2 = 0.001d0
      ftAS2 = 0.0d0
      ftAS = 0.d0
!      ftAS3 = -40.d0
      ! ftAS3 = 15.d0
      ftAS3 = 25.d0


#ifdef ELECTRO
      shiftSMZ=10.0D0
#else
      shiftSMZ=6.0D0
#endif
      shiftSMZ = shiftSMZ + ( dyssy/(TSTAR*1000.0d0) )

      ! SMZ = 1.d0 - (1.d0-SMZ_S)*( 0.5d0*( tanh( 1.d0*(timeBeat-  (shiftSMZ-6.0D0)        )) +1.d0 )  &
      !                     + 0.5d0*( tanh(-1.d0*(timeBeat-   (shiftSMZ+6.0D0)         )) +1.d0 )  & !qio
      !                     + 0.5d0*( tanh( 1.d0*(timeBeat- ((shiftSMZ-6.0D0)+Period))) +1.d0 ) -1.d0 )

      SMZ = 1.0d0
!->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
!->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
!Stress
pesoAV = dt/tframe
#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                                                                                                                    
#endif
      do i=1,nftot
         taufaceAV(i)   = taufaceAV(i)   + pesoAV*sqrt(tauface(i))
         pressfaceAV(i) = pressfaceAV(i) + pesoAV*sqrt(pressface(i))
      enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
! ELECTROPHYSIOLOGY
!->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
#ifdef ELECTRO
      if ((nbeat.GE.0).AND.(ECG_Ppro.LE.0)) then


         if (nat_pac.EQ.1) then !pacing naturale; bundles+Purk
#ifndef ELEGEO0
 !ALTRE UTILIZZATE PER EF
     !     inp=1
     !     vsi = vstartEF_2d(inp) ; vei = vendEF_2d(inp)
     !     esi = estartEF_2d(inp) ; eei = eendEF_2d(inp)
     !     fsi = fstartEF_2d(inp) ; fei = fendEF_2d(inp)
     !     call convert_geoEF_2d(vsi,vei,esi,eei,        &
     !        fsi,fei,xyzEF_2d(:,vsi:vei),            &
     !        vert_of_faceEF_2d(:,fsi:fei),           &
     !        tri_barEF_2d(:,fsi:fei))

     !    call calculate_normal(vsi,vei,fsi,             &
     !                         fei,xyzEF_2d(:,vsi:vei),      &
     !                         vert_of_faceEF_2d(:,fsi:fei),  &
     !                        tri_norEF_2d(:,fsi:fei))


     !    call calculate_distance(distEF_2d(esi:eei),vsi,vei,esi,eei,    &
     !                   xyzEF_2d(:,vsi:vei),vert_of_edgeEF_2d(:,esi:eei))

     !    call calculate_area(SurfaceEF_2d(inp),vsi,vei,fsi,fei,         &
     !                   xyzEF_2d(:,vsi:vei),vert_of_faceEF_2d(:,fsi:fei),    &
     !                   surEF_2d(fsi:fei))

     !    call calculate_ginterp_edge_facesEF_2d(vsi,vei,esi,eei,fsi,fei,         &
     !                   xyzEF_2d(:,vsi:vei),vert_of_edgeEF_2d(:,esi:eei),    &
     !                   face_of_edgeEF_2d(:,esi:eei),tri_barEF_2d(:,fsi:fei), &
     !                   versCFedgeEF_2d(:,esi:eei),distCFedgeEF_2d(esi:eei),g1interpedgeEF_2d(esi:eei))

     !    call calculate_normals_edge_facesEF_2d(vsi,vei,esi,eei,fsi,fei,xyzEF_2d(:,vsi:vei),  &
     !                   vert_of_edgeEF_2d(:,fsi:fei),edge_of_faceEF_2d(:,fsi:fei),   & 
     !                   tri_barEF_2d(:,fsi:fei),tri_norEF_2d(:,fsi:fei),normaledgeoffacesEF_2d(:,:,fsi:fei))



     !  do inp=1,NparticleEF_1d
     !     vsi = vstartEF_1d(inp) ; vei = vendEF_1d(inp)
     !     esi = estartEF_1d(inp) ; eei = eendEF_1d(inp)

     !     call convert_geoEF_1d(vsi,vei,esi,eei,        &
     !        xyzEF_1d(:,vsi:vei),            &
     !        vert_of_edgeEF_1d(:,esi:eei),           &
     !        edg_barEF_1d(:,esi:eei))

     !    call calculate_distance(distEF_1d(esi:eei),vsi,vei,esi,eei,    &
     !                   xyzEF_1d(:,vsi:vei),vert_of_edgeEF_1d(:,esi:eei))
     ! enddo !
!END ALTRE UTILIZZATE PER EF
#endif           
         ! do jjEFfast=1,stepEF_1d
         !    dtEFshift = real(1-jjEFfast)*dt/real(stepEF_1d)
         !    call ElectroRunEF_1d(timeBeat+dtEFshift)
         ! enddo
         ! do jjEFfast=1,stepEF_2d
         !    dtEFshift = real(1-jjEFfast)*dt/real(stepEF_2d)
         !    call ElectroRunEF_2d(timeBeat+dtEFshift)
         ! enddo
         endif !nat_pacing





         ! if (mod(ntime,stepEF_3d).EQ.0) then
         if ((mod(ntime,stepEF_3d).EQ.0).OR.(Ppro.GT.0)) then
#ifndef ELEGEO0
!ALTRE  UTILIZZATE PER EF
      call calculate_bar_cells(1,nvtot_3d,1,nctot_3d,xyz_3d(:,1:nvtot_3d),  &
           vert_of_cell_3d(:,1:nctot_3d),cell_bar(:,1:nctot_3d))

      do inp=1,1
         call calculate_area(Surface_3d(inp),vstart_3d(inp),vend_3d(inp),fstart_3d(inp),fend_3d(inp),         &
                    xyz_3d(:,vstart_3d(inp):vend_3d(inp)),vert_of_face_3d(:,fstart_3d(inp):fend_3d(inp)),    &
                           sur_3d(fstart_3d(inp):fend_3d(inp)))
      end do

      call calculate_ginterp_face_cells(1,nvtot_3d,1,nftot_3d,1,nctot_3d,         &
           xyz_3d(:,1:nvtot_3d),vert_of_face_3d(:,1:nftot_3d),    &
           cell_of_face_3d(:,1:nftot_3d),cell_bar(:,1:nctot_3d), &
           versCFface_3d(:,1:nftot_3d),distCFface_3d(1:nftot_3d),g1interpface_3d(1:nftot_3d))

      call calculate_normals_face_cells(1,nvtot_3d,1,nftot_3d,1,nctot_3d,xyz_3d(:,1:nftot_3d),  &
           vert_of_face_3d(:,1:nftot_3d),face_of_cell_3d(:,1:nctot_3d),   &
           cell_bar(:,1:nctot_3d),normalfaceofcells_3d(:,:,1:nctot_3d))
!END ALTRE UTILIZZATE PER EF
#endif         
            call ElectroRunEF_3d(timeBeat,nBeat) !solo con dt fiss         !TENSIONE ATTIVA
            call ElectroTensionEF_3d(timeBeat)
         endif


      endif !((nbeat.GE.0).AND.(ECG_Ppro.LE.0))


!->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
! ECG
!->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->    
115 format(11(2x,f12.4))
        !we calculate ecg at the same rate of saving ecg restarts through continua_ecg
      if ((mod(time,tpin).LT.dt).OR.(ECG_Ppro.GT.0)) then
!--------------------3D----------------------!
            call GradientEF_3d

            inp=1
            csi=cstart_3d(inp)
            cei=cend_3d(inp)
            !ATRI E VENTRICOLI 
            do j=1,ECG_nson
               call calculate_ECG(ECG_PHI(j),csi,cei,ECG_E(:,j),cell_bar(:,csi:cei), &
#ifdef ELEGEO0
                    vol0_3d(csi:cei),gradcell_3d(:,csi:cei), &
#else
                    vol_3d(csi:cei),gradcell_3d(:,csi:cei), &
#endif                    
                    LabelSmooth_cell(csi:cei),cell_to_chamb_3dBou(csi:cei),ECG_D)
            enddo
            if(myid.eq.0) write(196,115)timeMS,ECG_PHI(1:ECG_nson)
            if(myid.eq.0) flush(196)


            !Solo ATRI 
            do j=1,ECG_nson
               call calculate_ECG_Atr(ECG_PHI(j),csi,cei,ECG_E(:,j),cell_bar(:,csi:cei), &
#ifdef ELEGEO0
                    vol0_3d(csi:cei),gradcell_3d(:,csi:cei), &
#else
                    vol_3d(csi:cei),gradcell_3d(:,csi:cei), &
#endif                    
                    LabelSmooth_cell(csi:cei),cell_to_chamb_3dBou(csi:cei),ECG_D)
            enddo
            if(myid.eq.0) write(197,115)timeMS,ECG_PHI(1:ECG_nson)
            if(myid.eq.0) flush(197)

            !Solo VENTR
            do j=1,ECG_nson
               call calculate_ECG_Ventr(ECG_PHI(j),csi,cei,ECG_E(:,j),cell_bar(:,csi:cei), &
#ifndef ELEGEO0                    
                                  vol0_3d(csi:cei),gradcell_3d(:,csi:cei), &
#else                    
                                  vol_3d(csi:cei),gradcell_3d(:,csi:cei), &
#endif
                                  LabelSmooth_cell(csi:cei),cell_to_chamb_3dBou(csi:cei),ECG_D)
            enddo
            if(myid.eq.0) write(198,115)timeMS,ECG_PHI(1:ECG_nson)
            if(myid.eq.0) flush(198)


!--------------------2D----------------------!
            call GradientEF_2d
            fsi=1
            fei=nftotEF_2d
            !ATRI E VENTRICOLI 
            do j=1,ECG_nson
               call calculate_ECG_2D(ECG_PHI(j),fsi,fei,ECG_E(:,j),tri_barEF_2d(:,fsi:fei), &
                                  surEF_2d(fsi:fei),gradfaceEF_2d(:,fsi:fei), &
                                  ECG_D)
            enddo
            if(myid.eq.0) write(199,115)timeMS,ECG_PHI(1:ECG_nson)
            if(myid.eq.0) flush(199)

 
!         endif
         endif
         
#endif
!->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
! PRESSURE.OUT
!->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
!!!104 format(6(2x,f12.8),2x,2(1x,i5))
104 format(16(2x,f12.4))
      ! if (mod(ntime,100).EQ.0) then
      if (mod(time,tpin).LT.dt) then
         do inp=1,Nparticle
         vsi = vstart(inp) ; vei = vend(inp)
         esi = estart(inp) ; eei = eend(inp)
         fsi = fstart(inp) ; fei = fend(inp)
         call calculate_volume(Volume(inp),vsi,vei,fsi,fei,  &
                       xyz(:,vsi:vei),vert_of_face(:,fsi:fei))
         if (inp.EQ.1) then
         call calculate_volume_chamb(Volume_chamb(1:8),vsi,vei,fsi,fei,  &
                       xyz(:,vsi:vei),vert_of_face(:,fsi:fei),face_to_chamb4V(fsi:fei))
         endif
         enddo
      
         if(myid.eq.0) write(95,104)timeMS,pressLV,pressAO,pressLA,pressRV,pressAP,pressRA,Volume_chamb(1),Volume_chamb(2),Volume_chamb(3),Volume_chamb(4),&
              potEFnode(48475),& !Healthy
              potEFnode(43144),& !Scar core
              potEFnode(45278)   !Healtier Istmus
         if(myid.eq.0) flush(95)

!         stop
      endif


!->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
! ELLECTRO.OUT
!->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
#ifdef ELECTRO
      if (mod(time,tpin).LT.dt) then
      !calculate EF integrals
         integrale1=0.0D0
         do i=vstartEF_1d(3),vendEF_1d(3)
            integrale1=integrale1+potEFnode_1d(i)
         enddo
         integrale1=integrale1/nviEF_1d(3)

         integrale2=0.0D0
         do i=vstartEF_2d(1),vendEF_2d(1)
            integrale2=integrale2+potEFnode_2d(i)
         enddo
         integrale2=integrale2/nviEF_2d(1)

         !write(*,*) integrale1,integrale2
         integrale3=0.0D0
         do i=cstart_3d(1),cend_3d(1)
            integrale3=integrale3+potEFcell_3d(i)
         enddo
         integrale3=integrale3/nci_3d(1)

          if (myid .eq. 0) then
              write(*,*) "EF int",integrale1,integrale2,integrale3
!              write(*,*) "EF int",integrale3
          endif

         timeMS=time*(TSTAR*1000.d0)
105 format(12(2x,f12.8))
         !if(myid.eq.0) write(195,104)timeMS,integrale1,integrale2,integrale3!,potEFnode_1d(vAVmaster(1)),potEFnode_1d(vAVslave(1))  &
         if(myid.eq.0) write(195,104) integrale3!,potEFnode_1d(vAVmaster(1)),potEFnode_1d(vAVslave(1))  &

         !                                    ,potEFnode_1d(vAVmaster(2)),potEFnode_1d(vAVslave(2)),potEFnode_1d(vAVmaster(3)),potEFnode_1d(vAVslave(3)),potEFface_2d(fcheckEF_2d(1)),potEFface_2d(fcheckEF_2d(2))
         if(myid.eq.0) flush(195)

         !         stop

#ifdef EGM
204 format(1559(2x,f12.4)) 
!          !Compute EGMs for CARTO --> fileID==200
!          EGMvec(:)=0.0D0
!          !j=100
!          !$cuf kernel do(2)
!          EGMsum=0.0D0
!          do j=1,Negm
!             do i=1,nctot_3d
            
!                Delta_egm1 = (xyz_egm(1,j)-cell_bar(1,i))
!                Delta_egm2 = (xyz_egm(2,j)-cell_bar(2,i))
!                Delta_egm3 = (xyz_egm(3,j)-cell_bar(3,i))
!                invdist=1.0D0/sqrt(Delta_egm1*Delta_egm1+Delta_egm2*Delta_egm2+&
!                     Delta_egm3*Delta_egm3)
!                invdist3=invdist*invdist*invdist
! #ifdef ELEGEO0
!                volcell = vol0_3d(i)
! #else
!                volcell = vol_3d(i)
! #endif
!                EGMadd = 0.25D0*volcell*CARTO_Dcell3d(i)*(gradcell_3d(1,i)*Delta_egm1+&
!                     gradcell_3d(2,i)*Delta_egm2+gradcell_3d(3,i)*Delta_egm3)*invdist3/(pi*0.6667D0) !here 0.6667 is the blood magnetic conductivity
!                EGMvec(j) = EGMvec(j)+EGMadd
!                !aa=atomicadd(EGMvec(j),EGMadd)
!                !EGMvec(j)=EGMsum
!             enddo
!             !            !@cuf istat = cudaDeviceSynchronize !JDR TMP
           
!          enddo
!          !@cuf istat = cudaDeviceSynchronize !JDR TMP
         
!          if(myid.eq.0) write(200,204) (EGMvec(iii),iii=1,Negm) !TODO --> Sistemare scrittura EGM su file di testo, definire variabili per calcolo EGM
!          !if(myid.eq.0) write(200,*) EGMvec(100)
!          if(myid.eq.0) flush(200)
            inp=1
            csi=cstart_3d(inp)
            cei=cend_3d(inp)
            !ATRI E VENTRICOLI 
            do j=1,Negm
               call calculate_EGM(EGMvec(j),csi,cei,xyz_egm(:,j),cell_bar(:,csi:cei), &
#ifdef ELEGEO0
                    vol0_3d(csi:cei),gradcell_3d(:,csi:cei), &
#else
                    vol_3d(csi:cei),gradcell_3d(:,csi:cei), &
#endif                    
                    cell_to_chamb_3dBou(csi:cei),CARTO_Dcell3d)
            enddo
            if(myid.eq.0) write(200,204)timeMS,EGMvec(1:Negm)
            if(myid.eq.0) flush(200)

#endif
      endif
#endif

      return                                                            
      end                                                               
