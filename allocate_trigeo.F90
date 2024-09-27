      subroutine allocate_trigeo_3d
      use param
      use mls_param
      use mpih
      use mpi_param, only: kstart, kend
      implicit none

      integer::i

      write(*,*) "ALL 3d"
      allocate(xyz0_3d(3,nvtot_3d),xyz_3d(3,nvtot_3d),xyzold_3d(3,nvtot_3d))
      allocate(xyzv_3d(3,nvtot_3d),xyza_3d(3,nvtot_3d))
      allocate(xyzv0_3d(3,nvtot_3d),xyza0_3d(3,nvtot_3d))
      allocate(LabelSmooth(nvtot_3d))
      allocate(LabelSmooth_cell(nctot_3d))
      ! allocate(LabelDelay(nvtot_3d))
      allocate(LabelStenosi(nctot_3d))
      allocate(Segments_node_3d(nvtot_3d))
#ifdef SCALINGTEST
      allocate(timesforscaling(20))
#endif
!STRONG
#if defined(STRONG) || defined(STRONGRV)
      allocate(xyzp_3d(3,nvtot_3d),xyzvp_3d(3,nvtot_3d))
      allocate(xyzak_3d(3,nvtot_3d),xyzakm1_3d(3,nvtot_3d),xyzauk_3d(3,nvtot_3d),xyzakp1_3d(3,nvtot_3d))
      allocate(xyzvk_3d(3,nvtot_3d),xyzvkm1_3d(3,nvtot_3d),xyzvuk_3d(3,nvtot_3d),xyzvkp1_3d(3,nvtot_3d))
      allocate(xyzk_3d(3,nvtot_3d),xyzkp1_3d(3,nvtot_3d))
#endif
      allocate(vert_to_part_3d(nvtot_3d))
      allocate(edge_to_part_3d(netot_3d))
      allocate(face_to_part_3d(nftot_3d))
      allocate(cell_to_part_3d(nctot_3d))
      allocate(vert_to_chamb_3d(nvtot_3d))
      allocate(vert_to_chamb_3dBou(nvtot_3d))
      allocate(vert_to_chamb_3dBouOld(nvtot_3d))
      allocate(vert_to_chamb_3d4V(nvtot_3d))
      allocate(face_to_chamb_3d(nftot_3d))
      allocate(edge_to_chamb_3d(netot_3d))
      allocate(cell_to_chamb_3d(nctot_3d))
      allocate(cell_to_chamb_3dBou(nctot_3d))
      allocate(n_edge_of_vert_3d(nvtot_3d)) !n_vert_of_vert
      allocate(n_cell_of_vert_3d(nvtot_3d))
      allocate(vert_of_cell_3d(4,nctot_3d))
      allocate(vert_of_vert_3d(max_n_edge_of_vert_3d,nvtot_3d))
      allocate(cell_of_vert_3d(max_n_edge_of_vert_3d,nvtot_3d))
      allocate(vert_of_edge_3d(2,netot_3d))
      allocate(n_edge_of_cell_3d(nctot_3d)) 
      allocate(n_cell_of_edge_3d(netot_3d))
      allocate(edge_of_cell_3d(6,nctot_3d))
      allocate(cell_of_edge_3d(max_n_edge_of_vert_3d,netot_3d))
      allocate(vert_of_face_3d(3,nftot_3d))
      allocate(edge_of_vert_3d(max_n_edge_of_vert_3d,nvtot_3d))
      allocate(count_edge_of_vert_3d(nvtot_3d))
      allocate(face_of_cell_3d(4,nctot_3d))
      allocate(cell_of_face_3d(2,nftot_3d))
      allocate(versCFface_3d(3,nftot_3d))
      allocate(distCFface_3d(nftot_3d))
      allocate(g1interpface_3d(nftot_3d))
      allocate(Surface0_3d(Nparticle_3d),Surface_3d(Nparticle_3d))
!and more..
      allocate(AFung_3d(netot_3d))
      allocate(StiffV_3d(nvtot_3d))
      allocate(aalpha0_3d(3,nftot_3d),aalpha_3d(3,nftot_3d))
      allocate(dist0_3d(netot_3d),dist_3d(netot_3d))
      allocate(vol0_3d(nctot_3d),vol_3d(nctot_3d))
      allocate(sur0_3d(nftot_3d),sur_3d(nftot_3d))
      allocate(cell_bar(3,nctot_3d))
      allocate(normalfaceofcells_3d(3,4,nctot_3d))
      allocate(AmatrFibers_cell_3d(3,3,nctot_3d))
      allocate(AmatrFibers_node_3d(3,3,nvtot_3d))
      allocate(AmatrFibers_edge_3d(3,3,netot_3d))
      allocate(edgefiber_cosangle_3d(netot_3d))
      allocate(Mintcells_3d(3,3,nctot_3d))
      allocate(Mextcells_3d(3,3,nctot_3d))
      allocate(Mintfaces_3d(3,3,nftot_3d))
      allocate(Mextfaces_3d(3,3,nftot_3d))
      allocate(IstimEF_3d(nctot_3d))
      !S1S2 stimuli
      allocate(IstimEF_3dS1(nctot_3d))
      allocate(IstimEF_3dS2(nctot_3d))
      allocate(mass_of_vert_3d(nvtot_3d))
      allocate(boundary1_3d(2,nvtot_3d))


      allocate(potEFcell_3d(nctot_3d))
      allocate(potEFnode_3d(nvtot_3d))
      allocate(potEFface_3d(nftot_3d))
      allocate(gradcell_3d(3,nctot_3d))
      allocate(gradface_3d(3,nftot_3d))
      allocate(LabelSkipConn(nctot_3d))
#ifdef ELECTRO

#ifdef TP06
      allocate(XEF_3d(21,nctot_3d))
#endif
#ifdef COURTEMANCHE
      allocate(XEF_3d(21,nctot_3d))
#endif
#ifdef MINIMAL_MODEL      
      allocate(XEF_3d(4,nctot_3d))
#endif
      !-----------------------------
      !APD
      allocate(t_apd_3d(4,nctot_3d))
      allocate(old_1(nctot_3d))
      allocate(f_apd(nctot_3d))
      !-----------------------------                                                         
      !Laplace_tagging
      allocate(myotag(nctot_3d),meshquality_3d(nctot_3d))
      !----------------------------- 
      ! CARTO variables
      !----------------------------- 
#ifdef CARTO
      allocate(CARTO_Dnode3d(nvtot_3d),CARTO_Dcell3d(nctot_3d),CARTO_Dface3d(nftot_3d))
#endif

#ifdef EGM
      open(unit=15,file='./meshes/EGM_CARTO.txt',status='old')
      read(15,*) Negm
      allocate(xyz_egm(3,Negm),EGMvec(Negm))
      do i=1,Negm
         read(15,*)xyz_egm(1,i),xyz_egm(2,i),xyz_egm(3,i)
      enddo
      close(15)
#endif
 
      
      !-----------------------------
      !Scar 3D
      allocate(scar_cell(nctot_3d))
#ifdef BIDOMAIN
      allocate(potextEFcell_3d(nctot_3d))
      allocate(potextEFnode_3d(nvtot_3d))
      allocate(potextEFface_3d(nftot_3d))
      allocate(gradextcell_3d(3,nctot_3d))
      allocate(gradextface_3d(3,nftot_3d))
      allocate(QEFbido(nctot_3d,maxitEFbido+1))
      allocate(HEFbido(maxitEFbido+1,maxitEFbido))
      allocate(checkEFbido(nctot_3d,3))
      allocate(IPIVEFbido(maxitEFbido+1))
#endif
      
#endif
     ! allocate(boundary2_3d(2,nvtot_3d))
     ! allocate(BCoffset_3d(3,nvtot_3d))
      allocate(fpxyz_3d(3,nvtot_3d),fxyz_3d(3,nvtot_3d))

#ifdef STRONGRV
      allocate(fxyzm_3d(3,nvtot_3d))
#endif
      allocate(astressEFnode_3d(nvtot_3d))
      allocate(astressEFcell_3d(nctot_3d))
      allocate(astressEFedge_3d(netot_3d))
      allocate(EFtstart_3d(nvtot_3d))
      allocate(ke_3d(netot_3d),f_act_3d(netot_3d))
      return
      end
!=====================================================================
!=====================================================================
      subroutine allocate_trigeo
      use param
      use mls_param
      use mpih
      use mpi_param, only: kstart, kend
      implicit none 
      integer :: merr

      !-------------------------------------------------
      !-------------------------------------------------
      write(*,*) "ALL 2d"
      !-------------------------------------------------
      allocate(n_edge_of_vert(nvtot))
      allocate(vert_of_edge(2,netot))
      allocate(face_of_edge(2,netot))
      allocate(vert_of_face(3,nftot))
      allocate(edge_of_face(3,nftot))
      allocate(vert_of_vert(max_n_edge_of_vert,nvtot))
      allocate(edge_of_vert(max_n_edge_of_vert,nvtot))
      allocate(v1234(4,netot))
      allocate(boundary1(2,nvtot))
      ! allocate(boundary2(2,nvtot))
      ! allocate(BCoffset(3,nvtot))
      allocate(vert_to_part(nvtot))
      allocate(face_to_part(nftot))
      allocate(edge_to_part(netot))
      allocate(vert_to_chamb(nvtot))
      allocate(vert_to_chamb4V(nvtot))
      allocate(edge_to_chamb(netot))
      allocate(face_to_chamb(nftot))
      allocate(face_to_chamb4V(nftot))
!      allocate(LabelH(nftot))
      ! allocate(LabelBound(nvtot))
      !-------------------------------------------------
      allocate(bboxind(3,2,Nparticle))
      allocate(dum_for(nftot))
      allocate(dum_forv(nvtot))
      allocate(pind(6,nftot),pind_probeP(6,nftot),pind_probeN(6,nftot))
      allocate(pindr(6,nftot))
      allocate(pindv(6,nvtot))
      allocate(mvol(nftot))
      allocate(fcol(nftot))
      allocate(fcolv(nvtot))
      !-------------------------------------------------
      allocate(theta0(netot),theta(netot))
      allocate(dist00(netot),dist0(netot),dist(netot))
      allocate(sur0(nftot),sur(nftot))
      allocate(aalpha0(3,nftot),aalpha(3,nftot))
      !-------------------------------------------------
      allocate(dismax(maxnf,Nparticle),meanxyz(3,Nparticle),meanvel(3,Nparticle))

      allocate(Volume0(Nparticle),Volume(Nparticle))
      allocate(Surface0(Nparticle),Surface(Nparticle))
      allocate(Volume_chamb0(8),Volume_chamb(8),Volume_chambShift(8))
      allocate(Surface_chamb0(8),Surface_chamb(8),Surface_chambShift(8))
      ! allocate(Hboxx(0:n3),celvol(0:n3))
      ! allocate(Hboxxr(0:n3r),celvolr(0:n3r))
      allocate(shwtx(7),shwty(7),shwtz1(7),shwtz2(7))
      !-------------------------------------------------

      !-------------------------------------------------
      allocate(xyz0(3,nvtot),xyz(3,nvtot),xyzold(3,nvtot))
      allocate(xyzv(3,nvtot),xyza(3,nvtot))
      allocate(xyzv0(3,nvtot),xyza0(3,nvtot))

!STRONG
#if defined(STRONG) || defined(STRONGRV)
      allocate(xyzp(3,nvtot),xyzvp(3,nvtot))
      allocate(xyzak(3,nvtot),xyzakm1(3,nvtot),xyzauk(3,nvtot),xyzakp1(3,nvtot))
      allocate(xyzvk(3,nvtot),xyzvkm1(3,nvtot),xyzvuk(3,nvtot),xyzvkp1(3,nvtot))
      allocate(xyzk(3,nvtot),xyzkp1(3,nvtot))
#endif
      allocate(tri_ver(9,nftot),tri_vel(9,nftot))
      allocate(vel_tri(3,nftot),acc_tri(3,nftot))
      allocate(tri_bar(3,nftot),tri_nor(3,nftot))
      allocate(trcnt(1))
      allocate(mytr(nftot))
      !-------------------------------------------------

      !-------------------------------------------------
      allocate(fpxyz(3,nvtot))
      allocate(fxyz(3,nvtot))
      allocate(tauface(nftot),pressface(nftot))
      allocate(taufaceAV(nftot),pressfaceAV(nftot))
#ifdef STRONGRV
      allocate(fxyzm(3,nvtot))
#endif
      allocate(mass_of_vert(nvtot))
      allocate(AFung(netot))
      allocate(StiffV(nvtot))
      !-------------------------------------------------

      allocate(kb(netot), ke(netot))

 
      return
      end 


!=====================================================================                                                                               
      subroutine allocate_trigeoEF_2d
      use param
      use mls_param
      use mpih
      implicit none
      write(*,*) "ALL EF2d"
      !-------------------------------------------------                                                                                             
      allocate(n_edge_of_vertEF_2d(nvtotEF_2d))
      allocate(vert_of_edgeEF_2d(2,netotEF_2d))
      allocate(face_of_edgeEF_2d(2,netotEF_2d))
      allocate(vert_of_faceEF_2d(3,nftotEF_2d))
      allocate(edge_of_faceEF_2d(3,nftotEF_2d))
      allocate(vert_of_vertEF_2d(max_n_edge_of_vert,nvtotEF_2d))
      allocate(edge_of_vertEF_2d(max_n_edge_of_vert,nvtotEF_2d))
      allocate(face_of_vertEF_2d(max_n_edge_of_vert,nvtotEF_2d))
      allocate(n_face_of_vertEF_2d(nvtotEF_2d))
      write(*,*)"HELP1"
      allocate(dist0EF_2d(netotEF_2d),distEF_2d(netotEF_2d))
      allocate(sur0EF_2d(nftotEF_2d),surEF_2d(nftotEF_2d))
      !-------------------------------------------------                                                                                             
      allocate(xyz0EF_2d(3,nvtotEF_2d),xyzEF_2d(3,nvtotEF_2d),xyzSEF_2d(3,nvtotEF_2d),xyzS1EF_2d(3,nvtotEF_2d))
      allocate(tri_barEF_2d(3,nftotEF_2d),tri_norEF_2d(3,nftotEF_2d))
      allocate(Surface0EF_2d(NparticleEF_2d),SurfaceEF_2d(NparticleEF_2d))
      allocate(versCFedgeEF_2d(3,netotEF_2d))
      allocate(distCFedgeEF_2d(netotEF_2d))
      allocate(g1interpedgeEF_2d(netotEF_2d))
      allocate(normaledgeoffacesEF_2d(3,3,nftotEF_2d))
      allocate(AmatrFibersEF_2d(3,3,nftotEF_2d))
#ifdef STEWART09
      allocate(XEF_2d(20,nftotEF_2d))
#endif
      write(*,*)"HELP2"
#ifdef MV_PURK
      allocate(XEF_2d(4,nftotEF_2d))
#endif      
      ! allocate(Istim0EF_2d(nftotEF_2d))
      ! allocate(IstimEF_2d(nftotEF_2d))
      ! allocate(IstimLV0_2d(nftotEF_2d))
      ! allocate(IstimLV_2d(nftotEF_2d))
      ! allocate(IstimRV0_2d(nftotEF_2d))
      ! allocate(IstimRV_2d(nftotEF_2d))
      write(*,*)"HELP3"
      !S1S2 Stimuli
      !---------------------------------
      allocate(IstimEF_2dS1(nftotEF_2d))
      allocate(IstimEF_2dS2(nftotEF_2d))
      !--------------------------------- 
      ! CARTO variables
      !---------------------------------
      allocate(scar_tag_2d(nvtotEF_2d),CARTO_Dnode(nvtotEF_2d),CARTO_Dface(nftotEF_2d))
      write(*,*)"HELP4"
      !---------------------------------
      ! APD 2d
      !---------------------------------
      allocate(t_apd_2d(4,nvtotEF_2d),old_1_2d(nvtotEF_2d),f_apd_2d(nvtotEF_2d))
      allocate(CVface_EF_2d(nftotEF_2d))
      !---------------------------------
      allocate(MintfacesEF_2d(3,3,nftotEF_2d))
      allocate(MextfacesEF_2d(3,3,nftotEF_2d))
      allocate(MintedgesEF_2d(3,3,netotEF_2d))
      allocate(MextedgesEF_2d(3,3,netotEF_2d))
      allocate(gradfaceEF_2d(3,nftotEF_2d))
      allocate(gradedgeEF_2d(3,netotEF_2d))
      allocate(potEFnode_2d(nvtotEF_2d))
      allocate(potEFedge_2d(netotEF_2d))
      allocate(potEFface_2d(nftotEF_2d))
      allocate(LabelPurkSetto(nftotEF_2d))
      allocate(EFtstart_2d(nvtotEF_2d))
      write(*,*)"HELP5"
      return
      end
!=====================================================================                                                                               
!=====================================================================                                                                               
      subroutine allocate_trigeoEF_1d
      use param
      use mls_param
      use mpih
      implicit none
      write(*,*) "ALL EF1d"
      !-------------------------------------------------                                                                                             
      allocate(n_edge_of_vertEF_1d(nvtotEF_1d))
      allocate(vert_of_edgeEF_1d(2,netotEF_1d))
      allocate(vert_of_vertEF_1d(max_n_edge_of_vert_1d,nvtotEF_1d))
      allocate(edge_of_vertEF_1d(max_n_edge_of_vert_1d,nvtotEF_1d))
      allocate(dist0EF_1d(netotEF_1d),distEF_1d(netotEF_1d))
      allocate(edg_barEF_1d(3,netotEF_1d))
      !-------------------------------------------------                                                                                             
      allocate(xyz0EF_1d(3,nvtotEF_1d),xyzEF_1d(3,nvtotEF_1d),xyzSEF_1d(3,nvtotEF_1d))
#ifdef COURTEMANCHE
      allocate(XEF_1d(21,nvtot_3d))
#endif
#ifdef MINIMAL_MODEL      
      allocate(XEF_1d(4,nvtot_3d))
#endif
!      allocate(Istim0EF_1d(nvtotEF_1d))
      allocate(IstimEF_1d(nvtotEF_1d))
      allocate(MintnodesEF_1d(nvtotEF_1d))
!      allocate(MextnodesEF_1d(nvtotEF_1d))                                                                                                          
      allocate(MintedgesEF_1d(netotEF_1d))
!      allocate(MextedgesEF_1d(netotEF_1d))                                                                                                          
      allocate(potEFnode_1d(nvtotEF_1d))
      allocate(vert_to_partEF_1d(nvtotEF_1d))
      allocate(edge_to_partEF_1d(netotEF_1d))
      allocate(EFtstart_1d(nvtotEF_1d))
      return
      end
!=====================================================================  


!=====================================================================                                                                               
      subroutine allocate_trigeoCV_1d
      use param
      use mls_param
      use mpih
      implicit none

      !-------------------------------------------------                                     !CERCO DI DEFINIRNE UN PO MENO                                                        
      allocate(n_edge_of_vertCV_1d(nvtotCV_1d))
      allocate(vert_of_edgeCV_1d(2,netotCV_1d))
      allocate(vert_of_vertCV_1d(max_n_edge_of_vert_1d,nvtotCV_1d))
      allocate(edge_of_vertCV_1d(max_n_edge_of_vert_1d,nvtotCV_1d))
      ! allocate(dist0CV_1d(netotCV_1d),distCV_1d(netotCV_1d))
      ! allocate(edg_barCV_1d(3,netotCV_1d))
      !-------------------------------------------------                                                                                             
      allocate(xyz0CV_1d(3,nvtotCV_1d),xyzCV_1d(3,nvtotCV_1d),xyzSCV_1d(3,nvtotCV_1d))
      ! allocate(vert_to_partCV_1d(nvtotCV_1d))
      ! allocate(edge_to_partCV_1d(netotCV_1d))

      return
      end
!=====================================================================  
