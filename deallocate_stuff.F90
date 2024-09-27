      subroutine allocate_stuffPpro
      use param
      use mls_param
      use mpih
      use mpi_param, only: kstart, kend
      use ibm_param
      implicit none 

      mpunPpro = n1*n2*n3
      mpugeoPpro = n1*n2*n3

      allocate(indgeoPpro(4,mpunPpro,3),indgeoePPro(4,mpunPpro,3),indgeoeePpro(4,mpunPpro,3),distbPpro(4,mpunPpro),ntribfxPpro(4,mpunPpro))
      allocate(vbifxPpro(4,mpunPpro),vbfxPpro(4,mpunPpro))

      allocate(xyzbfxPpro(mpugeoPpro,9),areafxPpro(mpugeoPpro),barfxPpro(mpugeoPpro,3),nxyzfxPpro(mpugeoPpro,3))

      return
      end
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      subroutine deallocate_stuffPpro
      use param
      use mls_param
      use mpih
      use mpi_param, only: kstart, kend
      use ibm_param
      implicit none 

!      deallocate(indgeoPpro(4,mpunPpro,3),indgeoePPro(4,mpunPpro,3),indgeoeePpro(4,mpunPpro,3),distbPpro(4,mpunPpro))
      deallocate(indgeoPpro,indgeoePPro,indgeoeePpro,distbPpro)

      return
      end
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      subroutine deallocate_stuff
      use param
      use mls_param
      use mpih
      use mpi_param, only: kstart, kend
      use tile_arrays
      implicit none 

!da gcurv 3d
      deallocate(nvi_3d,nci_3d)
      deallocate(nei_3d,nfi_3d) !questi quando li legge? in read_geo_dim_3d
      deallocate(vstart_3d,vend_3d)
      deallocate(estart_3d,eend_3d)
      deallocate(fstart_3d,fend_3d)
      deallocate(cstart_3d,cend_3d)

!da gcurv 2d
      deallocate(nvi,nei,nfi)
      deallocate(vstart,vend)
      deallocate(estart,eend)
      deallocate(fstart,fend)

      deallocate(LabelSmooth)
      deallocate(LabelSmooth_cell)
      deallocate(LabelStenosi)
      deallocate(Segments_node_3d)
!3d
      deallocate(xyz0_3d,xyz_3d,xyzold_3d)
      deallocate(xyzv_3d,xyza_3d)
      deallocate(xyzv0_3d,xyza0_3d)
!STRONG
#if defined(STRONG) || defined(STRONGRV)
      deallocate(xyzp_3d,xyzvp_3d)
      deallocate(xyzak_3d,xyzakm1_3d,xyzauk_3d,xyzakp1_3d)
      deallocate(xyzvk_3d,xyzvkm1_3d,xyzvuk_3d,xyzvkp1_3d)
      deallocate(xyzk_3d3,nvtot_3d,xyzkp1_3d)
#endif
      deallocate(vert_to_part_3d)
      deallocate(face_to_part_3d)
      deallocate(edge_to_part_3d)
      deallocate(cell_to_part_3d)
      deallocate(vert_to_chamb_3d)
      deallocate(vert_to_chamb_3d4V)
      deallocate(vert_to_chamb_3dBou)
      deallocate(vert_to_chamb_3dBouOld)
      deallocate(cell_to_chamb_3dBou)     
      deallocate(edge_to_chamb_3d)
      deallocate(face_to_chamb_3d)
      deallocate(cell_to_chamb_3d)

      deallocate(n_edge_of_vert_3d) !n_vert_of_vert
      deallocate(n_cell_of_vert_3d)
      deallocate(vert_of_cell_3d)
      deallocate(vert_of_vert_3d)
      deallocate(cell_of_vert_3d)

      deallocate(vert_of_edge_3d)
      deallocate(n_edge_of_cell_3d)
      deallocate(n_cell_of_edge_3d)
      deallocate(edge_of_cell_3d)
      deallocate(cell_of_edge_3d)
      deallocate(vert_of_face_3d)
      deallocate(edge_of_vert_3d)
      deallocate(count_edge_of_vert_3d)
      deallocate(face_of_cell_3d)
      deallocate(cell_of_face_3d)
      deallocate(versCFface_3d)
      deallocate(distCFface_3d)
      deallocate(g1interpface_3d)
      deallocate(Surface0_3d)
      deallocate(Surface_3d)
!and more..
      deallocate(AFung_3d)
      deallocate(StiffV_3d)

      deallocate(aalpha0_3d,aalpha_3d)
      deallocate(dist0_3d,dist_3d)
      deallocate(vol0_3d,vol_3d)
      deallocate(sur0_3d,sur_3d)
      deallocate(cell_bar)
      deallocate(normalfaceofcells_3d)
      deallocate(AmatrFibers_cell_3d)
      deallocate(AmatrFibers_node_3d)
      deallocate(AmatrFibers_edge_3d)
      deallocate(edgefiber_cosangle_3d)
      deallocate(Mintcells_3d)
      deallocate(Mextcells_3d)
      deallocate(Mintfaces_3d)
      deallocate(Mextfaces_3d)
      deallocate(IstimEF_3d)
      deallocate(IstimEF_3dS1)
      deallocate(IstimEF_3dS2)
!      deallocate(IstimEF_2d)
      ! deallocate(IstimEF_3d)
      deallocate(mass_of_vert_3d)
!      deallocate(boundary1_3d)

     deallocate(boundary2_3d)
     deallocate(BCoffset_3d)


      deallocate(fpxyz_3d,fxyz_3d)
#ifdef STRONGRV
      deallocate(fxyzm_3d)
#endif
      deallocate(EFtstart_3d)
      deallocate(astressEFnode_3d)
      deallocate(astressEFcell_3d)
      deallocate(astressEFedge_3d)
      deallocate(f_act_3d)

      
      deallocate(ke_3d)
#ifdef ELECTRO
      deallocate(XEF_3d)
      deallocate(potEFcell_3d)
      deallocate(potEFnode_3d)
      deallocate(potEFface_3d)
      deallocate(gradcell_3d)
      deallocate(gradface_3d)
      deallocate(LabelSkipConn)

#ifdef BIDOMAIN      
      deallocate(potextEFcell_3d)
      deallocate(potextEFnode_3d)
      deallocate(potextEFface_3d)
      deallocate(gradextcell_3d)
      deallocate(gradextface_3d)
      deallocate(QEFbido)
      deallocate(HEFbido)
      deallocate(checkEFbido)
      deallocate(IPIVEFbido)
#endif
      
      deallocate(nviEF_1d,neiEF_1d)
      deallocate(vstartEF_1d,vendEF_1d)
      deallocate(estartEF_1d,eendEF_1d)

      deallocate(boundary2EF_1d)
      deallocate(BCoffsetEF_1d)
      deallocate(vAVmaster,vAVslave)

      deallocate(n_edge_of_vertEF_1d(nvtotEF_1d))
      deallocate(vert_of_edgeEF_1d(2,netotEF_1d))
      deallocate(vert_of_vertEF_1d(max_n_edge_of_vert_1d,nvtotEF_1d))
      deallocate(edge_of_vertEF_1d(max_n_edge_of_vert_1d,nvtotEF_1d))
      deallocate(dist0EF_1d(netotEF_1d),distEF_1d(netotEF_1d))
      deallocate(edg_barEF_1d(3,netotEF_1d))
      !-------------------------------------------------                                                                             
      deallocate(xyz0EF_1d,xyzEF_1d,xyzSEF_1d)
      deallocate(XEF_1d)
!      deallocate(Istim0EF_1d)
      deallocate(IstimEF_1d)
      deallocate(MintnodesEF_1d)
!      deallocate(MextnodesEF_1d)
      deallocate(MintedgesEF_1d)  
!      deallocate(MextedgesEF_1d)
      deallocate(potEFnode_1d)
      deallocate(vert_to_partEF_1d)
      deallocate(edge_to_partEF_1d)
      deallocate(EFtstart_1d)

      deallocate(nviEF_2d,neiEF_2d,nfiEF_2d)
      deallocate(vstartEF_2d,vendEF_2d)
      deallocate(estartEF_2d,eendEF_2d)
      deallocate(fstartEF_2d,fendEF_2d)

      deallocate(boundary2EF_2d)
      deallocate(BCoffsetEF_2d)
      deallocate(n_edge_of_vertEF_2d)
      deallocate(vert_of_edgeEF_2d)
      deallocate(face_of_edgeEF_2d)
      deallocate(vert_of_faceEF_2d)
      deallocate(edge_of_faceEF_2d)
      deallocate(vert_of_vertEF_2d)
      deallocate(edge_of_vertEF_2d)
      deallocate(face_of_vertEF_2d)
      deallocate(n_face_of_vertEF_2d)

      deallocate(dist0EF_2d,distEF_2d)
      deallocate(sur0EF_2d,surEF_2d)
      !-------------------------------------------------                                                                                             
      deallocate(xyz0EF_2d,xyzEF_2d,xyzSEF_2d,xyzS1EF_2d)
      deallocate(tri_barEF_2d,tri_norEF_2d)
      deallocate(Surface0EF_2d,SurfaceEF_2d)
      deallocate(versCFedgeEF_2d)
      deallocate(distCFedgeEF_2d)
      deallocate(g1interpedgeEF_2d)
      deallocate(normaledgeoffacesEF_2d)
      deallocate(AmatrFibersEF_2d)
      deallocate(XEF_2d)

      deallocate(MintfacesEF_2d)
      deallocate(MextfacesEF_2d)
      deallocate(MintedgesEF_2d)
      deallocate(MextedgesEF_2d)
      deallocate(gradfaceEF_2d)
      deallocate(gradedgeEF_2d)
      deallocate(potEFnode_2d)
      deallocate(potEFedge_2d)
      deallocate(potEFface_2d)
      deallocate(LabelPurkSetto)
      deallocate(EFtstart_2d)

      deallocate(vert_boundEF_2d)
      deallocate(vert_of_vert_boundEF_2d)
      deallocate(vPurkmaster)
      deallocate(fPurkslave)
      deallocate(vAtrimaster)
      deallocate(cAtrislave)
      deallocate(Bund_tstim)
      deallocate(fVentrmaster)
      deallocate(cVentrslave)
      deallocate(Purk_tstim)
      deallocate(fcheckEF_2d)


      deallocate(nviCV_1d,neiCV_1d)
      deallocate(vstartCV_1d,vendCV_1d)
      deallocate(estartCV_1d,eendCV_1d)

      deallocate(boundary2CV_1d)
      deallocate(BCoffsetCV_1d)
      deallocate(n_edge_of_vertCV_1d)
      deallocate(vert_of_edgeCV_1d)
      deallocate(vert_of_vertCV_1d)
      deallocate(edge_of_vertCV_1d)
      ! deallocate(dist0CV_1d,distCV_1d)
      ! deallocate(edg_barCV_1d)
      !-------------------------------------------------                               
      deallocate(xyz0CV_1d,xyzCV_1d,xyzSCV_1d)
      ! deallocate(vert_to_partCV_1d)
      ! deallocate(edge_to_partCV_1d)

#endif

!2D
      !-------------------------------------------------
      deallocate(n_edge_of_vert)
      deallocate(vert_of_edge)
      deallocate(face_of_edge)
      deallocate(vert_of_face)
      deallocate(edge_of_face)
      deallocate(vert_of_vert)
      deallocate(edge_of_vert)
      deallocate(v1234)
      deallocate(boundary2)
      deallocate(BCoffset)
      deallocate(vert_to_part)
      deallocate(face_to_part)
      deallocate(edge_to_part)
      deallocate(vert_to_chamb)
      deallocate(vert_to_chamb4V)
      deallocate(edge_to_chamb)
      deallocate(face_to_chamb)
      deallocate(face_to_chamb4V)

!      deallocate(LabelH)
      ! deallocate(body_nsl)
      ! deallocate(body_flux)
      ! deallocate(LabelBound)
      !-------------------------------------------------
      !-------------------------------------------------
      deallocate(bboxind)
      deallocate(dum_for)
      deallocate(dum_forv)
      deallocate(pind,pind_probeP,pind_probeN)
      deallocate(pindr)
      deallocate(pindv)
      deallocate(xyzMLSprobe,outMLSprobe,pindMLSprobe)
      deallocate(tpcnt,mypr)
      deallocate(mvol)
      deallocate(fcol)
      deallocate(fcolv)
      !-------------------------------------------------
      deallocate(theta0,theta)
      deallocate(dist00,dist0,dist)
      deallocate(sur0,sur)
      deallocate(aalpha0,aalpha)
      !-------------------------------------------------

      deallocate(dismax,meanxyz,meanvel)

      deallocate(Volume0,Volume)
      deallocate(Surface0,Surface)
      deallocate(Volume_chamb0,Volume_chamb,Volume_chambShift)
      deallocate(Surface_chamb0,Surface_chamb,Surface_chambShift)

      deallocate(Hboxx,celvol)
      deallocate(Hboxxr,celvolr)
      deallocate(shwtx,shwty,shwtz1,shwtz2)
      !-------------------------------------------------

      !-------------------------------------------------
      deallocate(xyz0,xyz,xyzold)
      deallocate(xyzv,xyza)
      deallocate(xyzv0,xyza0)
!STRONG
#if defined(STRONG) || defined(STRONGRV)
      deallocate(xyzp,xyzvp(3,nvtot))
      deallocate(xyzak,xyzakm1,xyzauk,xyzakp1)
      deallocate(xyzvk,xyzvkm1,xyzvuk,xyzvkp1)
      deallocate(xyzk,xyzkp1)
#endif

      deallocate(tri_ver,tri_vel)
      deallocate(vel_tri,acc_tri)
      deallocate(tri_bar,tri_nor)
      deallocate(trcnt)
      deallocate(tpcnt)
      deallocate(mytr)
      deallocate(mypr)
      !-------------------------------------------------
      !-------------------------------------------------
      deallocate(fpxyz)
      deallocate(fxyz)
      deallocate(tauface,pressface)
      deallocate(taufaceAV,pressfaceAV)
#ifdef STRONGRV
      deallocate(fxyzm)
#endif
      deallocate(mass_of_vert)
      deallocate(AFung)
      deallocate(StiffV)
      !-------------------------------------------------

      deallocate(kb, ke)
      deallocate(tag_2dwet)
!da SetTiling
      deallocate(tri_tiling)
      deallocate(albegaBar)
      deallocate(ntilei,tilestart,tileend)
      deallocate(faceid_t) ! face id of a specific tile      
      deallocate(pindt)
      deallocate(tstart)
      deallocate(tend)
      deallocate(tlcnt)
      deallocate(mytl)
#ifdef OFFSETBODY       
      deallocate(pindtOff)
#endif
      
#ifndef SOLOUNO
     deallocate(LabelVaortic)
     deallocate(LabelVmitral)
     deallocate(LabelVpulmo)
     deallocate(LabelVtricu)

#endif 


      return
      end 

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      subroutine deallocate_mls_local
      use param
      use mls_local
      use mpih
      use outflow_vars
      use mpi_param, only: kstart, kend, kstartr, kendr
      implicit none 
      integer :: merr
      !-------------------------------------------------
      deallocate(for_xc,stat=merr)
      
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory: forxc"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      !-------------------------------------------------
      deallocate(for_yc,stat=merr)

                                                                 
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory: foryc"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      !-------------------------------------------------
      deallocate(for_zc,stat=merr)
      
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory: forzc"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      !-------------------------------------------------
      deallocate(for_te,stat=merr)
                                                                 
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory: forte"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 

      !-------------------------------------------------
      deallocate(coll,stat=merr)
      
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory: coll"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      !-------------------------------------------------
      deallocate(for_sc,stat=merr)
      
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory: forsc"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 

!     Allocate outlfow variables
      deallocate(qb1s,qb1n)
      deallocate(qb2s,qb2n)
      deallocate(qb3s,qb3n)

      deallocate(dqb1s,dqb1n)
      deallocate(dqb2s,dqb2n)
      deallocate(dqb3s,dqb3n)

      deallocate(dq1x2o,dq2x2o,dq3x2o)

      return
      end 


