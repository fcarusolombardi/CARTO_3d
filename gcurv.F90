!************************************************************************
      subroutine gcurv
      use mpih
      use mpi_param,only: kstart,kend,kstartr,kendr
      use param
      use local_arrays
      use mls_param
      use mgrd_arrays
      use ibm_param
      use stat_arrays, only: timeint_cdsp
      use nvtx
      use probes

      implicit none
!@cuf integer :: istat
      integer :: ntstf,l
      real(DP)    :: cflm,cflmr
      real(DP)    :: dmaxc,dmaxr
      real(DP)    :: aengc,aengcQ1,aengcQ2,aengcQ3,aengr
      real(DP)    :: ti(2), tin(3)
      character*100 ipfi

      integer :: i,m,n,inp,ii,itime,chamb
      integer :: v1,v2,v3,contaste
      integer :: nvm,nem,nfm
      integer :: nvcum,nfcum,nBeatSwitch,nPPro
      character*50 :: outfileXML
      character*10 ::NumVer
      character*30 :: xCoord, yCoord, zCoord
      character*70 filname
      real(DP) :: facemass,totmass,nPProR
      real(DP) :: A1, A2,xV,yV,zV,xV1,yV1,zV1,rV1
      real(DP) :: OH,HP,OP,A,B,C,sogliola
      real(DP) :: x0,y0,z0,xH,yH,zH,xP,yP,zP

!
!     Code for the computation of three-dimensional incompressible flows    
!     in cartesian coordinates.                                             
!                                                                       
!     This code solves flow fields bounded in the x3 (zc). 
!     All boundaries can be no-slip or free-slip
!     by setting the appropriate indices in the input file.
!     The geometry  (a cylindrical can) is given in the subroutine cordi.
!     The discretization is uniform in the x1 (xc) and x2 (xc) 
!     directions, while can be non-uniform in the x3 (zc) direction.
!
!     The equations for the following variables
!                                                                       
!      q1=v(xc)    q2=v(yc)     q3=v(zc)  
!                                                                       
!     are discretized by finite-difference schemes.                    
!
!     All spatial derivatives are discretized by central second-order
!     accurate finite-difference schemes including the non linear terms.                                   
!
!     The nonlinear terms are treated explicitly, while the viscous terms
!     are computed implicitly. This would lead to the inversion of a
!     large banded matrix, that however is avoided by introducing
!     a factored scheme bringing to the solution of three tridiagonal
!     matrices for each velocity component (subroutine INVTR*).
!                              
!     In time a fractional-step procedure is used in the version of 
!     Nagi Mansour introducing the pressure in the first step.                         
!
!     The scalar quantity Phi, which projects the provisional velocity field
!     onto  a divergence free field, is solved by a direct method. 
!     For the x1 (xc) and x2 (yc) directions modified wave numbers coupled 
!     with trigonometric expansions (FFTs) are used. The equation is then
!     solved by simply inverting a tridiagonal matrix for the x3 (zc) direction.
!     No explicit boundary conditions are necessary for this Poisson equation.      
!                                                                       
!     Other details of the scheme are given in the introduction of the  
!     subroutine TSCHEM
!
!     timings
!
      tin(1) = MPI_WTIME()
!========================================================
!  matrix initialization
#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0)write(*,*)myid,'starting initia'
#endif
      if ((Ppro.LE.0).OR.(iread.EQ.1)) then
         call initia
      endif



if (imlsfor.EQ.1) then      
!========================================================
!  initialising tetrahedral geometries (3D)
!     initial allocations 
      allocate(nvi_3d(Nparticle_3d),nci_3d(Nparticle_3d))
      allocate(nei_3d(Nparticle_3d),nfi_3d(Nparticle_3d)) !questi quando li legge? in read_geo_dim_3d
      allocate(vstart_3d(Nparticle_3d),vend_3d(Nparticle_3d))
      allocate(estart_3d(Nparticle_3d),eend_3d(Nparticle_3d))
      allocate(fstart_3d(Nparticle_3d),fend_3d(Nparticle_3d))
      allocate(cstart_3d(Nparticle_3d),cend_3d(Nparticle_3d))

      nvtot_3d = 0 ; netot_3d =0; nftot_3d = 0; nctot_3d = 0;
    
      do inp=1,Nparticle_3d

      call read_geo_dim_3dVTK(nvi_3d(inp),nei_3d(inp),nfi_3d(inp),nci_3d(inp),geofile_3d(inp))

      if(myid.eq.0)then
        write(*,*)'Dim. Geo 3d: ',geofile_3d(inp)
        write(*,*)'Nodes,Edges,Faces,Cells: ',nvi_3d(inp),nei_3d(inp),nfi_3d(inp),nci_3d(inp)
      end if

      nvtot_3d = nvtot_3d + nvi_3d(inp)
      netot_3d = netot_3d + nei_3d(inp)
      nftot_3d = nftot_3d + nfi_3d(inp)
      nctot_3d = nctot_3d + nci_3d(inp)
      write(*,*)'NodesT,EdgesT,FacesT,CellsT: ',nvtot_3d,netot_3d,nftot_3d,nctot_3d
      end do

!     initialising vstart and vend for inp=1 and Nparticle
      vstart_3d(1) = 1  ;  estart_3d(1) = 1
      fstart_3d(1) = 1  ;  cstart_3d(1) = 1 
      vend_3d(1) = nvi_3d(1); eend_3d(1) = nei_3d(1)
      fend_3d(1)=nfi_3d(1)  ; cend_3d(1) = nci_3d(1)

      do inp=2,Nparticle_3d
        vstart_3d(inp) = vend_3d(inp-1)+1
        vend_3d(inp)   = vstart_3d(inp) + nvi_3d(inp)-1

        estart_3d(inp) = eend_3d(inp-1)+1
        eend_3d(inp)   = estart_3d(inp) + nei_3d(inp)-1

        fstart_3d(inp) = fend_3d(inp-1)+1
        fend_3d(inp)   = fstart_3d(inp) + nfi_3d(inp)-1

        cstart_3d(inp) = cend_3d(inp-1)+1
        cend_3d(inp)   = cstart_3d(inp) + nci_3d(inp)-1
      end do
!     end initialising v,e,fstart end  
      call read_geoSingleBody_3dVTK
      ! if(myid.eq.0)then
      ! write(*,*)'min x ',minval(xyz_3d(1,:)),maxval(xyz_3d(1,:)) 
      ! write(*,*)'min y ',minval(xyz_3d(2,:)),maxval(xyz_3d(2,:)) 
      ! write(*,*)'min z ',minval(xyz_3d(3,:)),maxval(xyz_3d(3,:)) 
      ! end if
!========================================================
!  initialising triangulated geometries (2D)
!     initial allocations 
      allocate(nvi(Nparticle),nei(Nparticle),nfi(Nparticle))
      allocate(vstart(Nparticle),vend(Nparticle))
      allocate(estart(Nparticle),eend(Nparticle))
      allocate(fstart(Nparticle),fend(Nparticle))


      nvtot = 0 ; netot = 0 ; nftot = 0    
      do inp=1,Nparticle
      call read_geo_dim(nvi(inp),nei(inp),nfi(inp),geofile(inp))

      if(myid.eq.0)then
        write(*,*)'Dim. Geo 2D: ',geofile(inp)
        write(*,*)'Nodes,Elements,Facets: ',nvi(inp),nei(inp),nfi(inp)
      end if

      nvtot = nvtot + nvi(inp)
      netot = netot + nei(inp)
      nftot = nftot + nfi(inp)

      end do

!     initialising vstart and vend for inp=1 and Nparticle
      vstart(1) = 1  ; estart(1) = 1 ; fstart(1) = 1
      vend(1) = nvi(1); eend(1)=nei(1) ; fend(1) = nfi(1)

      do inp=2,Nparticle
        vstart(inp) = vend(inp-1)+1
        vend(inp)   = vstart(inp) + nvi(inp)-1

        estart(inp) = eend(inp-1)+1
        eend(inp)   = estart(inp) + nei(inp)-1

        fstart(inp) = fend(inp-1)+1
        fend(inp)   = fstart(inp) + nfi(inp)-1
      end do
!     end initialising v,e,fstart end  

      !load 2D bodies
      call read_geoSingleBody
      ! if(myid.eq.0)then
      ! write(*,*)'min x ',minval(xyz(1,:)),maxval(xyz(1,:)) 
      ! write(*,*)'min y ',minval(xyz(2,:)),maxval(xyz(2,:)) 
      ! write(*,*)'min z ',minval(xyz(3,:)),maxval(xyz(3,:)) 
      ! end if

      call find_indices_2dstr
      write(*,*) "INDICES 2wet", inpstart_2dwet, inpend_2dwet
      write(*,*) nvstart_2dwet,nvend_2dwet,nestart_2dwet,neend_2dwet,nfstart_2dwet,nfend_2dwet
      write(*,*) "INDICES 2str", inpstart_2dstr, inpend_2dstr 
      write(*,*) nvstart_2dstr,nvend_2dstr,nestart_2dstr,neend_2dstr,nfstart_2dstr,nfend_2dstr
      write(*,*) nvtot,netot,nftot
      allocate(tag_2dwet(nvend_2dwet)) !again, we assume nvstart_2dwet=1 
      call tagging_2dwet
!=======================================================
      ! find the heart chambers
      !aggiunta recente: prepara volume shift per camere cardiache                                                                              
      call prepare_volume_chamb_shift
      call body1_to_chamb
!=======================================================
      ! hyperelastic material properties
#ifdef HYPERELASTIC
      call hyperelastic_properties
#endif
!========================================================
! boundary conditions for structures
!2D
      ! call find_boundaries(bcs,nvtot,netot,&
      !            face_of_edge, vert_of_edge, vert_to_part,&
      !            edge_to_part,xyz0,count2,boundary1,LabelBound)
      call find_boundaries(bcs,nvtot,netot,&
                 face_of_edge, vert_of_edge, vert_to_part,&
                 edge_to_part,xyz0,count2,boundary1)
      allocate(boundary2(2,count2))
      allocate(BCoffset(3,count2))  
      do ii=1,count2
         boundary2(1,ii)=boundary1(1,ii)
         boundary2(2,ii)=boundary1(2,ii)
      enddo
      deallocate(boundary1)

      call find_offset(nvtot,xyz0,vert_to_part,count2,boundary2,BCoffset)
!3D
      ! boundary1_3d(:,:) = 0 !initialise
      ! boundary2_3d(:,:) = 0 !initialise
      call find_boundaries_3d(nvtot_3d,vert_to_part_3d,vert_to_chamb_3d,&
                 xyz0_3d,count2_3d,boundary1_3d)
      allocate(boundary2_3d(2,count2_3d))
      allocate(BCoffset_3d(3,count2_3d))  
      do ii=1,count2_3d
         boundary2_3d(1,ii)=boundary1_3d(1,ii)
         boundary2_3d(2,ii)=boundary1_3d(2,ii)
      enddo
      deallocate(boundary1_3d)
      call find_offset_3d(nvtot_3d,xyz0,xyz0_3d,vert_to_part_3d,count2_3d,boundary2_3d,BCoffset_3d)

!========================================================
! concentrated mass 
      call find_mass
! #ifndef SOLOUNO
!       call find_movingprobes(nvtot,xyz0,vert_to_part)
! #endif
!========================================================
! pulmonary veins location
      ! call find_pulveins
!     interpolation routines from FV
!     ------------------------------
!          write(*,*) 'Read points outlet '
!          stri='InterpPoints.in'
!          call read_PInterp_dim(npInterp,stri)
!          call read_PInterp(xyzInterp(1:npInterp,1:3),npInterp,stri)
!          stri='QuadWeights.in'
!          call read_WQuad_dim(npQuad,npQuadIn,npQuadOut,stri)
!          call read_WQuad(WeiQuad(1:npInterp,1:3),npQuad,stri)
!          call search_nodicheckEF()
!     ------------------------------
   endif !ifmls
   call allocate_mls_local
   dt_o=dt
!========================================================
!     grid information  
      !
#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0)write(*,*)myid,'starting meshes'
#endif
      call meshes
#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0)write(*,*)myid,'starting cordin'
#endif
      call cordin
      
#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0)write(*,*)myid, 'starting coetar'
#endif
      call coetar
      !call coetarr
!m=======================================================
#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0)write(*,*)myid,'starting indic'
#endif
      call indic
      !call indicr 
!========================================================
      if(myid.eq.0) then
        write(*,754)n1,n2,n3
        write(*,756)mref1,mref2,mref3,lmax
        write(*,755)1.d0/dx1,1.d0/dx2,1.d0/dx3,dt,ntst
      endif
754   format(3x,'grid resolution: ','   n1=',i5, &
               '   n2=',i5,'   n3=',i5)
755   format(3x,'dx1=',e10.3,'   dx2=',e10.3,'   dx3=',e10.3, &
               '   dt=',e10.3,'   ntst=',i8,/)
756   format(3x,'refine scales:  ','    mref1=',i2,'    mref2=',i2, &
               '    mref3=',i2,'    lmax=',i2)

!     Locate the probes
      if (nson.ne.0) then
         call indson
         if (myid .eq. 0) then
          do i = 1,nson
            write(*,*)'Probe: ',i,(mason(i,l),l=1,3)
          end do
          do i = 1,nson
            write(*,'(3(2x,f12.8))') ,(coson(i,l),l=1,3)
          end do
          write(*,*)
        end if
      endif

! Locate the  aorta and ventricle
      cosonAO(1) = -0.55 
      cosonAO(2) = -0.9
      cosonAO(3) =  4.2
      cosonLV(1) =  0.8   
      cosonLV(2) = -0.3
      cosonLV(3) =  2.5
      cosonLA(1) =  0.2   
      cosonLA(2) =  0.23
      cosonLA(3) =  4.20

      cosonAP(1) =  0.3 
      cosonAP(2) = -1.0
      cosonAP(3) =  4.0
      cosonRV(1) = -0.7   
      cosonRV(2) = -0.6
      cosonRV(3) =  2.0
      cosonRA(1) = -1.8   
      cosonRA(2) = -0.4
      cosonRA(3) =  3.30

      call indsonAorta 


      if (myid .eq. 0 ) then
!       write(*,*) masonAO, masonLV
       write(*,'(3(2x,f12.8))') cosonAO
       write(*,'(3(2x,f12.8))') xc(masonAO(1)), yc(masonAO(2)), zc(masonAO(3))
       write(*,'(3(2x,f12.8))') cosonLV
       write(*,'(3(2x,f12.8))') xc(masonLV(1)), yc(masonLV(2)), zc(masonLV(3))
       write(*,'(3(2x,f12.8))') cosonLA
       write(*,'(3(2x,f12.8))') xc(masonLA(1)), yc(masonLA(2)), zc(masonLA(3))
      end if

 
#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0)write(*,*)myid,'starting write grid info'
#endif
      call write_grid_info

!=========================================================
#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0)write(*,*) myid,'starting phini'
#endif
      call phini

!     --------------------------------------------
!     Initialising Fadlun forcing
      !Initialize forclo array
      forclo = 1._DP

      if (infig.eq.1) then 
         nkix = 1
         nkfx = n3
         njix = 1
         njfx = n2
         niix = 1
         nifx = n1

!! raytracing
!     Solid body
!         stlfx='piastra23.stl'
!         if (myid.eq.0) write(*,*) 'Read geometry from ',stlfx
      !    call readgeo3ddim(nbfx,stlfx)
      !    call readgeo3d(mpugeo,xyzbfx,barfx,areafx,nbfx,stlfx)
      !    call readnormals(mpugeo,nxyzfx,nbfx,stlfx)

      !   call ibgrid3d(0,n1,n2,n3,xc,yc,zc &
      !                 ,mpugeo,nbfx,xyzbfx,nxyzfx,barfx,areafx &
      !                 ,mpun,mpun,npunifx,npunfx,indgeoee      &
      !                 ,indgeo,indgeoe,distb,0.0,0.0,0.0       &
      !                 ,vbfx,vbifx,ntribfx                     &
      !                 ,niix,nifx,njix,njfx,nkix,nkfx,PPro)

      ! if(myid.eq.0)then
      ! write(*,*) ' NPUNX = ',npunfx(1),npunifx(1)
      ! write(*,*) ' NPUNY = ',npunfx(2),npunifx(2)
      ! write(*,*) ' NPUNZ = ',npunfx(3),npunifx(3)
      ! write(*,*) ' NPUNP = ',npunfx(4),npunifx(4)
      ! end if

!-------------------------------------------------------
!LOAD CAP UNDER AORTA for Windkessel
!    Solid body AS
          stlfxAS = 'TappoPoint.stl'
         write(*,*) 'Read geometry from ',stlfxAS

         call readgeo3ddim(nbfxAS,stlfxAS)
         call readgeo3d(mpugeo,xyzbfxAS,barfxAS,areafxAS,nbfxAS,stlfxAS)
         call readnormals(mpugeo,nxyzfxAS,nbfxAS,stlfxAS)

        call ibgrid3d(0,n1,n2,n3,xc,yc,zc &
                      ,mpugeo,nbfxAS,xyzbfxAS,nxyzfxAS &
                      ,barfxAS,areafxAS                &
                      ,mpun,mpun,npunifxAS,npunfxAS    &
                      ,indgeoeeAS                      &
                      ,indgeoAS,indgeoeAS,distbAS      &
                      ,.0,0.0,0.0                      &
                      ,vbfxAS,vbifxAS,ntribfxAS        & 
                      ,niix,nifx,njix,njfx             &
                      ,nkix,nkfx,PPro)           

         write(*,*) ' NPUNX = ',npunfxAS(1),npunifxAS(1)
         write(*,*) ' NPUNY = ',npunfxAS(2),npunifxAS(2)
         write(*,*) ' NPUNZ = ',npunfxAS(3),npunifxAS(3)
         write(*,*) ' NPUNP = ',npunfxAS(4),npunifxAS(4)
         write(*,*)
 !-------------------------------------------------------


!     Logging direct-forcing IBM ray data exit infig
         if(myid.eq.0) then
         open(555,file='ray.data',form='unformatted')
         do i=1,4
            write(555)npunifx(i),npunfx(i)
            do n=1,npunfx(i)
               do m=1,3
                  write(555)indgeo(i,n,m)
                  write(555)indgeoe(i,n,m)
                  write(555)distb(i,n)
               enddo
            enddo
            do n=1,npunifx(i)
               do m=1,3
                  write(555)indgeoee(i,n,m)
               enddo
            enddo
         enddo
         close(555)
         end if


       if(numtasks.gt.1) call mpi_clean_ibm

!     debug file
      if(myid.eq.0)then

            open(412,file='indgeo.dat',status='unknown')
      do n=1,npunfx(1)
            write(412,*)indgeo(1,n,1:3)
      end do
            close(412)

      end if

      end if
!     --------------------------------------------

      do l=1,ndv                                                   
        vmax(l)=0.d0
      enddo                                   

#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0)write(*,*)myid,'starting densbo'
#endif
      call densbo
      call dsalbo
if (imlsfor.EQ.1) then            
!==============================================================
!     Initialize electrophysiology solver
!     XEF_3d will be overwritten by inirea      
#ifdef ELECTRO
   call ElectroStartEF_3d
   write(*,*) 'ElectroStart 3D done'
!  initialising triangulated geometries (1D)                                                                                                        
!     initial allocations                                                                                                                           
      allocate(nviEF_1d(NparticleEF_1d),neiEF_1d(NparticleEF_1d))
      allocate(vstartEF_1d(NparticleEF_1d),vendEF_1d(NparticleEF_1d))
      allocate(estartEF_1d(NparticleEF_1d),eendEF_1d(NparticleEF_1d))

      nvtotEF_1d = 0 ; netotEF_1d = 0
      do inp=1,NparticleEF_1d
      call read_geo_dim_1d(nviEF_1d(inp),neiEF_1d(inp),geofileEF_1d(inp))

      if(myid.eq.0)then
!        write(*,*)'Dim. Geo EF 1D: ',geofile(inp)                                                                                                  
        write(*,*)'Nodes,Elements EF1d: ',nviEF_1d(inp),neiEF_1d(inp)
      end if
      nvtotEF_1d = nvtotEF_1d + nviEF_1d(inp)
      netotEF_1d = netotEF_1d + neiEF_1d(inp)
      end do

!     initialising vstart and vend for inp=1 and Nparticle                                                                                          
      vstartEF_1d(1) = 1  ; estartEF_1d(1) = 1
      vendEF_1d(1) = nviEF_1d(1); eendEF_1d(1)=neiEF_1d(1)
      if (NparticleEF_1d.GT.1) then
      do inp=2,NparticleEF_1d
        vstartEF_1d(inp) = vendEF_1d(inp-1)+1
        vendEF_1d(inp)   = vstartEF_1d(inp) + nviEF_1d(inp)-1

        estartEF_1d(inp) = eendEF_1d(inp-1)+1
        eendEF_1d(inp)   = estartEF_1d(inp) + neiEF_1d(inp)-1
      end do
      endif
!     end initialising v,e,fstart end  

      write(*,*) vstartEF_1d(1),vendEF_1d(1)
      write(*,*) vstartEF_1d(2),vendEF_1d(2)
      write(*,*) vstartEF_1d(3),vendEF_1d(3)

      write(*,*) estartEF_1d(1),eendEF_1d(1)
      write(*,*) estartEF_1d(2),eendEF_1d(2)
      write(*,*) estartEF_1d(3),eendEF_1d(3)
      ! stop

      call read_geoSingleBodyEF_1d

      call ElectroStartEF_1d
      write(*,*) 'ElectroStart 1D done'

!1D
      allocate(boundary2EF_1d(nvtotEF_1d))
      allocate(BCoffsetEF_1d(3,nvtotEF_1d))  
      call find_tag_offsetEF_1d()


! !========================================================                                                                                      
!  initialising triangulated geometries (2D)                                                                                                        
!     initial allocations                                                                                                                           
      allocate(nviEF_2d(NparticleEF_2d),neiEF_2d(NparticleEF_2d),nfiEF_2d(NparticleEF_2d))
      allocate(vstartEF_2d(NparticleEF_2d),vendEF_2d(NparticleEF_2d))
      allocate(estartEF_2d(NparticleEF_2d),eendEF_2d(NparticleEF_2d))
      allocate(fstartEF_2d(NparticleEF_2d),fendEF_2d(NparticleEF_2d))

      nvtotEF_2d = 0 ; netotEF_2d = 0 ; nftotEF_2d = 0
      do inp=1,NparticleEF_2d
      call read_geo_dim(nviEF_2d(inp),neiEF_2d(inp),nfiEF_2d(inp),geofileEF_2d(inp))

      if(myid.eq.0)then
!        write(*,*)'Dim. Geo EF 2D: ',geofile(inp)                                                                                                  
        write(*,*)'Nodes,Elements,Facets EF2d: ',nviEF_2d(inp),neiEF_2d(inp),nfiEF_2d(inp)
      end if
      nvtotEF_2d = nvtotEF_2d + nviEF_2d(inp)
      netotEF_2d = netotEF_2d + neiEF_2d(inp)
      nftotEF_2d = nftotEF_2d + nfiEF_2d(inp)

      end do

!     initialising vstart and vend for inp=1 and Nparticle                                                                                          
      vstartEF_2d(1) = 1  ; estartEF_2d(1) = 1 ; fstartEF_2d(1) = 1
      vendEF_2d(1) = nviEF_2d(1); eendEF_2d(1)=neiEF_2d(1) ; fendEF_2d(1) = nfiEF_2d(1)
      if (NparticleEF_2d.GT.1) then
      do inp=2,NparticleEF_2d
        vstartEF_2d(inp) = vendEF_2d(inp-1)+1
        vendEF_2d(inp)   = vstartEF_2d(inp) + nviEF_2d(inp)-1

        estartEF_2d(inp) = eendEF_2d(inp-1)+1
        eendEF_2d(inp)   = estartEF_2d(inp) + neiEF_2d(inp)-1

        fstartEF_2d(inp) = fendEF_2d(inp-1)+1
        fendEF_2d(inp)   = fstartEF_2d(inp) + nfiEF_2d(inp)-1
      end do
      endif
!     end initialising v,e,fstart end                                                          

      call read_geoSingleBodyEF_2d
      call ElectroStartEF_2d
      write(*,*) 'ElectroStart 2D done'
!2D
      allocate(boundary2EF_2d(nvtotEF_2d))
      allocate(BCoffsetEF_2d(3,nvtotEF_2d))  
      call find_tag_offsetEF_2d()


!** Initialise EFtstart
      do i=1,nvtot_3d 
         EFtstart_3d(i) = -100000.d0
      enddo

      do i=1,nvtotEF_2d 
         EFtstart_2d(i) = -100000.d0
      enddo
      do i=1,nvtotEF_1d 
         EFtstart_1d(i) = -100000.d0
      enddo


#else
      tstartEF_LV = 0.d0/TSTAR
      tstartEF_LA = (PeriodDim-200.d0/1000.d0)/TSTAR
      do i=1,nvtot_3d
         potEFnode_3d(i) = 0.d0
         EFtstart_3d(i) = -10.d0
        
         chamb=vert_to_chamb_3d(i)
         if (chamb.EQ.1) then
            EFtstart_3d(i) = tstartEF_LV
         elseif (chamb.EQ.2) then
            EFtstart_3d(i) = tstartEF_LA
         elseif (chamb.EQ.3) then !-> paste LV
            EFtstart_3d(i) = tstartEF_LV
         elseif (chamb.EQ.4) then !-> paste LA
            EFtstart_3d(i) = tstartEF_LA
         endif
      enddo

#endif


!==============================================================
#ifdef CORONARYVEINS      
!     Initialize Coronary veins
      allocate(nviCV_1d(NparticleCV_1d),neiCV_1d(NparticleCV_1d))
      allocate(vstartCV_1d(NparticleCV_1d),vendCV_1d(NparticleCV_1d))
      allocate(estartCV_1d(NparticleCV_1d),eendCV_1d(NparticleCV_1d))

      nvtotCV_1d = 0 ; netotCV_1d = 0
      do inp=1,NparticleCV_1d
      call read_geo_dim_1d(nviCV_1d(inp),neiCV_1d(inp),geofileCV_1d(inp))

      if(myid.eq.0)then
!        write(*,*)'Dim. Geo EF 1D: ',geofile(inp)                                                                                                  
        write(*,*)'Nodes,Elements CV1d: ',nviCV_1d(inp),neiCV_1d(inp)
      end if
      nvtotCV_1d = nvtotCV_1d + nviCV_1d(inp)
      netotCV_1d = netotCV_1d + neiCV_1d(inp)
      end do

!    initialising vstart and vend for inp=1 and Nparticle                                                                                          
      vstartCV_1d(1) = 1  ; estartCV_1d(1) = 1
      vendCV_1d(1) = nviCV_1d(1); eendCV_1d(1)=neiCV_1d(1)
      if (NparticleCV_1d.GT.1) then
      do inp=2,NparticleCV_1d
        vstartCV_1d(inp) = vendCV_1d(inp-1)+1
        vendCV_1d(inp)   = vstartCV_1d(inp) + nviCV_1d(inp)-1

        estartCV_1d(inp) = eendCV_1d(inp-1)+1
        eendCV_1d(inp)   = estartCV_1d(inp) + neiCV_1d(inp)-1
      end do
      endif
!     end initialising v,e,fstart end  

      write(*,*) vstartCV_1d(1),vendCV_1d(1)
      write(*,*) estartCV_1d(1),eendCV_1d(1)

      call read_geoSingleBodyCV_1d

!1D
      allocate(boundary2CV_1d(nvtotCV_1d))
      allocate(BCoffsetCV_1d(3,nvtotCV_1d))  
      call find_tag_offsetCV_1d()
#endif

 endif !imlsfor

! !========================================================                                                                                      
if (ECG_Ppro.LE.0) then
!==============================================================
!     create the initial conditions
!
      if(nread.eq.0) then

        if(myid.eq.0) then
          write(*,'(6x,a,/)')'nread=0 ---> new initialization'
        endif

        ntime=0                                                         
        time=0.d0

        timeBeat=modulo(time,period)
        nBeat = floor(time/period)
        cflm=0.d0
        cflmr=0.d0
       
#ifdef DEBUG                               
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(myid.eq.0)write(*,*)myid,'starting inqpr'
#endif
        call inqpr

199    format(i8.8)
125    format(e24.16)

!============================================================= 
!
     else !nread
!
!=============================================================
        if(myid.eq.0) then
          write(*,'(6x,a,/)')'nread=1 ---> reading files'
        endif

#ifdef DEBUG
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
        if(myid.eq.0)write(*,*)myid, 'starting inirea' 
#endif

        call inirea
        call inistruc
!        tmax = tmax + time
      endif !nread
      nBeatSwitch=nBeat
      
      call update_both_ghosts(n1,n2,q1,kstart,kend)
      call update_both_ghosts(n1,n2,q2,kstart,kend)
      call update_both_ghosts(n1,n2,q3,kstart,kend)
      call update_both_ghosts(n1,n2,dens,kstart,kend)
      call update_both_ghosts(n1r,n2r,dsal,kstartr,kendr)
      call update_both_ghosts(n1,n2,pr,kstart,kend)


!============================================================
!    the boundary condition at upper and lower wall
!    for free-slip boundary 
#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0)write(*,*)myid, 'starting vel bc'
#endif
      call velbc
!============================================================
!   allocate memory for multigrid
#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0)write(*,*)myid, 'starting mem alloc multigrid'
#endif
      call mgrd_mem_alloc
      if(myid.eq.0)write(*,'(/,3x,a,/)')'multigrid mem allocated'
!=============================================================
!   for multi grid interpolation
#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0)write(*,*)myid, 'starting idc for mgrd'
#endif
      call mgrd_idc

!=============================================================
!    inital divergence check
#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0)write(*,*)myid, 'starting vel interpo'
#endif
      call mgrd_velitp
      call mgrd_dsalc
!============================================================
!   statistics initialization
#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0)write(*,*)myid,'starting initstst'
#endif
!      if (PPro.LE.0) call initstst
!=============================================================
!  for storing slices and mean profiles
#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0)write(*,*)myid,'starting slab_ini'
#endif
      call slab_ini
      call profiles_ini
!=============================================================
#ifdef MOVIE
      call inimov
#endif
if (imlsfor.EQ.1) then 
      !if (n3/2+1 .le. kend .and. n3/2+1 .ge. kstart) then
      ipfi = 'Start'
      call write_to_vtk_field(ipfi)
      if(myid.eq.0)then
        call write_to_vtk(ipfi,0)
        call write_to_vtk_3d(ipfi,0)
#ifdef ELECTRO
        !call write_to_vtkEF_1d(ipfi,0)
        call write_to_vtkEF_2d(ipfi,0)
#endif
        !call write_to_vtkCV_1d(ipfi,0)
        !call write_to_vtkCenSeg_1d(ipfi,0)
     endif
endif     
!=============================================================
#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
      if(myid.eq.0)write(*,*)myid, 'starting divgck'
#endif

      call divgck(dmaxc,dmaxr)
      call toteng(aengc,aengr)

!      call write_base
!      call write_fine

      if(myid.eq.0) then
       write(*,'(3x,a,es10.3,a,es10.3)') &
           'Max div base mesh:',dmaxc,  &
           '    Max div refined mesh:',dmaxr
       write(*,'(3x,a,es10.3,es10.3,es10.3,es10.3,a,es10.3,/)') &
           'tot KTE base mesha:',aengc,aengcQ1,aengcQ2,aengcQ3, &
           '    tot KTE refined mesh:',aengr
      endif

!==========================================================
      ntstf=ntst     
!===========================================================
      if(myid.eq.0) then
        write(*,711) trest,ntstf,tpin
      endif
711   format(3x,'check in cond : trest =',f10.1,  &
            3x,'ntstf =',i8,5x,'tpin =',es10.2,/)
#ifdef FLUID_STRUCTURE_SOLVER 
      call SetTiling
#endif
      if(myid.eq.0)print*,'Starting time iteration'
      if(myid.eq.0)print*,'-------------------------------'
!==============================================================
!  ********* starts the time dependent calculation ***

      tin(2) = MPI_WTIME()
      if(myid.eq.0) then
        write(*,'(3x,a,f10.3,a,/)') 'Initialization Time = ',  &
          tin(2)-tin(1), ' sec.'
      endif


! QUI METTI IF ELSE E LEGGI DA FILE

!Initialize SMZ
      SMZ = 1.d0
      !if (n3/2+1 .le. kend .and. n3/2+1 .ge. kstart) then

#ifdef USE_CUDA
      istat=CudaMemGetInfo(freeMem,totalMem)
      rTotalMem = totalMem/(1024.**2)
      rFreeMem =   freeMem/(1024.**2)
      rUsedMem = (totalMem-freeMem)/(1024**2)
      call MPI_ALLREDUCE(rUsedMem,maxUsedMem,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(rUsedMem,minUsedMem,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(rFreeMem,maxFreeMem,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(rFreeMem,minFreeMem,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      if(myid.eq.0)then
        write(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " GPU memory used: ",minUsedMem," - ",maxUsedMem," / ",rTotalMem," MBytes"
        write(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " GPU memory free: ",maxFreeMem," - ",minFreeMem," / ",rTotalMem," MBytes"
        print *," "
      endif
#endif


!m===============================================================
!m===============================================================
!m    MAIN LOOOOOP           
!m===============================================================
!m===============================================================
      DO ntime=1,ntstf             
        call nvtxStartRange("Timestep", 1)
!
!m================================================================
!     the calculation stops if 
!     the velocities are diverging for numerical
!     stability conditions (courant number restrictions)                
!
#ifdef DEBUG                                                               
        if(myid.eq.0)then
          write(*,*)' '
          write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$'
          write(*,*)'   NEW TIME STEP: ', ntime
          write(*,*)'starting cfl'
        endif 
#endif
        call nvtxStartRange("cfl", 2)
        call cfl(cflm)        
        call nvtxEndRange !cfl
!        call nvtxStartRange("cflr", 3)
!        call cflr(cflmr)
!        call nvtxEndRange

#ifdef DEBUG                                                               
        if(myid.eq.0)write(*,*)myid, 'finished cfl', cflm
#endif 

        if(idtv.eq.1.and.cflm.gt.1.0e-05) then
          if(ntime.gt.1) then
            dt=cflmax/cflm
            if(dt.gt.dtmax) dt=dtmax
          endif
            if(dt.lt.1.d-8) go to 166
        else
          cflm=cflm*dt
          if(cflm.gt.cfllim) go to 165
        endif

        beta=dt*nu*0.5d0
!=============================================================
!     TIME ADVANCE
!=============================================================
#ifdef DEBUG                                                               
        if(myid.eq.0)write(*,*)myid, 'starting tschem'
#endif

        ti(1) = MPI_WTIME()

        call nvtxStartRange("tschem", 4)

#ifdef LOOSE
        call tschem
#endif 
#ifdef STRONG
        call tschemStrong
#endif
#ifdef STRONGRV
        call tschemStrongRV
#endif
        call nvtxEndRange
115 format(1(2x,f12.4))
94 format(i2.2)
        time=time+dt           ! simulation time
        ti(2) = MPI_WTIME()    ! CPU time        
        timeBeat= modulo(time,period) !time within the heart beat
        nBeat = floor(time/period)

        if (nBeat.GT.nBeatSwitch) then
           do i=1,nPurkVentr
              Purk_tstim(i)=-10.d0
           enddo
           do i=1,nBundAtri
              Bund_tstim(i)=-10.d0
           enddo
           nBeatSwitch=nBeat
           write(*,*) "NUOVO CICLO", nBeatSwitch

           !f_apd_2d(:)=0.0D0
           !t_apd_2d(:,:)=0.0D0
           
           !salva tempi
           IF(myid.eq.0)then
              write(xCoord,94)nBeat-1
              filname = trim('t_act_3d_'//trim(xCoord)//'.uqout')
              open(117, file=filname, &
                   status='unknown',           &
                   access='sequential',        &
                   position='append')
              do i=1,nvtot_3d
                  write(117,115) EFtstart_3d(i)*(1000.D0*TSTAR)
              enddo
              close(117)

              filname = trim('t_act_2d_'//trim(xCoord)//'.uqout')
              open(117, file=filname, &
                   status='unknown',           &
                   access='sequential',        &
                   position='append')
              do i=1,nvtotEF_2d
                  write(117,115) EFtstart_2d(i)*(1000.D0*TSTAR)
              enddo

              close(117)
              filname = trim('t_act_1d_'//trim(xCoord)//'.uqout')
              open(117, file=filname, &
                   status='unknown',           &
                   access='sequential',        &
                   position='append')
              do i=1,nvtotEF_1d
                  write(117,115) EFtstart_1d(i)*(1000.D0*TSTAR)
              enddo
              close(117)
              
              filname = trim('t_act_AV_'//trim(xCoord)//'.uqout')
              open(117, file=filname, &
                   status='unknown',           &
                   access='sequential',        &
                   position='append')
                  write(117,115) EFtstart_1d(vAVslave(1))*(1000.D0*TSTAR)
                  write(117,115) EFtstart_1d(vAVslave(2))*(1000.D0*TSTAR)
                  write(117,115) EFtstart_1d(vAVmaster(3))*(1000.D0*TSTAR)
              close(117)

           endif

        endif

!=============================================================
!   computing statistics
!
          if (mod(time,tpin).LT.dt) then
          timeint_cdsp = timeint_cdsp + 1

!====  max velocity and nusselt numbers calculated by definition 
          call nvtxStartRange("vmaxv", 5)
          call vmaxv          
          call nvtxEndRange
          if(vmax(1).gt.1000.d0.and.vmax(2).gt.1000.d0) go to 266

!====  max divergence on base and refined meshes
          call nvtxStartRange("divgck", 9)
          call divgck(dmaxc,dmaxr)
          call nvtxEndRange
          if(dmaxc.gt.resid) go to 169

!====  total enegery on base and refined meshes
          call nvtxStartRange("toteng", 10)
          call toteng(aengc,aengr)
          call nvtxEndRange

          call nvtxStartRange("totengQ", 111)
          call totengQ(aengcQ1,aengcQ2,aengcQ3,aengr)
          call nvtxStartRange("totengQ", 111)

!====  display statistics
          if(myid.eq.0) then
            write(*,*) '---------------------------------------- '
            write(*,'(a,f12.4,a,i8,a,es10.2)') &
                 '  T = ',time,'   NTIME = ',ntime,'   DT = ',dt
            write(*,'(a,f12.4,a,es10.2)') &
                 '  Tms = ',time*(TSTAR*1000.d0),'  DTms = ',dt*(TSTAR*1000.d0)
            write(*,'(a,f8.4,a)') &
                 '  Iteration Time = ', ti(2)-ti(1), ' sec.'
            write(*,'(a,es10.3,a,es10.3)') &
                '  Max div base mesh:',dmaxc, &
                '  Max div refined mesh:',dmaxr
            write(*,'(a,es10.3,es10.3,es10.3,es10.3,a,es10.3)') &
                '  tot KTE base mesh:',aengc,aengcQ1,aengcQ2,aengcQ3, &
                '  tot KTE refined mesh:',aengr
#if defined(STRONG) || defined(STRONGRV)
            write(*,'(a,i8,es10.3)') &
                '  FSI iter and err: ',iter,max_error
#endif
#ifdef BIDOMAIN
            write(*,'(a,i8,es10.3,es10.3)') &
                '  GMRES: ',countjgmres,scarto1gmres,scarto2gmres
#endif
#ifdef SCALINGTEST
            write(*,'(a,es10.3)') 'expl prov vel : ',timesforscaling(1)
            write(*,'(a,es10.3)') 'update halo   : ',timesforscaling(2)
            write(*,'(a,es10.3)') 'force TiledOff: ',timesforscaling(3)
            write(*,'(a,es10.3)') 'force Tileed  : ',timesforscaling(4)
            write(*,'(a,es10.3)') 'apply force   : ',timesforscaling(5)
            write(*,'(a,es10.3)') 'compute divg  : ',timesforscaling(6)
            write(*,'(a,es10.3)') 'solve elliptic: ',timesforscaling(7)
            write(*,'(a,es10.3)') 'solenoidal vel: ',timesforscaling(8)
            write(*,'(a,es10.3)') 'compute press : ',timesforscaling(9)
            write(*,'(a,es10.3)') 'external force: ',timesforscaling(10)
            write(*,'(a,es10.3)') 'internal force: ',timesforscaling(11)
! #ifdef USE_CUDA
!       istat=CudaMemGetInfo(freeMem,totalMem)
!       rTotalMem = totalMem/(1024.**2)
!       rFreeMem =   freeMem/(1024.**2)
!       rUsedMem = (totalMem-freeMem)/(1024**2)
!       call MPI_ALLREDUCE(rUsedMem,maxUsedMem,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
!       call MPI_ALLREDUCE(rUsedMem,minUsedMem,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
!       call MPI_ALLREDUCE(rFreeMem,maxFreeMem,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
!       call MPI_ALLREDUCE(rFreeMem,minFreeMem,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
!       if(myid.eq.0)then !questo ridondante
!         write(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " GPU memory used: ",minUsedMem," - ",maxUsedMem," / ",rTotalMem," MBytes"
!         write(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " GPU memory free: ",maxFreeMem," - ",minFreeMem," / ",rTotalMem," MBytes"
!         print *," "
!       endif
! #endif

#endif
          endif !master thread
        endif !computing statistics

!   save slabs
!     RTO
        !if( dabs(time-dble(nint(time/tframe))*tframe) .lt. dt*0.5d0 )then
        !  call slab_rcd
        !  if(myid.eq.0) then
        !    write(*,*) '   ============================= '
        !    write(*,'(a,f12.4)') ' slab saved T = ',time
        !    write(*,*) '   ============================= '
        !  endif
        !endif
!     RTO
!===============================================================
! !   make movie
! #ifdef MOVIE
!         if( dabs(time-dble(nint(time/tframe))*tframe) .lt. 0.5d0*dt ) then
!           call mkmov_dsal
!         endif
! #endif

        !   save restart data and statistics
      ! if (nBeat.GE.1) trest=0.1D0 !ricambia
      ! if (nBeat.GE.1) tframe=0.01D0 !ricambia

!===============================================================
!   save visualization data for paraview
    if ( (modulo(time,tframe).LT.dt).AND.(Ppro.LE.0) ) then
!         write(ipfi,23) max(0,nint(time/tframe))
!         write(ipfi,23) ntime
        ! itime=nint(1000.*time)          !FV                                                  
        itime=nint(1000.*time*TSTAR) !FV cosi salva il ms                                                  
        itime=max(itime,0)
        write(ipfi,199) itime

  23    format(i5.5)
         !if (myid .eq. 0) call write_to_vtk(ipfi,1) !struc2d
         if (myid .eq. 0) call write_to_vtk_3d(ipfi,1) !struc3d
         !call write_to_vtk_field(ipfi) !2d slices
         ! if (myid.eq.0) then
         !    filname = 'apdcheck_endo_tmp.txt'
         !    open(121,file = filname)
         !    write(121,*) (-83.0D0+potEFnode_3d(1)*(85.7D0)),timeBeat*LSTAR*1000.0D0
         !    close(121)
         !    call system('cat apdcheck_endo_tmp.txt >> apdcheck_endo.txt')
         !    call system('rm apdcheck_endo_tmp.txt')
         ! endif
#ifdef ELECTRO
         !if (myid .eq. 0) call write_to_vtkEF_2d(ipfi,1) !2d purkinje
         !if (myid .eq. 0) call write_to_vtkEF_1d(ipfi,1) !1d bundle
#endif
#ifdef CORONARYVEINS
         !if (myid .eq. 0) call write_to_vtkCV_1d(ipfi,1) !coronarie
#endif         

     end if

!===============================================================
      
        if ( (modulo(time,trest).LT.dt).AND.(Ppro.LE.0) ) then
           tin(3) = MPI_WTIME()
           if(myid.eq.0)write(*,'(a,f9.2,a)') '  SAVE RESTART Time = ',(tin(3)-tin(2))/3600.0,'h'
           call mpi_write_continua
           call continua_str      
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)           
           call ststwr
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        endif

        if ( (modulo(time,tecg).LT.dt).AND.(Ppro.LE.0) ) then
           tin(3) = MPI_WTIME()
           if(myid.eq.0)write(*,'(a,f9.2,a)') '  SAVE ECG RESTART Time = ',(tin(3)-tin(2))/3600.0,'h'
           call continua_ecg
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        endif


!=============================================================
!   check if reaching time limitations
!=============================================================
!   exit when reach tmax
        if(time.ge.tmax) go to 333       !FV

!==========================================================
!   exit at onw days on marconi 
        if( ti(2)-tin(1) .gt. walltimemax ) go to 333

!============================================================
        call nvtxEndRange

      ENDDO   ! END OF TIME LOOP

333   continue
!==========================================================

      open(11,file='vtkfiles/data_cv_apd_12.txt')
      do i =cstart_3d(1),cend_3d(1)  
         write(11,*) t_apd_3d(1,i),t_apd_3d(2,i)
      enddo
      close(11)
      open(11,file='vtkfiles/data_cv_apd_34.txt')
      do i =cstart_3d(1),cend_3d(1)  
         write(11,*) t_apd_3d(3,i),t_apd_3d(4,i)
      enddo
      close(11)

      ! open(11,file='vtkfiles/data_cv_EF2d_apd_12.txt')
      ! do i =vstartEF_2d(1),vendEF_2d(1)  
      !    write(11,*) t_apd_2d(1,i),t_apd_2d(2,i)
      ! enddo
      ! close(11)
      ! open(11,file='vtkfiles/data_cv_EF2d_apd_34.txt')
      ! do i =vstartEF_2d(1),vendEF_2d(1)  
      !    write(11,*) t_apd_2d(3,i),t_apd_2d(4,i)
      ! enddo
      ! close(11)


      !call postprocessing routines
      if (PPro.GE.1) then
         call postpro 
         tin(3) = MPI_WTIME()
         if(myid.eq.0)write(*,'(a,f9.2,a)') '  SAVE PPRO RESTART Time = ',(tin(3)-tin(2))/3600.0,'h'
         call mpi_write_continua

         ! call continua_str      
         ! call MPI_BARRIER(MPI_COMM_WORLD,ierr)           
         ! call ststwr
         ! call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         !vtkfiles are saved anyway
      endif

!============================================================
!
      !data dump moved to 167 FV
      go to 167  

!==============================================================    
165   continue 
      if(myid.eq.0) then
        write(*,164) 
      endif 
164   format(10x,'cfl too large!  ') 
      go to 167 

!===============================================================       
166   continue
      if(myid.eq.0) then
        write(*,168) dt 
      endif
168   format(10x,'dt too small, DT= ',e14.7)
      go to 167 
                       
!================================================================          
266   continue
      if(myid.eq.0) then
        write(*,268)  
      endif 
268   format(10x,'velocities diverged!  ') 
      go to 167 
      
!================================================================          
169   continue
      if(myid.eq.0) then
        write(*,178) dmaxc
      endif  
178   format(10x,'too large local residue for mass conservation : ' ,e12.5,' at ')     
      call divgloc
 
!==============================================================                               
!
167   continue



else  !Post processing ECG
!==============================================================
!     create the initial conditions

   iread = 1
   nPProR = (nreadend-nreadstart+ECG_PPro)/ECG_PPro
   nPPro = nint( nPProR )
   do nwrit=nreadstart,nreadend,ECG_PPro            
      write(*,*) "ANALYZING", iread,"/",nPPro
      write(*,*) "NWRITE", nwrit

        ! if(myid.eq.0) then
        !   write(*,'(6x,a,/)')'nread=1 ---> reading files'
        ! endif

      call iniecg
!==========================================================
      ntstf=ntst     
!===========================================================
      if(myid.eq.0)print*,'Starting time iteration'
      if(myid.eq.0)print*,'-------------------------------'
!==============================================================
!  ********* starts the time dependent calculation ***


      ntime=1


!m===============================================================
!m===============================================================
!m    MAIN LOOOOOP           
!m===============================================================
!m===============================================================
!      DO ntime=1,ntstf             
!=============================================================
!     TIME ADVANCE
!=============================================================

        ! ti(1) = MPI_WTIME()
        call nvtxStartRange("tschem", 4)
        call preambolotsch
        call nvtxEndRange

        time=time+dt           ! simulation time
        timeBeat= modulo(time,period) !time within the heart beat
        nBeat = floor(time/period)
        iread = iread + 1
!      ENDDO   ! END OF TIME LOOP

     enddo


endif !end ECG postprocessing

      ! return  
      ! end 



!   save data before exit
      ! write(ipfi,23) max(0,nint(time/tframe))
      itime=nint(1000.*time)          !FV                                                  
      itime=nint(1000.*time*TSTAR) !cosi salva inMS
      itime=max(itime,0)
      write(ipfi,199) itime
      ! call write_to_vtk_field(ipfi)  !2d slices
      if(myid.eq.0)then !second input=2 in order to not save fields
        call write_to_vtk(ipfi,1) !struc2d
        call write_to_vtk_3d(ipfi,1) !struc3d
#ifdef ELECTRO
        call write_to_vtkEF_2d(ipfi,1) !2d
        call write_to_vtkEF_1d(ipfi,1) !1d bundle
#endif
        call write_to_vtkCV_1d(ipfi,1) !coronarie
      endif

      call phend


      return  
      end 

