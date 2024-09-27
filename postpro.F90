!************************************************************************
      subroutine postpro
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
      integer :: ntstf,l,itime
      real(DP)    :: cflm,cflmr
      real(DP)    :: dmaxc,dmaxr
      real(DP)    :: aengc,aengcQ1,aengcQ2,aengcQ3,aengr
      real(DP)    :: ti(2), tin(3)
      character*100 ipfi
      character*150 :: nomefile
      character *150,dimension(:),allocatable :: geofilePpro

      integer :: i,m,n,inp,ii,j,k
      integer :: ic,jc,kc
      integer :: v1,v2,v3
      integer :: nvm,nem,nfm
      integer :: NparticlePpro,kstartp
      character*50 :: outfileXML,namfi
      character*10 ::NumVer
      character*30 :: xCoord, yCoord, zCoord
      real(DP) :: facemass,totmass
      real(DP) :: A1, A2,xV,yV,zV,xV1,yV1,zV1,rV1
      real(DP), dimension(n1m,n2m,n3m) :: xg,yg,zg

      maskq1=0.0d0
      maskq2=0.0d0
      maskq3=0.0d0
      maskpr=0.0d0

      NparticlePpro = 6
      allocate(geofilePpro(NparticlePpro))
      geofilePpro(1)='./meshes/cuoreMask.gts'
      geofilePpro(2)='./meshes/tappoAOepi.gts'
      geofilePpro(3)='./meshes/tappoAOabd.gts'
      geofilePpro(4)='./meshes/tappoVP.gts'
      geofilePpro(5)='./meshes/tappoAP.gts'
      geofilePpro(6)='./meshes/tappoVC.gts'

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)           
      if(myid.eq.0) then
      nomefile='./meshes/cuoreMask.gts'
      call body2gts(1,nomefile)
      nomefile='./meshes/tappoAOepi.gts'
      call body2gts(6,nomefile)
      nomefile='./meshes/tappoAOabd.gts'
      call body2gts(7,nomefile)
      nomefile='./meshes/tappoVP.gts'
      call body2gts(8,nomefile)
      nomefile='./meshes/tappoAP.gts'
      call body2gts(9,nomefile)
      nomefile='./meshes/tappoVC.gts'
      call body2gts(10,nomefile)
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)           
      call raytracing_Mask(NparticlePpro,geofilePpro)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)           
      if(kstart.eq.1) then
         kstartp=2
      else
         kstartp=kstart
      end if

      do kc=kstartp,kend
         do jc=1,n2
            do ic=1,n1
 !                 q1(ic,jc,kc) = q1(ic,jc,kc)*maskq1(ic,jc,kc)
 !                 q2(ic,jc,kc) = q2(ic,jc,kc)*maskq2(ic,jc,kc)
 !                 q3(ic,jc,kc) = q3(ic,jc,kc)*maskq3(ic,jc,kc)
 !                 pr(ic,jc,kc) = pr(ic,jc,kc)*maskpr(ic,jc,kc)
                 q1(ic,jc,kc) = q1(ic,jc,kc)*maskpr(ic,jc,kc)
                 q2(ic,jc,kc) = q2(ic,jc,kc)*maskpr(ic,jc,kc)
                 q3(ic,jc,kc) = q3(ic,jc,kc)*maskpr(ic,jc,kc)
                 pr(ic,jc,kc) = pr(ic,jc,kc)*maskpr(ic,jc,kc)
            enddo
         enddo
      enddo
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)           

     ipfi='Mask'
     call write_to_vtk_Mask(ipfi)

!alternativa per salvare il 3D (h5 mi dava problemi con paraview su tesla)
 !     if (iread.EQ.1) then
 !      write(*,*)'-------- Grid file -----------'
 !      namfi='./cfield/grid.xyz'
 !      do k=1,n3m
 !         do j=1,n2m
 !            do i=1,n1m
 !               xg(i,j,k)=xm(i)
 !               yg(i,j,k)=ym(j)
 !               zg(i,j,k)=zm(k)
 !            enddo
 !         enddo
 !      enddo
 !     endif

 !      open(99,file=namfi,form='unformatted')
 !      write(99) 1
 !      write(99) n1m,n2m,n3m
 !      write(99) ((((xg(i,j,k)),   &
 !        i=1,n1m),j=1,n2m),k=1,n3m),  &
 !        ((((yg(i,j,k)), &
 !        i=1,n1m),j=1,n2m),k=1,n3m), &
 !        ((((zg(i,j,k)), &
 !        i=1,n1m),j=1,n2m),k=1,n3m) 
 !      close(99)

 !      write(*,*)'-------- ENd file -----------'

 !      itime=nint(1000.*time)          !FV                                                    
 !      itime=max(itime,0)
 !      write(ipfi,199) itime
 ! 199   format(i8.8)

 !       namfi='./cfield/field3d_'//trim(ipfi)//'.q'
 !      open(99,file=namfi,form='unformatted')
 !      write(99) 1
 !      write(99) n1m,n2m,n3m
 !      write(99) 0.0,0.0,0.0,0.0
 !      write(99)   &
 !        ((((pr(i,j,k)),  &
 !        i=1,n1m),j=1,n2m),k=1,n3m),  &
 !        ((((0.5*(q1(i,j,k)+q1(i+1,j,k))), &
 !        i=1,n1m),j=1,n2m),k=1,n3m),  &
 !        ((((0.5*(q2(i,j,k)+q2(i,j+1,k))), &
 !        i=1,n1m),j=1,n2m),k=1,n3m), &
 !        ((((0.5*(q3(i,j,k)+q3(i,j,k+1))), &
 !        i=1,n1m),j=1,n2m),k=1,n3m), &
 !        ((((maskpr(i,j,k)), &
 !        i=1,n1m),j=1,n2m),k=1,n3m)
 !      close(99)

      return  
      end 

!#####################################################################      
      subroutine body2gts(qualeB,nomegts)
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
      character*150 :: nomegts
      integer :: qualeB,vstartB,vendB,estartB,eendB,fstartB,fendB
      character*30 :: xCoord, yCoord, zCoord

      vstartB = vstart(qualeB)
      vendB   = vend(qualeB) 
      estartB = estart(qualeB)
      eendB   = eend(qualeB) 
      fstartB = fstart(qualeB)
      fendB   = fend(qualeB) 
 99   format(i8.8)
 125  format(e24.16)

if(myid.eq.0) then
      write(*,*) 'Write new .gts file'
      open(12,file=nomegts)
      write(xCoord,99) nvi(qualeB)
      write(yCoord,99) nei(qualeB)
      write(zCoord,99) nfi(qualeB)
      write(12,*)' '//trim(xCoord)//' '//trim(yCoord)//' '//trim(zCoord)
      do ii=vstartB,vendB
         write(xCoord,125) xyz(1,ii)
         write(yCoord,125) xyz(2,ii)
         write(zCoord,125) xyz(3,ii)
         write(12,*)' '//trim(xCoord)//' '//trim(yCoord)//' '//trim(zCoord)
      enddo
      do ii=estartB,eendB
         write(xCoord,99) vert_of_edge(1,ii)-(vstartB-1)
         write(yCoord,99) vert_of_edge(2,ii)-(vstartB-1)
         write(12,*)' '//trim(xCoord)//' '//trim(yCoord)
      enddo
      do ii=fstartB,fendB
         write(xCoord,99) edge_of_face(1,ii)-(estartB-1)
         write(yCoord,99) edge_of_face(2,ii)-(estartB-1)
         write(zCoord,99) edge_of_face(3,ii)-(estartB-1)
         write(12,*)' '//trim(xCoord)//' '//trim(yCoord)//' '//trim(zCoord)
      enddo
      close(12)
endif !myid

      return
      end subroutine body2gts
!#####################################################################      

      subroutine raytracing_Mask(NparticlePpro,geofilePpro)
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
      character*150,dimension(NparticlePpro) :: geofilePpro
      integer :: NparticlePpro
      integer, dimension (NparticlePpro) :: nviPpro,neiPpro,nfiPpro
      integer, dimension (NparticlePpro) :: vstartPpro,vendPpro,estartPpro,eendPpro,fstartPpro,fendPpro

      integer,allocatable, dimension(:) :: n_edge_of_vertPpro
      integer,allocatable, dimension(:,:) :: vert_of_edgePpro, face_of_edgePpro,vert_of_facePpro, edge_of_facePpro,vert_of_vertPpro, edge_of_vertPpro
      real(DP),allocatable, dimension(:) :: surPpro
     real(DP),allocatable, dimension(:,:) :: xyzPpro,xyzvPpro,xyzaPpro,tri_verPpro,tri_velPpro,tri_barPpro,tri_norPpro,vel_triPpro,acc_triPpro

      nvtotPpro = 0 ; netotPpro = 0 ; nftotPpro = 0    
      do inp=1,NparticlePpro
      call read_geo_dimPpro(nviPpro(inp),neiPpro(inp),nfiPpro(inp),geofilePpro(inp))

      nvtotPpro = nvtotPpro + nviPpro(inp)
      netotPpro = netotPpro + neiPpro(inp)
      nftotPpro = nftotPpro + nfiPpro(inp)
      enddo

      
      allocate(n_edge_of_vertPpro(nvtotPpro))
      allocate(vert_of_edgePpro(2,netotPpro), face_of_edgePpro(2,netotPpro))
      allocate(vert_of_facePpro(3,nftotPpro), edge_of_facePpro(3,nftotPpro))
      allocate(vert_of_vertPpro(max_n_edge_of_vert,nvtotPpro), edge_of_vertPpro(max_n_edge_of_vert,nvtotPpro))
      allocate(xyzPpro(3,nvtotPpro),xyzvPpro(3,nvtotPpro),xyzaPpro(3,nvtotPpro))
      allocate(surPpro(nftotPpro))
      allocate(tri_verPpro(9,nftotPpro),tri_velPpro(9,nftotPpro))
      allocate(tri_barPpro(3,nftotPpro),tri_norPpro(3,nftotPpro),vel_triPpro(3,nftotPpro),acc_triPpro(3,nftotPpro))



!     initialising vstart and vend for inp=1 and Nparticle
      vstartPpro(1) = 1  ; estartPpro(1) = 1 ; fstartPpro(1) = 1
      vendPpro(1) = nviPpro(1); eendPpro(1)=neiPpro(1) ; fendPpro(1) = nfiPpro(1)

      do inp=2,NparticlePpro
        vstartPpro(inp) = vendPpro(inp-1)+1
        vendPpro(inp)   = vstartPpro(inp) + nviPpro(inp)-1

        estartPpro(inp) = eendPpro(inp-1)+1
        eendPpro(inp)   = estartPpro(inp) + neiPpro(inp)-1

        fstartPpro(inp) = fendPpro(inp-1)+1
        fendPpro(inp)   = fstartPpro(inp) + nfiPpro(inp)-1
      end do
!     end initialising v,e,fstart end  


!      call read_geo_SingleBody(NparticleP,geofileP)
      call read_geoSingleBodyPpro(NparticlePpro,geofilePpro,nvtotPpro,netotPpro,nftotPpro,vstartPpro,vendPpro,estartPpro,eendPpro,fstartPpro,fendPpro,max_n_edge_of_vert, n_edge_of_vertPpro,vert_of_edgePpro, face_of_edgePpro,vert_of_facePpro, edge_of_facePpro ,vert_of_vertPpro, edge_of_vertPpro,xyzPpro,xyzvPpro,xyzaPpro,surPpro ,tri_verPpro,tri_velPpro,tri_barPpro,tri_norPpro,vel_triPpro,acc_triPpro)
      

      nbfxPpro=nftotPpro
      do l=1,nbfxPpro
          xyzbfxPpro(l,1:9) = tri_verPpro(1:9,l)
          barfxPpro(l,1:3) = tri_barPpro(1:3,l)
          areafxPpro(l) = surPpro(l)
          nxyzfxPpro(l,1:3)=tri_norPpro(1:3,l)
      enddo


      niix=1
      nifx=n1
      njix=1
      njfx=n2
      nkix=1
      nkfx=n3
        ! call ibgrid3dPpro(0,n1,n2,n3,xc,yc,zc &
        !               ,mpugeoPpro,nbfxPpro,xyzbfxPpro,nxyzfxPpro &
        !               ,barfxPpro,areafxPpro                &
        !               ,mpunPpro,mpunPpro,npunifxPpro,npunfxPpro    &
        !               ,indgeoeePpro                      &
        !               ,indgeoPpro,indgeoePpro,distbPpro      &
        !               ,.0,0.0,0.0                      &
        !               ,vbfxPpro,vbifxPpro,ntribfxPpro        & 
        !               ,niix,nifx,njix,njfx             & !noppro
        !               ,nkix,nkfx,PPro)           
!-------------------------------------------------------


      ! do n=1,npunifxPpro(1)
      !    i=indgeoeePpro(1,n,1)
      !    j=indgeoeePpro(1,n,2)
      !    k=indgeoeePpro(1,n,3)
      !    if(k .ge. kstart-1 .and. k .le. kend-1) then
      !       maskq1(i,j,k) = 1.
      !    endif
      ! enddo
      ! do n=1,npunifxPpro(2)
      !    i=indgeoeePpro(2,n,1)
      !    j=indgeoeePpro(2,n,2)
      !    k=indgeoeePpro(2,n,3)
      !    if(k .ge. kstart-1 .and. k .le. kend-1) then
      !       maskq2(i,j,k) = 1.
      !    endif
      ! enddo
      ! do n=1,npunifxPpro(3)
      !    i=indgeoeePpro(3,n,1)
      !    j=indgeoeePpro(3,n,2)
      !    k=indgeoeePpro(3,n,3)
      !    if(k .ge. kstart-1 .and. k .le. kend-1) then
      !       maskq3(i,j,k) = 1.
      !    endif
      ! enddo
      do n=1,npunifxPpro(4)
         i=indgeoeePpro(4,n,1)
         j=indgeoeePpro(4,n,2)
         k=indgeoeePpro(4,n,3)
         if(k .ge. kstart-1 .and. k .le. kend-1) then
            maskpr(i,j,k) = 1.
         endif
      enddo

call MPI_BARRIER(MPI_COMM_WORLD,ierr)           
      deallocate(n_edge_of_vertPpro)
      deallocate(vert_of_edgePpro, face_of_edgePpro)
      deallocate(vert_of_facePpro, edge_of_facePpro)
      deallocate(vert_of_vertPpro, edge_of_vertPpro)
      deallocate(xyzPpro,xyzvPpro,xyzaPpro)
      deallocate(surPpro)
      deallocate(tri_verPpro,tri_velPpro)
      deallocate(tri_barPpro,tri_norPpro,vel_triPpro,acc_triPpro)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)           
      return
      end subroutine  raytracing_Mask


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      subroutine write_to_vtk_Mask(ipfi)
      use param
      use mpi_param
      use mpih
      use local_arrays
      implicit none
     
      integer :: i,nslice
      character*70 filname
      character*100 ipfi
      character*2 id
      real(DP) :: sliceCoord, scarto



!trova slice                                                                                                                                                                 
      sliceCoord = 0.0D0
      scarto = 1000.
      do i=1,n2
         if (abs(yc(i)-sliceCoord).LT.scarto) then
            scarto=abs(yc(i)-sliceCoord)
            nslice = i
         endif
      enddo

!#ifdef USE_CUDA
!      q1(:,:,n3/2+1) = q1_d(:,:,n3/2+1)
!      q2(:,:,n3/2+1) = q2_d(:,:,n3/2+1)
!      q3(:,:,n3/2+1) = q3_d(:,:,n3/2+1)
!#endif
 
      write(id, '(i2.2)') myid
      filname = 'vtkfiles/field_01_'//id//'_'//trim(ipfi)//'.vtk'
     

!      filname = 'vtkfiles/field_01_'//id//'_'//ipfi//'.vtk'
      ! filname = 'vtkfiles/field_01_'//trim(ipfi)//'.vtk'

      open(121,file = filname)
!Header     
        write(121,'(a)')adjustl('# vtk DataFile Version 3.1')
        write(121,'(a)')adjustl('Stores full field data')
        write(121,'(a)')adjustl('ASCII')
        write(121,'(a)')adjustl('DATASET RECTILINEAR_GRID')
        write(121,'(a, 3I)')'DIMENSIONS ',n1, 1, kend-kstart+1
        write(121,*)'X_COORDINATES ',n1, ' FLOAT'
        write(121,*) xc
        write(121,*)'Y_COORDINATES ',1, ' FLOAT'
!        write(121,*) yc(n2/2 + 1)
        write(121,*) yc(nslice)
        write(121,*)'Z_COORDINATES ',kend-kstart+1, ' FLOAT'
        write(121,*) zc(kstart:kend)
        write(121,*)''
        
        write(121,*)''
        write(121,*)'POINT_DATA ', n1*1*(kend-kstart+1)
        write(121,*)'Scalars maskq1 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
!        write(121,*)q1(:,n2/2 + 1,kstart:kend)
        write(121,*)maskq1(:,nslice,kstart:kend)
        write(121,*)'Scalars maskq2 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
!        write(121,*)q2(:,n2/2 + 1,kstart:kend)
        write(121,*)maskq2(:,nslice,kstart:kend)
        write(121,*)'Scalars maskq3 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
        ! write(121,*)q3(:,n2/2 + 1,kstart:kend)
        write(121,*)maskq3(:,nslice,kstart:kend)
        write(121,*)'Scalars maskpr FLOAT'
        write(121,*)'LOOKUP_TABLE default'
        ! write(121,*)pr(:,n2/2 + 1,kstart:kend)
        write(121,*)maskpr(:,nslice,kstart:kend)

!        write(121,*)pr(:,n2/2 + 1,kstart:kend)

        ! write(121,*)''
        ! write(121,*)'POINT_DATA ', n1*1*(kend-kstart+1)
        ! write(121,*)'Scalars pr FLOAT'
        ! write(121,*)'LOOKUP_TABLE default'
        ! write(121,*)pr(:,n2/2 + 1,kstart:kend)

      close(121)



!Another one
!trova slice                                                                                                                                                                 
      sliceCoord = -1.60D0
      scarto = 1000.
      do i=1,n2
         if (abs(yc(i)-sliceCoord).LT.scarto) then
            scarto=abs(yc(i)-sliceCoord)
            nslice = i
         endif
      enddo

!#ifdef USE_CUDA
!      q1(:,:,n3/2+1) = q1_d(:,:,n3/2+1)
!      q2(:,:,n3/2+1) = q2_d(:,:,n3/2+1)
!      q3(:,:,n3/2+1) = q3_d(:,:,n3/2+1)
!#endif

!       write(id, '(i2.2)') myid
! !      filname = 'vtkfiles/field_02_'//id//'_'//ipfi//'.vtk'
!       filname = 'vtkfiles/field_02_'//trim(ipfi)//'.vtk'

      write(id, '(i2.2)') myid
      filname = 'vtkfiles/field_02_'//id//'_'//trim(ipfi)//'.vtk'

      open(121,file = filname)
!Header     
        write(121,'(a)')adjustl('# vtk DataFile Version 3.1')
        write(121,'(a)')adjustl('Stores full field data')
        write(121,'(a)')adjustl('ASCII')
        write(121,'(a)')adjustl('DATASET RECTILINEAR_GRID')
        write(121,'(a, 3I)')'DIMENSIONS ',n1, 1, kend-kstart+1
        write(121,*)'X_COORDINATES ',n1, ' FLOAT'
        write(121,*) xc
        write(121,*)'Y_COORDINATES ',1, ' FLOAT'
!        write(121,*) yc(n2/2 + 1)
        write(121,*) yc(nslice)
        write(121,*)'Z_COORDINATES ',kend-kstart+1, ' FLOAT'
        write(121,*) zc(kstart:kend)
        write(121,*)''
        
        write(121,*)''
        write(121,*)'POINT_DATA ', n1*1*(kend-kstart+1)
        write(121,*)'Scalars maskq1 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
!        write(121,*)q1(:,n2/2 + 1,kstart:kend)
        write(121,*)maskq1(:,nslice,kstart:kend)
        write(121,*)'Scalars maskq2 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
!        write(121,*)q2(:,n2/2 + 1,kstart:kend)
        write(121,*)maskq2(:,nslice,kstart:kend)
        write(121,*)'Scalars maskq3 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
        ! write(121,*)q3(:,n2/2 + 1,kstart:kend)
        write(121,*)maskq3(:,nslice,kstart:kend)
        write(121,*)'Scalars maskpr FLOAT'
        write(121,*)'LOOKUP_TABLE default'
        ! write(121,*)pr(:,n2/2 + 1,kstart:kend)
        write(121,*)maskpr(:,nslice,kstart:kend)

!        write(121,*)pr(:,n2/2 + 1,kstart:kend)

        ! write(121,*)''
        ! write(121,*)'POINT_DATA ', n1*1*(kend-kstart+1)
        ! write(121,*)'Scalars pr FLOAT'
        ! write(121,*)'LOOKUP_TABLE default'
        ! write(121,*)pr(:,n2/2 + 1,kstart:kend)

      close(121)


      return
      end

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
