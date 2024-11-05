!     -----------------------------------------------------------
!     Aux routines for reading in different geometries
!     -----------------------------------------------------------


!===================================================================
        subroutine read_geo_dim(nv,ne,nf,filename)
        
        implicit none
        character*50 :: filename
        integer :: nv,ne,nf

        open(11,file='meshes/'//trim(filename)//'.gts')
        read(11,*)nv,ne,nf
        close(11)

        
        return
        end subroutine read_geo_dim
!...................................................................
!===================================================================
        subroutine read_geo_dim_3dVTK(nv,ne,nf,nc,filename)
        use mls_param, only: max_n_edge_of_vert_3d

        implicit none

        character*50 :: filename, dammi1,dammi2,dammi3
        integer :: nv,ne,nf,nc, nfSur, numa,numb
        integer :: v1,v2,v3,v4,extsur
        integer :: i,cc,jj
        integer :: dummy1,dummy2,dummy3,dummy4,dummy5,dummy6
        integer :: vcc,vA,vB,vC,checkA,checkB,checkC
        integer :: countnf,countnfsur,ii,countCeqA
        integer :: jjA,jjB,jjC,cellA,cellB,cellC,cellFacing
        integer, dimension(:), allocatable:: n_edge_of_vert_loc
        integer, dimension(:,:), allocatable:: vert_of_vert_loc
        integer, dimension(:,:), allocatable:: vert_of_cell_loc
        integer, dimension(:), allocatable:: n_cell_of_vert_loc
        integer, dimension(:,:), allocatable:: cell_of_vert_loc

        open(109,file='meshes/'//trim(filename)//'.vtk')
        read(109,*)
        read(109,*)
        read(109,*)
        read(109,*)
        read(109,*) dammi1,nv,dammi2

        do i=1,nv
           read(109,*)
        end do

        read(109,*) dammi1,numa,numb
        close(109)
        nc = numb-4*numa
        nfSur = numa-nc
        write(*,*) nv,nc,nfSur


        if (nfSur.EQ.0) then
           extsur=0
        else
           extsur=1
           nf = nc*2 + nfSur/2; !number of faces
        endif


        allocate(vert_of_cell_loc(4,nc))
        allocate(vert_of_vert_loc(max_n_edge_of_vert_3d,nv))
        allocate(n_edge_of_vert_loc(nv))
        allocate(n_cell_of_vert_loc(nv))
        allocate(cell_of_vert_loc(max_n_edge_of_vert_3d,nv))

        n_edge_of_vert_loc = 0
        n_cell_of_vert_loc = 0

         !! read nodes coordinates
        open(109,file='meshes/'//trim(filename)//'.vtk')
        read(109,*)
        read(109,*)
        read(109,*)
        read(109,*)
        read(109,*)
        do i=1,nv
           read(109,*)
        end do        
        read(109,*) dammi1,dammi2,dammi3
         !! read surface triangles nodes (not used)
         if (extsur.EQ.1) then
         do i=1,nfSur
            read(109,*)
         end do
         endif

         !! read cell nodes 
         do i=1,nc
            read(109,*)dummy1,v1,v2,v3,v4
            v1=v1+1
            v2=v2+1
            v3=v3+1
            v4=v4+1
            vert_of_cell_loc(1,i)=v1
            vert_of_cell_loc(2,i)=v2
            vert_of_cell_loc(3,i)=v3
            vert_of_cell_loc(4,i)=v4

            n_cell_of_vert_loc(v1) = n_cell_of_vert_loc(v1)+1
            cell_of_vert_loc(n_cell_of_vert_loc(v1),v1)=i
            n_cell_of_vert_loc(v2) = n_cell_of_vert_loc(v2)+1
            cell_of_vert_loc(n_cell_of_vert_loc(v2),v2)=i
            n_cell_of_vert_loc(v3) = n_cell_of_vert_loc(v3)+1
            cell_of_vert_loc(n_cell_of_vert_loc(v3),v3)=i
            n_cell_of_vert_loc(v4) = n_cell_of_vert_loc(v4)+1
            cell_of_vert_loc(n_cell_of_vert_loc(v4),v4)=i

            do cc=1,4
               if (cc.EQ.1) then
                  vcc = v1; vA = v2; vB = v3; vC = v4;
               elseif (cc.EQ.2) then
                  vcc = v2; vA = v1; vB = v3; vC = v4;    
               elseif (cc.EQ.3) then
                  vcc = v3; vA = v1; vB = v2; vC = v4;    
               elseif (cc.EQ.4) then
                  vcc = v4; vA = v1; vB = v2; vC = v3;    
               endif

               checkA = 0; checkB = 0; checkC = 0;
               do jj = 1,n_edge_of_vert_loc(vcc)
                  if (vert_of_vert_loc(jj,vcc).EQ.vA) checkA = checkA + 1
                  if (vert_of_vert_loc(jj,vcc).EQ.vB) checkB = checkB + 1
                  if (vert_of_vert_loc(jj,vcc).EQ.vC) checkC = checkC + 1
               enddo
               if (checkA.EQ.0) then
                   n_edge_of_vert_loc(vcc) = n_edge_of_vert_loc(vcc)+1
                   vert_of_vert_loc(n_edge_of_vert_loc(vcc),vcc)=vA
               elseif (checkA.GT.1) then
                   write(*,*) "Error vert_of_vert";stop
               endif
               if (checkB.EQ.0) then
                   n_edge_of_vert_loc(vcc) = n_edge_of_vert_loc(vcc)+1
                   vert_of_vert_loc(n_edge_of_vert_loc(vcc),vcc)=vB
               elseif (checkB.GT.1) then
                  write(*,*) "Error vert_of_vert";stop
               endif
               if (checkC.EQ.0) then
                  n_edge_of_vert_loc(vcc) = n_edge_of_vert_loc(vcc)+1
                  vert_of_vert_loc(n_edge_of_vert_loc(vcc),vcc)=vC
               elseif (checkC.GT.1) then
                   write(*,*) "Error vert_of_vert";stop            
               endif

            enddo !cc

         end do ! cell nodes

         close(109)
! ! !   Completed reading the gts file
         ne = 0
         do i=1,nv
            ne = ne + n_edge_of_vert_loc(i)
         enddo
         ne = ne/2

!! conta facce per extsur=0
!! FACES
if (extsur.EQ.0) then
         countnf  =  0
         countnfSur = 0 !questo lo riazzero per ogni particle perche non lo uso
!         contafaccecell(:) = 0
         do ii=1,nc
            v1 = vert_of_cell_loc(1,ii)
            v2 = vert_of_cell_loc(2,ii)
            v3 = vert_of_cell_loc(3,ii)
            v4 = vert_of_cell_loc(4,ii) 

            do cc=1,4
               if (cc.EQ.1) then
               vA = v2; vB = v3; vC = v4; 
               elseif (cc.EQ.2) then
               vA = v1; vB = v3; vC = v4; 
               elseif (cc.EQ.3) then
               vA = v1; vB = v2; vC = v4; 
               elseif (cc.EQ.4) then
               vA = v1; vB = v2; vC = v3; 
               endif
    
               countCeqA = 0
               do jjA=1,n_cell_of_vert_loc(vA)
                  cellA = cell_of_vert_loc(jjA,vA)
                  do jjB=1,n_cell_of_vert_loc(vB)
                     cellB = cell_of_vert_loc(jjB,vB)
                     do jjC=1,n_cell_of_vert_loc(vC)
                        cellC = cell_of_vert_loc(jjC,vC)
                        if (cellC.EQ.cellA.AND.cellB.EQ.cellA.AND.cellA.NE.ii) then
                           countCeqA=countCeqA+1
                           cellFacing = cellA
                        endif
                     enddo
                  enddo
               enddo
    
               if (countCeqA.EQ.0) then !surface triangles
                  countnf = countnf+1
                  ! vert_of_face_3d(1,countnf) = vA
                  ! vert_of_face_3d(2,countnf) = vB
                  ! vert_of_face_3d(3,countnf) = vC
        
                  countnfSur = countnfSur+1;
                  !         vert_of_faceSur_3d(1,countnfSur) = vA;
                  !         vert_of_faceSur_3d(2,countnfSur) = vB;
                  !         vert_of_faceSur_3d(3,countnfSur) = vC;     
        
!                  contafaccecell(ii)=contafaccecell(ii)+1
                  ! face_of_cell_3d(contafaccecell(ii),ii)=countnf
                  ! cell_of_face_3d(1,countnf)=ii
                  ! cell_of_face_3d(2,countnf)=0

               elseif (countCeqA.EQ.1) then !internal triangle
                  if (cellFacing.GT.ii) then
                     countnf = countnf+1
!                      vert_of_face_3d(1,countnf) = vA
!                      vert_of_face_3d(2,countnf) = vB
!                      vert_of_face_3d(3,countnf) = vC        

!                      contafaccecell(ii)=contafaccecell(ii)+1
!                      face_of_cell_3d(contafaccecell(ii),ii)=countnf
!                      contafaccecell(cellFacing)=contafaccecell(cellFacing)+1
!                      face_of_cell_3d(contafaccecell(cellFacing),cellFacing)=countnf
!                      cell_of_face_3d(1,countnf)=ii
!                      cell_of_face_3d(2,countnf)=cellFacing
                  end if
               else
                  write(*,*) "Error faces 3d";stop
               endif
                          
                  
            enddo !  v1 v2 v3 v4  kk
         enddo ! cells ii

    nf=countnf
    nfSur=countnfSur
    write(*,*) "Tot. number Faces 3d",nf
    write(*,*) "Tot. number Faces 3d Sur",nfSur
endif        


        deallocate(vert_of_cell_loc)
        deallocate(vert_of_vert_loc)
        deallocate(n_edge_of_vert_loc)
        deallocate(n_cell_of_vert_loc)
        deallocate(cell_of_vert_loc)
        
        return
        end subroutine read_geo_dim_3dVTK
!...................................................................
!===================================================================                
      subroutine read_geoSingleBody_3dVTK
      use constants
      use param
      use mpih
!      use mpi_param, only: kstart, kend
      use mls_param
      !@cuf   use cudafor
      implicit none

      character*50 filename,strucfilename,dammi1,dammi2,dammi3
      character*100 ipfi,ipfip

      integer :: i,j,k,ii,jj,kk,v1,v2,v3,v4,inp,e1,e2,e3,count,countSTE,countSTElv,chamb
      real(DP) :: totmass,usVoldom,verchamb,s1,s2,s3,s4,supmin,supmax,ratios
      real(DP) :: xvDel,yvDel,zvDel,xDel,yDel,zDel,distvDel,tmaxdelay,distvardelay
      real(DP) :: fx,fy,fz,ex,ey,ez,cosAngFE
      real(DP) :: sum_sur,avg_sur
!      real(DP), dimension(Nparticle_3d) :: Surface0_3d
!      integer,dimension(:),allocatable :: nv_3d,ne_3d,nf_3d,nc_3d
      !commentare anche sotto
      integer :: v1d,v2d,v3d,v4d,extsur
      integer :: e1d,e2d,e3d,e4d
      integer :: f1,f2,f3,f4,seg1,seg2,seg3
      integer :: vsi,vei,esi,eei,fsi,fei,csi,cei
      integer :: nvT,nfsurT,neT,ncT,countne,numa,numb
      integer :: dummy,dummy1,dummy2,dummy3,dummy4,dummy5,dummy6
      integer :: cc,vcc,eA,vA,vB,vC,checkA,checkB,checkC,cellA,cellB,cellC
      integer :: jjA,jjB,jjC,countnf,countnfSur,countCeqA,cellFacing
      integer, dimension(nctot_3d) :: contafaccecell
      integer, dimension(17) :: count_seg
!@cuf   integer :: istat  
      !questo secondo me si puo commentare anche nel 2d
      ! do inp=1,Nparticle_3d
      !    open(109,file='meshes/'//geofile_3d(inp))
      !    read(109,*)nvi_3d(inp),nfsuri_3d(inp),nci_3d(inp)
      !    close(109)
      ! end do
!      stop
!     Allocations using vertices,edges and faces as parameters     
#ifdef USE_CUDA
      attributes(managed) :: contafaccecell
#endif

      call allocate_trigeo_3d
      LabelStenosi = 0
      n_edge_of_vert_3d = 0 !altri?
      n_cell_of_vert_3d = 0
      n_edge_of_cell_3d = 0
      n_cell_of_edge_3d = 0

!     Read in the vertices, edges and faces from the file
!     Also make connections between the three
      do inp=1,Nparticle_3d

         vsi = vstart_3d(inp) ; vei = vend_3d(inp)
         esi = estart_3d(inp) ; eei = eend_3d(inp)
         fsi = fstart_3d(inp) ; fei = fend_3d(inp)
         csi = cstart_3d(inp) ; cei = cend_3d(inp)

         ! Set particle labels
         vert_to_part_3d(vsi:vei) = inp
         edge_to_part_3d(esi:eei) = inp 
         face_to_part_3d(fsi:fei) = inp
         cell_to_part_3d(csi:cei) = inp

!open gts
        open(109,file='meshes/'//trim(geofile_3d(inp))//'.vtk')
        read(109,*)
        read(109,*)
        read(109,*)
        read(109,*)
        read(109,*) dammi1,nvT,dammi2

        do i=1,nvT
           read(109,*)
        end do

        read(109,*) dammi1,numa,numb
        close(109)
        ncT = numb-4*numa
        nfSurT = numa-ncT

        if (nfSurT.EQ.0) then
           extsur=0
        else
           extsur=1
!           nf = nc*2 + nfSur/2; !number of faces                                                                             
        endif

!         read nodes coordinates
        open(109,file='meshes/'//trim(geofile_3d(inp))//'.vtk')
        read(109,*)
        read(109,*)
        read(109,*)
        read(109,*)
        read(109,*)
        do i=vsi,vei
           read(109,*)xyz0_3d(1,i),xyz0_3d(2,i),xyz0_3d(3,i)
        end do
!        read(109,*)
        read(109,*) dammi1,dammi2,dammi3


!     position correction to bring in 0-6 !FV?
         xyzv0_3d = 0.0
         xyza0_3d = 0.0
         fxyz_3d  = 0.0      

         xyz_3d = xyz0_3d
         xyzv_3d = xyzv0_3d
         xyza_3d = xyza0_3d

         !! read surface triangles nodes (not used)
         if (extsur.EQ.1) then
         do i=1,nfsurT
            read(109,*)
         end do
         endif

!          !! read cell nodes 
         do i=csi,cei
            read(109,*)dummy1,v1d,v2d,v3d,v4d            

            v1 = v1d + vsi
            v2 = v2d + vsi
            v3 = v3d + vsi
            v4 = v4d + vsi

            vert_of_cell_3d(1,i)=v1
            vert_of_cell_3d(2,i)=v2
            vert_of_cell_3d(3,i)=v3
            vert_of_cell_3d(4,i)=v4

            n_cell_of_vert_3d(v1) = n_cell_of_vert_3d(v1)+1
            cell_of_vert_3d(n_cell_of_vert_3d(v1),v1)=i
            n_cell_of_vert_3d(v2) = n_cell_of_vert_3d(v2)+1
            cell_of_vert_3d(n_cell_of_vert_3d(v2),v2)=i
            n_cell_of_vert_3d(v3) = n_cell_of_vert_3d(v3)+1
            cell_of_vert_3d(n_cell_of_vert_3d(v3),v3)=i
            n_cell_of_vert_3d(v4) = n_cell_of_vert_3d(v4)+1
            cell_of_vert_3d(n_cell_of_vert_3d(v4),v4)=i

            do cc=1,4
               if (cc.EQ.1) then
                  vcc = v1; vA = v2; vB = v3; vC = v4;
               elseif (cc.EQ.2) then
                  vcc = v2; vA = v1; vB = v3; vC = v4;    
               elseif (cc.EQ.3) then
                  vcc = v3; vA = v1; vB = v2; vC = v4;    
               elseif (cc.EQ.4) then
                  vcc = v4; vA = v1; vB = v2; vC = v3;    
               endif

               checkA = 0; checkB = 0; checkC = 0;
               do jj = 1,n_edge_of_vert_3d(vcc)
                  if (vert_of_vert_3d(jj,vcc).EQ.vA) checkA = checkA + 1
                  if (vert_of_vert_3d(jj,vcc).EQ.vB) checkB = checkB + 1
                  if (vert_of_vert_3d(jj,vcc).EQ.vC) checkC = checkC + 1
               enddo
               if (checkA.EQ.0) then
                   n_edge_of_vert_3d(vcc) = n_edge_of_vert_3d(vcc)+1
                   vert_of_vert_3d(n_edge_of_vert_3d(vcc),vcc)=vA
               elseif (checkA.GT.1) then
                   write(*,*) "Error vert_of_vert";stop
               endif
               if (checkB.EQ.0) then
                   n_edge_of_vert_3d(vcc) = n_edge_of_vert_3d(vcc)+1
                   vert_of_vert_3d(n_edge_of_vert_3d(vcc),vcc)=vB
               elseif (checkB.GT.1) then
                  write(*,*) "Error vert_of_vert";stop
               endif
               if (checkC.EQ.0) then
                  n_edge_of_vert_3d(vcc) = n_edge_of_vert_3d(vcc)+1
                  vert_of_vert_3d(n_edge_of_vert_3d(vcc),vcc)=vC
               elseif (checkC.GT.1) then
                   write(*,*) "Error vert_of_vert";stop        
               endif


            enddo !cc

         end do ! cell nodes
#ifndef GEO_ANSA
         read(109,*) !2righe? !FV dammi1,dammi2,dammi3?
         read(109,*)
         do i=1,nfSurT
            read(109,*) dummy
!            write(*,*) dummy
         enddo
         do i=csi,cei
            read(109,*) dummy
!            write(*,*) dummy
         enddo
         read(109,*)
         read(109,*)
         read(109,*) 
         do i=vsi,vei
            read(109,*) verchamb
            vert_to_chamb_3d4V(i)=nint(verchamb)
         enddo
         read(109,*)
         read(109,*) 
         do i=vsi,vei
            read(109,*) verchamb
            ! vert_to_chamb_3d(i)=nint(verchamb) !uno dei due
         enddo
         read(109,*)
         read(109,*) 
         do i=vsi,vei
            read(109,*) verchamb
            vert_to_chamb_3d(i)=nint(verchamb) !uno dei due
         enddo
         read(109,*)
         read(109,*) 
         do i=vsi,vei
           read(109,*) LabelSmooth(i)
!            LabelSmooth(i)=1.0D0
        enddo
#else
        !open(209,file='meshes/septum_cell.txt')
        do i=csi,cei
           !read(209,*) verchamb
           !         cell_to_chamb_3d(i)=9!nint(verchamb) !uno dei due
           cell_to_chamb_3d(i)=1!nint(1)
        enddo
        !close(209)
#endif


#ifdef FIBERS
!#ifndef CARTO
        !Fibers nodes
        open(109,file='meshes/'//trim(geofile_3d(inp))//'_fibers.txt')
        read(109,*) !fiberx
        read(109,*) !fiberx
        do i=1,nvtot_3d
           read(109,*) AmatrFibers_node_3d(1,1,i)
           !AmatrFibers_node_3d(1,1,i) = 1.d0
        end do
        read(109,*) !fibery
        read(109,*) !fibery
        do i=1,nvtot_3d
           read(109,*) AmatrFibers_node_3d(2,1,i)
           !AmatrFibers_node_3d(2,1,i) = 0.d0
        end do
        read(109,*) !fiberz
        read(109,*) !fiberz
        do i=1,nvtot_3d
           read(109,*) AmatrFibers_node_3d(3,1,i)
           !AmatrFibers_node_3d(3,1,i) = 0.d0
        end do
        read(109,*) !normalsheetx
        read(109,*) !normalsheetx
        do i=1,nvtot_3d
           read(109,*) AmatrFibers_node_3d(1,3,i)
           !AmatrFibers_node_3d(1,3,i) = 0.d0
        end do
        read(109,*) !normalsheety
        read(109,*) !normalsheety
        do i=1,nvtot_3d
           read(109,*) AmatrFibers_node_3d(2,3,i)
           !AmatrFibers_node_3d(2,3,i) = 0.d0
        end do
        read(109,*) !normalsheetz
        read(109,*) !normalsheetz        
        do i=1,nvtot_3d
           read(109,*) AmatrFibers_node_3d(3,3,i)
           !AmatrFibers_node_3d(3,3,i) = 1.d0
        end do
        read(109,*) !sheetx
        read(109,*) !sheetx
        do i=1,nvtot_3d
           read(109,*) AmatrFibers_node_3d(1,2,i)
           !AmatrFibers_node_3d(1,2,i) = 0.d0
        end do
        read(109,*) !sheety
        read(109,*) !sheety
        do i=1,nvtot_3d
           read(109,*) AmatrFibers_node_3d(2,2,i)
           !AmatrFibers_node_3d(2,2,i) = 1.d0
        end do
        read(109,*) !sheetz
        read(109,*) !sheetz
        do i=1,nvtot_3d
           read(109,*) AmatrFibers_node_3d(3,2,i)
           !AmatrFibers_node_3d(3,2,i) = 0.d0
        end do
        ! read(109,*) !normalsheetx
        ! read(109,*) !normalsheetx
        ! do i=1,nvtot_3d
        !    read(109,*) AmatrFibers_node_3d(1,3,i)
        !    !AmatrFibers_node_3d(1,3,i) = 0.d0
        ! end do
        ! read(109,*) !normalsheety
        ! read(109,*) !normalsheety
        ! do i=1,nvtot_3d
        !    read(109,*) AmatrFibers_node_3d(2,3,i)
        !    !AmatrFibers_node_3d(2,3,i) = 0.d0
        ! end do
        ! read(109,*) !normalsheetz
        ! read(109,*) !normalsheetz        
        ! do i=1,nvtot_3d
        !    read(109,*) AmatrFibers_node_3d(3,3,i)
        !    !AmatrFibers_node_3d(3,3,i) = 1.d0
        ! end do

        close(109)
! #else
!         !$cuf kernel do (1) 
!         do i=1,nvtot_3d
! !           read(109,*) AmatrFibers_node_3d(1,1,i)
!            AmatrFibers_node_3d(1,1,i) = 1.d0
!         end do
!         !@cuf istat = cudaDeviceSynchronize !JDR TMP
!         !$cuf kernel do (1) 
!         do i=1,nvtot_3d
!            ! read(109,*) AmatrFibers_node_3d(2,1,i)
!            AmatrFibers_node_3d(2,1,i) = 0.d0
!         end do
!         !@cuf istat = cudaDeviceSynchronize !JDR TMP
!         !$cuf kernel do (1) 
!         do i=1,nvtot_3d
!            ! read(109,*) AmatrFibers_node_3d(3,1,i)
!            AmatrFibers_node_3d(3,1,i) = 0.d0
!         end do
!         !@cuf istat = cudaDeviceSynchronize !JDR TMP
!         !$cuf kernel do (1) 
!         do i=1,nvtot_3d
!            ! read(109,*) AmatrFibers_node_3d(1,2,i)
!            AmatrFibers_node_3d(1,2,i) = 0.d0
!         end do
!         !@cuf istat = cudaDeviceSynchronize !JDR TMP
!         !$cuf kernel do (1) 
!         do i=1,nvtot_3d
!            ! read(109,*) AmatrFibers_node_3d(2,2,i)
!            AmatrFibers_node_3d(2,2,i) = 1.d0
!         end do
!         !@cuf istat = cudaDeviceSynchronize !JDR TMP
!         !$cuf kernel do (1) 
!         do i=1,nvtot_3d
!            ! read(109,*) AmatrFibers_node_3d(3,2,i)
!            AmatrFibers_node_3d(3,2,i) = 0.d0
!         end do
!         !@cuf istat = cudaDeviceSynchronize !JDR TMP
!         !$cuf kernel do (1) 
!         do i=1,nvtot_3d
!            ! read(109,*) AmatrFibers_node_3d(1,3,i)
!            AmatrFibers_node_3d(1,3,i) = 0.d0
!         end do
!         !@cuf istat = cudaDeviceSynchronize !JDR TMP
!         !$cuf kernel do (1) 
!         do i=1,nvtot_3d
!            ! read(109,*) AmatrFibers_node_3d(2,3,i)
!            AmatrFibers_node_3d(2,3,i) = 0.d0
!         end do
!         !@cuf istat = cudaDeviceSynchronize !JDR TMP
!         !$cuf kernel do (1) 
!         do i=1,nvtot_3d
!            ! read(109,*) AmatrFibers_node_3d(3,3,i)
!            AmatrFibers_node_3d(3,3,i) = 1.d0
!         end do
!         !@cuf istat = cudaDeviceSynchronize !JDR TMP
!#endif
        !Fibers cells
        !$cuf kernel do (1) 
        do i=1,nctot_3d
           v1=vert_of_cell_3d(1,i)
           v2=vert_of_cell_3d(2,i)
           v3=vert_of_cell_3d(3,i)
           v4=vert_of_cell_3d(4,i)

           do j=1,3
              do k=1,3
                 AmatrFibers_cell_3d(j,k,i)=0.25D0*(AmatrFibers_node_3d(j,k,v1)+AmatrFibers_node_3d(j,k,v2)+AmatrFibers_node_3d(j,k,v3)+AmatrFibers_node_3d(j,k,v4))        
              enddo
           enddo

        enddo
        !@cuf istat = cudaDeviceSynchronize !JDR TMP
#else
        call calculate_fiberdirection(vsi,vei,csi,cei,xyz0_3d(:,vsi:vei),  &
                       cell_bar(:,csi:cei),AmatrFibers_cell_3d(:,:,csi:cei))
#endif


        
         close(109)

!! EDGES
         countne  =  esi -1
         do ii=vsi,vei
            do jj=1,n_edge_of_vert_3d(ii)
               v1 = vert_of_vert_3d(jj,ii)
               if (v1.GT.ii) then
               countne = countne + 1;

               vert_of_edge_3d(1,countne)=ii 
               vert_of_edge_3d(2,countne)=v1

               do kk = 1,n_cell_of_vert_3d(ii)
                  cellA = cell_of_vert_3d(kk,ii)
                  do cc=1,4
                     v2 = vert_of_cell_3d(cc,cellA)
                     if(v1.EQ.v2) then
                     n_edge_of_cell_3d(cellA)=n_edge_of_cell_3d(cellA)+1
                     edge_of_cell_3d(n_edge_of_cell_3d(cellA),cellA) = countne 
                     n_cell_of_edge_3d(countne)=n_cell_of_edge_3d(countne)+1
                     cell_of_edge_3d(n_cell_of_edge_3d(countne),countne)=cellA
                     endif
                  enddo
               enddo
               endif
            end do
         enddo

         if (countne.NE.eei) then
            write(*,*) "Error edges of particle" 
            stop
         endif
         !$cuf kernel do (1)
         do i=csi,cei
            if (n_edge_of_cell_3d(i).NE.6) then
               write(*,*) "Error edges of cell",n_edge_of_cell_3d(i)
               stop
            endif
         enddo
         !@cuf istat = cudaDeviceSynchronize !JDR TMP 




!! EDGE_OF_VERT
         count_edge_of_vert_3d=0
         do ii=esi,eei
            v1 = vert_of_edge_3d(1,ii) 
            v2 = vert_of_edge_3d(2,ii) 
            count = 0
            do jj=1,n_edge_of_vert_3d(v1)
               vA = vert_of_vert_3d(jj,v1)               
               if (vA.EQ.v2) then
                  edge_of_vert_3d(jj,v1)=ii
                  count_edge_of_vert_3d(v1)=count_edge_of_vert_3d(v1)+1
                  count = count +1
               endif
               if (count.GT.1) then
                  write(*,*) "Error edge_of_vert";stop
               endif
            enddo

            count = 0
            do jj=1,n_edge_of_vert_3d(v2)
               vB = vert_of_vert_3d(jj,v2)               
               if (vB.EQ.v1) then
                  edge_of_vert_3d(jj,v2)=ii
                  count_edge_of_vert_3d(v2)=count_edge_of_vert_3d(v2)+1
                  count = count + 1
               endif
               if (count.GT.1) then
                  write(*,*) "Error edge_of_vert";stop
               endif
            enddo
         enddo
         !%check
         !$cuf kernel do (1)
         do ii=vsi,vei
            if (count_edge_of_vert_3d(ii).NE.n_edge_of_vert_3d(ii)) then
               write(*,*) "Error edge_of_vert";stop
            endif
         enddo
         !@cuf istat = cudaDeviceSynchronize !JDR TMP
         
         !check edge_of_vert vert_of_vert
!         !$cuf kernel do (1) 
         do ii=vsi,vei
            do jj=1,n_edge_of_vert_3d(ii)
               vA = vert_of_vert_3d(jj,ii)
               eA = edge_of_vert_3d(jj,ii)
               
               v1 = vert_of_edge_3d(1,eA) 
               v2 = vert_of_edge_3d(2,eA) 

               if (((v1.EQ.ii).AND.(v2.EQ.vA)).OR.((v1.EQ.vA).AND.(v2.EQ.ii))) then
                  !test passed
               else
                  write(*,*) "Error edge_of_vert vert_of_vert";stop
               endif
               
            enddo
         enddo
!         !@cuf istat = cudaDeviceSynchronize !JDR TMP 

         
#ifdef GEO_ANSA
        
        !$cuf kernel do (1)
        do i=csi,cei
           
           v1=vert_of_cell_3d(1,i)
           v2=vert_of_cell_3d(2,i)
           v3=vert_of_cell_3d(3,i)
           v4=vert_of_cell_3d(4,i)
           LabelSmooth(v1)=1
           LabelSmooth(v2)=1
           LabelSmooth(v3)=1
           LabelSmooth(v4)=1
           if (cell_to_chamb_3d(i).NE.9)then
              LabelSmooth(v1)=1
              LabelSmooth(v2)=1
              LabelSmooth(v3)=1
              LabelSmooth(v4)=1
           endif

           if (cell_to_chamb_3d(i).EQ.1)then
              vert_to_chamb_3d(v1)=1
              vert_to_chamb_3d(v2)=1
              vert_to_chamb_3d(v3)=1
              vert_to_chamb_3d(v4)=1
              vert_to_chamb_3d4V(v1)=1
              vert_to_chamb_3d4V(v2)=1
              vert_to_chamb_3d4V(v3)=1
              vert_to_chamb_3d4V(v4)=1
           elseif (cell_to_chamb_3d(i).EQ.2)then
              vert_to_chamb_3d(v1)=2
              vert_to_chamb_3d(v2)=2
              vert_to_chamb_3d(v3)=2
              vert_to_chamb_3d(v4)=2
              vert_to_chamb_3d4V(v1)=2
              vert_to_chamb_3d4V(v2)=2
              vert_to_chamb_3d4V(v3)=2
              vert_to_chamb_3d4V(v4)=2
           elseif (cell_to_chamb_3d(i).EQ.3)then
              vert_to_chamb_3d(v1)=3
              vert_to_chamb_3d(v2)=3
              vert_to_chamb_3d(v3)=3
              vert_to_chamb_3d(v4)=3
              vert_to_chamb_3d4V(v1)=3
              vert_to_chamb_3d4V(v2)=3
              vert_to_chamb_3d4V(v3)=3
              vert_to_chamb_3d4V(v4)=3
           elseif (cell_to_chamb_3d(i).EQ.4)then
              vert_to_chamb_3d(v1)=4
              vert_to_chamb_3d(v2)=4
              vert_to_chamb_3d(v3)=4
              vert_to_chamb_3d(v4)=4
              vert_to_chamb_3d4V(v1)=4
              vert_to_chamb_3d4V(v2)=4
              vert_to_chamb_3d4V(v3)=4
              vert_to_chamb_3d4V(v4)=4
           elseif (cell_to_chamb_3d(i).EQ.5)then
              vert_to_chamb_3d(v1)=5
              vert_to_chamb_3d(v2)=5
              vert_to_chamb_3d(v3)=5
              vert_to_chamb_3d(v4)=5
              vert_to_chamb_3d4V(v1)=5
              vert_to_chamb_3d4V(v2)=5
              vert_to_chamb_3d4V(v3)=5
              vert_to_chamb_3d4V(v4)=5
           elseif (cell_to_chamb_3d(i).EQ.6)then
              vert_to_chamb_3d(v1)=6
              vert_to_chamb_3d(v2)=6
              vert_to_chamb_3d(v3)=6
              vert_to_chamb_3d(v4)=6
              vert_to_chamb_3d4V(v1)=6
              vert_to_chamb_3d4V(v2)=6
              vert_to_chamb_3d4V(v3)=6
              vert_to_chamb_3d4V(v4)=6
           elseif (cell_to_chamb_3d(i).EQ.7)then
              vert_to_chamb_3d(v1)=7
              vert_to_chamb_3d(v2)=7
              vert_to_chamb_3d(v3)=7
              vert_to_chamb_3d(v4)=7
              vert_to_chamb_3d4V(v1)=7
              vert_to_chamb_3d4V(v2)=7
              vert_to_chamb_3d4V(v3)=7
              vert_to_chamb_3d4V(v4)=7
           elseif (cell_to_chamb_3d(i).EQ.8)then
              vert_to_chamb_3d(v1)=8
              vert_to_chamb_3d(v2)=8
              vert_to_chamb_3d(v3)=8
              vert_to_chamb_3d(v4)=8
              vert_to_chamb_3d4V(v1)=8
              vert_to_chamb_3d4V(v2)=8
              vert_to_chamb_3d4V(v3)=8
              vert_to_chamb_3d4V(v4)=8
           elseif (cell_to_chamb_3d(i).EQ.9)then
              vert_to_chamb_3d(v1)=9
              vert_to_chamb_3d(v2)=9
              vert_to_chamb_3d(v3)=9
              vert_to_chamb_3d(v4)=9
              vert_to_chamb_3d4V(v1)=1
              vert_to_chamb_3d4V(v2)=1
              vert_to_chamb_3d4V(v3)=1
              vert_to_chamb_3d4V(v4)=1
              
           endif
           
        enddo
        !@cuf istat = cudaDeviceSynchronize !JDR TMP
        
#endif

!! FACES
         countnf  =  fsi -1
         countnfSur = 0 !questo lo riazzero per ogni particle perche non lo uso
         contafaccecell(:) = 0
         do ii=csi,cei
            v1 = vert_of_cell_3d(1,ii)
            v2 = vert_of_cell_3d(2,ii)
            v3 = vert_of_cell_3d(3,ii)
            v4 = vert_of_cell_3d(4,ii)    

            LabelSmooth_cell(ii)=0.25*(LabelSmooth(v1)+LabelSmooth(v2)+LabelSmooth(v3)+LabelSmooth(v4))   
            do cc=1,4
               if (cc.EQ.1) then
               vA = v2; vB = v3; vC = v4; 
               elseif (cc.EQ.2) then
               vA = v1; vB = v3; vC = v4; 
               elseif (cc.EQ.3) then
               vA = v1; vB = v2; vC = v4; 
               elseif (cc.EQ.4) then
               vA = v1; vB = v2; vC = v3; 
               endif
    
               countCeqA = 0
               do jjA=1,n_cell_of_vert_3d(vA)
                  cellA = cell_of_vert_3d(jjA,vA)
                  do jjB=1,n_cell_of_vert_3d(vB)
                     cellB = cell_of_vert_3d(jjB,vB)
                     do jjC=1,n_cell_of_vert_3d(vC)
                        cellC = cell_of_vert_3d(jjC,vC)
                        if (cellC.EQ.cellA.AND.cellB.EQ.cellA.AND.cellA.NE.ii) then
                           countCeqA=countCeqA+1
                           cellFacing = cellA
                        endif
                     enddo
                  enddo
               enddo
    
               if (countCeqA.EQ.0) then !surface triangles
                  countnf = countnf+1
                  vert_of_face_3d(1,countnf) = vA
                  vert_of_face_3d(2,countnf) = vB
                  vert_of_face_3d(3,countnf) = vC
        
                  countnfSur = countnfSur+1;
                  !         vert_of_faceSur_3d(1,countnfSur) = vA;
                  !         vert_of_faceSur_3d(2,countnfSur) = vB;
                  !         vert_of_faceSur_3d(3,countnfSur) = vC;     
        
                  contafaccecell(ii)=contafaccecell(ii)+1
                  face_of_cell_3d(contafaccecell(ii),ii)=countnf
                  cell_of_face_3d(1,countnf)=ii
                  cell_of_face_3d(2,countnf)=0

               elseif (countCeqA.EQ.1) then !internal triangle
                  if (cellFacing.GT.ii) then
                     countnf = countnf+1
                     vert_of_face_3d(1,countnf) = vA
                     vert_of_face_3d(2,countnf) = vB
                     vert_of_face_3d(3,countnf) = vC        

                     contafaccecell(ii)=contafaccecell(ii)+1
                     face_of_cell_3d(contafaccecell(ii),ii)=countnf
                     contafaccecell(cellFacing)=contafaccecell(cellFacing)+1
                     face_of_cell_3d(contafaccecell(cellFacing),cellFacing)=countnf
                     cell_of_face_3d(1,countnf)=ii
                     cell_of_face_3d(2,countnf)=cellFacing
                  end if
               else
                  write(*,*) "Error faces 3d";stop
               endif
                          
                  
            enddo !  v1 v2 v3 v4  kk
         enddo ! cells ii
        
!%check 
         if (countnf.NE.(fei)) then; write(*,*) "Error number of faces 3d";stop;endif
         if (nfSurT.NE.countnfSur) write(*,*) "Error number of surface faces 3d"
         write(*,*) "Check nfSurT",nfSurT
         write(*,*) "Check countnfSur",countnfSur
         
         !$cuf kernel do (1) 
         do ii=csi,cei
            if (contafaccecell(ii).NE.4) then
               write(*,*) "Error number of faces 3d", contafaccecell(ii) !Here add something for shell writing in parallel
               stop
            endif
         enddo
         !@cuf istat = cudaDeviceSynchronize !JDR TMP 
         


!initialize stuff
        call calculate_anglealp(aalpha0_3d(:,fsi:fei),vsi,vei,fsi,fei,  &
                       xyz0_3d(:,vsi:vei),                            &
                       vert_of_face_3d(:,fsi:fei))

        call calculate_distance(dist0_3d(esi:eei),vsi,vei,esi,eei,    &
                       xyz0_3d(:,vsi:vei),vert_of_edge_3d(:,esi:eei))

        call calculate_volume_cells(vsi,vei,csi,cei,xyz0_3d(:,vsi:vei),  &
                       vert_of_cell_3d(:,csi:cei),vol0_3d(csi:cei))
        vol_3d=vol0_3d

        call calculate_bar_cells(vsi,vei,csi,cei,xyz0_3d(:,vsi:vei),  &
                       vert_of_cell_3d(:,csi:cei),cell_bar(:,csi:cei))


        call calculate_area(Surface0_3d(inp),vsi,vei,fsi,fei,         &
                       xyz0_3d(:,vsi:vei),vert_of_face_3d(:,fsi:fei),    &
                       sur0_3d(fsi:fei))
        sur_3d=sur0_3d

        call calculate_ginterp_face_cells(vsi,vei,fsi,fei,csi,cei,         &
                       xyz0_3d(:,vsi:vei),vert_of_face_3d(:,fsi:fei),    &
                       cell_of_face_3d(:,fsi:fei),cell_bar(:,csi:cei), &
                       versCFface_3d(:,fsi:fei),distCFface_3d(fsi:fei),g1interpface_3d(fsi:fei))

        call calculate_normals_face_cells(vsi,vei,fsi,fei,csi,cei,xyz0_3d(:,vsi:vei),  &
                       vert_of_face_3d(:,fsi:fei),face_of_cell_3d(:,csi:cei),   & 
                       cell_bar(:,csi:cei),normalfaceofcells_3d(:,:,csi:cei))




        !Edges cells
        !$cuf kernel do (1)  
        do i=1,netot_3d
           v1=vert_of_edge_3d(1,i)
           v2=vert_of_edge_3d(2,i)
           do j=1,3
              do k=1,3
                 AmatrFibers_edge_3d(j,k,i)=0.5D0*(AmatrFibers_node_3d(j,k,v1)+AmatrFibers_node_3d(j,k,v2))
              enddo
           enddo
        enddo
        !@cuf istat = cudaDeviceSynchronize !JDR TMP 

        !edge-fiber angle                                                                                                                      !$cuf kernel do (1)               
        do i = 1,netot_3d !add angle dependence
           
           fx = AmatrFibers_edge_3d(1,1,i)
           fy = AmatrFibers_edge_3d(2,1,i)
           fz = AmatrFibers_edge_3d(3,1,i)
           
           v1 = vert_of_edge_3d(1,i)
           v2 = vert_of_edge_3d(2,i)
           ex = (xyz0_3d(1,v1)-xyz0_3d(1,v2))/dist0_3d(i)
           ey = (xyz0_3d(2,v1)-xyz0_3d(2,v2))/dist0_3d(i)
           ez = (xyz0_3d(3,v1)-xyz0_3d(3,v2))/dist0_3d(i)
           
           cosAngFE = abs(ex*fx + ey*fy + ez*fz)
           edgefiber_cosangle_3d(i) = cosAngFE
        enddo
        !@cuf istat = cudaDeviceSynchronize !JDR TMP

      !SEGMENTS HEART
      ! xyz_seg(:,:) = 0.D0
      ! count_seg(:)=0
      !   open(109,file='meshes/'//trim(geofile_3d(inp))//'_segments.txt')
      !   do i=1,nvtot_3d
      !      read(109,*) verchamb
      !      seg1=nint(verchamb)
      !      Segments_node_3d(i)=seg1
      !      if (seg1.GE.1) then
      !         count_seg(seg1)=count_seg(seg1)+1
      !         xyz_seg(1:3,seg1)=xyz_seg(1:3,seg1)+xyz_3d(1:3,i)
      !      endif
      !   end do
      !   close(109)
      !   do seg1=1,17
      !      xyz_seg(1:3,seg1)=xyz_seg(1:3,seg1)/real(count_seg(seg1))
      !      write(*,*) "Center seg",seg1,xyz_seg(1:3,seg1)
      !   enddo

       ! stop
!face
      ! do f1=fstart_3d(1),fend_3d(1)
      !    v1 = vert_of_face_3d(1,f1)
      !    v2 = vert_of_face_3d(2,f1)
      !    v3 = vert_of_face_3d(3,f1)
      !    seg1=Segments_node_3d(v1)
      !    seg2=Segments_node_3d(v2)
      !    seg3=Segments_node_3d(v3)
      !    if ((seg1.EQ.seg2).AND.(seg2.EQ.seg3)) then
      !       Segments_face_3d(f1)=seg1
      !    elseif ((seg1.EQ.seg2).AND.(seg2.NE.seg3)) then
      !       Segments_face_3d(f1)=seg1
      !    elseif ((seg1.EQ.seg3).AND.(seg2.NE.seg3)) then
      !       Segments_face_3d(f1)=seg1
      !    elseif ((seg2.EQ.seg3).AND.(seg2.NE.seg1)) then
      !       Segments_face_3d(f1)=seg2
      !    elseif ((seg2.NE.seg3).AND.(seg2.NE.seg1)) then
      !       Segments_face_3d(f1)=MIN(seg1,seg2,seg3)
      !    endif
      ! enddo


      ! !AVG SUR
      ! sum_sur=0.0D0
      ! countSTElv=0
      ! do i=fsi,fei
      !    chamb=face_to_chamb_3d(i)
      !    if ((chamb.EQ.1).or.(chamb.EQ.3)) then
      !       sum_sur=sum_sur+sur_3d(i)
      !       countSTElv=countSTElv+1
      !    endif
      ! enddo
      ! avg_sur = sum_sur/countSTElv

      ! open(109,file='meshes/c2c_3d_check.txt')
      ! do i=csi,cei
      !    write(109,*)cell_to_chamb_3d(i)
      ! enddo
      ! close(109)

      ! !LABEL SKEWED TRET
      ! countSTE=0
      ! countSTElv=0
      !   do i=csi,cei
      !      f1=face_of_cell_3d(1,i)
      !      f2=face_of_cell_3d(2,i)
      !      f3=face_of_cell_3d(3,i)
      !      f4=face_of_cell_3d(4,i)

      !      chamb=cell_to_chamb_3d(i)
      !      s1=sur_3d(f1)
      !      s2=sur_3d(f2)
      !      s3=sur_3d(f3)
      !      s4=sur_3d(f4)

      !      supmin=min(s1,s2,s3,s4)
      !      supmax=max(s1,s2,s3,s4)

      !      ratios=supmin/supmax
      !     ! if ((ratios.LT.0.10).OR.(supmin.LT.3E-3)) then                                 
      !      if ((supmin.LT.1.0E-3)) then                                 
      !         LabelStenosi(i)=1
      !         countSTE=countSTE+1
      !         if ((cell_to_chamb_3d(i).EQ.1).or.(cell_to_chamb_3d(i).EQ.3)) then
      !            countSTElv=countSTElv+1
      !         endif
      !      endif
      !   enddo
      !   write(*,*) "countSTE = ", countSTE, countSTElv, nctot_3d
        
      end do !inp
    
      return
      end subroutine read_geoSingleBody_3dVTK
!...................................................................
!===================================================================                         
      subroutine body1_to_chamb
      use constants
      use param
      use mpih
!      use mpi_param, only: kstart, kend
      use mls_param
      !@cuf   use cudafor 
      implicit none
!@cuf integer :: istat
      integer :: i,j,k,v1,v2,v3,v4,e1,e2,f1,c1,jj
      integer :: vsi,vei,esi,eei,fsi,fei,loop
      integer :: chamb1,chamb2,chamb3,chamb4
      real(DP) :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xB,yB,zB
      real(DP) :: x0,y0,z0,alp,xsec,ysec,zsec,dblocca,dblocca2
      real(DP) :: coeff1,coeff2,s1,s2,s3,s4,supmin,supmax,ratios
      integer :: checkLV,checkLA,countSTE,chambcheck
      integer :: f2,f3,f4,chamb
      integer, dimension(9) :: countSTEvett

      checkLV=0
      checkLA=0
!     Read in the vertices, edges and faces from the file
!     Also make connections between the three

!vert_to_chamb (2d e 3d) ci sono gia
!mi serve face_to_chamb_2d e cell_to_chamb_3d

!------------2D-------------------
!#ifndef GEO_ANSA !HERE GENERATE INPUT_FILE fROM ANSA 
      !edges
      !$cuf kernel do (1)
      do e1=estart(1),eend(1)
         v1 = vert_of_edge(1,e1)
         v2 = vert_of_edge(2,e1)
         chamb1=vert_to_chamb(v1)
         chamb2=vert_to_chamb(v2)
         if (chamb1.EQ.chamb2) then
            edge_to_chamb(e1)=chamb1 !setto ventricolo
         else
            edge_to_chamb(e1)=MIN(chamb1,chamb2)
         endif
      enddo
      !@cuf istat = cudaDeviceSynchronize !JDR TMP
      
      !face
      !$cuf kernel do (1) 
      do f1=fstart(1),fend(1)
         v1 = vert_of_face(1,f1)
         v2 = vert_of_face(2,f1)
         v3 = vert_of_face(3,f1)
         chamb1=vert_to_chamb(v1)
         chamb2=vert_to_chamb(v2)
         chamb3=vert_to_chamb(v3)
         if ((chamb1.EQ.chamb2).AND.(chamb2.EQ.chamb3)) then
            face_to_chamb(f1)=chamb1 
         elseif ((chamb1.EQ.chamb2).AND.(chamb2.NE.chamb3)) then
            face_to_chamb(f1)=chamb1
         elseif ((chamb1.EQ.chamb3).AND.(chamb2.NE.chamb3)) then
            face_to_chamb(f1)=chamb1
         elseif ((chamb2.EQ.chamb3).AND.(chamb2.NE.chamb1)) then
            face_to_chamb(f1)=chamb2
         elseif ((chamb2.NE.chamb3).AND.(chamb2.NE.chamb1)) then
            face_to_chamb(f1)=MIN(chamb1,chamb2,chamb3)
         endif
      enddo
      !@cuf istat = cudaDeviceSynchronize !JDR TMP
      
      !face4V
      !$cuf kernel do (1) 
      do f1=fstart(1),fend(1)
         v1 = vert_of_face(1,f1)
         v2 = vert_of_face(2,f1)
         v3 = vert_of_face(3,f1)
         chamb1=vert_to_chamb4V(v1)
         chamb2=vert_to_chamb4V(v2)
         chamb3=vert_to_chamb4V(v3)
         if ((chamb1.EQ.chamb2).AND.(chamb2.EQ.chamb3)) then
            face_to_chamb4V(f1)=chamb1 
         elseif ((chamb1.EQ.chamb2).AND.(chamb2.NE.chamb3)) then
            face_to_chamb4V(f1)=chamb1
         elseif ((chamb1.EQ.chamb3).AND.(chamb2.NE.chamb3)) then
            face_to_chamb4V(f1)=chamb1
         elseif ((chamb2.EQ.chamb3).AND.(chamb2.NE.chamb1)) then
            face_to_chamb4V(f1)=chamb2
         elseif ((chamb2.NE.chamb3).AND.(chamb2.NE.chamb1)) then
            face_to_chamb4V(f1)=MIN(chamb1,chamb2,chamb3)
         endif
      enddo
      !@cuf istat = cudaDeviceSynchronize !JDR TMP 
!#endif
      
!-----------VOLUMEchamb---------------
      vsi = vstart(1) ; vei = vend(1)
      esi = estart(1) ; eei = eend(1)
      fsi = fstart(1) ; fei = fend(1)
      ! call calculate_volume(Volume(1),vsi,vei,fsi,fei,  &
      !                  xyz0(:,vsi:vei),vert_of_face(:,fsi:fei))
      call calculate_volume_chamb(Volume_chamb0(1:8),vsi,vei,fsi,fei,  &
                       xyz(:,vsi:vei),vert_of_face(:,fsi:fei),face_to_chamb4V(fsi:fei))
      Volume_chamb = Volume_chamb0
      write(*,*) "Volume0 LV : ",Volume_chamb0(1)
      write(*,*) "Volume0 LA : ",Volume_chamb0(2)
      write(*,*) "Volume0 RV : ",Volume_chamb0(3)
      write(*,*) "Volume0 RA : ",Volume_chamb0(4)
      write(*,*) "Volume0 AO : ",Volume_chamb0(5)
      write(*,*) "Volume0 VP : ",Volume_chamb0(6)
      write(*,*) "Volume0 AP : ",Volume_chamb0(7)
      write(*,*) "Volume0 VC : ",Volume_chamb0(8)
      
      call calculate_area(Surface(1),vsi,vei,fsi,fei,         &
            xyz0(:,vsi:vei),vert_of_face(:,fsi:fei),    &
            sur0(fsi:fei))
      call calculate_area_chamb(Surface_chamb0(1:8),vsi,vei,fsi,fei,  &
                       xyz0(:,vsi:vei),vert_of_face(:,fsi:fei),face_to_chamb4V(fsi:fei))
      Surface_chamb = Surface_chamb0
      write(*,*) "Surface0 LV : ",Surface_chamb0(1)
      write(*,*) "Surface0 LA : ",Surface_chamb0(2)
      write(*,*) "Surface0 RV : ",Surface_chamb0(3)
      write(*,*) "Surface0 RA : ",Surface_chamb0(4)
      write(*,*) "Surface0 AO : ",Surface_chamb0(5)
      write(*,*) "Surface0 VP : ",Surface_chamb0(6)
      write(*,*) "Surface0 AP : ",Surface_chamb0(7)
      write(*,*) "Surface0 VC : ",Surface_chamb0(8)
      
!------------3D-------------------
#ifndef GEO_ANSA
      !edges
      !$cuf kernel do (1) 
      do e1=estart_3d(1),eend_3d(1)
         v1 = vert_of_edge_3d(1,e1)
         v2 = vert_of_edge_3d(2,e1)
         chamb1=vert_to_chamb_3d(v1)
         chamb2=vert_to_chamb_3d(v2)
         if (chamb1.EQ.chamb2) then
            edge_to_chamb_3d(e1)=chamb1 !setto ventricolo
         else
            edge_to_chamb_3d(e1)=MIN(chamb1,chamb2)
         endif
      enddo
      !@cuf istat = cudaDeviceSynchronize !JDR TMP
      
      !face
      !$cuf kernel do (1) 
      do f1=fstart_3d(1),fend_3d(1)
         v1 = vert_of_face_3d(1,f1)
         v2 = vert_of_face_3d(2,f1)
         v3 = vert_of_face_3d(3,f1)
         chamb1=vert_to_chamb_3d(v1)
         chamb2=vert_to_chamb_3d(v2)
         chamb3=vert_to_chamb_3d(v3)
         if ((chamb1.EQ.chamb2).AND.(chamb2.EQ.chamb3)) then
            face_to_chamb_3d(f1)=chamb1 
         elseif ((chamb1.EQ.chamb2).AND.(chamb2.NE.chamb3)) then
            face_to_chamb_3d(f1)=chamb1
         elseif ((chamb1.EQ.chamb3).AND.(chamb2.NE.chamb3)) then
            face_to_chamb_3d(f1)=chamb1
         elseif ((chamb2.EQ.chamb3).AND.(chamb2.NE.chamb1)) then
            face_to_chamb_3d(f1)=chamb2
         elseif ((chamb2.NE.chamb3).AND.(chamb2.NE.chamb1)) then
            face_to_chamb_3d(f1)=MIN(chamb1,chamb2,chamb3)
         endif
      enddo
      !@cuf istat = cudaDeviceSynchronize !JDR TMP
      
      !cells
      !$cuf kernel do (1) 
      do c1=cstart_3d(1),cend_3d(1)
         v1 = vert_of_cell_3d(1,c1)
         v2 = vert_of_cell_3d(2,c1)
         v3 = vert_of_cell_3d(3,c1)
         v4 = vert_of_cell_3d(4,c1)
         chamb1=vert_to_chamb_3d(v1)
         chamb2=vert_to_chamb_3d(v2)
         chamb3=vert_to_chamb_3d(v3)
         chamb4=vert_to_chamb_3d(v4)
         if ((chamb1.EQ.chamb2).AND.(chamb2.EQ.chamb3).AND.(chamb3.EQ.chamb4)) then
            cell_to_chamb_3d(c1)=chamb1 
         elseif ((chamb1.EQ.chamb2).AND.(chamb2.EQ.chamb3).AND.(chamb3.NE.chamb4)) then
            cell_to_chamb_3d(c1)=chamb1
         elseif ((chamb1.EQ.chamb2).AND.(chamb2.EQ.chamb4).AND.(chamb3.NE.chamb2)) then
            cell_to_chamb_3d(c1)=chamb1
         elseif ((chamb1.EQ.chamb3).AND.(chamb3.EQ.chamb4).AND.(chamb1.NE.chamb2)) then
            cell_to_chamb_3d(c1)=chamb1
         elseif ((chamb2.EQ.chamb3).AND.(chamb3.EQ.chamb4).AND.(chamb1.NE.chamb2)) then
            cell_to_chamb_3d(c1)=chamb2
         else
            cell_to_chamb_3d(c1)=MIN(chamb1,chamb2,chamb3,chamb4)
         endif
      enddo
      !@cuf istat = cudaDeviceSynchronize !JDR TMP 
#endif

      !vertLV vertL
      !$cuf kernel do (1)
      do v1=vstart_3d(1),vend_3d(1)
         chamb1=vert_to_chamb_3d(v1)
         if ((checkLV.EQ.0).AND.(chamb1.EQ.1)) then
            checkLV=1
            vertLV=v1
         endif
         if ((checkLA.EQ.0).AND.(chamb1.EQ.2)) then
            checkLA=1
            vertLA=v1
         endif         
      enddo
      !@cuf istat = cudaDeviceSynchronize !JDR TMP 
      write(*,*) "vertLV", vertLV
      write(*,*) "vertLA", vertLA


!       countSTE = 0
!       countSTEvett = 0
!         do i=cstart_3d(1),cend_3d(1)
!            f1=face_of_cell_3d(1,i)
!            f2=face_of_cell_3d(2,i)
!            f3=face_of_cell_3d(3,i)
!            f4=face_of_cell_3d(4,i)

!            chamb=cell_to_chamb_3d(i)

!            s1=sur_3d(f1)
!            s2=sur_3d(f2)
!            s3=sur_3d(f3)
!            s4=sur_3d(f4)

!            supmin=min(s1,s2,s3,s4)
!            supmax=max(s1,s2,s3,s4)

!            ratios=supmin/supmax
!            if ((supmin.LT.1.5E-3)) then                                 
! !              LabelStenosi(i)=1
!               countSTE=countSTE+1
!               countSTEvett(chamb)=countSTEvett(chamb)+1
!            endif
!         enddo
!         write(*,*) "countSTE = ", countSTE, nctot_3d
!         write(*,*) "countSTEchamb = ", countSTEvett

!################ PER DETTAGLI ECG ##################
!escludi zone di transizione tra una camera e l'altra 
!perche creano gradienti spuri
!BOU new chamb array ecg
      !vert
      !$cuf kernel do (1)
      do v1=vstart_3d(1),vend_3d(1)
         vert_to_chamb_3dBou(v1)=vert_to_chamb_3d(v1)
      enddo
      !@cuf istat = cudaDeviceSynchronize !JDR TMP

do loop=1,1
         
      !$cuf kernel do (1) 
      do v1=vstart_3d(1),vend_3d(1)
         vert_to_chamb_3dBouOld(v1)=vert_to_chamb_3dBou(v1)
      enddo
      !@cuf istat = cudaDeviceSynchronize !JDR TMP

      !$cuf kernel do (1) 
      do v1=vstart_3d(1),vend_3d(1)
         chamb=vert_to_chamb_3dBouOld(v1)
         if ( (chamb.EQ.2).OR.(chamb.EQ.4) ) then

            do jj=1,n_edge_of_vert_3d(v1)
                v2 = vert_of_vert_3d(jj,v1)
                chambcheck=vert_to_chamb_3dBouOld(v2)
                if (chamb.NE.chambcheck) then
                   vert_to_chamb_3dBou(v1)=chambcheck
                endif
            enddo

         endif ! if atri
      enddo
      !@cuf istat = cudaDeviceSynchronize !JDR TMP 
 enddo


 !cellsBou
 !$cuf kernel do (1) 
 do c1=cstart_3d(1),cend_3d(1)
    v1 = vert_of_cell_3d(1,c1)
    v2 = vert_of_cell_3d(2,c1)
    v3 = vert_of_cell_3d(3,c1)
    v4 = vert_of_cell_3d(4,c1)
    chamb1=vert_to_chamb_3dBou(v1)
    chamb2=vert_to_chamb_3dBou(v2)
    chamb3=vert_to_chamb_3dBou(v3)
    chamb4=vert_to_chamb_3dBou(v4)
    if ((chamb1.EQ.chamb2).AND.(chamb2.EQ.chamb3).AND.(chamb3.EQ.chamb4)) then
       cell_to_chamb_3dBou(c1)=chamb1 
    elseif ((chamb1.EQ.chamb2).AND.(chamb2.EQ.chamb3).AND.(chamb3.NE.chamb4)) then
       cell_to_chamb_3dBou(c1)=chamb1
    elseif ((chamb1.EQ.chamb2).AND.(chamb2.EQ.chamb4).AND.(chamb3.NE.chamb2)) then
       cell_to_chamb_3dBou(c1)=chamb1
    elseif ((chamb1.EQ.chamb3).AND.(chamb3.EQ.chamb4).AND.(chamb1.NE.chamb2)) then
       cell_to_chamb_3dBou(c1)=chamb1
    elseif ((chamb2.EQ.chamb3).AND.(chamb3.EQ.chamb4).AND.(chamb1.NE.chamb2)) then
       cell_to_chamb_3dBou(c1)=chamb2
    else
       cell_to_chamb_3dBou(c1)=MIN(chamb1,chamb2,chamb3,chamb4)
    endif
 enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP 

      ! !unify vert_to_chamb_3d at interfaces (atria-septum-vessels/ventricles-septum-vessels)
      ! do c1=cstart_3d(1),cend_3d(1)

      !    if (cell_to_chamb_3d(c1).eq.3) then !RV
      !       v1 = vert_of_cell_3d(1,c1)
      !       v2 = vert_of_cell_3d(2,c1)
      !       v3 = vert_of_cell_3d(3,c1)
      !       v4 = vert_of_cell_3d(4,c1)

      !       vert_to_chamb_3d(v1)=3
      !       vert_to_chamb_3d(v2)=3
      !       vert_to_chamb_3d(v3)=3
      !       vert_to_chamb_3d(v4)=3
      !    elseif (cell_to_chamb_3d(c1).eq.1) then !LV
      !       v1 = vert_of_cell_3d(1,c1)
      !       v2 = vert_of_cell_3d(2,c1)
      !       v3 = vert_of_cell_3d(3,c1)
      !       v4 = vert_of_cell_3d(4,c1)

      !       vert_to_chamb_3d(v1)=1
      !       vert_to_chamb_3d(v2)=1
      !       vert_to_chamb_3d(v3)=1
      !       vert_to_chamb_3d(v4)=1
      !    elseif (cell_to_chamb_3d(c1).eq.4) then !RA
      !       v1 = vert_of_cell_3d(1,c1)
      !       v2 = vert_of_cell_3d(2,c1)
      !       v3 = vert_of_cell_3d(3,c1)
      !       v4 = vert_of_cell_3d(4,c1)

      !       vert_to_chamb_3d(v1)=4
      !       vert_to_chamb_3d(v2)=4
      !       vert_to_chamb_3d(v3)=4
      !       vert_to_chamb_3d(v4)=4
      !    elseif (cell_to_chamb_3d(c1).eq.2) then !LA
      !       v1 = vert_of_cell_3d(1,c1)
      !       v2 = vert_of_cell_3d(2,c1)
      !       v3 = vert_of_cell_3d(3,c1)
      !       v4 = vert_of_cell_3d(4,c1)

      !       vert_to_chamb_3d(v1)=2
      !       vert_to_chamb_3d(v2)=2
      !       vert_to_chamb_3d(v3)=2
      !       vert_to_chamb_3d(v4)=2
      !    elseif (cell_to_chamb_3d(c1).eq.7) then !PA
      !       v1 = vert_of_cell_3d(1,c1)
      !       v2 = vert_of_cell_3d(2,c1)
      !       v3 = vert_of_cell_3d(3,c1)
      !       v4 = vert_of_cell_3d(4,c1)

      !       vert_to_chamb_3d(v1)=7
      !       vert_to_chamb_3d(v2)=7
      !       vert_to_chamb_3d(v3)=7
      !       vert_to_chamb_3d(v4)=7
      !    elseif (cell_to_chamb_3d(c1).eq.5) then !AO
      !       v1 = vert_of_cell_3d(1,c1)
      !       v2 = vert_of_cell_3d(2,c1)
      !       v3 = vert_of_cell_3d(3,c1)
      !       v4 = vert_of_cell_3d(4,c1)

      !       vert_to_chamb_3d(v1)=5
      !       vert_to_chamb_3d(v2)=5
      !       vert_to_chamb_3d(v3)=5
      !       vert_to_chamb_3d(v4)=5
      !    endif
        
      ! enddo



      !check 

      return
      end subroutine body1_to_chamb
!...................................................................
!===================================================================                         
                                                                                  
      subroutine read_geoSingleBody
      use constants
      use param
      use mpih
      use mpi_param, only: kstart, kend
      use mls_param
      !@cuf use cudafor
      implicit none
      !@cuf integer :: istat
      character*50 filename,strucfilename
      character*100 ipfi,ipfip

      integer :: i,j,k,v1,v2,v3,inp,e1,e2,e3,count
      real(DP) :: totmass,usVoldom
      integer,dimension(:),allocatable :: nv,ne,nf
      real(DP) :: xmax,xmin,ymax,ymin,zmax,zmin

      integer :: v1d,v2d
      integer :: e1d,e2d,e3d
      integer :: vsi,vei,esi,eei,fsi,fei,csi,cei

      do inp=1,Nparticle

      open(109,file='meshes/'//trim(geofile(inp))//'.gts')
      read(109,*)nvi(inp),nei(inp),nfi(inp)
      !maxnv=max(maxnv,nv(inp))
      !maxne=max(maxne,ne(inp))
      !maxnf=max(maxnf,nf(inp))
      close(109)
      end do
!     Allocations using vertices,edges and faces as parameters     

      call allocate_trigeo
!     Read in the vertices, edges and faces from the file
!     Also make connections between the three
      n_edge_of_vert=0
      face_of_edge=0

      do inp=1,Nparticle

      vsi = vstart(inp) ; vei = vend(inp)
      esi = estart(inp) ; eei = eend(inp)
      fsi = fstart(inp) ; fei = fend(inp)

      ! Set particle labels
      vert_to_part(vsi:vei) = inp
      edge_to_part(esi:eei) = inp
      face_to_part(fsi:fei) = inp

        open(109,file='meshes/'//trim(geofile(inp))//'.gts')
        read(109,*)nvi(inp),nei(inp),nfi(inp)
        do i=vsi,vei
          read(109,*)xyz0(1,i),xyz0(2,i),xyz0(3,i)
            xyz0(1,i) = xyz0(1,i)
            xyz0(2,i) = xyz0(2,i)
            xyz0(3,i) = xyz0(3,i)
        end do
!     position correction to bring in 0-6
 
      xyzv0 = 0.0
      xyza0 = 0.0
      fxyz  = 0.0      

      xyz = xyz0
      xyzv = xyzv0
      xyza = xyza0

        

        do i=esi,eei
          read(109,*)v1,v2

          v1d = v1 + vsi-1
          v2d = v2 + vsi-1
          vert_of_edge(1,i)=v1d
          vert_of_edge(2,i)=v2d 

          n_edge_of_vert(v1d)=n_edge_of_vert(v1d)+1
          n_edge_of_vert(v2d)=n_edge_of_vert(v2d)+1

          vert_of_vert(n_edge_of_vert(v1d),v1d)=v2d           
          vert_of_vert(n_edge_of_vert(v2d),v2d)=v1d           

          edge_of_vert(n_edge_of_vert(v1d),v1d)=i
          edge_of_vert(n_edge_of_vert(v2d),v2d)=i
        enddo

        do i=fsi,fei
          read(109,*) e1d, e2d, e3d
                    
          edge_of_face(1,i) = e1d + esi - 1
          edge_of_face(2,i) = e2d + esi - 1
          edge_of_face(3,i) = e3d + esi - 1
           
        end do
!        !$cuf kernel do (1)
        do i=fsi,fei
           e1=edge_of_face(1,i)
           e2=edge_of_face(2,i)

           if (vert_of_edge(2,e1).eq.vert_of_edge(1,e2)) then
              v1=vert_of_edge(1,e1)
              v2=vert_of_edge(2,e1)
              v3=vert_of_edge(2,e2)
           elseif(vert_of_edge(2,e1).eq.vert_of_edge(2,e2)) then
              v1=vert_of_edge(1,e1)
              v2=vert_of_edge(2,e1)
              v3=vert_of_edge(1,e2)
           elseif(vert_of_edge(1,e1).eq.vert_of_edge(1,e2)) then
              v1=vert_of_edge(2,e1)
              v2=vert_of_edge(1,e1)
              v3=vert_of_edge(2,e2)
           else 
              v1=vert_of_edge(2,e1)
              v2=vert_of_edge(1,e1)
              v3=vert_of_edge(1,e2)
           endif 

              vert_of_face(1,i)=v1
              vert_of_face(2,i)=v2
              vert_of_face(3,i)=v3
           enddo
!           !@cuf istat = cudaDeviceSynchronize !JDR TMP

          close(109)
!   Completed reading the gts file
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Check - vertex cannot be connected to itself
          !$cuf kernel do (1)
          do i=vsi,vei
              do j=1,n_edge_of_vert(i)
                 if (vert_of_vert(j,i).eq.i)  then
                       write(*,*)'Error ',vert_of_vert(j,i),i;stop
                 endif
              enddo
           enddo
           !@cuf istat = cudaDeviceSynchronize !JDR TMP
           
           !Check - for faces and edges
           !$cuf kernel do (1)
           do i=fsi,fei
              e1=edge_of_face(1,i)
              e2=edge_of_face(2,i)
              e3=edge_of_face(3,i)
              if (face_of_edge(1,e1).eq.0) then
                 face_of_edge(1,e1)=i
              elseif (face_of_edge(2,e1).eq.0) then
                 face_of_edge(2,e1)=i
              else
                 write(*,*)'Edge error1', i,e1,e2,e3;stop
              endif
              if (face_of_edge(1,e2).eq.0) then
                 face_of_edge(1,e2)=i
              elseif (face_of_edge(2,e2).eq.0) then
                 face_of_edge(2,e2)=i
              else
                 write(*,*)'Edge error2';stop
              endif
              if (face_of_edge(1,e3).eq.0) then
                 face_of_edge(1,e3)=i
              elseif (face_of_edge(2,e3).eq.0) then
                 face_of_edge(2,e3)=i
              else
                 write(*,*)'Edge error3';stop
              endif
           enddo 
           !@cuf istat = cudaDeviceSynchronize !JDR TMP
           
           !Check
           count=0
           do i=esi,eei
              if (face_of_edge(1,i).eq.face_of_edge(2,i)) then
                 write(*,*)'Error on edges ';stop
              endif
              if (face_of_edge(1,i).eq.0.or.face_of_edge(2,i).eq.0) then
                 count=count+1
              endif
           enddo
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      call convert_geo(vsi,vei,esi,eei,        &
            fsi,fei,xyz(:,vsi:vei),            &
            xyzv(:,vsi:vei),xyza(:,vsi:vei),   &
            vert_of_face(:,fsi:fei),           &
            tri_ver(:,fsi:fei),                &
            tri_vel(:,fsi:fei),                &
            tri_bar(:,fsi:fei),                &
            vel_tri(:,fsi:fei),                &
            acc_tri(:,fsi:fei))

        call calculate_normal(vsi,vei,fsi,             &
                             fei,xyz0(:,vsi:vei),      &
                             vert_of_face(:,fsi:fei),  &
                            tri_nor(:,fsi:fei))

        call calculate_angle(theta0(esi:eei),esi,eei,fsi,fei,      &
                       face_of_edge(:,esi:eei),tri_nor(:,fsi:fei))

        call calculate_anglealp(aalpha0(:,fsi:fei),vsi,vei,fsi,fei,  &
                       xyz0(:,vsi:vei),                            &
                       vert_of_face(:,fsi:fei))

        call calculate_distance(dist0(esi:eei),vsi,vei,esi,eei,    &
                       xyz0(:,vsi:vei),vert_of_edge(:,esi:eei))

        call calculate_area(Surface0(inp),vsi,vei,fsi,fei,         &
                       xyz0(:,vsi:vei),vert_of_face(:,fsi:fei),    &
                       sur0(fsi:fei))
        sur=sur0
        write(*,*) "HEEEEEEEERE"
        call calculate_volume(Volume0(inp),vsi,vei,fsi,fei,  &
                       xyz0(:,vsi:vei),vert_of_face(:,fsi:fei))

        call find_quartet(v1234(:,esi:eei),vsi,vei,esi,eei,fsi,fei,       &
                      xyz0(:,vsi:vei),                                    &
                      face_of_edge(:,esi:eei),edge_of_face(:,fsi:fei),    &
                      vert_of_edge(:,esi:eei),tri_nor(:,fsi:fei))

        if(myid.eq.0)write(*,*)'Finished reading geo. no. ',inp


      totmass = 1.0
      !uspm = 1.0

      end do
    
!     save for aortic pre-tension 
      dist00 = dist0 


#ifndef SOLOUNO
      !CARICA LABEL LEAFLETS
      inp=2 !Valvola aortica
      vsi = vstart(inp) ; vei = vend(inp)
      esi = estart(inp) ; eei = eend(inp)
      fsi = fstart(inp) ; fei = fend(inp)
      allocate(LabelVaortic(vsi:vei))
      
      open(109,file='meshes/'//trim(geofile(inp))//'_leaflets.txt')
!      open(109,file='meshes/AorticValve9V1_leaflets.txt')
      do i=vsi,vei
         read(109,*) LabelVaortic(i)
      end do
      close(109)

      inp=3 !Valvola mitrale
      vsi = vstart(inp) ; vei = vend(inp)
      esi = estart(inp) ; eei = eend(inp)
      fsi = fstart(inp) ; fei = fend(inp)
      allocate(LabelVmitral(vsi:vei))
      
      open(109,file='meshes/'//trim(geofile(inp))//'_leaflets.txt')
!      open(109,file='meshes/MitralValve6V2_leaflets.txt')
      do i=vsi,vei
         read(109,*) LabelVmitral(i)
      end do
      close(109)

      inp=4 !Valvola polmonare
      vsi = vstart(inp) ; vei = vend(inp)
      esi = estart(inp) ; eei = eend(inp)
      fsi = fstart(inp) ; fei = fend(inp)
      allocate(LabelVpulmo(vsi:vei))
      
      open(109,file='meshes/'//trim(geofile(inp))//'_leaflets.txt')
      do i=vsi,vei
         read(109,*) LabelVpulmo(i)
      end do
      close(109)

      inp=5 !Valvola tricuspide
      vsi = vstart(inp) ; vei = vend(inp)
      esi = estart(inp) ; eei = eend(inp)
      fsi = fstart(inp) ; fei = fend(inp)
      allocate(LabelVtricu(vsi:vei))
      
      open(109,file='meshes/'//trim(geofile(inp))//'_leaflets.txt')
      do i=vsi,vei
         read(109,*) LabelVtricu(i)
      end do
      close(109)
#endif
      
        return
        end subroutine read_geoSingleBody
!===================================================================
        ! subroutine find_boundaries(bcs,nvtot,netot,face_of_edge, &
        !                vert_of_edge,vert_to_part,edge_to_part,xyz,count2,boundary1,LabelBound)
        subroutine find_boundaries(bcs,nvtot,netot,face_of_edge, &
                       vert_of_edge,vert_to_part,edge_to_part,xyz,count2,boundary1)
        use constants
        use mls_param

        implicit none
        integer :: bcs,nvtot,netot,count2,nbound
        integer :: i, j,v1, v2,inp,vsi,vei,qualeoff
        integer, dimension (*) :: vert_to_part
        integer, dimension (*) :: edge_to_part
        integer, dimension (2,*) :: face_of_edge, vert_of_edge
        integer, dimension (2,nvtot) :: boundary1
        ! integer, dimension (nvtot) :: LabelBound
        real(DP), dimension (3,*) :: xyz
        real(DP) pig, Rsec, zsecmin, off0, offmin,off
        real(DP),dimension(:,:),allocatable ::xyz_bound
        pig=acos(-1.0)

        ! LabelBound=0.0D0

        off0 = 100000.
        count2 = 0

#ifndef SOLOUNO
! SOSTITUIRE CON UN LOOP 2,3,4,5
!VALVOLA AORTICA
!carica 
      inp=2
      vsi=vstart(inp)
      vei=vend(inp)
      open(109,file='meshes/'//trim(geofile(inp))//'.txt')
      read(109,*)nbound
      allocate(xyz_bound(3,nbound))
      do i=1,nbound
         read(109,*) xyz_bound(1,i),xyz_bound(2,i),xyz_bound(3,i)
         xyz_bound(1,i) = xyz_bound(1,i)
         xyz_bound(2,i) = xyz_bound(2,i)
         xyz_bound(3,i) = xyz_bound(3,i)
      end do
      close(109)
!trova      
      do j=1,nbound
         offmin=off0
         
         do i=vsi,vei
              off = sqrt((xyz(1,i)-xyz_bound(1,j))**2+(xyz(2,i)-xyz_bound(2,j))**2+(xyz(3,i)-xyz_bound(3,j))**2)
              if (off.LT.offmin) then                 
                 qualeoff = i
                 offmin = off
              endif
           enddo

           count2 = count2 +1
           boundary1(1,count2) = qualeoff
           boundary1(2,count2) = 1
           ! LabelBound(qualeoff)=1.0D0
       enddo
       deallocate(xyz_bound)
       

!VALVOLA MITRALE
!carica 
      inp=3
      vsi=vstart(inp)
      vei=vend(inp)
      open(109,file='meshes/'//trim(geofile(inp))//'.txt')      
      read(109,*)nbound
      allocate(xyz_bound(3,nbound))
      do i=1,nbound
         read(109,*) xyz_bound(1,i),xyz_bound(2,i),xyz_bound(3,i)
         xyz_bound(1,i) = xyz_bound(1,i)
         xyz_bound(2,i) = xyz_bound(2,i)
         xyz_bound(3,i) = xyz_bound(3,i)
      end do
      close(109)
!trova      
      do j=1,nbound
         offmin=off0

         do i=vsi,vei
              off = sqrt((xyz(1,i)-xyz_bound(1,j))**2+(xyz(2,i)-xyz_bound(2,j))**2+(xyz(3,i)-xyz_bound(3,j))**2)
              if (off.LT.offmin) then                 
                 qualeoff = i
                 offmin = off
              endif
           enddo

           count2 = count2 +1
           boundary1(1,count2) = qualeoff
           boundary1(2,count2) = 1
           ! LabelBound(qualeoff)=1.0D0
       enddo
       deallocate(xyz_bound)


!VALVOLA POLMONARE
!carica 
      inp=4
      vsi=vstart(inp)
      vei=vend(inp)
      open(109,file='meshes/'//trim(geofile(inp))//'.txt')
      ! open(109,file='meshes/'//trim(geofile(inp))//'r.txt')      
      read(109,*)nbound
      allocate(xyz_bound(3,nbound))
      do i=1,nbound
         read(109,*) xyz_bound(1,i),xyz_bound(2,i),xyz_bound(3,i)
         xyz_bound(1,i) = xyz_bound(1,i)
         xyz_bound(2,i) = xyz_bound(2,i)
         xyz_bound(3,i) = xyz_bound(3,i)
      end do
      close(109)
!trova      
      do j=1,nbound
         offmin=off0

         do i=vsi,vei
              off = sqrt((xyz(1,i)-xyz_bound(1,j))**2+(xyz(2,i)-xyz_bound(2,j))**2+(xyz(3,i)-xyz_bound(3,j))**2)
              if (off.LT.offmin) then                 
                 qualeoff = i
                 offmin = off
              endif
           enddo

           count2 = count2 +1
           boundary1(1,count2) = qualeoff
           boundary1(2,count2) = 1
           ! LabelBound(qualeoff)=1.0D0
       enddo
       deallocate(xyz_bound)

!VALVOLA TRICUSPIDE
!carica 
      inp=5
      vsi=vstart(inp)
      vei=vend(inp)
      open(109,file='meshes/'//trim(geofile(inp))//'.txt')      
      read(109,*)nbound
      allocate(xyz_bound(3,nbound))
      do i=1,nbound
         read(109,*) xyz_bound(1,i),xyz_bound(2,i),xyz_bound(3,i)
         xyz_bound(1,i) = xyz_bound(1,i)
         xyz_bound(2,i) = xyz_bound(2,i)
         xyz_bound(3,i) = xyz_bound(3,i)
      end do
      close(109)
!trova      
      do j=1,nbound
         offmin=off0

         do i=vsi,vei
              off = sqrt((xyz(1,i)-xyz_bound(1,j))**2+(xyz(2,i)-xyz_bound(2,j))**2+(xyz(3,i)-xyz_bound(3,j))**2)
              if (off.LT.offmin) then                 
                 qualeoff = i
                 offmin = off
              endif
           enddo

           count2 = count2 +1
           boundary1(1,count2) = qualeoff
           boundary1(2,count2) = 1
           ! LabelBound(qualeoff)=1.0D0
       enddo
       deallocate(xyz_bound)
       

       
       !TAPPI
       !carica da 6 a 10
       vsi=vstart(6)
       vei=vend(10)
       do i=vsi,vei
          count2 = count2 +1
          boundary1(1,count2) = i
          boundary1(2,count2) = 1
          ! LabelBound(qualeoff)=1.0D0
       enddo

#endif
       
        end subroutine find_boundaries
!===================================================================  
!francesco
        subroutine find_offset(nvtot, xyz, vert_to_part,count2,boundary2,BCoffset)
        use constants
        use mls_param, only: vstart,vend  
        implicit none
        integer :: nvtot,count2,i,j,v1,v2,inp
        integer nvstart,nvend,qualeoff
        integer, dimension (*) :: vert_to_part
!        integer, dimension (2,*) :: face_of_edge, vert_of_edge
        integer, dimension (2,count2) :: boundary2
        real, dimension(3,count2):: BCoffset
        real, dimension(3):: ttt
        real(DP), dimension (3,*) :: xyz
        real(DP) pig,off0, offmin, off
        pig=acos(-1.0)


        off0 = 100000.

        do i=1,count2
           v1  = boundary2(1,i)
           inp = boundary2(2,i) 

           if (inp.GE.1) then
              nvstart = vstart(inp) !attenzione ora punta di nuovo al 2d
              nvend = vend(inp)      

              offmin = off0
              do j=nvstart,nvend             
                 off = sqrt((xyz(1,v1)-xyz(1,j))**2+(xyz(2,v1)-xyz(2,j))**2+(xyz(3,v1)-xyz(3,j))**2)
                 if (off.LT.offmin) then                 
                    qualeoff = j
                    offmin = off
                 endif
              enddo
!           write(*,*) i,qualeoff,vert_to_part(qualeoff)
              boundary2(2,i) = qualeoff           
              BCoffset(1:3,i) = xyz(1:3,v1)-xyz(1:3,qualeoff)              
           endif
        enddo

        end subroutine find_offset
!===================================================================
        subroutine find_boundaries_3d(nvtot, vert_to_part,vert_to_chamb,&
                                   xyz,count2,boundary1)
        use constants
        use mls_param, only: vstart_3d,vend_3d

        implicit none
        integer :: nvtot,count2
        integer :: i, v1, v2,inp,chamb,vsi,vei,nbound,j,qualeoff
        integer, dimension (*) :: vert_to_part
        integer, dimension (*) :: vert_to_chamb
        integer, dimension (2,nvtot) :: boundary1
        real(DP), dimension (3,*) :: xyz
        real(DP) dblocca, xP1, zP1, xP2, zP2, alpha, xsparti !ventricle
        real(DP) dbloccaA, x0, y0, z0, alp, xsec, ysec, zsec !aorta
        real(DP) pig, Rsec,xV,yV,zV,xV1,yV1,zV1,off0, offmin,off
        real(DP),dimension(:,:),allocatable ::xyz_bound
        pig=acos(-1.0)

        off0 = 100000.
        count2 = 0

        inp=1
        vsi=vstart_3d(inp)
        vei=vend_3d(inp)
        do i=vsi,vei
           chamb=vert_to_chamb(i)
           if (chamb.EQ.5) then !aorta
              if (xyz(3,i).LT.2.0) then
                 count2 = count2 +1
                 boundary1(1,count2) = i
                 boundary1(2,count2) = -1                 
              endif
           endif

           if (chamb.EQ.7) then !arteria polmonare
              if (xyz(1,i).GT.1.2) then
                 count2 = count2 +1
                 boundary1(1,count2) = i
                 boundary1(2,count2) = -1                 
              endif

              if (xyz(1,i).LT.-1.2) then
                 count2 = count2 +1
                 boundary1(1,count2) = i
                 boundary1(2,count2) = -1                 
              endif
           endif

        enddo

        


        
! !AORTA
! !carica 
!       inp=1
!       vsi=vstart_3d(inp)
!       vei=vend_3d(inp)
!       open(109,file='meshes/AortaAdd.txt')
!       read(109,*)nbound
!       allocate(xyz_bound(3,nbound))
!       do i=1,nbound
!          read(109,*) xyz_bound(1,i),xyz_bound(2,i),xyz_bound(3,i)
!       end do
!       close(109)
! !trova      
!       do j=1,nbound
!          offmin=off0

!          do i=vsi,vei
!               off = sqrt((xyz(1,i)-xyz_bound(1,j))**2+(xyz(2,i)-xyz_bound(2,j))**2+(xyz(3,i)-xyz_bound(3,j))**2)
!               if (off.LT.offmin) then                 
!                  qualeoff = i
!                  offmin = off
!               endif
!            enddo

!            count2 = count2 +1
!            boundary1(1,count2) = qualeoff
!            boundary1(2,count2) = -1
!        enddo
!        deallocate(xyz_bound)


       
        end subroutine find_boundaries_3d
!===================================================================  
!francesco
        subroutine find_offset_3d(nvtot, xyz, xyz_3d,vert_to_part,count2,boundary2,BCoffset)
        use constants
!        use mls_param, only: vstart,vend  
        implicit none
        integer :: nvtot,count2,i,j,v1,v2,inp
        integer nvstart,nvend,nvstart_3d,nvend_3d, qualeoff
        integer, dimension (*) :: vert_to_part
!        integer, dimension (2,*) :: face_of_edge, vert_of_edge
        integer, dimension (2,count2) :: boundary2
        real, dimension(3,count2):: BCoffset
        real, dimension(3):: ttt
        real(DP), dimension (3,*) :: xyz
        real(DP), dimension (3,*) :: xyz_3d
        real(DP) pig,off0, offmin, off
        pig=acos(-1.0)


        off0 = 100000.
        do i=1,count2
           v1  = boundary2(1,i)
           inp = boundary2(2,i) 
!            if (inp.GE.1) then
!               write(*,*) "AAAAA"
!               nvstart_3d = vstart_3d(inp) !attenzione ora punta al 3d 
!               nvend_3d = vend_3d(inp)      

!               offmin = off0
!               do j=nvstart_3d,nvend_3d             
! !                 off = sqrt((xyz(1,v1)-xyz_3d(1,j))**2+(xyz(2,v1)-xyz_3d(2,j))**2+(xyz(3,v1)-xyz_3d(3,j))**2)
!                  off = sqrt((xyz_3d(1,v1)-xyz_3d(1,j))**2+(xyz_3d(2,v1)-xyz_3d(2,j))**2+(xyz_3d(3,v1)-xyz_3d(3,j))**2)
!                  if (off.LT.offmin) then                 
!                     qualeoff = j
!                     offmin = off
!                  endif
!               enddo
! !           write(*,*) i,qualeoff,vert_to_part(qualeoff)
!               boundary2(2,i) = qualeoff           
!               BCoffset(1:3,i) = xyz_3d(1:3,v1)-xyz_3d(1:3,qualeoff)              
!            endif
        enddo

        end subroutine find_offset_3d
!===================================================================

!===================================================================
!francesco
        subroutine find_movingprobes(nvtot,xyz,vert_to_part)
        use constants
        use mls_param, only: vstart,vend 
        use probes 
        use param

        implicit none
        integer :: i, v1, v2, nvtot,dummN,nvstart,nvend,inp
        integer, dimension (*) :: vert_to_part
        real, dimension (3,nvtot) :: xyz
        real x0, y0, z0, alp, xsec, ysec, zsec, zsecAA !aorta
        real pig, dist1,dist2
        pig=acos(-1.0)


!# BLOCCA VERTICI SELEZIONATI QUI
        !utili 
        x0 = -0.437
        y0 =  0.
        z0 = +0.251
        alp = 41.*pig/180.

        !ventricolo            
        inp = 1
        nvstart = vstart(inp) 
        nvend = vend(inp)      
        dist1 = -100
        dist2 =  100
        do i=nvstart,nvend
        if (xyz(3,i).GT.1.0.AND.xyz(3,i).LT.1.5) then                                                                                                                                                         
            if (xyz(1,i).GT.dist1) then
               dist1 = xyz(1,i)
               nodesMovingProbes(1,1) = i
            endif
            if (xyz(1,i).LT.dist2) then
               dist2 = xyz(1,i)
               nodesMovingProbes(1,2) = i
            endif
         endif
        enddo

        !trova zmin aletta aortica
        inp = 3
        nvstart = vstart(inp) 
        inp = 5
        nvend = vend(inp)      
        zsecAA =  100
        do i=nvstart,nvend
            xsec = (xyz(1,i)-x0)*cos(alp)-(xyz(3,i)-z0)*sin(alp)
            zsec = (xyz(1,i)-x0)*sin(alp)+(xyz(3,i)-z0)*cos(alp)
            ysec =  xyz(2,i)
            if (zsec.LT.zsecAA) then 
               zsecAA = zsec
            endif
        enddo

        !aorta
        inp = 2
        nvstart = vstart(inp) 
        nvend = vend(inp)      
        dist1 = -100
        dist2 =  100
        do i=nvstart,nvend
           xsec = (xyz(1,i)-x0)*cos(alp)-(xyz(3,i)-z0)*sin(alp)
           zsec = (xyz(1,i)-x0)*sin(alp)+(xyz(3,i)-z0)*cos(alp)
           if ((xyz(2,i).GT.-0.5).AND.(zsec.LT.(zsecAA-0.15)).AND.(zsec.GT.(zsecAA-0.18))) then                                                                                                                               
             if (xsec.GT.dist1) then
               dist1 = xsec
               nodesMovingProbes(2,1) = i
              endif
             if (xsec.LT.dist2) then
               dist2 = xsec
               nodesMovingProbes(2,2) = i
             endif
            endif
         enddo


        !atrio            
        inp = 1
        nvstart = vstart(inp) 
        nvend = vend(inp)      
        dist1 = -100
        dist2 =  100
        do i=nvstart,nvend
        if (xyz(3,i).LT.-0.3.AND.xyz(3,i).GT.-0.5) then                                      
            if (xyz(1,i).GT.dist1) then
               dist1 = xyz(1,i)
               nodesMovingProbes(3,1) = i
            endif
            if (xyz(1,i).LT.dist2) then
               dist2 = xyz(1,i)
               nodesMovingProbes(3,2) = i
            endif
         endif
        enddo

        write(*,*) "MovingProbe Ventr1", xyz(1:3,nodesMovingProbes(1,1)) 
        write(*,*) "MovingProbe Ventr2", xyz(1:3,nodesMovingProbes(1,2)) 
        write(*,*) "MovingProbe Aorta1", xyz(1:3,nodesMovingProbes(2,1)) 
        write(*,*) "MovingProbe Aorta2", xyz(1:3,nodesMovingProbes(2,2)) 
        write(*,*) "MovingProbe Atrio1", xyz(1:3,nodesMovingProbes(3,1)) 
        write(*,*) "MovingProbe Atrio2", xyz(1:3,nodesMovingProbes(3,2)) 
        end subroutine find_movingprobes
!===================================================================

        subroutine read_PInterp_dim(np,filename)

        implicit none
        character*50 :: filename
        integer :: np

        open(11,file=filename)
        read(11,*)np
        close(11)

        write(*,*)'Npoints for Interpolation: ',np

        return
        end subroutine read_PInterp_dim
!...................................................................  
!===================================================================  
        subroutine read_PInterp(xyzInterp,np,filename)
        use constants

        implicit none
        character*50 :: filename
        integer :: np, i
        real(DP), dimension (np,3) :: xyzInterp

        open(11,file=filename)
        read(11,*)
        do i=1,np
           read(11,*)xyzInterp(i,1),xyzInterp(i,2),xyzInterp(i,3)
        enddo
        close(11)

        return
        end subroutine read_PInterp
!................................................................... 
!===================================================================
        subroutine read_WQuad_dim(np,npIn,npOut,filename)

        implicit none
        character*50 :: filename
        integer :: np, npIn, npOut

        open(11,file=filename)
        read(11,*)np, npIn, npOut
        close(11)
        
        if (np.ne.npIn+npOut) then
        write(*,*) 'Npoints Inlet + Npoints Outlet not equal to NpointsQuad'
        stop
        endif

        write(*,*)'Npoints for Inlet  quadrature:',npIn
        write(*,*)'Npoints for Outlet quadrature:',npOut

        return
        end subroutine read_WQuad_dim
!...................................................................                                                                                                                                                                         
!===================================================================                                                                                                                                                                         
        subroutine read_WQuad(WeiQuad,np,filename)
        use constants

        implicit none
        character*50 :: filename
        integer :: np, i
        real(DP), dimension (np,3) :: WeiQuad

        open(11,file=filename)
        read(11,*)
        do i=1,np
           read(11,*)WeiQuad(i,1),WeiQuad(i,2),WeiQuad(i,3)
        enddo
        close(11)

        return
        end subroutine read_WQuad
!...................................................................         
!................................................................... 
!===================================================================
        subroutine hyperelastic_properties
        use param
        use local_arrays
        use mls_param

        implicit none
        character*50 :: filename
        integer :: i,inp,chamb
        real(DP) :: A1, A2,xV,yV,zV,xV1,yV1,zV1,rV1,cosangFE

!FV new routine to be checked
        
! 2D !TO DO: add angle dependence for leaflets
      do inp = 1,Nparticle !TO DO: set here rkc values of each body
         if (inp.EQ.2) then  
            A1 = 49.558
            A2 = 5.2871
         elseif (inp.EQ.3) then  
            A1 = 49.558
            A2 = 5.2871
         elseif (inp.EQ.4) then  
            A1 = 49.558
            A2 = 5.2871
         elseif (inp.EQ.5) then  
            A1 = 49.558
            A2 = 5.2871
         endif
           
         rkc(inp) = rke(inp)/A1
 ! (Note, rck=rke/A1 is set here only to have a relation between the linear elastic modulus and the nonlinaer one at low strain,
 ! but not nessary)             
      enddo
      
      do i = 1,netot  !TO DO: set here the values of A1,A2 for each body
         inp = edge_to_part(i)
         if (inp.EQ.2) then  
            A1 = 49.558
            A2 = 5.2871
         elseif (inp.EQ.3) then  
            A1 = 49.558
            A2 = 5.2871
         elseif (inp.EQ.4) then  
            A1 = 49.558
            A2 = 5.2871
         elseif (inp.EQ.5) then  
            A1 = 49.558
            A2 = 5.2871
         endif         
         AFung(i) = A1
      enddo



! 3D
!TO DO: add here rkc for each heart chamber, put correct values
      A1 = 49.558
      rkc_3d(1) = rke_3d(1)/A1    !LV
      rkc_3d(2) = rke_3d(1)/A1    !LA
      rkc_3d(3) = rke_3d(1)/A1    !RV
      rkc_3d(4) = rke_3d(1)/A1    !RA
      rkc_3d(5) = rke_3d(1)/A1    !AO
      rkc_3d(6) = rke_3d(1)/A1    !PV
      rkc_3d(7) = rke_3d(1)/A1    !PA
      rkc_3d(8) = rke_3d(1)/A1    !VC
      rkc_3d(9) = rke_3d(1)/A1    !SP
      

     
      do i = 1,netot_3d 
         chamb = edge_to_chamb_3d(i)
          if (chamb.EQ.1) then !TO DO: set here the values for each chamber
             A1 = 49.558
             A2 = 5.2871
          elseif (chamb.EQ.2) then
             A1 = 49.558
             A2 = 5.2871
          elseif (chamb.EQ.3) then
             A1 = 49.558
             A2 = 5.2871
          elseif (chamb.EQ.4) then
             A1 = 49.558
             A2 = 5.2871
          elseif (chamb.EQ.5) then
             A1 = 49.558
             A2 = 5.2871
          elseif (chamb.EQ.6) then
             A1 = 49.558
             A2 = 5.2871
          elseif (chamb.EQ.7) then
             A1 = 49.558
             A2 = 5.2871
          elseif (chamb.EQ.8) then
             A1 = 49.558
             A2 = 5.2871
          elseif (chamb.EQ.9) then
             A1 = 49.558
             A2 = 5.2871
          endif

          cosAngFE   = edgefiber_cosangle_3d(i)
          AFung_3d(i) = sqrt(A1**2*cosAngFE**2+A2**2*(1-cosAngFE**2))
      enddo


      StiffV(:)=0.
      StiffV_3d(:)=0.

        return
        end subroutine hyperelastic_properties
!...................................................................                                                                                                                                                                         
!===================================================================                                                          
!................................................................... 
        subroutine find_mass
        use param
        use local_arrays
        use mls_param

        implicit none
        character*50 :: filename
        integer :: i,inp,v1,v2,v3,v4
        real(DP) :: facemass,totmass,cellmass

!2D
        mass_of_vert(:) = 0.0
        do inp=1,Nparticle
           totmass = Surface0(inp)*thck(inp)*rhos(inp)
           rpm(inp) = totmass/real(nvi(inp))

           do i=fstart(inp),fend(inp)
              v1 = vert_of_face(1,i)
              v2 = vert_of_face(2,i)
              v3 = vert_of_face(3,i)

              facemass = sur0(i)*thck(inp)*rhos(inp)
              
              mass_of_vert(v1) = mass_of_vert(v1) + facemass/3.D0
              mass_of_vert(v2) = mass_of_vert(v2) + facemass/3.D0
              mass_of_vert(v3) = mass_of_vert(v3) + facemass/3.D0
           end do

        end do

!3D
        mass_of_vert_3d(:) = 0.0
        do inp=1,Nparticle_3d
           do i=cstart_3d(inp),cend_3d(inp)
              v1 = vert_of_cell_3d(1,i)
              v2 = vert_of_cell_3d(2,i)
              v3 = vert_of_cell_3d(3,i)
              v4 = vert_of_cell_3d(4,i)

              if (LabelStenosi(i).EQ.1) then
                 cellmass = vol0_3d(i)*rhos_3d(inp)*2.0D0
              else
                 cellmass = vol0_3d(i)*rhos_3d(inp)
              endif

              mass_of_vert_3d(v1) = mass_of_vert_3d(v1) + cellmass/4.D0
              mass_of_vert_3d(v2) = mass_of_vert_3d(v2) + cellmass/4.D0
              mass_of_vert_3d(v3) = mass_of_vert_3d(v3) + cellmass/4.D0
              mass_of_vert_3d(v4) = mass_of_vert_3d(v4) + cellmass/4.D0
           end do

        end do



        return
        end subroutine find_mass
!...................................................................                                                          
!===================================================================                                                          
      subroutine find_indices_2dstr
      use constants
      use param
      use mpih
      use mls_param
      implicit none
      
      integer :: wetold, wetnew, inp, count

      if (Nparticle.GT.1) then
!check wet_3d surfaces are at the beginning of gtsfiles
!(they could be also at the end but avoided for now)
      count = 0 
      wetold = wet_3d(1)
      do inp=2,Nparticle
         wetnew = wet_3d(inp)
        if ((wet_3d(1).EQ.0).AND.(wetnew.NE.0)) then
           write(*,*) "Error wet_3d"
           stop
        endif
        if ((wetnew.EQ.0).AND.(wetold.NE.0)) count=count+1
        if ((wetold.EQ.0).AND.(wetnew.NE.0)) count=count+1
         wetold = wetnew
      enddo
      if (count.GT.1) then 
         write(*,*) "Error wet_3d"
         stop
      endif
      if ((count.GE.1).AND.(wet_3d(1).EQ.0)) then
         write(*,*) "Error wet_3d"
         stop
      endif


!find indices wet part and structural      
      if (wet_3d(1).GT.0) then
         inpstart_2dwet = 1
         nvstart_2dwet = vstart(inpstart_2dwet)
         nestart_2dwet = estart(inpstart_2dwet)
         nfstart_2dwet = fstart(inpstart_2dwet)

         count = 0
         wetold = wet_3d(1)        
         do inp=2,Nparticle
            wetnew = wet_3d(inp)
            if ((wetnew.EQ.0).AND.(wetold.NE.0)) then
               count=count+1

               inpend_2dwet = inp-1
               nvend_2dwet = vend(inpend_2dwet)
               neend_2dwet = eend(inpend_2dwet)
               nfend_2dwet = fend(inpend_2dwet)

               inpstart_2dstr = inp
            endif
            wetold = wetnew
         enddo
         if ( count.NE.1 ) then 
            write(*,*) "Error wet_3d"
            stop
         endif
      else
         inpstart_2dwet = 0
         inpend_2dwet = 0
         nvstart_2dwet = 0
         nvend_2dwet = 0
         nestart_2dwet = 0
         neend_2dwet = 0
         nfstart_2dwet = 0
         nfend_2dwet = 0

         inpstart_2dstr = 1
      endif

      inpend_2dstr = Nparticle !eventual 3d_wet are only at the beginning

      nvstart_2dstr = vstart(inpstart_2dstr)
      nvend_2dstr = vend(inpend_2dstr)
      nestart_2dstr = estart(inpstart_2dstr)
      neend_2dstr = eend(inpend_2dstr)
      nfstart_2dstr = fstart(inpstart_2dstr)
      nfend_2dstr = fend(inpend_2dstr)


!additional check
      if (wet_3d(1).GT.0) then                
         if ((nvstart_2dstr-nvend_2dwet).NE.1) then; write(*,*) "Error str wet";stop;endif
         if ((nestart_2dstr-neend_2dwet).NE.1) then; write(*,*) "Error str wet";stop;endif
         if ((nfstart_2dstr-nfend_2dwet).NE.1) then; write(*,*) "Error str wet";stop;endif

         if (nvstart_2dwet.NE.1) then; write(*,*) "Error str wet";stop;endif
         if (nestart_2dwet.NE.1) then; write(*,*) "Error str wet";stop;endif
         if (nfstart_2dwet.NE.1) then; write(*,*) "Error str wet";stop;endif

         if (nvend_2dstr.NE.nvtot) then; write(*,*) "Error str wet";stop;endif
         if (neend_2dstr.NE.netot) then; write(*,*) "Error str wet";stop;endif
         if (nfend_2dstr.NE.nftot) then; write(*,*) "Error str wet";stop;endif

         if (inpstart_2dwet.NE.1) then; write(*,*) "Error str wet";stop;endif
         if ((inpstart_2dstr-inpend_2dwet).NE.1) then; write(*,*) "Error str wet";stop;endif
         if (inpend_2dstr.NE.Nparticle) then; write(*,*) "Error str wet";stop;endif

      else
         if ((nvstart_2dwet.NE.0).AND.(nvend_2dwet.NE.0)) then; write(*,*) "Error str wet";stop;endif
         if ((nestart_2dwet.NE.0).AND.(neend_2dwet.NE.0)) then; write(*,*) "Error str wet";stop;endif
         if ((nfstart_2dwet.NE.0).AND.(nfend_2dwet.NE.0)) then; write(*,*) "Error str wet";stop;endif

         if ((nvstart_2dstr.NE.1).AND.(nvend_2dstr.NE.nvtot)) then; write(*,*) "Error str wet";stop;endif
         if ((nestart_2dstr.NE.1).AND.(neend_2dstr.NE.netot)) then; write(*,*) "Error str wet";stop;endif
         if ((nfstart_2dstr.NE.1).AND.(nfend_2dstr.NE.nftot)) then; write(*,*) "Error str wet";stop;endif

         if ((inpstart_2dwet.NE.0).AND.(inpend_2dwet.NE.0)) then; write(*,*) "Error str wet";stop;endif
         if ((inpstart_2dstr.NE.1).AND.(inpend_2dstr.NE.Nparticle)) then; write(*,*) "Error str wet";stop;endif
      endif


      else
         inpstart_2dwet = 1
         inpend_2dwet = 1
         nvstart_2dwet = 1
         nvend_2dwet = vend(inpend_2dwet)
         nestart_2dwet = 1
         neend_2dwet = eend(inpend_2dwet)
         nfstart_2dwet = 1
         nfend_2dwet = fend(inpend_2dwet)

         inpstart_2dstr = 0
         inpend_2dstr = 0
         nvstart_2dstr = 0
         nvend_2dstr = 0
         nestart_2dstr = 0
         neend_2dstr = 0
         nfstart_2dstr = 0
         nfend_2dstr = 0
      endif




      return
      end subroutine find_indices_2dstr

!===================================================================  
! !francesco
        subroutine tagging_2dwet
        use constants
        use param
        use mls_param
        implicit none

        integer :: i,j,k,v1,v2,inp, inp_wet3d,qualeoff
        integer :: nvstart_3d,nvend_3d
        real(DP), dimension(nvend_2dwet):: scarti
        real(DP) pig,off0, offmin, off
        pig=acos(-1.0)


        off0 = 100000.
        do i=nvstart_2dwet,nvend_2dwet
           inp = vert_to_part(i)
           inp_wet3d = wet_3d(inp)

           nvstart_3d = vstart_3d(inp_wet3d) !attenzione ora punta al 3d 
           nvend_3d   = vend_3d(inp_wet3d)      

           offmin = off0
           do j=nvstart_3d,nvend_3d             
              off = sqrt((xyz(1,i)-xyz_3d(1,j))**2+(xyz(2,i)-xyz_3d(2,j))**2+(xyz(3,i)-xyz_3d(3,j))**2)
              if (off.LT.offmin) then                 
                 qualeoff = j
                 offmin = off
              endif
           enddo
           tag_2dwet(i) = qualeoff
           scarti(i) = offmin
           vert_to_chamb(i)=vert_to_chamb_3d(qualeoff)
           vert_to_chamb4V(i)=vert_to_chamb_3d4V(qualeoff)
       enddo


        do i=nvstart_2dwet,nvend_2dwet
           !write(*,*) "YYY", tag_2dwet(i) ,vert_to_part_3d(qualeoff), scarti(i)
          if (scarti(i).GE.1E-3) then; write(*,*) "Error 2dwet tagging"; stop; endif 
        enddo
        end subroutine tagging_2dwet
! !===================================================================
!===================================================================
        subroutine read_geo_dimPpro(nv,ne,nf,filename)
        
        implicit none
        character*50 :: filename
        integer :: nv,ne,nf
        open(11,file=filename)
        read(11,*)nv,ne,nf
        close(11)

        
        return
        end subroutine read_geo_dimPpro
!...................................................................
!===================================================================                                                                                                       
      subroutine read_geoSingleBodyPpro(NparticleP,geofileP,nvtotP,netotP,nftotP,vstartP,vendP,estartP,eendP,fstartP,fendP,max_n_edge_of_vert, n_edge_of_vertP,vert_of_edgeP, face_of_edgeP,vert_of_faceP, edge_of_faceP ,vert_of_vertP, edge_of_vertP,xyzP,xyzvP,xyzaP,surP ,tri_verP,tri_velP,tri_barP,tri_norP,vel_triP,acc_triP)
      use constants
      use param
      ! use mpih
      ! use mpi_param, only: kstart, kend
      ! use mls_param
      implicit none

      character*50 filename,strucfilename
      character*100 ipfi,ipfip
      character*150,dimension(NparticleP) :: geofileP  

      integer :: i,j,k,v1,v2,v3,inp,e1,e2,e3,count
      real(DP) :: totmass,usVoldom
      integer,dimension(:),allocatable :: nv,ne,nf
      real(DP) :: xmax,xmin,ymax,ymin,zmax,zmin

      integer :: v1d,v2d
      integer :: e1d,e2d,e3d
      integer :: vsi,vei,esi,eei,fsi,fei,csi,cei
      integer :: NparticleP,nvtotP,netotP,nftotP,max_n_edge_of_vert
      integer, dimension (NparticleP) :: nviP,neiP,nfiP,vstartP,vendP,estartP,eendP,fstartP,fendP
      integer, dimension (nvtotP) :: n_edge_of_vertP
      integer, dimension (2,netotP) :: vert_of_edgeP, face_of_edgeP
      integer, dimension (3,nftotP) :: vert_of_faceP, edge_of_faceP 
      integer, dimension (max_n_edge_of_vert,nvtotP) :: vert_of_vertP, edge_of_vertP
      real(DP), dimension (3,nvtotP) :: xyzP,xyzvP,xyzaP
      real(DP), dimension (nftotP) :: surP 
      real(DP), dimension (NparticleP) :: SurfaceP
      real(DP), dimension (9,nftotP) :: tri_verP,tri_velP
      real(DP), dimension (3,nftotP) :: tri_barP,tri_norP,vel_triP,acc_triP

      do inp=1,NparticleP
!      open(109,file='meshes/'//geofileP(inp))
      open(109,file=geofileP(inp))
      read(109,*)nviP(inp),neiP(inp),nfiP(inp)
      close(109)
      end do

!     Read in the vertices, edges and faces from the file
!     Also make connections between the three
      n_edge_of_vertP=0
      face_of_edgeP=0

      do inp=1,NparticleP

      vsi = vstartP(inp) ; vei = vendP(inp)
      esi = estartP(inp) ; eei = eendP(inp)
      fsi = fstartP(inp) ; fei = fendP(inp)

      ! Set particle labels
      ! vert_to_part(vsi:vei) = inp
      ! edge_to_part(esi:eei) = inp
      ! face_to_part(fsi:fei) = inp

        open(109,file=geofileP(inp))
        read(109,*)nviP(inp),neiP(inp),nfiP(inp)
        do i=vsi,vei
          read(109,*)xyzP(1,i),xyzp(2,i),xyzP(3,i)
            xyzP(1,i) = xyzP(1,i)
            xyzP(2,i) = xyzP(2,i)
            xyzP(3,i) = xyzP(3,i)
        end do

        do i=esi,eei
          read(109,*)v1,v2

          v1d = v1 + vsi-1
          v2d = v2 + vsi-1
          vert_of_edgeP(1,i)=v1d
          vert_of_edgeP(2,i)=v2d 

          n_edge_of_vertP(v1d)=n_edge_of_vertP(v1d)+1
          n_edge_of_vertP(v2d)=n_edge_of_vertP(v2d)+1

          vert_of_vertP(n_edge_of_vertP(v1d),v1d)=v2d           
          vert_of_vertP(n_edge_of_vertP(v2d),v2d)=v1d           

          edge_of_vertP(n_edge_of_vertP(v1d),v1d)=i
          edge_of_vertP(n_edge_of_vertP(v2d),v2d)=i
        enddo

        do i=fsi,fei
          read(109,*) e1d, e2d, e3d
                    
          edge_of_faceP(1,i) = e1d + esi - 1
          edge_of_faceP(2,i) = e2d + esi - 1
          edge_of_faceP(3,i) = e3d + esi - 1
           
        end do
 
        do i=fsi,fei
           e1=edge_of_faceP(1,i)
           e2=edge_of_faceP(2,i)

           if (vert_of_edgeP(2,e1).eq.vert_of_edgeP(1,e2)) then
              v1=vert_of_edgeP(1,e1)
              v2=vert_of_edgeP(2,e1)
              v3=vert_of_edgeP(2,e2)
           elseif(vert_of_edgeP(2,e1).eq.vert_of_edgeP(2,e2)) then
              v1=vert_of_edgeP(1,e1)
              v2=vert_of_edgeP(2,e1)
              v3=vert_of_edgeP(1,e2)
           elseif(vert_of_edgeP(1,e1).eq.vert_of_edgeP(1,e2)) then
              v1=vert_of_edgeP(2,e1)
              v2=vert_of_edgeP(1,e1)
              v3=vert_of_edgeP(2,e2)
           else 
              v1=vert_of_edgeP(2,e1)
              v2=vert_of_edgeP(1,e1)
              v3=vert_of_edgeP(1,e2)
           endif 

              vert_of_faceP(1,i)=v1
              vert_of_faceP(2,i)=v2
              vert_of_faceP(3,i)=v3
           enddo

          close(109)
!   Completed reading the gts file

!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Check - vertex cannot be connected to itself
           do i=vsi,vei
              do j=1,n_edge_of_vertP(i)
                 if (vert_of_vertP(j,i).eq.i)  then
                       write(*,*)'Error ',vert_of_vertP(j,i),i
                       stop
                    endif
              enddo
           enddo
!     Check - for faces and edges
           do i=fsi,fei
              e1=edge_of_faceP(1,i)
              e2=edge_of_faceP(2,i)
              e3=edge_of_faceP(3,i)
              if (face_of_edgeP(1,e1).eq.0) then
                 face_of_edgeP(1,e1)=i
              elseif (face_of_edgeP(2,e1).eq.0) then
                 face_of_edgeP(2,e1)=i
              else
                 write(*,*)'Edge error1', i,e1,e2,e3;stop
              endif
              if (face_of_edgeP(1,e2).eq.0) then
                 face_of_edgeP(1,e2)=i
              elseif (face_of_edgeP(2,e2).eq.0) then
                 face_of_edgeP(2,e2)=i
              else
                 write(*,*)'Edge error2';stop
              endif
              if (face_of_edgeP(1,e3).eq.0) then
                 face_of_edgeP(1,e3)=i
              elseif (face_of_edgeP(2,e3).eq.0) then
                 face_of_edgeP(2,e3)=i
              else
                 write(*,*)'Edge error3';stop
              endif
           enddo 

           !Check
           count=0
           do i=esi,eei
              if (face_of_edgeP(1,i).eq.face_of_edgeP(2,i)) then
                 write(*,*)'Error on edges ';stop
              endif
              if (face_of_edgeP(1,i).eq.0.or.face_of_edgeP(2,i).eq.0) then
                 count=count+1
              endif
           enddo
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!non so perche mi da problemi cuda...
      call convert_geoCPU(vsi,vei,esi,eei,        &
            fsi,fei,xyzP(:,vsi:vei),            &
            xyzvP(:,vsi:vei),xyzaP(:,vsi:vei),   &
            vert_of_faceP(:,fsi:fei),           &
            tri_verP(:,fsi:fei),                &
            tri_velP(:,fsi:fei),                &
            tri_barP(:,fsi:fei),                &
            vel_triP(:,fsi:fei),                &
            acc_triP(:,fsi:fei))

      call calculate_normalCPU(vsi,vei,fsi,             &
                             fei,xyzP(:,vsi:vei),      &
                             vert_of_faceP(:,fsi:fei),  &
                             tri_norP(:,fsi:fei))

!         call calculate_angle(theta0(esi:eei),esi,eei,fsi,fei,      &
!                        face_of_edge(:,esi:eei),tri_nor(:,fsi:fei))

!         call calculate_anglealp(aalpha0(:,fsi:fei),vsi,vei,fsi,fei,  &
!                        xyz0(:,vsi:vei),                            &
!                        vert_of_face(:,fsi:fei))

!         call calculate_distance(dist0(esi:eei),vsi,vei,esi,eei,    &
!                        xyz0(:,vsi:vei),vert_of_edge(:,esi:eei))

        call calculate_areaCPU(SurfaceP(inp),vsi,vei,fsi,fei,         &
                       xyzP(:,vsi:vei),vert_of_faceP(:,fsi:fei),    &
                       surP(fsi:fei))
!         sur=sur0

!         call calculate_volume(Volume0(inp),vsi,vei,fsi,fei,  &
!                        xyz0(:,vsi:vei),vert_of_face(:,fsi:fei))

!         call find_quartet(v1234(:,esi:eei),vsi,vei,esi,eei,fsi,fei,       &
!                       xyz0(:,vsi:vei),                                    &
!                       face_of_edge(:,esi:eei),edge_of_face(:,fsi:fei),    &
!                       vert_of_edge(:,esi:eei),tri_nor(:,fsi:fei))

!         if(myid.eq.0)write(*,*)'Finished reading geo. no. ',inp


!       totmass = 1.0
!       !uspm = 1.0

      end do

        return
        end subroutine read_geoSingleBodyPpro
!===================================================================
      subroutine prepare_volume_chamb_shift
      use constants
      use param
      ! use mpih
      ! use mpi_param, only: kstart, kend
      use mls_param
      implicit none

      character*150 filename,strucfilename
      character*100 ipfi,ipfip
      character*150  geofileP  
      integer :: nvtotP,netotP,nftotP,inp
      real(DP) :: SurfaceP,VolumeP
      real(DP) :: SurTap_AO,VolTap_AO
      real(DP) :: SurTap_MI,VolTap_MI
      real(DP) :: SurTap_PU,VolTap_PU
      real(DP) :: SurTap_TR,VolTap_TR

      Volume_chambShift=0.0

#ifndef SOLOUNO
      write(*,*) "Load vessels Caps"
      !LEFT HEART
      !valvola aortica
      inp=2
      filename='./meshes/'//trim(geofile(inp))//'_tappo.gts'
      call read_geo_dimPpro(nvtotP,netotP,nftotP,filename)
      call read_geoSingleBodyTappiChamb(filename,nvtotP,netotP,nftotP,max_n_edge_of_vert,VolumeP,SurfaceP)
      SurTap_AO = abs(SurfaceP)
      VolTap_AO = abs(VolumeP)
      !valvola mitrale
      inp=3
      filename='./meshes/'//trim(geofile(inp))//'_tappo.gts'
      call read_geo_dimPpro(nvtotP,netotP,nftotP,filename)
      call read_geoSingleBodyTappiChamb(filename,nvtotP,netotP,nftotP,max_n_edge_of_vert,VolumeP,SurfaceP)
      SurTap_MI = abs(SurfaceP)
      VolTap_MI = abs(VolumeP)
      
      !left ventricle
      Surface_chambShift(1) = SurTap_AO + SurTap_MI
      Volume_chambShift(1) = VolTap_AO + VolTap_MI

      !left atrium
      Surface_chambShift(2) = SurTap_MI
      Volume_chambShift(2) = VolTap_MI


      !RIGHT HEART
      !valvola polmonare
      inp=4
      filename='./meshes/'//trim(geofile(inp))//'_tappo.gts'
      call read_geo_dimPpro(nvtotP,netotP,nftotP,filename)
      call read_geoSingleBodyTappiChamb(filename,nvtotP,netotP,nftotP,max_n_edge_of_vert,VolumeP,SurfaceP)
      SurTap_PU = abs(SurfaceP)
      VolTap_PU = abs(VolumeP)
      !valvola tricuspide
      inp=5
      filename='./meshes/'//trim(geofile(inp))//'_tappo.gts'
      call read_geo_dimPpro(nvtotP,netotP,nftotP,filename)
      call read_geoSingleBodyTappiChamb(filename,nvtotP,netotP,nftotP,max_n_edge_of_vert,VolumeP,SurfaceP)
      SurTap_TR = abs(SurfaceP)
      VolTap_TR = abs(VolumeP)

      !right ventricle
      Surface_chambShift(3) = SurTap_PU + SurTap_TR
      Volume_chambShift(3) = VolTap_PU + VolTap_TR

      !right atrium
      Surface_chambShift(4) = SurTap_TR
      Volume_chambShift(4) = VolTap_TR
#endif

      
      !veins/arteries no matters
      write(*,*) "Surface Shift", Surface_chambShift
      write(*,*) "Volume Shift", Volume_chambShift

      return
      end subroutine prepare_volume_chamb_shift
!===================================================================

      
      subroutine read_geoSingleBodyTappiChamb(geofileP,nvtotP,netotP,nftotP,max_n_edge_of_vert,VolumeP,SurfaceP)
      use constants
      use param
      ! use mpih
      ! use mpi_param, only: kstart, kend
      ! use mls_param
      implicit none

      character*50 filename,strucfilename
      character*100 ipfi,ipfip
      character*150  geofileP  

      integer :: i,j,k,v1,v2,v3,e1,e2,e3,count
      real(DP) :: totmass,usVoldom
      integer,dimension(:),allocatable :: nv,ne,nf
      real(DP) :: xmax,xmin,ymax,ymin,zmax,zmin

      integer :: v1d,v2d
      integer :: e1d,e2d,e3d
      integer :: vsi,vei,esi,eei,fsi,fei,csi,cei
      integer :: NparticleP,nvtotP,netotP,nftotP,max_n_edge_of_vert
      integer :: nviP,neiP,nfiP,vstartP,vendP,estartP,eendP,fstartP,fendP
      integer, dimension (nvtotP) :: n_edge_of_vertP
      integer, dimension (2,netotP) :: vert_of_edgeP, face_of_edgeP
      integer, dimension (3,nftotP) :: vert_of_faceP, edge_of_faceP 
      integer, dimension (max_n_edge_of_vert,nvtotP) :: vert_of_vertP, edge_of_vertP
      real(DP), dimension (3,nvtotP) :: xyzP,xyzvP,xyzaP
      real(DP), dimension (nftotP) :: surP 
      real(DP) :: SurfaceP,VolumeP
      real(DP), dimension (9,nftotP) :: tri_verP,tri_velP
      real(DP), dimension (3,nftotP) :: tri_barP,tri_norP,vel_triP,acc_triP


!     Read in the vertices, edges and faces from the file
!     Also make connections between the three
      n_edge_of_vertP=0
      face_of_edgeP=0

      vsi=1
      vei=nvtotP
      esi=1
      eei=netotP
      fsi=1
      fei=nftotP


      ! Set particle labels
      ! vert_to_part(vsi:vei) = inp
      ! edge_to_part(esi:eei) = inp
      ! face_to_part(fsi:fei) = inp

        open(109,file=geofileP)
        read(109,*)nviP,neiP,nfiP
        do i=vsi,vei
          read(109,*)xyzP(1,i),xyzp(2,i),xyzP(3,i)
            xyzP(1,i) = xyzP(1,i)
            xyzP(2,i) = xyzP(2,i)
            xyzP(3,i) = xyzP(3,i)
        end do

        do i=esi,eei
          read(109,*)v1,v2

          v1d = v1 + vsi-1
          v2d = v2 + vsi-1
          vert_of_edgeP(1,i)=v1d
          vert_of_edgeP(2,i)=v2d 

          n_edge_of_vertP(v1d)=n_edge_of_vertP(v1d)+1
          n_edge_of_vertP(v2d)=n_edge_of_vertP(v2d)+1

          vert_of_vertP(n_edge_of_vertP(v1d),v1d)=v2d           
          vert_of_vertP(n_edge_of_vertP(v2d),v2d)=v1d           

          edge_of_vertP(n_edge_of_vertP(v1d),v1d)=i
          edge_of_vertP(n_edge_of_vertP(v2d),v2d)=i
        enddo

        do i=fsi,fei
          read(109,*) e1d, e2d, e3d
                    
          edge_of_faceP(1,i) = e1d + esi - 1
          edge_of_faceP(2,i) = e2d + esi - 1
          edge_of_faceP(3,i) = e3d + esi - 1
           
        end do
 
        do i=fsi,fei
           e1=edge_of_faceP(1,i)
           e2=edge_of_faceP(2,i)

           if (vert_of_edgeP(2,e1).eq.vert_of_edgeP(1,e2)) then
              v1=vert_of_edgeP(1,e1)
              v2=vert_of_edgeP(2,e1)
              v3=vert_of_edgeP(2,e2)
           elseif(vert_of_edgeP(2,e1).eq.vert_of_edgeP(2,e2)) then
              v1=vert_of_edgeP(1,e1)
              v2=vert_of_edgeP(2,e1)
              v3=vert_of_edgeP(1,e2)
           elseif(vert_of_edgeP(1,e1).eq.vert_of_edgeP(1,e2)) then
              v1=vert_of_edgeP(2,e1)
              v2=vert_of_edgeP(1,e1)
              v3=vert_of_edgeP(2,e2)
           else 
              v1=vert_of_edgeP(2,e1)
              v2=vert_of_edgeP(1,e1)
              v3=vert_of_edgeP(1,e2)
           endif 

              vert_of_faceP(1,i)=v1
              vert_of_faceP(2,i)=v2
              vert_of_faceP(3,i)=v3
           enddo

          close(109)
!   Completed reading the gts file

!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Check - vertex cannot be connected to itself
           do i=vsi,vei
              do j=1,n_edge_of_vertP(i)
                 if (vert_of_vertP(j,i).eq.i)  then
                       write(*,*)'Error ',vert_of_vertP(j,i),i
                       stop
                    endif
              enddo
           enddo
!     Check - for faces and edges
           do i=fsi,fei
              e1=edge_of_faceP(1,i)
              e2=edge_of_faceP(2,i)
              e3=edge_of_faceP(3,i)
              if (face_of_edgeP(1,e1).eq.0) then
                 face_of_edgeP(1,e1)=i
              elseif (face_of_edgeP(2,e1).eq.0) then
                 face_of_edgeP(2,e1)=i
              else
                 write(*,*)'Edge error1', i,e1,e2,e3;stop
              endif
              if (face_of_edgeP(1,e2).eq.0) then
                 face_of_edgeP(1,e2)=i
              elseif (face_of_edgeP(2,e2).eq.0) then
                 face_of_edgeP(2,e2)=i
              else
                 write(*,*)'Edge error2';stop
              endif
              if (face_of_edgeP(1,e3).eq.0) then
                 face_of_edgeP(1,e3)=i
              elseif (face_of_edgeP(2,e3).eq.0) then
                 face_of_edgeP(2,e3)=i
              else
                 write(*,*)'Edge error3';stop
              endif
           enddo 

           !Check
           count=0
           do i=esi,eei
              if (face_of_edgeP(1,i).eq.face_of_edgeP(2,i)) then
                 write(*,*)'Error on edges ';stop
              endif
              if (face_of_edgeP(1,i).eq.0.or.face_of_edgeP(2,i).eq.0) then
                 count=count+1
              endif
           enddo
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!non so perche mi da problemi cuda...
      call convert_geoCPU(vsi,vei,esi,eei,        &
            fsi,fei,xyzP(:,vsi:vei),            &
            xyzvP(:,vsi:vei),xyzaP(:,vsi:vei),   &
            vert_of_faceP(:,fsi:fei),           &
            tri_verP(:,fsi:fei),                &
            tri_velP(:,fsi:fei),                &
            tri_barP(:,fsi:fei),                &
            vel_triP(:,fsi:fei),                &
            acc_triP(:,fsi:fei))

      call calculate_normalCPU(vsi,vei,fsi,             &
                             fei,xyzP(:,vsi:vei),      &
                             vert_of_faceP(:,fsi:fei),  &
                             tri_norP(:,fsi:fei))

!         call calculate_angle(theta0(esi:eei),esi,eei,fsi,fei,      &
!                        face_of_edge(:,esi:eei),tri_nor(:,fsi:fei))

!         call calculate_anglealp(aalpha0(:,fsi:fei),vsi,vei,fsi,fei,  &
!                        xyz0(:,vsi:vei),                            &
!                        vert_of_face(:,fsi:fei))

!         call calculate_distance(dist0(esi:eei),vsi,vei,esi,eei,    &
!                        xyz0(:,vsi:vei),vert_of_edge(:,esi:eei))

        call calculate_areaCPU(SurfaceP,vsi,vei,fsi,fei,         &
                       xyzP(:,vsi:vei),vert_of_faceP(:,fsi:fei),    &
                       surP(fsi:fei))
!         sur=sur0

        call calculate_volumeCPU(VolumeP,vsi,vei,fsi,fei,  &
                       xyzP(:,vsi:vei),vert_of_faceP(:,fsi:fei))

!         call find_quartet(v1234(:,esi:eei),vsi,vei,esi,eei,fsi,fei,       &
!                       xyz0(:,vsi:vei),                                    &
!                       face_of_edge(:,esi:eei),edge_of_face(:,fsi:fei),    &
!                       vert_of_edge(:,esi:eei),tri_nor(:,fsi:fei))

!         if(myid.eq.0)write(*,*)'Finished reading geo. no. ',inp


!       totmass = 1.0
!       !uspm = 1.0


        return
        end subroutine read_geoSingleBodyTappiChamb
!===================================================================
        







!===================================================================                         
                                                                                  
      subroutine read_geoSingleBodyEF_2d
      use constants
      use param
      use mpih
      use mpi_param, only: kstart, kend
      use mls_param
      implicit none

      character*50 filename,strucfilename
      character*100 ipfi,ipfip

      integer :: i,j,k,v1,v2,v3,inp,e1,e2,e3,count
      integer :: f1,f2,checkF1,checkF2,chamb,fminsetto
      real(DP) :: totmass,usVoldom
      integer,dimension(:),allocatable :: nv,ne,nf
      real(DP),dimension(:,:),allocatable :: xyz_setto
      integer,dimension(:),allocatable :: n_edge_of_vert_setto
      integer,dimension(:,:),allocatable :: edge_of_face_setto
      integer,dimension(:,:),allocatable :: edge_of_vert_setto
      integer,dimension(:,:),allocatable :: vert_of_edge_setto
      integer,dimension(:,:),allocatable :: vert_of_vert_setto
      integer,dimension(:,:),allocatable :: vert_of_face_setto

      real(DP) :: xmax,xmin,ymax,ymin,zmax,zmin
      real(DP) :: xbar,ybar,zbar,offmin,off,offmin0
      real(DP) :: xbarslave,ybarslave,zbarslave

      integer :: v1d,v2d,vmaster,fslave,fmaster,cslave
      integer :: e1d,e2d,e3d,c1,c2
      integer :: densCoupBundPurk,densCoupBundAtri,densCoupBundVentr
      integer :: vsi,vei,esi,eei,fsi,fei,csi,cei,vv,vnext
      integer :: nBundPurk1,nBundAtri1,nPurkVentr1,count3EF_2d
      real(DP) :: A1,A2,A3,A4,A5,A6,Lvm,Lhm,SIndex,Laxis
      real(DP) :: xAp,yAp,zAp,xApTr,yApTr,zAptr,xNORM,yNORM,zNORM,sca
      integer :: nv_setto,ne_setto,nf_setto

      do inp=1,NparticleEF_2d

      open(109,file='meshes/'//trim(geofileEF_2d(inp))//'.gts')
      read(109,*)nviEF_2d(inp),neiEF_2d(inp),nfiEF_2d(inp)
      !maxnv=max(maxnv,nv(inp))
      !maxne=max(maxne,ne(inp))
      !maxnf=max(maxnf,nf(inp))
      close(109)
      end do

!     Allocations using vertices,edges and faces as parameters     

      call allocate_trigeoEF_2d

!     Read in the vertices, edges and faces from the file
!     Also make connections between the three
      n_edge_of_vertEF_2d=0
      face_of_edgeEF_2d=0

      do inp=1,NparticleEF_2d

      vsi = vstartEF_2d(inp) ; vei = vendEF_2d(inp)
      esi = estartEF_2d(inp) ; eei = eendEF_2d(inp)
      fsi = fstartEF_2d(inp) ; fei = fendEF_2d(inp)

        open(109,file='meshes/'//trim(geofileEF_2d(inp))//'.gts')
        read(109,*)nviEF_2d(inp),neiEF_2d(inp),nfiEF_2d(inp)
        do i=vsi,vei
          read(109,*)xyz0EF_2d(1,i),xyz0EF_2d(2,i),xyz0EF_2d(3,i)
          ! xyz0EF_2d(1,i)=xyz0EF_2d(1,i)
          ! xyz0EF_2d(2,i)=xyz0EF_2d(2,i)
          ! xyz0EF_2d(3,i)=xyz0EF_2d(3,i)
        end do

!     position correction to bring in 0-6
 
      xyzEF_2d = xyz0EF_2d
        do i=esi,eei
          read(109,*)v1,v2

          v1d = v1 + vsi-1
          v2d = v2 + vsi-1
          vert_of_edgeEF_2d(1,i)=v1d
          vert_of_edgeEF_2d(2,i)=v2d 

          n_edge_of_vertEF_2d(v1d)=n_edge_of_vertEF_2d(v1d)+1
          n_edge_of_vertEF_2d(v2d)=n_edge_of_vertEF_2d(v2d)+1

          vert_of_vertEF_2d(n_edge_of_vertEF_2d(v1d),v1d)=v2d           
          vert_of_vertEF_2d(n_edge_of_vertEF_2d(v2d),v2d)=v1d           

          edge_of_vertEF_2d(n_edge_of_vertEF_2d(v1d),v1d)=i
          edge_of_vertEF_2d(n_edge_of_vertEF_2d(v2d),v2d)=i
        enddo

        do i=fsi,fei
          read(109,*) e1d, e2d, e3d
                    
          edge_of_faceEF_2d(1,i) = e1d + esi - 1
          edge_of_faceEF_2d(2,i) = e2d + esi - 1
          edge_of_faceEF_2d(3,i) = e3d + esi - 1
           
        end do
 
        do i=fsi,fei
           e1=edge_of_faceEF_2d(1,i)
           e2=edge_of_faceEF_2d(2,i)

           if (vert_of_edgeEF_2d(2,e1).eq.vert_of_edgeEF_2d(1,e2)) then
              v1=vert_of_edgeEF_2d(1,e1)
              v2=vert_of_edgeEF_2d(2,e1)
              v3=vert_of_edgeEF_2d(2,e2)
           elseif(vert_of_edgeEF_2d(2,e1).eq.vert_of_edgeEF_2d(2,e2)) then
              v1=vert_of_edgeEF_2d(1,e1)
              v2=vert_of_edgeEF_2d(2,e1)
              v3=vert_of_edgeEF_2d(1,e2)
           elseif(vert_of_edgeEF_2d(1,e1).eq.vert_of_edgeEF_2d(1,e2)) then
              v1=vert_of_edgeEF_2d(2,e1)
              v2=vert_of_edgeEF_2d(1,e1)
              v3=vert_of_edgeEF_2d(2,e2)
           else 
              v1=vert_of_edgeEF_2d(2,e1)
              v2=vert_of_edgeEF_2d(1,e1)
              v3=vert_of_edgeEF_2d(1,e2)
           endif 

              vert_of_faceEF_2d(1,i)=v1
              vert_of_faceEF_2d(2,i)=v2
              vert_of_faceEF_2d(3,i)=v3
           enddo

          close(109)

!   Completed reading the gts file
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Check - vertex cannot be connected to itself
           do i=vsi,vei
              do j=1,n_edge_of_vertEF_2d(i)
                 if (vert_of_vertEF_2d(j,i).eq.i)  &
                       write(*,*)'Error ',vert_of_vertEF_2d(j,i),i
              enddo
           enddo
!     Check - for faces and edges
           do i=fsi,fei
              e1=edge_of_faceEF_2d(1,i)
              e2=edge_of_faceEF_2d(2,i)
              e3=edge_of_faceEF_2d(3,i)
              if (face_of_edgeEF_2d(1,e1).eq.0) then
                 face_of_edgeEF_2d(1,e1)=i
              elseif (face_of_edgeEF_2d(2,e1).eq.0) then
                 face_of_edgeEF_2d(2,e1)=i
              else
                 write(*,*)'Edge error1', i,e1,e2,e3;stop
              endif
              if (face_of_edgeEF_2d(1,e2).eq.0) then
                 face_of_edgeEF_2d(1,e2)=i
              elseif (face_of_edgeEF_2d(2,e2).eq.0) then
                 face_of_edgeEF_2d(2,e2)=i
              else
                 write(*,*)'Edge error2';stop
              endif
              if (face_of_edgeEF_2d(1,e3).eq.0) then
                 face_of_edgeEF_2d(1,e3)=i
              elseif (face_of_edgeEF_2d(2,e3).eq.0) then
                 face_of_edgeEF_2d(2,e3)=i
              else
                 write(*,*)'Edge error3';stop
              endif
           enddo 

           !Check
           count=0
           do i=esi,eei
              if (face_of_edgeEF_2d(1,i).eq.face_of_edgeEF_2d(2,i)) then
                 write(*,*)'Error on edges '
              endif
              if (face_of_edgeEF_2d(1,i).eq.0.or.face_of_edgeEF_2d(2,i).eq.0) then
                 count=count+1
              endif
           enddo
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         call convert_geoEF_2d(vsi,vei,esi,eei,        &
            fsi,fei,xyzEF_2d(:,vsi:vei),            &
            vert_of_faceEF_2d(:,fsi:fei),           &
            tri_barEF_2d(:,fsi:fei))

        call calculate_normal(vsi,vei,fsi,             &
                             fei,xyz0EF_2d(:,vsi:vei),      &
                             vert_of_faceEF_2d(:,fsi:fei),  &
                            tri_norEF_2d(:,fsi:fei))


        call calculate_distance(dist0EF_2d(esi:eei),vsi,vei,esi,eei,    &
                       xyz0EF_2d(:,vsi:vei),vert_of_edgeEF_2d(:,esi:eei))

        call calculate_area(Surface0EF_2d(inp),vsi,vei,fsi,fei,         &
                       xyz0EF_2d(:,vsi:vei),vert_of_faceEF_2d(:,fsi:fei),    &
                       sur0EF_2d(fsi:fei))
        surEF_2d=sur0EF_2d
        distEF_2d=dist0EF_2d


        if(myid.eq.0)write(*,*)'Finished reading geo. no. ',inp,esi,eei


        call calculate_ginterp_edge_facesEF_2d(vsi,vei,esi,eei,fsi,fei,         &
                       xyz0EF_2d(:,vsi:vei),vert_of_edgeEF_2d(:,esi:eei),    &
                       face_of_edgeEF_2d(:,esi:eei),tri_barEF_2d(:,fsi:fei), &
                       versCFedgeEF_2d(:,esi:eei),distCFedgeEF_2d(esi:eei),g1interpedgeEF_2d(esi:eei))

        call calculate_normals_edge_facesEF_2d(vsi,vei,esi,eei,fsi,fei,xyz0EF_2d(:,vsi:vei),  &
                       vert_of_edgeEF_2d(:,esi:eei),edge_of_faceEF_2d(:,fsi:fei),   & 
                       tri_barEF_2d(:,fsi:fei),tri_norEF_2d(:,fsi:fei),normaledgeoffacesEF_2d(:,:,fsi:fei))

        call calculate_fiberdirectionEF_2d(vsi,vei,fsi,fei,xyz0EF_2d(:,vsi:vei),  &
                       tri_barEF_2d(:,fsi:fei),tri_norEF_2d(:,fsi:fei),AmatrFibersEF_2d(:,:,fsi:fei))

        face_of_vertEF_2d(:,:)=0
        n_face_of_vertEF_2d(:)=0
        do i=1,nvtotEF_2d
           count=0
           do j=1,n_edge_of_vertEF_2d(i)
              e1=edge_of_vertEF_2d(j,i);
              f1=face_of_edgeEF_2d(1,e1);
              f2=face_of_edgeEF_2d(2,e1);

              checkF1=0
              do k=1,max_n_edge_of_vert
                 if (f1.EQ.face_of_vertEF_2d(k,i)) checkF1=checkF1+1
              enddo
              if ((checkF1.EQ.0).AND.(f1.NE.0)) then
                 count=count+1;
                 face_of_vertEF_2d(count,i)=f1;
              endif

              checkF2=0
              do k=1,max_n_edge_of_vert
                 if (f2.EQ.face_of_vertEF_2d(k,i)) checkF2=checkF2+1
              enddo
              if ((checkF2.EQ.0).AND.(f2.NE.0)) then
                 count=count+1;
                 face_of_vertEF_2d(count,i)=f2;
              endif
           enddo !j
           n_face_of_vertEF_2d(i)=count;
        enddo

      end do
      
!##################################################################
!more connettivity for smoothing while advected
      count2EF_2d = 0
      do i=1,netotEF_2d
         f1=face_of_edgeEF_2d(1,i)
         f2=face_of_edgeEF_2d(2,i)
         if ((f1.EQ.0).OR.(f2.EQ.0)) then !cosi li conto tutti 2volte
            v1=vert_of_edgeEF_2d(1,i);   
            v2=vert_of_edgeEF_2d(2,i);
            count2EF_2d = count2EF_2d + 1;
            count2EF_2d = count2EF_2d + 1;
         endif
      enddo

      count3EF_2d = 0;
      allocate(vert_boundEF_2d(count2EF_2d))
      do i=1,netotEF_2d
         f1=face_of_edgeEF_2d(1,i)
         f2=face_of_edgeEF_2d(2,i)
         if ((f1.EQ.0).OR.(f2.EQ.0)) then !cosi li conto tutti 2volte
            v1=vert_of_edgeEF_2d(1,i) ;    
            v2=vert_of_edgeEF_2d(2,i);
            count3EF_2d = count3EF_2d + 1;
            vert_boundEF_2d(count3EF_2d) = v1;
            count3EF_2d = count3EF_2d + 1;
            vert_boundEF_2d(count3EF_2d) = v2;
         endif
      enddo

      if (count3EF_2d.NE.count2EF_2d) then
         write(*,*) "error1 EF_2d", count2EF_2d,count3EF_2d; stop; 
      endif

      ! do i=1,count2EF_2d
      !    write(*,*) "A",i,vert_boundEF_2d(i);
      ! enddo
      ! stop


      allocate(vert_of_vert_boundEF_2d(2,count2EF_2d))
      vert_of_vert_boundEF_2d=0.0D0
      do i=1,count2EF_2d
         vv=vert_boundEF_2d(i);
         do j=1,n_edge_of_vertEF_2d(vv)
            e1=edge_of_vertEF_2d(j,vv);
            f1=face_of_edgeEF_2d(1,e1);
            f2=face_of_edgeEF_2d(2,e1);
!            write(*,*) "S",i,vv,j,e1,f1,f2
            if ((f1.EQ.0).OR.(f2.EQ.0)) then
               v1=vert_of_edgeEF_2d(1,e1);
               v2=vert_of_edgeEF_2d(2,e1);
               if (v1.EQ.vv) then 
                  vnext=v2;
               elseif (v2.EQ.vv) then
                  vnext=v1;
               else
                   write(*,*) "error2 EF_2d"; stop
               endif

         
               if (vert_of_vert_boundEF_2d(1,i).EQ.0) then
                  vert_of_vert_boundEF_2d(1,i)=vnext;
               else
                  if (vert_of_vert_boundEF_2d(2,i).EQ.0) then
                     vert_of_vert_boundEF_2d(2,i)=vnext;
                  else
                     write(*,*) "error3 EF_2d",i,j; stop
                  endif
               endif

            endif
         enddo
      
      enddo

!##################################################################
!COUPLING boundles-->Purkinje
!conta nodi coupling 
      nBundPurk=0
      do i=vstartEF_1d(3),vendEF_1d(3)
         if ((n_edge_of_vertEF_1d(i).EQ.1).AND.(i.NE.vAVslave(3))) then
            nBundPurk=nBundPurk+1
         endif
      enddo
      write(*,*) "nBundPurk = ", nBundPurk

      allocate(vPurkmaster(nBundPurk))
      allocate(fPurkslave(nBundPurk))

      nBundPurk1=0
      do i=vstartEF_1d(3),vendEF_1d(3)
         if ((n_edge_of_vertEF_1d(i).EQ.1).AND.(i.NE.vAVslave(3))) then
            nBundPurk1=nBundPurk1+1
            vPurkmaster(nBundPurk1)=i               
         endif
      enddo
      if (nBundPurk1.NE.nBundPurk)then;write(*,*)"err BundPurk";stop;endif


!trova nodi coupling
      offmin0 = 1E5
      do i=1,nBundPurk
         vmaster=vPurkmaster(i)
         offmin=offmin0
         do j=1,nftotEF_2d
            off = sqrt((xyzEF_1d(1,vmaster)-tri_barEF_2d(1,j))**2+(xyzEF_1d(2,vmaster)-tri_barEF_2d(2,j))**2+(xyzEF_1d(3,vmaster)-tri_barEF_2d(3,j))**2)
            if (off.LT.offmin) then
               fslave = j
               offmin = off
            endif
         enddo

         fPurkslave(i)=fslave         
         write(*,*) "vmaster", xyzEF_1d(1:3,vmaster)
         write(*,*) "fslave",  tri_barEF_2d(1:3,fslave)
      enddo

!COUPLING boundles-->Atri
!conta nodi coupling 
      count=0
      nBundAtri=0
      ! densCoupBundAtri=10
      densCoupBundAtri=2
      do i=vstartEF_1d(1),vendEF_1d(1)
!         if (xyzEF_1d(1,i).GT.-11.) then !ricambia
            count=count+1
            if (mod(count,densCoupBundAtri).EQ.0) then
               nBundAtri=nBundAtri+1
            endif
!         endif
      enddo
      write(*,*) "nBundAtri = ", nBundAtri

      allocate(vAtrimaster(nBundAtri))
      allocate(cAtrislave(nBundAtri))

      allocate(Bund_tstim(nBundAtri))
      do i=1,nBundAtri
         Bund_tstim(i)=-10.d0
      enddo

      count=0
      nBundAtri1=0
      do i=vstartEF_1d(1),vendEF_1d(1)
!         if (xyzEF_1d(1,i).GT.-11.) then !ricambia
            count=count+1
            if (mod(count,densCoupBundAtri).EQ.0) then
               nBundAtri1=nBundAtri1+1
               vAtrimaster(nBundAtri1)=i               
            endif
!         endif
      enddo
      if (nBundAtri1.NE.nBundAtri)then;write(*,*)"err BundAtri";stop;endif

!trova nodi coupling
      offmin0 = 1E5
      csi = cstart_3d(1) ; cei = cend_3d(1)
      do i=1,nBundAtri
         vmaster=vAtrimaster(i)
         offmin=offmin0
         do j=csi,cei
            chamb=cell_to_chamb_3d(j)
            if ((chamb.EQ.2).OR.(chamb.EQ.4)) then
               off = sqrt((xyzEF_1d(1,vmaster)-cell_bar(1,j))**2+(xyzEF_1d(2,vmaster)-cell_bar(2,j))**2+(xyzEF_1d(3,vmaster)-cell_bar(3,j))**2)
               if (off.LT.offmin) then
                  cslave = j
                  offmin = off
               endif
            endif
         enddo

         cAtrislave(i)=cslave         

         ! write(*,*) "vmaster", xyzEF_1d(1:3,vmaster)
         ! write(*,*) "cslave",  cell_bar(1:3,cslave)

      enddo



!CARICA SETTO PER ESCLUDERLO
      LabelPurkSetto(:)=0
      if (EscSettoPurk.EQ.1) then
         inp=1
         open(109,file='meshes/'//trim(geofileEF_2d(inp))//'_setto.gts')
         read(109,*)nv_setto,ne_setto,nf_setto
         close(109)

         allocate(xyz_setto(3,nv_setto))
         allocate(n_edge_of_vert_setto(nv_setto))
         allocate(edge_of_face_setto(3,nf_setto))
         allocate(edge_of_vert_setto(max_n_edge_of_vert,nv_setto))
         allocate(vert_of_edge_setto(2,ne_setto))
         allocate(vert_of_face_setto(3,nf_setto))
         allocate(vert_of_vert_setto(max_n_edge_of_vert,nv_setto))

         n_edge_of_vert_setto=0
         !ok usare queste variabili perche siamo fuori loop corpi                                                                                                                          
         vsi=1
         vei=nv_setto
         esi=1
         eei=ne_setto
         fsi=1
         fei=nf_setto

         open(109,file='meshes/'//trim(geofileEF_2d(inp))//'_setto.gts')
         read(109,*)nv_setto,ne_setto,nf_setto
         do i=vsi,vei
            read(109,*)xyz_setto(1,i),xyz_setto(2,i),xyz_setto(3,i)
            xyz_setto(1,i)=xyz_setto(1,i)
            xyz_setto(2,i)=xyz_setto(2,i)
            xyz_setto(3,i)=xyz_setto(3,i)
         end do
         do i=esi,eei
          read(109,*)v1,v2

          v1d = v1 + vsi-1
          v2d = v2 + vsi-1
          vert_of_edge_setto(1,i)=v1d
          vert_of_edge_setto(2,i)=v2d

          n_edge_of_vert_setto(v1d)=n_edge_of_vert_setto(v1d)+1
          n_edge_of_vert_setto(v2d)=n_edge_of_vert_setto(v2d)+1

          vert_of_vert_setto(n_edge_of_vert_setto(v1d),v1d)=v2d
          vert_of_vert_setto(n_edge_of_vert_setto(v2d),v2d)=v1d

          edge_of_vert_setto(n_edge_of_vert_setto(v1d),v1d)=i
          edge_of_vert_setto(n_edge_of_vert_setto(v2d),v2d)=i
        enddo

        do i=fsi,fei
          read(109,*) e1d, e2d, e3d

          edge_of_face_setto(1,i) = e1d + esi - 1
          edge_of_face_setto(2,i) = e2d + esi - 1
          edge_of_face_setto(3,i) = e3d + esi - 1

        end do
        do i=fsi,fei
           e1=edge_of_face_setto(1,i)
           e2=edge_of_face_setto(2,i)

           if (vert_of_edge_setto(2,e1).eq.vert_of_edge_setto(1,e2)) then
              v1=vert_of_edge_setto(1,e1)
              v2=vert_of_edge_setto(2,e1)
              v3=vert_of_edge_setto(2,e2)
           elseif(vert_of_edge_setto(2,e1).eq.vert_of_edge_setto(2,e2)) then
              v1=vert_of_edge_setto(1,e1)
              v2=vert_of_edge_setto(2,e1)
              v3=vert_of_edge_setto(1,e2)
           elseif(vert_of_edge_setto(1,e1).eq.vert_of_edge_setto(1,e2)) then
              v1=vert_of_edge_setto(2,e1)
              v2=vert_of_edge_setto(1,e1)
              v3=vert_of_edge_setto(2,e2)
           else
              v1=vert_of_edge_setto(2,e1)
              v2=vert_of_edge_setto(1,e1)
              v3=vert_of_edge_setto(1,e2)
           endif

              vert_of_face_setto(1,i)=v1
              vert_of_face_setto(2,i)=v2
              vert_of_face_setto(3,i)=v3
           enddo

           close(109)
           offmin0 = 1E5
        do i=1,nf_setto
           v1=vert_of_face_setto(1,i)
           v2=vert_of_face_setto(2,i)
           v3=vert_of_face_setto(3,i)
           xbar=(xyz_setto(1,v1)+xyz_setto(1,v2)+xyz_setto(1,v3))/3.0D0
           ybar=(xyz_setto(2,v1)+xyz_setto(2,v2)+xyz_setto(2,v3))/3.0D0
           zbar=(xyz_setto(3,v1)+xyz_setto(3,v2)+xyz_setto(3,v3))/3.0D0

           offmin=offmin0
           do j=1,nftotEF_2d
                 off = sqrt((tri_barEF_2d(1,j)-xbar)**2+(tri_barEF_2d(2,j)-ybar)**2+(tri_barEF_2d(3,j)-zbar)**2)
                 if (off.LT.offmin) then
                    fminsetto = j
                    offmin = off
                 endif
           enddo

           LabelPurkSetto(fminsetto)=1

        enddo


      endif !end LABELSETTO 
      if (EscSettoPurk.EQ.1) then
         deallocate(xyz_setto)
         deallocate(n_edge_of_vert_setto)
         deallocate(edge_of_face_setto)
         deallocate(edge_of_vert_setto)
         deallocate(vert_of_edge_setto)
         deallocate(vert_of_face_setto)
         deallocate(vert_of_vert_setto)
      endif


      
      
!COUPLING Purk-->Ventr
!conta nodi coupling 
      count=00
      nPurkVentr=0
      densCoupBundVentr=1
      do i=1,nftotEF_2d
         if (LabelPurkSetto(i).EQ.0) then !evita setto
            count=count+1
            if (mod(count,densCoupBundVentr).EQ.0) then
               nPurkVentr=nPurkVentr+1
            endif
         endif
      enddo
      write(*,*) "nPurkVentr", nPurkVentr

      allocate(fVentrmaster(nPurkVentr))
      allocate(cVentrslave(nPurkVentr))

      allocate(Purk_tstim(nPurkVentr))
      do i=1,nPurkVentr
         Purk_tstim(i)=-10.d0
      enddo

      count=0
      nPurkVentr1=0
!      densCoup=30
      do i=1,nftotEF_2d
         if (LabelPurkSetto(i).EQ.0) then !evita setto
            count=count+1
            if (mod(count,densCoupBundVentr).EQ.0) then
               nPurkVentr1=nPurkVentr1+1
               fVentrmaster(nPurkVentr1)=i               
            endif
         endif
      enddo
      if (nPurkVentr1.NE.nPurkVentr)then;write(*,*)"err BundVentr";stop;endif

!trova nodi coupling
      offmin0 = 1E5
      csi = cstart_3d(1) ; cei = cend_3d(1)
      do i=1,nPurkVentr
         fmaster=fVentrmaster(i)
         offmin=offmin0
         do j=csi,cei
            chamb=cell_to_chamb_3d(j)
            if ((chamb.EQ.1).OR.(chamb.EQ.3)) then
               off = sqrt((tri_barEF_2d(1,fmaster)-cell_bar(1,j))**2+(tri_barEF_2d(2,fmaster)-cell_bar(2,j))**2+(tri_barEF_2d(3,fmaster)-cell_bar(3,j))**2)
               if (off.LT.offmin) then
                  cslave = j
                  offmin = off
               endif
            endif
         enddo

         cVentrslave(i)=cslave         

!          write(*,*) "vmaster", tri_barEF_2d(1:3,fmaster)
! !         write(*,*) "fslave",  xbarslave,ybarslave,zbarslave
!          write(*,*) "cslate",  cell_bar(1:3,cslave)

      enddo
      

! !TROVA NODI CHECK EF
!       allocate(fcheckEF_2d(2))
!       zmin = 1E5
!       fsi = fstartEF_2d(1) ; fei = fendEF_2d(1)
!       do i=fsi,fei
!          if (tri_barEF_2d(3,i).LT.zmin) then
!             zmin = tri_barEF_2d(3,i)
!             fcheckEF_2d(1)=i
!          endif
!       enddo

!       zmax = -1E5
!       fsi = fstartEF_2d(1) ; fei = fendEF_2d(1)
!       do i=fsi,fei
!          if ((tri_barEF_2d(3,i).GT.zmax).AND.(tri_barEF_2d(1,i).GT.1.5)) then
!             zmax = tri_barEF_2d(3,i)
!             fcheckEF_2d(2)=i
!          endif
!       enddo
!       write(*,*) "checkEF_2d",fcheckEF_2d(1:2)
!       write(*,*) "checkEF_2d",tri_barEF_2d(1:3,fcheckEF_2d(1))
!       write(*,*) "checkEF_2d",tri_barEF_2d(1:3,fcheckEF_2d(2))


!####################################################
!################ PER DETTAGLI ECG ##################
!####################################################
      LabelSkipConn(:)=0
      if (EscCellConn.EQ.1) then    
         do i=1,nBundAtri
            cslave=cAtrislave(i)
            LabelSkipConn(cslave)=1
            do j=1,4 
               f1=face_of_cell_3d(j,cslave)
               c1=cell_of_face_3d(1,f1)
               if (c1.NE.0) LabelSkipConn(c1)=1
               c2=cell_of_face_3d(2,f1)
               if (c2.NE.0) LabelSkipConn(c2)=1
            enddo
         enddo

         do i=1,nPurkVentr 
            cslave=cVentrslave(i)
            LabelSkipConn(cslave)=1
            do j=1,4 
               f1=face_of_cell_3d(j,cslave)
               c1=cell_of_face_3d(1,f1)
               if (c1.NE.0) LabelSkipConn(c1)=1
               c2=cell_of_face_3d(2,f1)
               if (c2.NE.0) LabelSkipConn(c2)=1
            enddo
         enddo
      endif

      EscZonSet=EscZonSet/(LSTAR*1000.d0)
      if (EscZonSet.GE.0.D0) then
         xAp=39.12/(LSTAR*1000.d0)
         yAp=-10.1/(LSTAR*1000.d0)
         zAp=-4.96/(LSTAR*1000.d0)
         xNORM=-0.55 !componenti normale (non necessariamente normalizzate)
         yNORM=0.06
         zNORM=0.83
         sca=sqrt(xNORM**2+yNORM**2+zNORM**2)
         xNORM=xNORM/sca
         yNORM=yNORM/sca
         zNORM=zNORM/sca
         
         xApTr=xAp+EscZonSet*xNORM !coordinate apice traslate (punto sul piano)
         yApTr=yAp+EscZonSet*yNORM
         zApTr=zAp+EscZonSet*zNORM
         
         do i=1,nctot_3d
            chamb=cell_to_chamb_3d(i)
            if ((chamb.EQ.2).OR.(chamb.EQ.4)) then !escludi atri
            else
               xbar=cell_bar(1,i) !coordinate baricentro cella
               ybar=cell_bar(2,i)
               zbar=cell_bar(3,i)
               sca= (xApTr-xbar)*xNORM+(yApTr-ybar)*yNORM+(zApTr-zbar)*zNORM
               if (sca.LT.0.0D0) then !controllo prodotto scalare
                  LabelSkipConn(i)=1
               endif
            endif !escludi atri
         enddo
         

      endif !EscZonSet

         



    
        return
        end subroutine read_geoSingleBodyEF_2d
!===================================================================
        subroutine read_geo_dim_1d(nv,ne,filename)
        
        implicit none
        character*50 :: filename
        integer :: nv,ne

        open(11,file='meshes/'//trim(filename)//'.gts')        
        read(11,*)nv,ne
        close(11)
       
        return
        end subroutine read_geo_dim_1d
!...................................................................
!===================================================================
      subroutine read_geoSingleBodyEF_1d
      use constants
      use param
      use mpih
      use mpi_param, only: kstart, kend
      use mls_param
      implicit none

      character*50 filename,strucfilename
      character*100 ipfi,ipfip

      integer :: i,j,k,v1,v2,v3,inp,e1,e2,e3,count
      integer :: f1,f2,checkF1,checkF2,nBundAV1
      integer :: vsi1,vei1,vsi2,vei2
      real(DP) :: totmass,usVoldom
      real(DP) :: offmin,off
      integer,dimension(:),allocatable :: nv,ne,nf
      real(DP) :: xmax,xmin,ymax,ymin,zmax,zmin
      integer,dimension(3) :: tipsAVnode
      real(DP),dimension(3) :: distAVnode

      integer :: v1d,v2d
      integer :: e1d,e2d,e3d
      integer :: vsi,vei,esi,eei,fsi,fei,csi,cei


      do inp=1,NparticleEF_1d
      write(*,*) trim(geofileEF_1d(inp) )
      open(109,file='meshes/'//trim(geofileEF_1d(inp))//'.gts')
      read(109,*)nviEF_1d(inp),neiEF_1d(inp)
      close(109)
      end do

!     Allocations using vertices,edges and faces as parameters     

      call allocate_trigeoEF_1d

!     Read in the vertices, edges and faces from the file
!     Also make connections between the three
      n_edge_of_vertEF_1d(:)=0


      do inp=1,NparticleEF_1d

         vsi = vstartEF_1d(inp) ; vei = vendEF_1d(inp)
         esi = estartEF_1d(inp) ; eei = eendEF_1d(inp)
         vert_to_partEF_1d(vsi:vei) = inp
         edge_to_partEF_1d(esi:eei) = inp

         write(*,*) nviEF_1d(inp),neiEF_1d(inp)
         open(109,file='meshes/'//trim(geofileEF_1d(inp))//'.gts')
!        open(109,file='meshes/'//geofileEF_1d(inp))
        read(109,*)nviEF_1d(inp),neiEF_1d(inp)
        do i=vsi,vei
          read(109,*)xyz0EF_1d(1,i),xyz0EF_1d(2,i),xyz0EF_1d(3,i)
          xyz0EF_1d(1,i)=xyz0EF_1d(1,i)
          xyz0EF_1d(2,i)=xyz0EF_1d(2,i)
          xyz0EF_1d(3,i)=xyz0EF_1d(3,i)
        end do

        do i=esi,eei
          read(109,*)v1,v2
          v1d = v1 + vsi-1
          v2d = v2 + vsi-1
          vert_of_edgeEF_1d(1,i)=v1d
          vert_of_edgeEF_1d(2,i)=v2d 

          n_edge_of_vertEF_1d(v1d)=n_edge_of_vertEF_1d(v1d)+1
          n_edge_of_vertEF_1d(v2d)=n_edge_of_vertEF_1d(v2d)+1

          edge_of_vertEF_1d(n_edge_of_vertEF_1d(v1d),v1d)=i
          edge_of_vertEF_1d(n_edge_of_vertEF_1d(v2d),v2d)=i
        enddo

        close(109)

      enddo


!     position correction to bring in 0-6 
      xyzEF_1d = xyz0EF_1d
      ! n_edge_of_vertEF_1d(:)=0
      !   do i=esi,eei
      !     read(109,*)v1,v2
      !     v1d = v1 + vsi-1
      !     v2d = v2 + vsi-1
      !     vert_of_edgeEF_1d(1,i)=v1d
      !     vert_of_edgeEF_1d(2,i)=v2d 

      !     n_edge_of_vertEF_1d(v1d)=n_edge_of_vertEF_1d(v1d)+1
      !     n_edge_of_vertEF_1d(v2d)=n_edge_of_vertEF_1d(v2d)+1

      !     edge_of_vertEF_1d(n_edge_of_vertEF_1d(v1d),v1d)=i
      !     edge_of_vertEF_1d(n_edge_of_vertEF_1d(v2d),v2d)=i
      !   enddo



      do inp=1,NparticleEF_1d
         vsi = vstartEF_1d(inp) ; vei = vendEF_1d(inp)
         esi = estartEF_1d(inp) ; eei = eendEF_1d(inp)

         call convert_geoEF_1d(vsi,vei,esi,eei,        &
            xyzEF_1d(:,vsi:vei),            &
            vert_of_edgeEF_1d(:,esi:eei),           &
            edg_barEF_1d(:,esi:eei))

        call calculate_distance(dist0EF_1d(esi:eei),vsi,vei,esi,eei,    &
                       xyz0EF_1d(:,vsi:vei),vert_of_edgeEF_1d(:,esi:eei))


        do i=vsi,vei
           do j=1,n_edge_of_vertEF_1d(i)
              e1=edge_of_vertEF_1d(j,i);
              v1=vert_of_edgeEF_1d(1,e1)
              v2=vert_of_edgeEF_1d(2,e1)
              if ((v1.NE.i).AND.(v2.EQ.i)) then
                 vert_of_vertEF_1d(j,i)=v1
              elseif ((v1.EQ.i).AND.(v2.NE.i)) then
                 vert_of_vertEF_1d(j,i)=v2
              else
                 write(*,*) "ErrorEF_1d"                 
              endif
           enddo
        enddo

      end do
      distEF_1d=dist0EF_1d

!TROVA 3NODI ESTREMITA NODI AV
      nBundAV=0
      vsi = vstartEF_1d(2) ; vei = vendEF_1d(2)
      do i=vsi,vei
         if (n_edge_of_vertEF_1d(i).EQ.1) then
            nBundAV=nBundAV+1        
         endif
      enddo
      if (nBundAV.NE.3) then; write(*,*) "error AV node"; stop; endif
      allocate(vAVmaster(nBundAV),vAVslave(nBundAV))
      nBundAV1=0
      do i=vsi,vei
         if (n_edge_of_vertEF_1d(i).EQ.1) then
            nBundAV1=nBundAV1+1
            tipsAVnode(nBundAV1)=i
         endif
      enddo
      if (nBundAV1.NE.nBundAV) then; write(*,*) "error AV1 node"; stop; endif

!DETERMINA NODO MASTER DELL'AV CONNESSO AI BUNDLES VENTRICOLARI (E I 2SLAVES CONNESSI AI BUNDLES ATRIALI SENZA ORDINE)
      do i=1,3
         v1=tipsAVnode(i)
         offmin = 1E5
         vsi = vstartEF_1d(3) ; vei = vendEF_1d(3)
         do v2=vsi,vei
               off = sqrt((xyzEF_1d(1,v1)-xyzEF_1d(1,v2))**2+(xyzEF_1d(2,v1)-xyzEF_1d(2,v2))**2+(xyzEF_1d(3,v1)-xyzEF_1d(3,v2))**2)
            if (off.LT.offmin) then
               offmin = off
            endif
         enddo
         distAVnode(i)=offmin
      enddo

      if ((distAVnode(1).LT.distAVnode(2)).AND.(distAVnode(1).LT.distAVnode(3))) then
         vAVmaster(3)=tipsAVnode(1)
         vAVslave(1)=tipsAVnode(2)
         vAVslave(2)=tipsAVnode(3)
      elseif ((distAVnode(2).LT.distAVnode(1)).AND.(distAVnode(2).LT.distAVnode(3))) then
         vAVmaster(3)=tipsAVnode(2)
         vAVslave(1)=tipsAVnode(1)
         vAVslave(2)=tipsAVnode(3)
      elseif ((distAVnode(3).LT.distAVnode(1)).AND.(distAVnode(3).LT.distAVnode(2))) then
         vAVmaster(3)=tipsAVnode(3)
         vAVslave(1)=tipsAVnode(1)
         vAVslave(2)=tipsAVnode(2)
      else
         write(*,*) "error AV2 node", tipsAVnode(1:3);stop
      endif

!TROVA MASTERS AND SLAVES MANCANTI
!trova master1
         offmin = 1E5
         v1=vAVslave(1)
         vsi = vstartEF_1d(1) ; vei = vendEF_1d(1)         
         do v2=vsi,vei
               off = sqrt((xyzEF_1d(1,v1)-xyzEF_1d(1,v2))**2+(xyzEF_1d(2,v1)-xyzEF_1d(2,v2))**2+(xyzEF_1d(3,v1)-xyzEF_1d(3,v2))**2)
            if (off.LT.offmin) then
               offmin = off
               vAVmaster(1)=v2
            endif
         enddo
!trova master2
         offmin = 1E5
         v1=vAVslave(2)
         vsi = vstartEF_1d(1) ; vei = vendEF_1d(1)         
         do v2=vsi,vei
               off = sqrt((xyzEF_1d(1,v1)-xyzEF_1d(1,v2))**2+(xyzEF_1d(2,v1)-xyzEF_1d(2,v2))**2+(xyzEF_1d(3,v1)-xyzEF_1d(3,v2))**2)
            if (off.LT.offmin) then
               offmin = off
               vAVmaster(2)=v2
            endif
         enddo
!trova slave3
         offmin = 1E5
         v1=vAVmaster(3)
         vsi = vstartEF_1d(3) ; vei = vendEF_1d(3)         
         do v2=vsi,vei
               off = sqrt((xyzEF_1d(1,v1)-xyzEF_1d(1,v2))**2+(xyzEF_1d(2,v1)-xyzEF_1d(2,v2))**2+(xyzEF_1d(3,v1)-xyzEF_1d(3,v2))**2)
            if (off.LT.offmin) then
               offmin = off
               vAVslave(3)=v2
            endif
         enddo

!PRINT
         do i=1,nBundAV
            v1=vAVmaster(i)
            v2=vAVslave(i)
            ! write(*,*) i, " vAVmaster ", xyzEF_1d(1:3,v1)
            ! write(*,*) i, " vAVslace  ", xyzEF_1d(1:3,v2)
         enddo


        return
        end subroutine read_geoSingleBodyEF_1d
!===================================================================

!===================================================================
      subroutine read_geoSingleBodyCV_1d
      use constants
      use param
      use mpih
      use mpi_param, only: kstart, kend
      use mls_param
      implicit none

      character*50 filename,strucfilename
      character*100 ipfi,ipfip

      integer :: i,j,k,v1,v2,v3,inp,e1,e2,e3,count
      integer :: f1,f2,checkF1,checkF2,nBundAV1
      integer :: vsi1,vei1,vsi2,vei2
      real(DP) :: totmass,usVoldom
      real(DP) :: offmin,off
      integer,dimension(:),allocatable :: nv,ne,nf
      real(DP) :: xmax,xmin,ymax,ymin,zmax,zmin
      integer,dimension(3) :: tipsAVnode
      real(DP),dimension(3) :: distAVnode

      integer :: v1d,v2d
      integer :: e1d,e2d,e3d
      integer :: vsi,vei,esi,eei,fsi,fei,csi,cei

      do inp=1,NparticleCV_1d

      open(109,file='meshes/'//geofileCV_1d(inp))
      read(109,*)nviCV_1d(inp),neiCV_1d(inp)
      close(109)
      end do

!     Allocations using vertices,edges and faces as parameters     

      call allocate_trigeoCV_1d

!     Read in the vertices, edges and faces from the file
!     Also make connections between the three
      ! n_edge_of_vertCV_1d(:)=0


      n_edge_of_vertCV_1d(:)=0.0D0
      do inp=1,NparticleCV_1d

         vsi = vstartCV_1d(inp) ; vei = vendCV_1d(inp)
         esi = estartCV_1d(inp) ; eei = eendCV_1d(inp)
         ! vert_to_partCV_1d(vsi:vei) = inp
         ! edge_to_partCV_1d(esi:eei) = inp

        open(109,file='meshes/'//geofileCV_1d(inp))
        read(109,*)nviCV_1d(inp),neiCV_1d(inp)
        do i=vsi,vei
          read(109,*)xyz0CV_1d(1,i),xyz0CV_1d(2,i),xyz0CV_1d(3,i)
          xyz0CV_1d(1,i)=xyz0CV_1d(1,i)
          xyz0CV_1d(2,i)=xyz0CV_1d(2,i)
          xyz0CV_1d(3,i)=xyz0CV_1d(3,i)
        end do

        do i=esi,eei
          read(109,*)v1,v2
          v1d = v1 + vsi-1
          v2d = v2 + vsi-1
          vert_of_edgeCV_1d(1,i)=v1d
          vert_of_edgeCV_1d(2,i)=v2d 

          n_edge_of_vertCV_1d(v1d)=n_edge_of_vertCV_1d(v1d)+1
          n_edge_of_vertCV_1d(v2d)=n_edge_of_vertCV_1d(v2d)+1

          edge_of_vertCV_1d(n_edge_of_vertCV_1d(v1d),v1d)=i
          edge_of_vertCV_1d(n_edge_of_vertCV_1d(v2d),v2d)=i
        enddo

        close(109)

      enddo


!     position correction to bring in 0-6 
      xyzCV_1d = xyz0CV_1d

      do inp=1,NparticleCV_1d
         vsi = vstartCV_1d(inp) ; vei = vendCV_1d(inp)
         esi = estartCV_1d(inp) ; eei = eendCV_1d(inp)

         ! call convert_geoEF_1d(vsi,vei,esi,eei,        &
         !    xyzCV_1d(:,vsi:vei),            &
         !    vert_of_edgeECV_1d(:,esi:eei),           &
         !    edg_barCV_1d(:,esi:eei))

      !   call calculate_distance(dist0EF_1d(esi:eei),vsi,vei,esi,eei,    &
      !                  xyz0EF_1d(:,vsi:vei),vert_of_edgeEF_1d(:,esi:eei))


        do i=vsi,vei
           do j=1,n_edge_of_vertCV_1d(i)
              e1=edge_of_vertCV_1d(j,i);
              v1=vert_of_edgeCV_1d(1,e1)
              v2=vert_of_edgeCV_1d(2,e1)
              if ((v1.NE.i).AND.(v2.EQ.i)) then
                 vert_of_vertCV_1d(j,i)=v1
              elseif ((v1.EQ.i).AND.(v2.NE.i)) then
                 vert_of_vertCV_1d(j,i)=v2
              else
                 write(*,*) "ErrorCV_1d"                 
              endif
           enddo
        enddo

      end do
      ! distEF_1d=dist0EF_1d



        return
        end subroutine read_geoSingleBodyCV_1d
!===================================================================



!===================================================================
        subroutine find_tag_offsetEF_2d
        use constants
        use mls_param

        implicit none
        integer :: i, j,v1, v2,inp,vsi,vei,qualeoff
        real(DP) pig, Rsec, zsecmin, off0, offmin,off
        pig=acos(-1.0)


        off0 = 100000.

        inp = 1
        vsi = vstart(inp) 
        vei = vend(inp)

        do i=1,nvtotEF_2d
           offmin=off0
           do j=vsi,vei
              off = sqrt((xyz0EF_2d(1,i)-xyz0(1,j))**2+(xyz0EF_2d(2,i)-xyz0(2,j))**2+(xyz0EF_2d(3,i)-xyz0(3,j))**2)
              if (off.LT.offmin) then                 
                 qualeoff = j
                 offmin = off
              endif
           enddo
           boundary2EF_2d(i) = qualeoff
           BCoffsetEF_2d(1:3,i) = xyz0EF_2d(1:3,i)-xyz0(1:3,qualeoff)              
       enddo

       

        ! do i=1,nvtotEF_2d
        !    vsi=boundary2EF_2d(i) 
        !    write(*,*) "A",xyz0EF_2d(1:3,i)
        !    write(*,*) "B",xyz0(1:3,vsi)
        ! enddo

        end subroutine find_tag_offsetEF_2d
!===================================================================  
!===================================================================
        subroutine find_tag_offsetEF_1d
        use constants
        use mls_param

        implicit none
        integer :: i, j,v1, v2,inp,vsi,vei,qualeoff
        real(DP) pig, Rsec, zsecmin, off0, offmin,off
        pig=acos(-1.0)


        off0 = 100000.0

        inp = 1
        vsi = vstart(inp) 
        vei = vend(inp)

        do i=1,nvtotEF_1d
           offmin=off0
           do j=vsi,vei
              off = sqrt((xyz0EF_1d(1,i)-xyz0(1,j))**2+(xyz0EF_1d(2,i)-xyz0(2,j))**2+(xyz0EF_1d(3,i)-xyz0(3,j))**2)
              if (off.LT.offmin) then                 
                 qualeoff = j
                 offmin = off
              endif
           enddo
           boundary2EF_1d(i) = qualeoff
           BCoffsetEF_1d(1:3,i) = xyz0EF_1d(1:3,i)-xyz0(1:3,qualeoff)              
       enddo

       

        ! do i=1,nvtotEF_1d
        !    vsi=boundary2EF_1d(i) 
        !    write(*,*) "A",xyz0EF_1d(1:3,i)
        !    write(*,*) "B",xyz0(1:3,vsi)
        ! enddo

        end subroutine find_tag_offsetEF_1d
!===================================================================  
!===================================================================
        subroutine find_tag_offsetCV_1d
        use constants
        use mls_param

        implicit none
        integer :: i, j,v1, v2,inp,vsi,vei,qualeoff,fsi,fei,k,c1,c2,v3
        real(DP) pig, Rsec, zsecmin, off0, offmin,off
        real(DP),dimension(nvtotCV_1d) ::vett_offmin
        pig=acos(-1.0)


        off0 = 100000.0

        inp = 1
        vsi = vstart_3d(inp) 
        vei = vend_3d(inp)
        fsi = fstart_3d(inp) 
        fei = fend_3d(inp)

       !  do i=1,nvtotCV_1d
       !     offmin=off0
       !     do j=vsi,vei
       !        c1=

       !        off = sqrt((xyz0CV_1d(1,i)-xyz0_3d(1,j))**2+(xyz0CV_1d(2,i)-xyz0_3d(2,j))**2+(xyz0CV_1d(3,i)-xyz0_3d(3,j))**2)
       !        if (off.LT.offmin) then                 
       !           qualeoff = j
       !           offmin = off
       !        endif
       !     enddo
       !     boundary2CV_1d(i) = qualeoff
       !     BCoffsetCV_1d(1:3,i) = xyz0CV_1d(1:3,i)-xyz0_3d(1:3,qualeoff)              
       ! enddo

        ! do i=1,nvtotCV_1d
        !    offmin=off0
        !    do k=fsi,fei
        !       c1=cell_of_face_3d(1,k)
        !       c2=cell_of_face_3d(1,k)

              
        vett_offmin(:) = off0
        do j=fsi,fei
           c1=cell_of_face_3d(1,j)
           c2=cell_of_face_3d(2,j)
           
           if ((c1.EQ.0).OR.(c2.EQ.0)) then
              v1=vert_of_face_3d(1,j)
              v2=vert_of_face_3d(2,j)
              v3=vert_of_face_3d(3,j)
              do i=1,nvtotCV_1d                 
                 !v1
                 off = sqrt((xyz0CV_1d(1,i)-xyz0_3d(1,v1))**2+(xyz0CV_1d(2,i)-xyz0_3d(2,v1))**2+(xyz0CV_1d(3,i)-xyz0_3d(3,v1))**2)
                 if (off.LT.vett_offmin(i)) then
                    boundary2CV_1d(i) = v1
                    BCoffsetCV_1d(1:3,i) = xyz0CV_1d(1:3,i)-xyz0_3d(1:3,v1)
                    vett_offmin(i)=off
                 endif
                 !v2
                 off = sqrt((xyz0CV_1d(1,i)-xyz0_3d(1,v2))**2+(xyz0CV_1d(2,i)-xyz0_3d(2,v2))**2+(xyz0CV_1d(3,i)-xyz0_3d(3,v2))**2)
                 if (off.LT.vett_offmin(i)) then
                    boundary2CV_1d(i) = v2
                    BCoffsetCV_1d(1:3,i) = xyz0CV_1d(1:3,i)-xyz0_3d(1:3,v2)
                    vett_offmin(i)=off
                 endif
                 !v3
                 off = sqrt((xyz0CV_1d(1,i)-xyz0_3d(1,v3))**2+(xyz0CV_1d(2,i)-xyz0_3d(2,v3))**2+(xyz0CV_1d(3,i)-xyz0_3d(3,v3))**2)
                 if (off.LT.vett_offmin(i)) then
                    boundary2CV_1d(i) = v3
                    BCoffsetCV_1d(1:3,i) = xyz0CV_1d(1:3,i)-xyz0_3d(1:3,v3)
                    vett_offmin(i)=off
                 endif
              enddo
           endif !if c1 c2

        enddo !loop facce 3d

        end subroutine find_tag_offsetCV_1d
!===================================================================  
