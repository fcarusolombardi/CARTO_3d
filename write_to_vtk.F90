!-------------------------------------------------
!     write structure data to vtk file
!-------------------------------------------------

      subroutine write_to_vtk(ipfi,optionvtk)
      use mpi_param
      use mls_param
      use mpih
      use param

      implicit none
     
      integer :: i,inp,itime
      integer :: tnv, tnf,optionvtk
      integer, dimension(nftot) :: cell_type      
      real(DP), dimension(maxnf) :: scalar_face
      integer, dimension(3,nftot) :: vert_face_dum
      
      real(DP) :: tprfi
      character*70 filname
      character*100 ipfi,ipfip

      cell_type(:) = 5 !5 is for triangular elements
      

      do inp=1,Nparticle
       scalar_face(1:maxnf) = 0.0
       vert_face_dum(1:3,fstart(inp):fend(inp))=vert_of_face(1:3,fstart(inp):fend(inp))-vstart(inp)
      end do

do inp=1,Nparticle
      tnv = vend(inp) - vstart(inp) + 1
      tnf = fend(inp) - fstart(inp) + 1
      write(ipfip,94)inp
   94 format(i2.2)

      filname = 'vtkfiles/struc_'//trim(ipfip)//'_'//trim(ipfi)//'.vtk'

      open(121,file = filname)
!Header     
        write(121,'(a)')adjustl('# vtk DataFile Version 3.1')
        write(121,'(a)')adjustl('Stores triangular gts mesh')
        write(121,'(a)')adjustl('ASCII')
        write(121,'(a)')adjustl('DATASET UNSTRUCTURED_GRID')
        write(121,*)''
        write(121,*)'POINTS ',tnv,' FLOAT'
      do i=vstart(inp),vend(inp)
        write(121,*)xyz(1:3,i)
      end do
        write(121,*)''
        write(121,*)'CELLS ',tnf, 4*tnf
      do i=fstart(inp),fend(inp)
        write(121,*)'3 ',vert_face_dum(1:3,i)
      end do
        write(121,*)''
        write(121,*)'CELL_TYPES ',tnf
        write(121,*)cell_type(fstart(inp):fend(inp))
! 
!        write(121,*)''
!        write(121,*)'CELL_DATA ',fend
!        write(121,*)scalar_face(:)

        if (optionvtk.EQ.0) then
           write(121,*)'' 
           write(121,*)'POINT_DATA ',nvi(inp)
           write(121,*)'Scalars vertchamb FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart(inp),vend(inp)
              write(121,*) vert_to_chamb(i)
           end do

           write(121,*)''
           write(121,*)'CELL_DATA ',nfi(inp)
           write(121,*)'Scalars facechamb FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=fstart(inp),fend(inp)
              write(121,*) face_to_chamb(i)
           end do

           write(121,*)'' 
           write(121,*)'POINT_DATA ',nvi(inp)
!           write(121,*)'Scalars LabelBound FLOAT'
           write(121,*)'Scalars vertchamb4V FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart(inp),vend(inp)
              write(121,*) vert_to_chamb4V(i)
           enddo

           write(121,*)''
           write(121,*)'CELL_DATA ',nfi(inp)
           write(121,*)'Scalars facechamb4V FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=fstart(inp),fend(inp)
              write(121,*) face_to_chamb4V(i)
           end do


           if (inp.EQ.2) then
           write(121,*)''
           write(121,*)'POINT_DATA ',nvi(inp)
           write(121,*)'Scalars LabelLeaflet FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart(inp),vend(inp)
              write(121,*) LabelVaortic(i)
           enddo
           endif

           if (inp.EQ.3) then
           write(121,*)''
           write(121,*)'POINT_DATA ',nvi(inp)
           write(121,*)'Scalars LabelLeaflet FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart(inp),vend(inp)
              write(121,*) LabelVmitral(i)
           enddo
           endif

           if (inp.EQ.4) then
           write(121,*)''
           write(121,*)'POINT_DATA ',nvi(inp)
           write(121,*)'Scalars LabelLeaflet FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart(inp),vend(inp)
              write(121,*) LabelVpulmo(i)
           enddo
           endif

           if (inp.EQ.5) then
           write(121,*)''
           write(121,*)'POINT_DATA ',nvi(inp)
           write(121,*)'Scalars LabelLeaflet FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart(inp),vend(inp)
              write(121,*) LabelVtricu(i)
           enddo
           endif
        
           
        elseif (optionvtk.EQ.1) then
           ! write(121,*)'' 
           ! write(121,*)'POINT_DATA ',nvi(inp)
           ! write(121,*)'Scalars LabelBound FLOAT'
           ! write(121,*)'LOOKUP_TABLE default'
           ! do i=vstart(inp),vend(inp)
           !    write(121,*) LabelBound(i)
           ! enddo


           write(121,*)''
           write(121,*)'CELL_DATA ',nfi(inp)
           write(121,*)'Scalars tauface FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=fstart(inp),fend(inp)
              write(121,*) sqrt(tauface(i)) !sqrt here
           end do

           write(121,*)''
!           write(121,*)'CELL_DATA ',nfi(inp)
           write(121,*)'Scalars pressface FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=fstart(inp),fend(inp)
              write(121,*) sqrt(pressface(i)) !sqrt here
           end do

           write(121,*)''
!           write(121,*)'CELL_DATA ',nfi(inp)
           write(121,*)'Scalars taufaceAV FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=fstart(inp),fend(inp)
              write(121,*) taufaceAV(i)
              taufaceAV(i) = 0.0d0 !put it to zero for new average
           end do

           write(121,*)''
!           write(121,*)'CELL_DATA ',nfi(inp)
           write(121,*)'Scalars pressfaceAV FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=fstart(inp),fend(inp)
              write(121,*) pressfaceAV(i)
              pressfaceAV(i) = 0.0d0 !put it to zero for new average
           end do


        endif !optionvtk
      
      close(121)
     
end do

!stop
      return
      end
      

!-------------------------------------------------
!     write structure data 3d to vtk file
!-------------------------------------------------

      subroutine write_to_vtk_3d(ipfi,optionvtk)
      use mpi_param
      use mls_param
      use mpih
      use param

      implicit none
     
      integer :: i,inp,itime,chamb,j,ie
      integer :: tnv, tnc,optionvtk,cstart,cend
      integer, dimension(nctot_3d) :: cell_type      
!      real(DP), dimension(maxnf) :: scalar_face
      integer, dimension(4,nctot_3d) :: vert_cell_dum
      real(DP) :: x_inf_loc,y_inf_loc,z_inf_loc,size_inf
      real(DP) :: campo_inf,num,den,d,d0,kei,epsG
      real(DP) :: tprfi,stimleads
      character*70 filname
      character*100 ipfi,ipfip


      

104 format((2x,f12.8))

      cell_type(:) = 10 !10 is for tretrahedral elements
      

      do inp=1,Nparticle_3d
!       scalar_face(1:maxnf) = 0.0
       vert_cell_dum(1:4,cstart_3d(inp):cend_3d(inp))=vert_of_cell_3d(1:4,cstart_3d(inp):cend_3d(inp))-vstart_3d(inp)
      end do

      do inp=1,Nparticle_3d
      tnv = vend_3d(inp) - vstart_3d(inp) + 1
      tnc = cend_3d(inp) - cstart_3d(inp) + 1
      write(ipfip,94)inp
   94 format(i2.2)


      filname = 'vtkfiles/struc3d_'//trim(ipfip)//'_'//trim(ipfi)//'.vtk'

      open(121,file = filname)
!Header     
        write(121,'(a)')adjustl('# vtk DataFile Version 3.1')
        write(121,'(a)')adjustl('Stores tetrahdral mesh')
        write(121,'(a)')adjustl('ASCII')
        write(121,'(a)')adjustl('DATASET UNSTRUCTURED_GRID')
!        write(121,*)''
        write(121,*)'POINTS ',tnv,' FLOAT'
      do i=vstart_3d(inp),vend_3d(inp)
        write(121,*)xyz_3d(1:3,i)
      end do
        write(121,*)''
        write(121,*)'CELLS ',tnc, 5*tnc
      do i=cstart_3d(inp),cend_3d(inp)
        write(121,*)'4 ',vert_cell_dum(1:4,i)
      end do
        write(121,*)''
        write(121,*)'CELL_TYPES ',tnc
        write(121,*)cell_type(cstart_3d(inp):cend_3d(inp))
! 
!        write(121,*)''
!        write(121,*)'CELL_DATA ',fend
!        write(121,*)scalar_face(:)

        if (optionvtk.EQ.0) then

           if (inp.EQ.1) then
           !questo per dati vertici FV
           write(121,*)'' 
           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars vertchamb FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) vert_to_chamb_3d(i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars vert_to_chamb_3dBou FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) vert_to_chamb_3dBou(i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars LabelSmooth FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,104) LabelSmooth(i)
           end do

!            write(121,*)'' 
! !           write(121,*)'POINT_DATA ',nvi_3d(inp)
!            write(121,*)'Scalars LabelDelay FLOAT'
!            write(121,*)'LOOKUP_TABLE default'
!            do i=vstart_3d(inp),vend_3d(inp)
!               write(121,104) LabelDelay(i)
!            end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars vertchamb4V FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) vert_to_chamb_3d4V(i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars fiber_x FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) AmatrFibers_node_3d(1,1,i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars fiber_y FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) AmatrFibers_node_3d(2,1,i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars fiber_z FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) AmatrFibers_node_3d(3,1,i)
           end do
           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars sheetfiber_x FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) AmatrFibers_node_3d(1,2,i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars sheetfiber_y FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) AmatrFibers_node_3d(2,2,i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars sheetfiber_z FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) AmatrFibers_node_3d(3,2,i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars crossfiber_x FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) AmatrFibers_node_3d(1,3,i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars crossfiber_y FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) AmatrFibers_node_3d(2,3,i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars crossfiber_z FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) AmatrFibers_node_3d(3,3,i)
           end do


           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars Segments FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) Segments_node_3d(i)
           end do

           if (minf.GE.1) then 
           size_inf=xyz_inf(4)
           x_inf_loc=xyz_inf(1)
           y_inf_loc=xyz_inf(2)
           z_inf_loc=xyz_inf(3)           
           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars Inf FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              campo_inf=1.0D0 - exp( -( ( (xyz0_3d(1,i)-x_inf_loc)**2+(xyz0_3d(2,i)-y_inf_loc)**2+(xyz0_3d(3,i)-z_inf_loc)**2)**4/size_inf**8 ) )
              write(121,*) campo_inf
           end do
           endif

           write(121,*)''
           write(121,*)'CELL_DATA ',nci_3d(inp)
           write(121,*)'Scalars chambersCell FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) cell_to_chamb_3d(i)
           end do

           write(121,*)''
!           write(121,*)'CELL_DATA ',nci_3d(inp)                                                                                                          
           write(121,*)'Scalars chambersCellBou FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) cell_to_chamb_3dBou(i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars fiber_x FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) AmatrFibers_cell_3d(1,1,i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars fiber_y FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) AmatrFibers_cell_3d(2,1,i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars fiber_z FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) AmatrFibers_cell_3d(3,1,i)
           end do
           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars sheetfiber_x FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) AmatrFibers_cell_3d(1,2,i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars sheetfiber_y FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) AmatrFibers_cell_3d(2,2,i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars sheetfiber_z FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) AmatrFibers_cell_3d(3,2,i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars crossfiber_x FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) AmatrFibers_cell_3d(1,3,i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars crossfiber_y FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) AmatrFibers_cell_3d(2,3,i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars crossfiber_z FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) AmatrFibers_cell_3d(3,3,i)
           end do


           write(121,*)'' 
!           write(121,*)'CELL_DATA ',nci_3d(inp)
           write(121,*)'Scalars LabelSkipConn  FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) LabelSkipConn(i)
           end do
           
           write(121,*)'' 
           !write(121,*)'CELL_DATA ',nci_3d(inp)
           write(121,*)'Scalars Istim FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) IstimEF_3dS1(i)
           end do
           
        endif


        elseif (optionvtk.EQ.1) then

!            write(121,*)'' 
!            write(121,*)'POINT_DATA ',nvi_3d(inp)
!            write(121,*)'Scalars LabelSmooth FLOAT'
!            write(121,*)'LOOKUP_TABLE default'
!            do i=vstart_3d(inp),vend_3d(inp)
!               write(121,*) LabelSmooth(i)
!            end do


!            write(121,*)'' 
! !           write(121,*)'POINT_DATA ',nvi_3d(inp)
!            write(121,*)'Scalars myocardium FLOAT'
!            write(121,*)'LOOKUP_TABLE default'
!            do i=vstart_3d(inp),vend_3d(inp)
!               chamb=vert_to_chamb_3d4V(i)
!               if (chamb.LE.4) then
!                  write(121,*) 1.0D0
!               else
!                  write(121,*) 0.0D0
!               endif
!            end do

#ifdef ELECTRO
           write(121,*)'' 
          write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars LAT_node FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) LATnode_3d(i)
           end do

           write(121,*)'' 
          !write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars vtrans FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) potEFnode_3d(i)
           end do

           write(121,*)'' 
           write(121,*)'Scalars CVx FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) CVnode(1,i)
           end do
           write(121,*)'' 
           write(121,*)'Scalars CVy FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) CVnode(2,i)
           end do
           write(121,*)'' 
           write(121,*)'Scalars CVz FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) CVnode(3,i)
           end do
           
           write(121,*)'' 
           write(121,*)'Scalars CV_div FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) CVdiv(i)
           end do

           write(121,*)'' 
           write(121,*)'Scalars CV_rotx FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) CVrot(1,i)
           end do

           write(121,*)'' 
           write(121,*)'Scalars CV_roty FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) CVrot(2,i)
           end do

           write(121,*)'' 
           write(121,*)'Scalars CV_rotz FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) CVrot(3,i)
           end do
           
           write(121,*)''
           write(121,*)'Scalars vtrans_dt FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) XEF_3ddt(i)
           end do
           
           write(121,*)'' 
          ! write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars vert2chamb FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) vert_to_chamb_3d(i)
           end do

          !  write(121,*)'' 
          ! ! write(121,*)'POINT_DATA ',nvi_3d(inp)
          !  write(121,*)'Scalars atens FLOAT'
          !  write(121,*)'LOOKUP_TABLE default'
          !  do i=vstart_3d(inp),vend_3d(inp)
          !     write(121,*) astressEFnode_3d(i)
          !  end do

           write(121,*)'' 
           write(121,*)'CELL_DATA ',nci_3d(inp)
           ! !write(121,*)'Scalars myotag FLOAT'
           ! write(121,*)'Scalars label_smooth FLOAT'
           ! write(121,*)'LOOKUP_TABLE default'
           ! do i=cstart_3d(inp),cend_3d(inp)
           !    write(121,*) LabelSmooth_cell(i)
           ! end do

           ! write(121,*)'' 
           !write(121,*)'Scalars myotag FLOAT'
           write(121,*)'Scalars CARTO_Dcell FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) CARTO_Dcell3d(i)
           end do
           write(121,*)'' 
           write(121,*)'Scalars CVcellx FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) CVcell(1,i)
           end do
           write(121,*)'' 
           write(121,*)'Scalars CVcelly FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) CVcell(2,i)
           end do
           write(121,*)'' 
           write(121,*)'Scalars CVcellz FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) CVcell(3,i)
           end do
           ! write(121,*)'' 
           ! !write(121,*)'Scalars myotag FLOAT'
           ! write(121,*)'Scalars grad_y FLOAT'
           ! write(121,*)'LOOKUP_TABLE default'
           ! do i=cstart_3d(inp),cend_3d(inp)
           !    write(121,*) gradcell_3d(2,i)
           ! end do
           
           ! write(121,*)'' 
           ! !write(121,*)'Scalars myotag FLOAT'
           ! write(121,*)'Scalars pot_cell FLOAT'
           ! write(121,*)'LOOKUP_TABLE default'
           ! do i=cstart_3d(inp),cend_3d(inp)
           !    write(121,*) potEFcell_3d(i)
           ! end do
           
           ! write(121,*)'' 
           ! !write(121,*)'CELL_DATA ',nci_3d(inp)
           ! write(121,*)'Scalars c2chamb FLOAT'
           ! write(121,*)'LOOKUP_TABLE default'
           ! do i=cstart_3d(inp),cend_3d(inp)
           !    write(121,*) cell_to_chamb_3d(i)
           ! end do           
           write(121,*)'' 
           !write(121,*)'CELL_DATA ',nci_3d(inp)
           write(121,*)'Scalars APD FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) t_apd_3d(2,i)-t_apd_3d(1,i)
           end do
 
           write(121,*)'' 
           !write(121,*)'CELL_DATA ',nci_3d(inp)
           write(121,*)'Scalars LAT FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) t_apd_3d(1,i)!-t_apd_3d(2,i)
           end do

           ! write(121,*)'' 
           ! !write(121,*)'CELL_DATA ',nci_3d(inp)
           ! write(121,*)'Scalars APD_S2 FLOAT'
           ! write(121,*)'LOOKUP_TABLE default'
           ! do i=cstart_3d(inp),cend_3d(inp)
           !    write(121,*) t_apd_3d(4,i)-t_apd_3d(3,i)
           ! end do

           ! write(121,*)'' 
           ! !write(121,*)'CELL_DATA ',nci_3d(inp)
           ! write(121,*)'Scalars DI FLOAT'
           ! write(121,*)'LOOKUP_TABLE default'
           ! do i=cstart_3d(inp),cend_3d(inp)
           !    write(121,*) t_apd_3d(3,i)-t_apd_3d(2,i)
           ! end do

           
           write(121,*)'' 
           !write(121,*)'CELL_DATA ',nci_3d(inp)
           write(121,*)'Scalars Stenosi FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=cstart_3d(inp),cend_3d(inp)
              write(121,*) LabelStenosi(i)
           end do

           ! write(121,*)'' 
           ! !write(121,*)'CELL_DATA ',nci_3d(inp)
           ! write(121,*)'Scalars mesh_skewness FLOAT'
           ! write(121,*)'LOOKUP_TABLE default'
           ! do i=cstart_3d(inp),cend_3d(inp)
           !    write(121,*) meshquality_3d(i)
           ! end do

           ! write(121,*)'' 
           ! !write(121,*)'CELL_DATA ',nci_3d(inp)
           ! write(121,*)'Scalars cell_volume FLOAT'
           ! write(121,*)'LOOKUP_TABLE default'
           ! do i=cstart_3d(inp),cend_3d(inp)
           !    write(121,*) vol0_3d(i)
           ! end do

           ! write(121,*)'' 
           ! !write(121,*)'CELL_DATA ',nci_3d(inp)
           ! write(121,*)'Scalars Istim FLOAT'
           ! write(121,*)'LOOKUP_TABLE default'
           ! do i=cstart_3d(inp),cend_3d(inp)
           !    write(121,*) IstimEF_3dS1(i)
           ! end do
           
           
!            write(121,*)'' 
!            write(121,*)'CELL_DATA ',nci_3d(inp)
!            write(121,*)'Scalars Istim FLOAT'
!            write(121,*)'LOOKUP_TABLE default'
!            do i=cstart_3d(inp),cend_3d(inp)
!               stimleads = IstimEF_3d(i) 
! !              stimleads=IstimRV_3d(i)
!               write(121,*) stimleads
!            end do


#endif

if (Ppro.GT.0) then

   if (ntst.GT.0) then !passice stress
      write(121,*)'' 
      ! write(121,*)'POINT_DATA ',nvi_3d(inp)
      write(121,*)'Scalars pstress FLOAT'
      write(121,*)'LOOKUP_TABLE default'      
      do i = vstart_3d(inp),vend_3d(inp)
         num=0.D0
         den=0.D0
         do j = 1, n_edge_of_vert_3d(i)
            ie = edge_of_vert_3d(j,i)
            d = dist_3d(ie)
            d0 = dist0_3d(ie)
            kei = ke_3d(ie)
            epsG = (d-d0) / d
            ! !stress for visualization
            num = num + kei*epsG/(d/2)
            den = den + 1.D0/(d/2)
         enddo
!         stress_node_3d(i)=num/den
         write(121,*) num/den
      enddo
   endif
      
           

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars fiber_x FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) AmatrFibers_node_3d(1,1,i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars fiber_y FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) AmatrFibers_node_3d(2,1,i)
           end do

           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars fiber_z FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              write(121,*) AmatrFibers_node_3d(3,1,i)
           end do

           if (minf.GE.1) then 
           size_inf=xyz_inf(4)
           x_inf_loc=xyz_inf(1)
           y_inf_loc=xyz_inf(2)
           z_inf_loc=xyz_inf(3)           
           write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
           write(121,*)'Scalars Inf FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstart_3d(inp),vend_3d(inp)
              campo_inf=1.0D0 - exp( -( ( (xyz0_3d(1,i)-x_inf_loc)**2+(xyz0_3d(2,i)-y_inf_loc)**2+(xyz0_3d(3,i)-z_inf_loc)**2)**4/size_inf**8 ) )
              write(121,*) campo_inf
           end do
        endif !minf

        write(121,*)'' 
!           write(121,*)'POINT_DATA ',nvi_3d(inp)
        write(121,*)'Scalars vertchamb4V FLOAT'
        write(121,*)'LOOKUP_TABLE default'
        do i=vstart_3d(inp),vend_3d(inp)
           write(121,*) vert_to_chamb_3d4V(i)
        end do

endif !Ppro


      endif !optionvtk
      close(121)
      
      end do 
      return
      end


!-------------------------------------------------
!     write field data to vtk file
!-------------------------------------------------

      subroutine write_to_vtk_field(ipfi)
      use param
      use mpi_param
      use mpih
      use local_arrays
      implicit none
     
      integer :: i,nslice
      character*150 filname
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
        write(121,*)'Scalars q1 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
!        write(121,*)q1(:,n2/2 + 1,kstart:kend)
        write(121,*)q1(:,nslice,kstart:kend)
        write(121,*)'Scalars q2 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
!        write(121,*)q2(:,n2/2 + 1,kstart:kend)
        write(121,*)q2(:,nslice,kstart:kend)
        write(121,*)'Scalars q3 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
        ! write(121,*)q3(:,n2/2 + 1,kstart:kend)
        write(121,*)q3(:,nslice,kstart:kend)
        write(121,*)'Scalars pr FLOAT'
        write(121,*)'LOOKUP_TABLE default'
        ! write(121,*)pr(:,n2/2 + 1,kstart:kend)
        write(121,*)pr(:,nslice,kstart:kend)

!        write(121,*)pr(:,n2/2 + 1,kstart:kend)

        ! write(121,*)''
        ! write(121,*)'POINT_DATA ', n1*1*(kend-kstart+1)
        ! write(121,*)'Scalars pr FLOAT'
        ! write(121,*)'LOOKUP_TABLE default'
        ! write(121,*)pr(:,n2/2 + 1,kstart:kend)

      close(121)

!Another one
!trova slice                                                                                                                                                                 
      sliceCoord = 1.5180D0
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
      !filname = 'vtkfiles/field_02_'//trim(ipfi)//'.vtk'
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
        write(121,*)'Scalars q1 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
!        write(121,*)q1(:,n2/2 + 1,kstart:kend)
        write(121,*)q1(:,nslice,kstart:kend)
        write(121,*)'Scalars q2 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
!        write(121,*)q2(:,n2/2 + 1,kstart:kend)
        write(121,*)q2(:,nslice,kstart:kend)
        write(121,*)'Scalars q3 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
        ! write(121,*)q3(:,n2/2 + 1,kstart:kend)
        write(121,*)q3(:,nslice,kstart:kend)
        write(121,*)'Scalars pr FLOAT'
        write(121,*)'LOOKUP_TABLE default'
        ! write(121,*)pr(:,n2/2 + 1,kstart:kend)
        write(121,*)pr(:,nslice,kstart:kend)

!        write(121,*)pr(:,n2/2 + 1,kstart:kend)

        ! write(121,*)''
        ! write(121,*)'POINT_DATA ', n1*1*(kend-kstart+1)
        ! write(121,*)'Scalars pr FLOAT'
        ! write(121,*)'LOOKUP_TABLE default'
        ! write(121,*)pr(:,n2/2 + 1,kstart:kend)

      close(121)


!Another one
!trova slice                                                                                                                                                                 
      sliceCoord = 0.00D0
      scarto = 1000.
      do i=1,n1
         if (abs(xc(i)-sliceCoord).LT.scarto) then
            scarto=abs(xc(i)-sliceCoord)
            nslice = i
         endif
      enddo

!#ifdef USE_CUDA
!      q1(:,:,n3/2+1) = q1_d(:,:,n3/2+1)
!      q2(:,:,n3/2+1) = q2_d(:,:,n3/2+1)
!      q3(:,:,n3/2+1) = q3_d(:,:,n3/2+1)
!#endif

      write(id, '(i2.2)') myid
!      filname = 'vtkfiles/field_03_'//trim(ipfi)//'.vtk'
      filname = 'vtkfiles/field_03_'//id//'_'//trim(ipfi)//'.vtk'
      open(121,file = filname)
!Header     
        write(121,'(a)')adjustl('# vtk DataFile Version 3.1')
        write(121,'(a)')adjustl('Stores full field data')
        write(121,'(a)')adjustl('ASCII')
        write(121,'(a)')adjustl('DATASET RECTILINEAR_GRID')
        write(121,'(a, 3I)')'DIMENSIONS ',1, n2, kend-kstart+1
        write(121,*)'X_COORDINATES ',1, ' FLOAT'
        write(121,*) xc(nslice)
        write(121,*)'Y_COORDINATES ',n2, ' FLOAT'
!        write(121,*) yc(n2/2 + 1)
        write(121,*) yc
        write(121,*)'Z_COORDINATES ',kend-kstart+1, ' FLOAT'
        write(121,*) zc(kstart:kend)
        write(121,*)''
        
        write(121,*)''
        write(121,*)'POINT_DATA ', 1*n2*(kend-kstart+1)
        write(121,*)'Scalars q1 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
!        write(121,*)q1(:,n2/2 + 1,kstart:kend)
        write(121,*)q1(nslice,:,kstart:kend)
        write(121,*)'Scalars q2 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
!        write(121,*)q2(:,n2/2 + 1,kstart:kend)
        write(121,*)q2(nslice,:,kstart:kend)
        write(121,*)'Scalars q3 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
        ! write(121,*)q3(:,n2/2 + 1,kstart:kend)
        write(121,*)q3(nslice,:,kstart:kend)
        write(121,*)'Scalars pr FLOAT'
        write(121,*)'LOOKUP_TABLE default'
        ! write(121,*)pr(:,n2/2 + 1,kstart:kend)
        write(121,*)pr(nslice,:,kstart:kend)

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
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      subroutine write_to_vtkEF_2d(ipfi,optionvtk)
      use mpi_param
      use mls_param
      use mpih
      use param
!@cuf   use cudafor 
      implicit none
     
      integer :: i,inp,itime,j,v1,v2,v3,vv,fslave
      integer :: tnv, tnf,optionvtk,kk
      integer, dimension(nftotEF_2d) :: cell_type      
      real(DP), dimension(maxnf) :: scalar_face,stimleads
      integer, dimension(3,nftotEF_2d) :: vert_face_dum
      
      real(DP) :: tprfi,vettox,vettoy,vettoz,smufa
      character*70 filname
      character*100 ipfi,ipfip
!@cuf   integer :: istat 


!move EF_2d
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do i = 1,nvtotEF_2d
        v2 = boundary2EF_2d(i)
           xyzEF_2d(1,i) = xyz(1,v2)+BCoffsetEF_2d(1,i)   
           xyzEF_2d(2,i) = xyz(2,v2)+BCoffsetEF_2d(2,i)   
           xyzEF_2d(3,i) = xyz(3,v2)+BCoffsetEF_2d(3,i)   
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP



!smooth EF_2d
#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                                                                                                                     
#endif
      do i = 1,nvtotEF_2d
           xyzSEF_2d(1,i) = xyzEF_2d(1,i)
           xyzSEF_2d(2,i) = xyzEF_2d(2,i)
           xyzSEF_2d(3,i) = xyzEF_2d(3,i)
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP   

      smufa=0.01D0
      do kk=1,1

!smooth bulk
#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                                                                                                                     
#endif
         do i=1,nvtotEF_2d
            vettox = 0.0D0
            vettoy = 0.0D0
            vettoz = 0.0D0
            do j=1,n_edge_of_vertEF_2d(i)
               v1=vert_of_vertEF_2d(j,i);
               vettox=vettox+xyzEF_2d(1,v1);        
               vettoy=vettoy+xyzEF_2d(2,v1);        
               vettoz=vettoz+xyzEF_2d(3,v1);        
            enddo
            vettox=vettox-n_edge_of_vertEF_2d(i)*xyzEF_2d(1,i);
            vettoy=vettoy-n_edge_of_vertEF_2d(i)*xyzEF_2d(2,i);
            vettoz=vettoz-n_edge_of_vertEF_2d(i)*xyzEF_2d(3,i);

            xyzSEF_2d(1,i)=xyzEF_2d(1,i)+smufa*vettox  
            xyzSEF_2d(2,i)=xyzEF_2d(2,i)+smufa*vettoy  
            xyzSEF_2d(3,i)=xyzEF_2d(3,i)+smufa*vettoz  
               
         enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP   


         do i=1,nBundPurk !sono solo 2
            fslave=fPurkslave(i)
            v1=vert_of_faceEF_2d(1,fslave)
            v2=vert_of_faceEF_2d(2,fslave)
            v3=vert_of_faceEF_2d(3,fslave)
            xyzSEF_2d(1:3,v1)=xyzEF_2d(1:3,v1)
            xyzSEF_2d(1:3,v2)=xyzEF_2d(1:3,v2)
            xyzSEF_2d(1:3,v3)=xyzEF_2d(1:3,v3)
         enddo

#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif       
         do i=1,nvtotEF_2d
            xyzS1EF_2d(1,i) = xyzSEF_2d(1,i) 
            xyzS1EF_2d(2,i) = xyzSEF_2d(2,i) 
            xyzS1EF_2d(3,i) = xyzSEF_2d(3,i) 
         enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP   

!smooth bordo
#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                                                                                                                     
#endif       
         do i=1,count2EF_2d
            vv=vert_boundEF_2d(i)
            v1=vert_of_vert_boundEF_2d(1,i)
            v2=vert_of_vert_boundEF_2d(2,i)
      
            vettox=xyzSEF_2d(1,v1)+xyzSEF_2d(1,v2)-2.0d0*xyzSEF_2d(1,vv)
            vettoy=xyzSEF_2d(2,v1)+xyzSEF_2d(2,v2)-2.0d0*xyzSEF_2d(2,vv)
            vettoz=xyzSEF_2d(3,v1)+xyzSEF_2d(3,v2)-2.0d0*xyzSEF_2d(3,vv)
            
            xyzS1EF_2d(1,vv)=xyzS1EF_2d(1,vv)+smufa*vettox    
            xyzS1EF_2d(2,vv)=xyzS1EF_2d(2,vv)+smufa*vettoy    
            xyzS1EF_2d(3,vv)=xyzS1EF_2d(3,vv)+smufa*vettoz    
         enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP   



#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                                                                                                                     
#endif       
         do i=1,nvtotEF_2d
            xyzEF_2d(1,i) = xyzS1EF_2d(1,i) 
            xyzEF_2d(2,i) = xyzS1EF_2d(2,i) 
            xyzEF_2d(3,i) = xyzS1EF_2d(3,i) 
         enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP   

      enddo


      
      cell_type(:) = 5 !5 is for triangular elements

      ! do inp=1,NparticleEF_2d
      !  scalar_face(1:maxnf) = 0.0
      !  vert_face_dum(1:3,fstartEF_2d(inp):fendEF_2d(inp))=vert_of_faceEF_2d(1:3,fstartEF_2d(inp):fendEF_2d(inp))-vstartEF_2d(inp)
      ! end do

! do inp=1,NparticleEF_2d
!       tnv = vendEF_2d(inp) - vstartEF_2d(inp) + 1
      !       tnf = fendEF_2d(inp) - fstartEF_2d(inp) + 1
      tnv=nvtotEF_2d
      tnf=nftotEF_2d
!       write(ipfip,94)inp
!    94 format(i2.2)

      ! filname = 'vtkfiles/strucEF_2d_'//trim(ipfip)//'_'//trim(ipfi)//'.vtk'
      filname = 'vtkfiles/strucEF_2d_'//trim(ipfi)//'.vtk'

      open(121,file = filname)
!Header     
        write(121,'(a)')adjustl('# vtk DataFile Version 3.1')
        write(121,'(a)')adjustl('Stores triangular gts mesh')
        write(121,'(a)')adjustl('ASCII')
        write(121,'(a)')adjustl('DATASET UNSTRUCTURED_GRID')
        write(121,*)''
        write(121,*)'POINTS ',tnv,' FLOAT'
!        do i=vstartEF_2d(inp),vendEF_2d(inp)
      do i=1,nvtotEF_2d
        write(121,*)xyzEF_2d(1:3,i)
      end do
        write(121,*)''
        write(121,*)'CELLS ',tnf, 4*tnf
!        do i=fstartEF_2d(inp),fendEF_2d(inp)
      do i=1,nftotEF_2d
!         write(121,*)'3 ',vert_face_dum(1:3,i)
         write(121,*)'3 ',vert_of_faceEF_2d(1:3,i)-1
      end do
        write(121,*)''
        write(121,*)'CELL_TYPES ',tnf
        ! write(121,*)cell_type(fstartEF_2d(inp):fendEF_2d(inp))
        write(121,*)cell_type(1:nftotEF_2d)
! 
!        write(121,*)''
!        write(121,*)'CELL_DATA ',fend
!        write(121,*)scalar_face(:)

        if (optionvtk.EQ.0) then
           write(121,*)''
           ! write(121,*)'CELL_DATA ',nfiEF_2d(inp)
           write(121,*)'CELL_DATA ',nftotEF_2d
           write(121,*)'Scalars LabelPurkSetto FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           ! do i=fstartEF_2d(inp),fendEF_2d(inp)
           do i=1,nftotEF_2d
              write(121,*) LabelPurkSetto(i)
           end do
           
           write(121,*)'Scalars Istim_S1 FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           ! do i=fstartEF_2d(inp),fendEF_2d(inp)
           do i=1,nftotEF_2d
              write(121,*) IstimEF_2dS1(i)
           end do
           
           write(121,*)'Scalars Istim_S2 FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           ! do i=fstartEF_2d(inp),fendEF_2d(inp)
           do i=1,nftotEF_2d
              write(121,*) IstimEF_2dS2(i)
           end do
           
           write(121,*)'' 
           ! write(121,*)'POINT_DATA ',nviEF_2d(inp)
           write(121,*)'POINT_DATA ',nvtotEF_2d
           write(121,*)'Scalars CARTO_D FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           ! do i=vstartEF_2d(inp),vendEF_2d(inp)
           do i=1,nvtotEF_2d
              write(121,*) CARTO_Dnode(i)
           end do
           
        elseif (optionvtk.EQ.1) then
           write(121,*)'' 
           ! write(121,*)'POINT_DATA ',nviEF_2d(inp)
           write(121,*)'POINT_DATA ',nvtotEF_2d
           write(121,*)'Scalars potEF_2d FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           ! do i=vstartEF_2d(inp),vendEF_2d(inp)
           do i=1,nvtotEF_2d
              !write(121,*) potEFnode_2d(i)
              write(121,*) potEFnode_2d(i)
           end do
           write(121,*)'' 
           ! write(121,*)'POINT_DATA ',nviEF_2d(inp)
           !write(121,*)'POINT_DATA ',nvtotEF_2d
           write(121,*)'Scalars Diff_node FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           ! do i=vstartEF_2d(inp),vendEF_2d(inp)
           do i=1,nvtotEF_2d
              !write(121,*) potEFnode_2d(i)
              write(121,*) CARTO_Dnode(i)
           end do
           write(121,*)''
           ! write(121,*)'POINT_DATA ',nviEF_2d(inp)
           !write(121,*)'POINT_DATA ',nvtotEF_2d
           write(121,*)'Scalars LAT FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           ! do i=vstartEF_2d(inp),vendEF_2d(inp)
           do i=1,nvtotEF_2d
              !write(121,*) potEFnode_2d(i)
              write(121,*) t_apd_2d(1,i)
           end do
           write(121,*)'' 
           ! write(121,*)'POINT_DATA ',nviEF_2d(inp)
           write(121,*)'CELL_DATA ',nftotEF_2d
           write(121,*)'Scalars Diffusivity_face FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           ! do i=vstartEF_2d(inp),vendEF_2d(inp)
           do i=1,nftotEF_2d
              write(121,*) CARTO_Dface(i)
           end do
           write(121,*)'' 
           ! write(121,*)'POINT_DATA ',nviEF_2d(inp)
           !write(121,*)'CELL_DATA ',nftotEF_2d
           write(121,*)'Scalars CV_face FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           ! do i=vstartEF_2d(inp),vendEF_2d(inp)
           do i=1,nftotEF_2d
              write(121,*) CVface_EF_2d(i)
           end do
           write(121,*)'' 
           ! write(121,*)'POINT_DATA ',nviEF_2d(inp)
           !write(121,*)'CELL_DATA ',nftotEF_2d
           write(121,*)'Scalars Diff_CV_face FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           ! do i=vstartEF_2d(inp),vendEF_2d(inp)
           do i=1,nftotEF_2d
              write(121,*) sqrt(CVface_EF_2d(i))
           end do
        endif !optionvtk
      
      close(121)
     
! end do
      return
      end
      
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      subroutine write_to_vtkEF_1d(ipfi,optionvtk)
      use mpi_param
      use mls_param
      use mpih
      use param
!@cuf   use cudafor
      implicit none
     
      integer :: i,inp,itime,v2,v1,kk
      integer :: tnv, tne,optionvtk,vmaster,vslave
      integer, dimension(netotEF_1d) :: cell_type      
      integer, dimension(2,netotEF_1d) :: vert_edge_dum
      
      real(DP) :: tprfi,vettox,vettoy,vettoz,smufa
      character*70 filname
      character*100 ipfi
!@cuf   integer :: istat 

!move EF_1d
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do i = 1,nvtotEF_1d
        v2 = boundary2EF_1d(i)
           xyzEF_1d(1,i) = xyz(1,v2)+BCoffsetEF_1d(1,i)   
           xyzEF_1d(2,i) = xyz(2,v2)+BCoffsetEF_1d(2,i)   
           xyzEF_1d(3,i) = xyz(3,v2)+BCoffsetEF_1d(3,i)   
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP


!garantisci collegamento master/slav
       !fix master/slave nodes                                                                                                                                               
      do i=1,nBundAV
         vmaster=vAVmaster(i)
         vslave=vAVslave(i)
         xyzEF_1d(1:3,vslave)=xyzEF_1d(1:3,vmaster)
      enddo


!smooth EF_1d
#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                                                                                                                     
#endif
      do i = 1,nvtotEF_1d
           xyzSEF_1d(1,i) = xyzEF_1d(1,i)
           xyzSEF_1d(2,i) = xyzEF_1d(2,i)
           xyzSEF_1d(3,i) = xyzEF_1d(3,i)
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP   
      
      smufa=0.5D0
      do kk=1,500
#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                                                                                                                     
#endif
         do i=1,nvtotEF_1d
            if (n_edge_of_vertEF_1d(i).EQ.2) then   
               v1=vert_of_vertEF_1d(1,i);
               v2=vert_of_vertEF_1d(2,i);
        
               vettox = xyzEF_1d(1,v1)+xyzEF_1d(1,v2)-2.0d0*xyzEF_1d(1,i);
               vettoy = xyzEF_1d(2,v1)+xyzEF_1d(2,v2)-2.0d0*xyzEF_1d(2,i);
               vettoz = xyzEF_1d(3,v1)+xyzEF_1d(3,v2)-2.0d0*xyzEF_1d(3,i);

               xyzSEF_1d(1,i) = xyzEF_1d(1,i) + smufa*vettox
               xyzSEF_1d(2,i) = xyzEF_1d(2,i) + smufa*vettoy
               xyzSEF_1d(3,i) = xyzEF_1d(3,i) + smufa*vettoz
            endif
         enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP   

       !fix master/slave nodes, non si muovono durante lo smoothing                                                                                                                           
         do i=1,nBundAV
            vmaster=vAVmaster(i)
            vslave=vAVslave(i)
            xyzSEF_1d(1:3,vmaster)=xyzEF_1d(1:3,vmaster)
            xyzSEF_1d(1:3,vslave)=xyzEF_1d(1:3,vslave)
         enddo

#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                                                                                                                     
#endif       
         do i=1,nvtotEF_1d
            xyzEF_1d(1,i) = xyzSEF_1d(1,i) 
            xyzEF_1d(2,i) = xyzSEF_1d(2,i) 
            xyzEF_1d(3,i) = xyzSEF_1d(3,i) 
         enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP   

      enddo

      cell_type(:) = 3 !3 is for linear elements
      do i=1, netotEF_1d
        vert_edge_dum(1,i) = vert_of_edgeEF_1d(1,i)-1
        vert_edge_dum(2,i) = vert_of_edgeEF_1d(2,i)-1
      end do
      
104 format((2x,f12.8))
   94 format(i2.2)

      filname = 'vtkfiles/strucEF_1d_'//trim(ipfi)//'.vtk'

      open(121,file = filname)
!Header     
        write(121,'(a)')adjustl('# vtk DataFile Version 3.1')
        write(121,'(a)')adjustl('Stores triangular gts mesh')
        write(121,'(a)')adjustl('ASCII')
        write(121,'(a)')adjustl('DATASET UNSTRUCTURED_GRID')
        write(121,*)''
        write(121,*)'POINTS ',nvtotEF_1d,' FLOAT'
      do i=1,nvtotEF_1d
        write(121,*)xyzEF_1d(1:3,i)
      end do
        write(121,*)''
        write(121,*)'CELLS ',netotEF_1d, 3*netotEF_1d
      do i=1, netotEF_1d
        write(121,*)'2 ',vert_edge_dum(1:2,i)
      end do
        write(121,*)''
        write(121,*)'CELL_TYPES ',netotEF_1d
        write(121,*)cell_type(1:netotEF_1d)
! 
!        write(121,*)''
!        write(121,*)'CELL_DATA ',fend
!        write(121,*)scalar_face(:)

        if (optionvtk.EQ.0) then
!           write(121,*)''
           write(121,*)'POINT_DATA ',nvtotEF_1d
           write(121,*)'Scalars Istim FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=1,nvtotEF_1d
              write(121,104) IstimEF_1d(i)
           end do
        elseif (optionvtk.EQ.1) then
!           write(121,*)'' 
           write(121,*)'POINT_DATA ',nvtotEF_1d
           write(121,*)'Scalars potEF_1d FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=1,nvtotEF_1d
              write(121,104) potEFnode_1d(i)
           end do

        endif !optionvtk
      
      close(121)
     

      return
      end

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      subroutine write_to_vtkEF_1dSplit(ipfi,optionvtk)
      use mpi_param
      use mls_param
      use mpih
      use param

      implicit none
     
      integer :: i,inp,itime
      integer :: tnv, tne,optionvtk
      integer, dimension(netotEF_1d) :: cell_type      
      real(DP), dimension(maxnf) :: scalar_face
      integer, dimension(2,netotEF_1d) :: vert_edge_dum
      
      real(DP) :: tprfi
      character*70 filname
      character*100 ipfi,ipfip

      cell_type(:) = 3 !3 is for linear elements
      

      do inp=1,NparticleEF_1d
       scalar_face(1:maxnf) = 0.0
       vert_edge_dum(1:2,estartEF_1d(inp):eendEF_1d(inp))=vert_of_edgeEF_1d(1:2,estartEF_1d(inp):eendEF_1d(inp))-vstartEF_1d(inp)
      end do

do inp=1,NparticleEF_1d
      tnv = vendEF_1d(inp) - vstartEF_1d(inp) + 1
      tne = eendEF_1d(inp) - estartEF_1d(inp) + 1
      write(ipfip,94)inp
   94 format(i2.2)

      filname = 'vtkfiles/strucEF_1d_'//trim(ipfip)//'_'//trim(ipfi)//'.vtk'

      open(121,file = filname)
!Header     
        write(121,'(a)')adjustl('# vtk DataFile Version 3.1')
        write(121,'(a)')adjustl('Stores triangular gts mesh')
        write(121,'(a)')adjustl('ASCII')
        write(121,'(a)')adjustl('DATASET UNSTRUCTURED_GRID')
        write(121,*)''
        write(121,*)'POINTS ',tnv,' FLOAT'
      do i=vstartEF_1d(inp),vendEF_1d(inp)
        write(121,*)xyzEF_1d(1:3,i)
      end do
        write(121,*)''
        write(121,*)'CELLS ',tne, 3*tne
      do i=estartEF_1d(inp),eendEF_1d(inp)
        write(121,*)'2 ',vert_edge_dum(1:2,i)
      end do
        write(121,*)''
        write(121,*)'CELL_TYPES ',tne
        write(121,*)cell_type(estartEF_1d(inp):eendEF_1d(inp))
! 
!        write(121,*)''
!        write(121,*)'CELL_DATA ',fend
!        write(121,*)scalar_face(:)

        if (optionvtk.EQ.0) then
!           write(121,*)''
           write(121,*)'POINT_DATA ',nviEF_1d(inp)
           write(121,*)'Scalars Istim FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstartEF_1d(inp),vendEF_1d(inp)
              write(121,*) IstimEF_1d(i)
           end do
        elseif (optionvtk.EQ.1) then
!           write(121,*)'' 
           write(121,*)'POINT_DATA ',nviEF_1d(inp)
           write(121,*)'Scalars potEF_1d FLOAT'
           write(121,*)'LOOKUP_TABLE default'
           do i=vstartEF_1d(inp),vendEF_1d(inp)
              write(121,*) potEFnode_1d(i)
           end do

        endif !optionvtk
      
      close(121)
     
end do
      return
      end
      


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      subroutine write_to_vtkCV_1d(ipfi,optionvtk)
      use mpi_param
      use mls_param
      use mpih
      use param
!@cuf   use cudafor
      implicit none
     
      integer :: i,inp,itime,v2,v1,kk
      integer :: tnv, tne,optionvtk,vmaster,vslave
      integer, dimension(netotCV_1d) :: cell_type      
      integer, dimension(2,netotCV_1d) :: vert_edge_dum
      
      real(DP) :: tprfi,vettox,vettoy,vettoz,smufa
      character*70 filname
      character*100 ipfi
!@cuf   integer :: istat 

!move CV_1d
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do i = 1,nvtotCV_1d
        v2 = boundary2CV_1d(i)
           xyzCV_1d(1,i) = xyz_3d(1,v2)+BCoffsetCV_1d(1,i)   
           xyzCV_1d(2,i) = xyz_3d(2,v2)+BCoffsetCV_1d(2,i)   
           xyzCV_1d(3,i) = xyz_3d(3,v2)+BCoffsetCV_1d(3,i)   
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP



!smooth CV_1d
#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                                                                                                                     
#endif
      do i = 1,nvtotCV_1d
           xyzSCV_1d(1,i) = xyzCV_1d(1,i)
           xyzSCV_1d(2,i) = xyzCV_1d(2,i)
           xyzSCV_1d(3,i) = xyzCV_1d(3,i)
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP   
      
      smufa=0.5D0
      do kk=1,500
#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                                                                                                                     
#endif
         do i=1,nvtotCV_1d
            if (n_edge_of_vertCV_1d(i).EQ.2) then   
               v1=vert_of_vertCV_1d(1,i);
               v2=vert_of_vertCV_1d(2,i);
        
               vettox = xyzCV_1d(1,v1)+xyzCV_1d(1,v2)-2.0d0*xyzCV_1d(1,i);
               vettoy = xyzCV_1d(2,v1)+xyzCV_1d(2,v2)-2.0d0*xyzCV_1d(2,i);
               vettoz = xyzCV_1d(3,v1)+xyzCV_1d(3,v2)-2.0d0*xyzCV_1d(3,i);

               xyzSCV_1d(1,i) = xyzCV_1d(1,i) + smufa*vettox
               xyzSCV_1d(2,i) = xyzCV_1d(2,i) + smufa*vettoy
               xyzSCV_1d(3,i) = xyzCV_1d(3,i) + smufa*vettoz
            endif
         enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP   

#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                                                                                                                     
#endif       
         do i=1,nvtotCV_1d
            xyzCV_1d(1,i) = xyzSCV_1d(1,i) 
            xyzCV_1d(2,i) = xyzSCV_1d(2,i) 
            xyzCV_1d(3,i) = xyzSCV_1d(3,i) 
         enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP   

      enddo





      cell_type(:) = 3 !3 is for linear elements

      do i=1, netotCV_1d
        vert_edge_dum(1,i) = vert_of_edgeCV_1d(1,i)-1
        vert_edge_dum(2,i) = vert_of_edgeCV_1d(2,i)-1
      end do
      
104 format((2x,f12.8))
   94 format(i2.2)

      filname = 'vtkfiles/strucCV_1d_'//trim(ipfi)//'.vtk'

      open(121,file = filname)
!Header     
        write(121,'(a)')adjustl('# vtk DataFile Version 3.1')
        write(121,'(a)')adjustl('Stores triangular gts mesh')
        write(121,'(a)')adjustl('ASCII')
        write(121,'(a)')adjustl('DATASET UNSTRUCTURED_GRID')
        write(121,*)''
        write(121,*)'POINTS ',nvtotCV_1d,' FLOAT'
      do i=1,nvtotCV_1d
        write(121,*)xyzCV_1d(1:3,i)
      end do
        write(121,*)''
        write(121,*)'CELLS ',netotCV_1d, 3*netotCV_1d
      do i=1, netotCV_1d
        write(121,*)'2 ',vert_edge_dum(1:2,i)
      end do
        write(121,*)''
        write(121,*)'CELL_TYPES ',netotCV_1d
        write(121,*)cell_type(1:netotCV_1d)
! 
!        write(121,*)''
!        write(121,*)'CELL_DATA ',fend
!        write(121,*)scalar_face(:)

!         if (optionvtk.EQ.0) then
! !           write(121,*)''
!            write(121,*)'POINT_DATA ',nvtotCV_1d
!            write(121,*)'Scalars Istim FLOAT'
!            write(121,*)'LOOKUP_TABLE default'
!            do i=1,nvtotCV_1d
!               write(121,104) IstimCV_1d(i)
!            end do
!         elseif (optionvtk.EQ.1) then
! !           write(121,*)'' 
!            write(121,*)'POINT_DATA ',nvtotCV_1d
!            write(121,*)'Scalars potCV_1d FLOAT'
!            write(121,*)'LOOKUP_TABLE default'
!            do i=1,nvtotCV_1d
!               write(121,104) potCVnode_1d(i)
!            end do

!         endif !optionvtk
      
      close(121)

! !leads
!       filname = 'vtkfiles/struc_leads_'//trim(ipfi)//'.csv'

!       open(121,file = filname)
! !Header     
! !        write(121,'(a)')adjustl('x coord, y coord, z coord, scalar')
!         write(121,'(a)')adjustl('x coord, y coord, z coord')
!         ! write(121,*)cell_bar(1:3,lead_RA),1
!         ! write(121,*)cell_bar(1:3,lead_RV),2
!         ! write(121,*)cell_bar(1:3,lead_LV),3
!         write(121,*)cell_bar(1,lead_RA),cell_bar(2,lead_RA),cell_bar(3,lead_RA)
!         write(121,*)cell_bar(1,lead_RV),cell_bar(2,lead_RV),cell_bar(3,lead_RV)
!         write(121,*)cell_bar(1,lead_LV),cell_bar(2,lead_LV),cell_bar(3,lead_LV)
!       close(121)     




!leads
!LV
!       filname = 'vtkfiles/struc_leadLV_'//trim(ipfi)//'.vtk'
!       open(121,file = filname)
!         write(121,'(a)')adjustl('# vtk DataFile Version 3.1')
!         write(121,'(a)')adjustl('Stores triangular gts mesh')
!         write(121,'(a)')adjustl('ASCII')
!         write(121,'(a)')adjustl('DATASET UNSTRUCTURED_GRID')
!         write(121,*)''
!         write(121,*)'POINTS ',2,' FLOAT'
!         write(121,*)cell_bar(1:3,lead_LV)
!         write(121,*)cell_bar(1:3,lead_LV)
!         write(121,*)''
!         write(121,*)'CELLS    1    3'
!         write(121,*)'2    0     1'
!         write(121,*)''
!         write(121,*)'CELL_TYPES   1'
!         write(121,*)'     3'
!       close(121)


! !RV
!       filname = 'vtkfiles/struc_leadRV_'//trim(ipfi)//'.vtk'
!       open(121,file = filname)
!         write(121,'(a)')adjustl('# vtk DataFile Version 3.1')
!         write(121,'(a)')adjustl('Stores triangular gts mesh')
!         write(121,'(a)')adjustl('ASCII')
!         write(121,'(a)')adjustl('DATASET UNSTRUCTURED_GRID')
!         write(121,*)''
!         write(121,*)'POINTS ',2,' FLOAT'
!         write(121,*)cell_bar(1:3,lead_RV)
!         write(121,*)cell_bar(1:3,lead_RV)
!         write(121,*)''
!         write(121,*)'CELLS    1    3'
!         write(121,*)'2    0     1'
!         write(121,*)''
!         write(121,*)'CELL_TYPES   1'
!         write(121,*)'     3'
!       close(121)


! !RA
!       filname = 'vtkfiles/struc_leadRA_'//trim(ipfi)//'.vtk'
!       open(121,file = filname)
!         write(121,'(a)')adjustl('# vtk DataFile Version 3.1')
!         write(121,'(a)')adjustl('Stores triangular gts mesh')
!         write(121,'(a)')adjustl('ASCII')
!         write(121,'(a)')adjustl('DATASET UNSTRUCTURED_GRID')
!         write(121,*)''
!         write(121,*)'POINTS ',2,' FLOAT'
!         write(121,*)cell_bar(1:3,lead_RA)
!         write(121,*)cell_bar(1:3,lead_RA)
!         write(121,*)''
!         write(121,*)'CELLS    1    3'
!         write(121,*)'2    0     1'
!         write(121,*)''
!         write(121,*)'CELL_TYPES   1'
!         write(121,*)'     3'
!       close(121)







        
      return
      end

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      subroutine write_to_vtkCenSeg_1d(ipfi,optionvtk)
      use mpi_param
      use mls_param
      use mpih
      use param
!@cuf   use cudafor
      implicit none
     
      integer :: i,inp,itime,v2,v1,kk,netotseg
      integer :: tnv, tne,optionvtk,vmaster,vslave
      integer, dimension(netotCV_1d) :: cell_type      
      integer, dimension(2,netotCV_1d) :: vert_edge_dum
      
      real(DP) :: tprfi,vettox,vettoy,vettoz,smufa
      character*70 filname
      character*100 ipfi
!@cuf   integer :: istat 

      netotseg=16
      cell_type(:) = 3 !3 is for linear elements

104 format((2x,f12.8))
   94 format(i2.2)

      filname = 'vtkfiles/strucCenSeg_1d_'//trim(ipfi)//'.vtk'

      open(121,file = filname)
!Header     
        write(121,'(a)')adjustl('# vtk DataFile Version 3.1')
        write(121,'(a)')adjustl('Stores triangular gts mesh')
        write(121,'(a)')adjustl('ASCII')
        write(121,'(a)')adjustl('DATASET UNSTRUCTURED_GRID')
        write(121,*)''
        write(121,*)'POINTS ',17,' FLOAT'
      do i=1,17
        write(121,*)xyz_seg(1:3,i)
      end do
        write(121,*)''
        write(121,*)'CELLS ',netotseg, 3*netotseg
      do i=1, netotseg
        write(121,*)'2 ',i,i+1
      end do
        write(121,*)''
        write(121,*)'CELL_TYPES ',netotseg
        write(121,*)cell_type(1:netotseg)
        close(121)
      return
      end
