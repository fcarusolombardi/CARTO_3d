!================================================
      subroutine mpi_write_continua
      use param
      use mpih
      use mpi_param, only: kstart,kend,kstartr,kendr
      use local_arrays, only: dens,q2,q3,q1,dsal,pr
      use hdf5
      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_q1
      integer(HID_T) :: dset_q2
      integer(HID_T) :: dset_q3
      integer(HID_T) :: dset_dens
      integer(HID_T) :: dset_dsal
      integer(HID_T) :: dset_pr

      integer(HID_T) :: plist_id

      integer(HSIZE_T), dimension(3) :: dims
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims, itime

      character filnam1*250
      character filnam2*250
      character filnam3*250
      character filnam4*250
      character filnam5*250
      character filnam6*250
      character*100 ipfi

      call h5open_f(hdf_error)

!   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!   Form the name of the file

      ! itime=nint(1000.*time)          !FV                                                      
      itime=nint(1000.*time*TSTAR) !cosi salva in ms
      itime=max(itime,0)
      write(ipfi,199) itime
 199   format(i8.8)

      if (PPro.LE.0) then
         filnam1 = 'restart/continua_dens'//trim(ipfi)//'.h5'
         filnam2 = 'restart/continua_q1'//trim(ipfi)//'.h5'
         filnam3 = 'restart/continua_q2'//trim(ipfi)//'.h5'
         filnam4 = 'restart/continua_q3'//trim(ipfi)//'.h5'
         filnam5 = 'restart/continua_dsal'//trim(ipfi)//'.h5'
         filnam6 = 'restart/continua_pr'//trim(ipfi)//'.h5'
      else
         filnam1 = 'cfield/continua_dens'//trim(ipfi)//'.h5'
         filnam2 = 'cfield/continua_q1'//trim(ipfi)//'.h5'
         filnam3 = 'cfield/continua_q2'//trim(ipfi)//'.h5'
         filnam4 = 'cfield/continua_q3'//trim(ipfi)//'.h5'
         filnam5 = 'cfield/continua_dsal'//trim(ipfi)//'.h5'
         filnam6 = 'cfield/continua_pr'//trim(ipfi)//'.h5'
      endif

!   Set offsets and element counts
      ndims = 3

      dims(1)=n1
      dims(2)=n2
      dims(3)=n3m

      data_count(1) = n1
      data_count(2) = n2
      data_count(3) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstart-1

!   dens
#ifdef CALTEMP
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id, &
                      hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'dens', H5T_NATIVE_DOUBLE, &
                     filespace, dset_dens, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset_dens, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)
      call h5dwrite_f(dset_dens, H5T_NATIVE_DOUBLE, &
                     dens(1:n1,1:n2,kstart:kend), dims, & 
                     hdf_error, file_space_id = slabspace, &
                     mem_space_id = memspace, &
                     xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_dens, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)
#endif

!   q1

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id, &
                      hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vth', H5T_NATIVE_DOUBLE, &
                     filespace, dset_q1, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset_q1, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)
      call h5dwrite_f(dset_q1, H5T_NATIVE_DOUBLE, &
                     q1(1:n1,1:n2,kstart:kend), dims, &
                     hdf_error, file_space_id = slabspace, &
                     mem_space_id = memspace,  &
                     xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q1, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!   q2

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam3, H5F_ACC_TRUNC_F, file_id, &
                      hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vr', H5T_NATIVE_DOUBLE, &
                     filespace, dset_q2, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset_q2, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)
      call h5dwrite_f(dset_q2, H5T_NATIVE_DOUBLE, &
                      q2(1:n1,1:n2,kstart:kend), dims, &
                      hdf_error, file_space_id = slabspace, &
                      mem_space_id = memspace, &
                      xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q2, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!   q3

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam4, H5F_ACC_TRUNC_F, file_id, &
                      hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE, &
                     filespace, dset_q3, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset_q3, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)
      call h5dwrite_f(dset_q3, H5T_NATIVE_DOUBLE, &
                     q3(1:n1,1:n2,kstart:kend), dims,  &
                     hdf_error, file_space_id = slabspace, &
                     mem_space_id = memspace, &
                     xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q3, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!   pressure

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam6, H5F_ACC_TRUNC_F, file_id, &
                      hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Pr', H5T_NATIVE_DOUBLE, &
                     filespace, dset_pr, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_pr, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)
      call h5dwrite_f(dset_pr, H5T_NATIVE_DOUBLE, &
                     pr(1:n1,1:n2,kstart:kend), dims, &
                     hdf_error, file_space_id = slabspace, &
                     mem_space_id = memspace, &
                     xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_pr, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!   Set offsets and element counts

      ndims = 3
      dims(1)=n1r
      dims(2)=n2r
      dims(3)=n3mr

      data_count(1) = n1r
      data_count(2) = n2r
      data_count(3) = kendr-kstartr+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstartr-1

!   dsal
#ifdef CALSCAL
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam5, H5F_ACC_TRUNC_F, file_id, &
                      hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'dsal', H5T_NATIVE_DOUBLE, &
                      filespace, dset_dsal, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_dsal, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)
      call h5dwrite_f(dset_dsal, H5T_NATIVE_DOUBLE, &
                     dsal(1:n1r,1:n2r,kstartr:kendr), dims, &
                     hdf_error, file_space_id = slabspace, &
                     mem_space_id = memspace, &
                     xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_dsal, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)
#endif

      call h5close_f(hdf_error)

      if (myid .eq. 0) then
         if (PPro.LE.0) then
            open(13,file='restart/continua_grid'//trim(ipfi)//'.dat',status='unknown')
            rewind(13)                                                      
            write(13,'(3i8)') n1,n2,n3
            write(13,'(3f25.15)') rext1,rext2,time
            write(13,'(i8,f25.15)') istr3,str3
            write(13,'(3i8)') mref1, mref2, mref3
            close(13)

!  xmf file for continua data files.
            open(45,file='restart/field3d_'//trim(ipfi)//'.xmf',status='unknown')
         else
            open(45,file='cfield/field3d_'//trim(ipfi)//'.xmf',status='unknown')
         endif !Ppro
         rewind(45)
         write(45,'("<?xml version=""1.0"" ?>")')
         write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
         write(45,'("<Xdmf Version=""2.0"">")')
         write(45,'("<Domain>")')
         
         write(45,'("<Grid Name=""DDC"" GridType=""Uniform"">")')
         write(45,'("<Topology TopologyType=""3DRectMesh""  &
              NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2,n1
         write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
         write(45,'("<DataItem Dimensions=""",i4,""" &
              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1
         write(45,'("field_gridc.h5:/xc")')
         write(45,'("</DataItem>")')
         write(45,'("<DataItem Dimensions=""",i4,""" &
              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2
         write(45,'("field_gridc.h5:/yc")')
         write(45,'("</DataItem>")')
         write(45,'("<DataItem Dimensions=""",i4,""" &
              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m
         write(45,'("field_gridc.h5:/zc")')
         write(45,'("</DataItem>")')
         write(45,'("</Geometry>")')
         !u velocity
         write(45,'("<Attribute Name=""u"" & 
              AttributeType=""Scalar"" Center=""Node"">")')
         write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" &
              NumberType=""Float"" Precision=""4"" Format=""HDF"">")') n3m,n2,n1
         write(45,'("continua_q1'//trim(ipfi)//'.h5:/Vth")')
         write(45,'("</DataItem>")')
         write(45,'("</Attribute>")')
         !v velocity
         write(45,'("<Attribute Name=""v"" &
              AttributeType=""Scalar"" Center=""Node"">")')
         write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" &
              NumberType=""Float"" Precision=""4"" Format=""HDF"">")') n3m,n2,n1
         write(45,'("continua_q2'//trim(ipfi)//'.h5:/Vr")')
         write(45,'("</DataItem>")')
         write(45,'("</Attribute>")')
         !w velocity
         write(45,'("<Attribute Name=""w"" &
              AttributeType=""Scalar"" Center=""Node"">")')
         write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" &
              NumberType=""Float"" Precision=""4"" Format=""HDF"">")') n3m,n2,n1
         write(45,'("continua_q3'//trim(ipfi)//'.h5:/Vz")')
         write(45,'("</DataItem>")')
         write(45,'("</Attribute>")')
         !pressure
         write(45,'("<Attribute Name=""pr"" &
              AttributeType=""Scalar"" Center=""Node"">")')
         write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" &
              NumberType=""Float"" Precision=""4"" Format=""HDF"">")') n3m,n2,n1
         write(45,'("continua_pr'//trim(ipfi)//'.h5:/Pr")')
         write(45,'("</DataItem>")')
         write(45,'("</Attribute>")')
#ifdef CALTEMP        
         write(45,'("<Attribute Name=""T"" &
              AttributeType=""Scalar"" Center=""Node"">")')
         write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" &
              NumberType=""Float"" Precision=""4"" Format=""HDF"">")') n3m,n2,n1
         write(45,'("continua_dens'//trim(ipfi)//'.h5:/dens")')
         write(45,'("</DataItem>")')
         write(45,'("</Attribute>")')
#endif        
#ifdef CALSCAL
         write(45,'("<Attribute Name=""S"" &
              AttributeType=""Scalar"" Center=""Node"">")')
         write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" &
              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')  &
              n3mr,n2r,n1r
         write(45,'("continua_dsal'//trim(ipfi)//'.h5:/dsal")')
         write(45,'("</DataItem>")')
         write(45,'("</Attribute>")')
#endif        
         write(45,'("<Time Value=""",e12.5,""" />")') time
         write(45,'("</Grid>")')
        
         write(45,'("</Domain>")')
         write(45,'("</Xdmf>")')
         close(45) 
         
      endif !myid

      return
      end subroutine mpi_write_continua

!================================================      
      subroutine mpi_read_continua(n1o,n2o,n3o,ks,ke,intvar,qua)
      use mpih
      use param
      use hdf5
      implicit none
      integer, intent(in) :: ks,ke,n2o,n1o,n3o
      real(DP), dimension(1:n1o,1:n2o,ks-lvlhalo:ke+lvlhalo)::qua
      integer k,j,i

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_qua

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims,itime

      integer, intent(in) :: intvar
      character(250)  :: filnam1
      character(100) :: ipfi
      character(10)  :: dsetname

      call h5open_f(hdf_error)

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!   Select file and dataset based on intvar
      write(ipfi,199) nwrit !FV
 199   format(i8.8)

      select case (intvar)
        case (1)
          dsetname = trim('Vth')
          filnam1 = trim(trim(folderload)//'continua_q1'//trim(ipfi)//'.h5')
        case (2)
          dsetname = trim('Vr')
          filnam1 = trim(trim(folderload)//'continua_q2'//trim(ipfi)//'.h5')
        case (3)
          dsetname = trim('Vz')
          filnam1 = trim(trim(folderload)//'continua_q3'//trim(ipfi)//'.h5')
        case (4)
          dsetname = trim('dens')
          filnam1 = trim(trim(folderload)//'continua_dens'//trim(ipfi)//'.h5')
        case (5)
          dsetname = trim('dsal')
          filnam1 = trim(trim(folderload)//'continua_dsal'//trim(ipfi)//'.h5')
        case (6)
          dsetname = trim('Pr')
          filnam1 = trim(trim(folderload)//'continua_pr'//trim(ipfi)//'.h5')

      end select

      do k=ks,ke
!$OMP  PARALLEL DO &
!$OMP PRIVATE(i,j)
        do j=1,n2o
          do i=1,n1o
            qua(i,j,k)=0.d0
          enddo
        enddo
!$OMP  END PARALLEL DO
      enddo

!   Set offsets and element counts
   
      ndims = 3

      dims(1)=n1o
      dims(2)=n2o
      dims(3)=n3o-1


      data_count(1) = n1o
      data_count(2) = n2o
      data_count(3) = ke-ks+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = ks-1

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fopen_f(filnam1, H5F_ACC_RDONLY_F, file_id, &
                    hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id,hdf_error)

      call h5dopen_f(file_id, dsetname, dset_qua, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset_qua, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)

      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)

      call h5dread_f(dset_qua, H5T_NATIVE_DOUBLE, & 
                    qua(1:n1o,1:n2o,ks:ke), dims, & 
                    hdf_error, file_space_id = slabspace, &
                    mem_space_id = memspace, & 
                    xfer_prp = plist_id)

      call h5dclose_f(dset_qua, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5pclose_f(plist_id, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      call h5close_f(hdf_error)

      if(myid.eq.0)write(*,'(5x,a)')'reading complete: '//filnam1

      return
      end subroutine mpi_read_continua

!============================================================
      subroutine write_grid_info
      use mpih
      use param
      use hdf5

      IMPLICIT none

      integer hdf_error

      integer(HID_T) :: file_grid
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid

      integer(HSIZE_T) :: dims_grid(1)

      character namfile*70

      call h5open_f(hdf_error)

      IF(myid.eq.0)THEN      

      if (PPro.LE.0) then
         namfile='./restart/field_gridc.h5'
      else
         namfile='./cfield/field_gridc.h5'
      endif

      call h5fcreate_f(namfile,H5F_ACC_TRUNC_F,file_grid,hdf_error)

      dims_grid(1)=n1

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'xf', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xc(1:n1), &
             dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'xc', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xm(1:n1), &
             dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n2
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'yf', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, yc(1:n2), &
            dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'yc', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ym(1:n2), &
            dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n3m
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'zf', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zc(1:n3m), &
             dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'zc', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zm(1:n3m), &
             dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5fclose_f(file_grid, hdf_error)

#ifdef CALSCAL      
      if (PPro.LE.0) then
         namfile='./restart/field_gridd.h5'
      else
         namfile='./cfield/field_gridd.h5'
      endif

      call h5fcreate_f(namfile,H5F_ACC_TRUNC_F,file_grid,hdf_error)

      dims_grid(1)=n1r
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'xf', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xcr(1:n1r), &
             dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'xc', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xmr(1:n1r), &
             dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n2r
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'yf', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ycr(1:n2r), &
            dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'yc', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ymr(1:n2r), &
            dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n3mr
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'zf', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zcr(1:n3mr), &
             dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'zc', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zmr(1:n3mr), &
             dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5fclose_f(file_grid, hdf_error)
#endif

      ENDIF

      call h5close_f(hdf_error)

      return
      end subroutine write_grid_info
!============================================================
