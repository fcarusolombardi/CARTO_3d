      subroutine slab_rcd
      use param
      use local_arrays, only: dens,q1,q2,q3,dsal
#ifdef USE_BUDA
      use local_arrays, only: dens_d,q1_d,q2_d,q3_d,dsal_d
#endif
      use mpih
      use mpi_param
      use slab_param
      USE hdf5

      implicit none
      integer i,j,is
      real(DP) :: q3cc(n1m,n2m)
      real(DP), allocatable, dimension(:,:) :: q1xz,q2xz,q3xz,texz,saxz
      real(DP), allocatable, dimension(:,:) :: q1yz,q2yz,q3yz,teyz,sayz
      character filename*100
      character dataname*50

      integer :: hdf_error
      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: dset_var
      integer(HSIZE_T) :: dims(2)
      integer :: ndims

      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace
      integer(HID_T) :: plist_id

      integer(HSIZE_T), dimension(2) :: data_count
      integer(HSSIZE_T), dimension(2) :: data_offset

      integer :: comm, info

#ifdef USE_BUDA
      q1 = q1_d
      q2 = q2_d
      q3 = q3_d
      dens = dens_d
      dsal = dsal_d
#endif

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      call h5open_f(hdf_error)

      do is=1,nslab
      if(myid.eq.idslab(is)) then

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(j,i)
       do j=1,n2m
        do i=1,n1m
         q3cc(i,j) = (q3(i,j,kslab(is))+q3(i,j,kslab(is)+1))*0.5d0
        enddo
       enddo
!$OMP  END PARALLEL DO

       write(filename,'(a,i2.2,a)')'data/slab_z',is,'.h5'
       call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, &
                      file_id, hdf_error)

       ndims = 2
       dims(1) = n1m
       dims(2) = n2m

       write(dataname,'(a,i6.6)')'q1/T',nint(time)
       call h5screate_simple_f(ndims, dims, &
                              filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname),    &
                       H5T_NATIVE_DOUBLE,           &
                       filespace,                   &
                       dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE,    & 
                      q1(1:n1m,1:n2m,kslab(is)), dims, &
                      hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'q2/T',nint(time)
       call h5screate_simple_f(ndims, dims,            &
                              filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname),       &
                       H5T_NATIVE_DOUBLE,              &
                       filespace,                      &
                       dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE,    &
                      q2(1:n1m,1:n2m,kslab(is)), dims, &
                      hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'q3/T',nint(time)
       call h5screate_simple_f(ndims, dims,            &
                              filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname),       &
                       H5T_NATIVE_DOUBLE,              &
                       filespace,                      &
                       dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE,    &
                      q3cc(1:n1m,1:n2m), dims,         &
                      hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'te/T',nint(time)
       call h5screate_simple_f(ndims, dims,             &
                              filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname),        &
                       H5T_NATIVE_DOUBLE,               &
                       filespace,                       &
                       dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE,        &
                      dens(1:n1m,1:n2m,kslab(is)), dims,   &
                      hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       call h5fclose_f(file_id, hdf_error)

      endif
      enddo

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do is=1,nslab
      if(myid.eq.idslabr(is)) then

        write(filename,'(a,i2.2,a)')'data/slab_z',is,'.h5'
        call h5fopen_f (trim(filename), H5F_ACC_RDWR_F,       &
                      file_id, hdf_error)

        ndims = 2
        dims(1) = n1mr
        dims(2) = n2mr
        write(dataname,'(a,i6.6)')'sa/T',nint(time)
        call h5screate_simple_f(ndims, dims,            &
                              filespace, hdf_error)
        call h5dcreate_f(file_id, trim(dataname),       &
                       H5T_NATIVE_DOUBLE,               &
                       filespace,                       &
                       dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE,         &
                      dsal(1:n1mr,1:n2mr,kslabr(is)), dims,  &
                      hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5fclose_f(file_id, hdf_error)

      endif
      enddo

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      allocate(q1yz(1:n2m,kstart:kend))
      allocate(q2yz(1:n2m,kstart:kend))
      allocate(q3yz(1:n2m,kstart:kend))
      allocate(teyz(1:n2m,kstart:kend))
      allocate(sayz(1:n2mr,kstartr:kendr))

      write(filename,'(a)')'data/slab_xmid.h5'
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fopen_f(trim(filename), H5F_ACC_RDWR_F,      &
            file_id, hdf_error,access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      ndims = 2
      dims(1) = n2m
      dims(2) = n3m
      data_count(1) = n2m
      data_count(2) = kend-kstart+1
      data_offset(1) = 0
      data_offset(2) = kstart-1

      q1yz(1:n2m,kstart:kend) = q1(n1m/2,1:n2m,kstart:kend)
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      write(dataname,'(a,i6.6)')'q1/T',nint(time)
      call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE,   &
                      filespace, dset_var, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_var, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F,        &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,      &
                             hdf_error)
      call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE,                   &
        q1yz(1:n2m,kstart:kend), dims, hdf_error,                    &
        file_space_id = slabspace, mem_space_id = memspace,          &
        xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_var, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)

      q2yz(1:n2m,kstart:kend) = q2(n1m/2,1:n2m,kstart:kend)
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      write(dataname,'(a,i6.6)')'q2/T',nint(time)
      call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE,   &
                      filespace, dset_var, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_var, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F,       &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,     &
                             hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE,               &
        q2yz(1:n2m,kstart:kend), dims, hdf_error,                 &
        file_space_id = slabspace, mem_space_id = memspace,       &
        xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_var, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)

      q3yz(1:n2m,kstart:kend) = q3(n1m/2,1:n2m,kstart:kend)
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      write(dataname,'(a,i6.6)')'q3/T',nint(time)
      call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE,   &
                      filespace, dset_var, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_var, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F,        &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,       &
                             hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE,                  &
        q3yz(1:n2m,kstart:kend), dims, hdf_error,                    &
        file_space_id = slabspace, mem_space_id = memspace,          &
        xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_var, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)

      teyz(1:n2m,kstart:kend) = dens(n1m/2,1:n2m,kstart:kend)
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      write(dataname,'(a,i6.6)')'te/T',nint(time)
      call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE,    &
                      filespace, dset_var, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_var, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F,         &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,       &
                             hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE,                    &
        teyz(1:n2m,kstart:kend), dims, hdf_error,                      &
        file_space_id = slabspace, mem_space_id = memspace,            &
        xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_var, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)

      dims(1) = n2mr
      dims(2) = n3mr
      data_count(1) = n2mr
      data_count(2) = kendr-kstartr+1
      data_offset(1) = 0
      data_offset(2) = kstartr-1

      sayz(1:n2mr,kstartr:kendr) = dsal(n1mr/2,1:n2mr,kstartr:kendr)
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      write(dataname,'(a,i6.6)')'sa/T',nint(time)
      call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE, &
                      filespace, dset_var, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_var, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F,      &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,    &
                             hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE,                 &
        sayz(1:n2mr,kstartr:kendr), dims, hdf_error,                &
        file_space_id = slabspace, mem_space_id = memspace,         &
        xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_var, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)

      call h5fclose_f(file_id, hdf_error)

      if(allocated(q1yz)) deallocate(q1yz)
      if(allocated(q2yz)) deallocate(q2yz)
      if(allocated(q3yz)) deallocate(q3yz)
      if(allocated(teyz)) deallocate(teyz)
      if(allocated(sayz)) deallocate(sayz)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      allocate(q1xz(1:n1m,kstart:kend))
      allocate(q2xz(1:n1m,kstart:kend))
      allocate(q3xz(1:n1m,kstart:kend))
      allocate(texz(1:n1m,kstart:kend))
      allocate(saxz(1:n1mr,kstartr:kendr))

      write(filename,'(a)')'data/slab_ymid.h5'
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fopen_f(trim(filename), H5F_ACC_RDWR_F,       &
            file_id, hdf_error,access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      ndims = 2
      dims(1) = n1m
      dims(2) = n3m
      data_count(1) = n1m
      data_count(2) = kend-kstart+1
      data_offset(1) = 0
      data_offset(2) = kstart-1

      q1xz(1:n1m,kstart:kend) = q1(1:n1m,n2m/2,kstart:kend)
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      write(dataname,'(a,i6.6)')'q1/T',nint(time)
      call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE,     &
                      filespace, dset_var, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_var, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F,          &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,        &
                             hdf_error)
      call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE,                &
        q1xz(1:n1m,kstart:kend), dims, hdf_error,                 &
        file_space_id = slabspace, mem_space_id = memspace,       &
        xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_var, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)

      q2xz(1:n1m,kstart:kend) = q2(1:n1m,n2m/2,kstart:kend)
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      write(dataname,'(a,i6.6)')'q2/T',nint(time)
      call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE,   &
                      filespace, dset_var, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_var, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F,       &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,     &
                             hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE,                 &
        q2xz(1:n1m,kstart:kend), dims, hdf_error,                   &
        file_space_id = slabspace, mem_space_id = memspace,         &
        xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_var, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)

      q3xz(1:n1m,kstart:kend) = q3(1:n1m,n2m/2,kstart:kend)
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      write(dataname,'(a,i6.6)')'q3/T',nint(time)
      call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE,    &
                      filespace, dset_var, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_var, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F,         &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,       &
                             hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE,                   &
        q3xz(1:n1m,kstart:kend), dims, hdf_error,                     &
        file_space_id = slabspace, mem_space_id = memspace,           &
        xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_var, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)

      texz(1:n1m,kstart:kend) = dens(1:n1m,n2m/2,kstart:kend)
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      write(dataname,'(a,i6.6)')'te/T',nint(time)
      call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE,    &
                      filespace, dset_var, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_var, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F,         &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,       &
                             hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE,                   &
        texz(1:n1m,kstart:kend), dims, hdf_error,                     &
        file_space_id = slabspace, mem_space_id = memspace,           &
        xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_var, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)

      dims(1) = n1mr
      dims(2) = n3mr
      data_count(1) = n1mr
      data_count(2) = kendr-kstartr+1
      data_offset(1) = 0
      data_offset(2) = kstartr-1

      saxz(1:n1mr,kstartr:kendr) = dsal(1:n1mr,n2mr/2,kstartr:kendr)
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      write(dataname,'(a,i6.6)')'sa/T',nint(time)
      call h5dcreate_f(file_id, trim(dataname), H5T_NATIVE_DOUBLE,    &
                      filespace, dset_var, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_var, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F,        &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,        &
                             hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE,                    &
        saxz(1:n1mr,kstartr:kendr), dims, hdf_error,                   &
        file_space_id = slabspace, mem_space_id = memspace,            &
        xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_var, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)

      call h5fclose_f(file_id, hdf_error)

      if(allocated(q1xz)) deallocate(q1xz)
      if(allocated(q2xz)) deallocate(q2xz)
      if(allocated(q3xz)) deallocate(q3xz)
      if(allocated(texz)) deallocate(texz)
      if(allocated(saxz)) deallocate(saxz)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      call h5close_f(hdf_error)
      return
      end subroutine slab_rcd

!==============================================================
      subroutine slab_ini
      use param
      use mpih
      use mpi_param
      use slab_param
      use hdf5

      implicit none
      integer,dimension(0:numtasks-1) :: ksd,ked,ksdr,kedr,cod
      integer i,j,k
      integer kdums,kdume

      real(DP) dst,dstmin

      integer :: hdf_error
      integer(HID_T) :: file_id
      integer(HID_T) :: group_id
      integer(HID_T) :: filespace
      integer(HID_T) :: dset_var
      integer(HSIZE_T) :: dims(1)
      integer :: ndims
      character filename*100
      logical :: tag_exist

      call h5open_f(hdf_error)

!=================================================
!m  reading z locations of slab

      zslab(1) = sbl/2.d0
      zslab(2) = sbl
      zslab(3) = sbl*2.d0
      zslab(4) = alx3*0.1d0
      zslab(5) = alx3*0.2d0
      zslab(6) = alx3*0.3d0
      zslab(7) = alx3*0.4d0
      zslab(8) = alx3*0.5d0
      zslab(9) = alx3*0.6d0
      zslab(10) = alx3*0.7d0
      zslab(11) = alx3*0.8d0
      zslab(12) = alx3*0.9d0
      zslab(13) = alx3-zslab(3)
      zslab(14) = alx3-zslab(2)
      zslab(15) = alx3-zslab(1)

!=================================================
!m  determine the processor output eachslab
      kdums =1 ; kdume = 1
      do i=0,numtasks-1
        call block(n3m, numtasks, i, ksd(i), ked(i), kdums,kdume, cod)
      enddo
      cod = cod * mref3
      do i=0,numtasks-1
        ksdr(i) = sum(cod(0:i))-cod(i)+1
        kedr(i) = ksdr(i) + cod(i) - 1
      enddo

      idslab = -1
      kslab = 1
      do j=1,nslab
        dstmin = 1.d0
        do i=0,numtasks-1
          do k=ksd(i),ked(i)
            dst = dabs(zslab(j)-zm(k))
            if(dst.lt.dstmin) then
              idslab(j) = i
              kslab(j) = k
              dstmin = dst
            endif
          enddo
        enddo
      enddo

      idslabr = -1
      kslabr = 1
      do j=1,nslab
        dstmin = 1.d0
        do i=0,numtasks-1
          do k=ksdr(i),kedr(i)
            dst = dabs(zslab(j)-zmr(k))
            if(dst.lt.dstmin) then
              idslabr(j) = i
              kslabr(j) = k
              dstmin = dst
            endif
          enddo
        enddo
      enddo

      if(myid.eq.0)then
        open(51,file='fact/zloc_slab.out',status='unknown')
        do i=1,nslab
         write(51,'(f16.8,4i6)')zslab(i),kslab(i), &
               kslabr(i),idslab(i),idslabr(i)
        enddo
        close(51)
      endif

!m====================================================
!m  create hdf5 file if reset or new simulation 

      if(myid.eq.0)then

!m  reset or new simulation

      if(ireset.eq.1 .or. nread.eq.0) then

        ndims = 1
        dims(1) = 1
        do i=1,nslab
          write(filename,'(a,i2.2,a)')'data/slab_z',i,'.h5'
          call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, &
                          file_id, hdf_error)

          call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error)
          call h5dcreate_f(file_id, 'zloc', H5T_NATIVE_DOUBLE, &
                          filespace, dset_var, hdf_error)
          call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, zslab(i), &
                          dims,hdf_error)
          call h5dclose_f(dset_var, hdf_error)
          call h5sclose_f(filespace, hdf_error)


          call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error)
          call h5dcreate_f(file_id, 'Lx', H5T_NATIVE_DOUBLE, &
                          filespace, dset_var, hdf_error)
          call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, rext1, &
                          dims,hdf_error)
          call h5dclose_f(dset_var, hdf_error)
          call h5sclose_f(filespace, hdf_error)

          call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error)
          call h5dcreate_f(file_id, 'n1m', H5T_NATIVE_INTEGER, &
                          filespace, dset_var, hdf_error)
          call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n1m, &
                          dims,hdf_error)
          call h5dclose_f(dset_var, hdf_error)
          call h5sclose_f(filespace, hdf_error)

          call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error)
          call h5dcreate_f(file_id, 'n1mr', H5T_NATIVE_INTEGER, &
                          filespace, dset_var, hdf_error)
          call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n1mr, &
                          dims,hdf_error)
          call h5dclose_f(dset_var, hdf_error)
          call h5sclose_f(filespace, hdf_error)


          call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error) 
          call h5dcreate_f(file_id, 'Ly', H5T_NATIVE_DOUBLE, &
                          filespace, dset_var, hdf_error)
          call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, rext2, &
                          dims,hdf_error)
          call h5dclose_f(dset_var, hdf_error)
          call h5sclose_f(filespace, hdf_error)

          call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error)
          call h5dcreate_f(file_id, 'n2m', H5T_NATIVE_INTEGER, &
                          filespace, dset_var, hdf_error)
          call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n2m, &
                          dims,hdf_error)
          call h5dclose_f(dset_var, hdf_error)
          call h5sclose_f(filespace, hdf_error)

          call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error)
          call h5dcreate_f(file_id, 'n2mr', H5T_NATIVE_INTEGER, &
                          filespace, dset_var, hdf_error)
          call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n2mr, &
                          dims,hdf_error)
          call h5dclose_f(dset_var, hdf_error)
          call h5sclose_f(filespace, hdf_error)

          call h5gcreate_f(file_id, 'q1', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)

          call h5gcreate_f(file_id, 'q2', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)

          call h5gcreate_f(file_id, 'q3', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)

          call h5gcreate_f(file_id, 'te', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)

          call h5gcreate_f(file_id, 'sa', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)

          call h5fclose_f(file_id, hdf_error)
        enddo


        ndims = 1
        dims(1) = 1
        write(filename,'(a)')'data/slab_xmid.h5'
        call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, &
                        file_id, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'Ly', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, rext2, &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error)
        call h5dcreate_f(file_id, 'n2m', H5T_NATIVE_INTEGER, &
                          filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n2m, &
                          dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'n2mr', H5T_NATIVE_INTEGER, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n2mr, &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'n3m', H5T_NATIVE_INTEGER, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n3m, &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'n3mr', H5T_NATIVE_INTEGER, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n3mr, &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dims(1) = n2m
        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'ym', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, ym(1:n2m), &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dims(1) = n2mr
        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'ymr', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, ymr(1:n2mr), &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dims(1) = n3m
        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'zm', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, zm(1:n3m), &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dims(1) = n3mr
        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'zmr', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, zmr(1:n3mr), &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5gcreate_f(file_id, 'q1', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'q2', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'q3', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'te', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'sa', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5fclose_f(file_id, hdf_error)


        ndims = 1
        dims(1) = 1
        write(filename,'(a)')'data/slab_ymid.h5'
        call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, &
                        file_id, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'Lx', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, rext1, &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error)
        call h5dcreate_f(file_id, 'n1m', H5T_NATIVE_INTEGER, &
                          filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n1m, &
                          dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'n1mr', H5T_NATIVE_INTEGER, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n1mr, &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'n3m', H5T_NATIVE_INTEGER, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n3m, &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'n3mr', H5T_NATIVE_INTEGER, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n3mr, &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dims(1) = n1m
        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'xm', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, xm(1:n1m), &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dims(1) = n1mr
        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'xmr', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, xmr(1:n1mr), &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dims(1) = n3m
        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'zm', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, zm(1:n3m), &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dims(1) = n3mr
        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'zmr', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, zmr(1:n3mr), &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5gcreate_f(file_id, 'q1', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'q2', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'q3', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'te', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'sa', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5fclose_f(file_id, hdf_error)

      endif

!m if continue, check if file exist
      if(ireset.eq.0 .and. nread.eq.1) then

        ndims = 1
        dims(1) = 1
        do i=1,nslab
          write(filename,'(a,i2.2,a)')'data/slab_z',i,'.h5'
          inquire(file=filename,exist=tag_exist)

          if(.not.tag_exist)then
          call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, &
                          file_id, hdf_error)

          call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error)
          call h5dcreate_f(file_id, 'zloc', H5T_NATIVE_DOUBLE, &
                          filespace, dset_var, hdf_error)
          call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, zslab(i), &
                          dims,hdf_error)
          call h5dclose_f(dset_var, hdf_error)
          call h5sclose_f(filespace, hdf_error)


          call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error)
          call h5dcreate_f(file_id, 'Lx', H5T_NATIVE_DOUBLE, &
                          filespace, dset_var, hdf_error)
          call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, rext1, &
                          dims,hdf_error)
          call h5dclose_f(dset_var, hdf_error)
          call h5sclose_f(filespace, hdf_error)

          call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error)
          call h5dcreate_f(file_id, 'n1m', H5T_NATIVE_INTEGER, &
                          filespace, dset_var, hdf_error)
          call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n1m, &
                          dims,hdf_error)
          call h5dclose_f(dset_var, hdf_error)
          call h5sclose_f(filespace, hdf_error)

          call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error)
          call h5dcreate_f(file_id, 'n1mr', H5T_NATIVE_INTEGER, &
                          filespace, dset_var, hdf_error)
          call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n1mr, &
                          dims,hdf_error)
          call h5dclose_f(dset_var, hdf_error)
          call h5sclose_f(filespace, hdf_error)


          call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error)
          call h5dcreate_f(file_id, 'Ly', H5T_NATIVE_DOUBLE, &
                          filespace, dset_var, hdf_error)
          call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, rext2, &
                          dims,hdf_error)
          call h5dclose_f(dset_var, hdf_error)
          call h5sclose_f(filespace, hdf_error)

          call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error)
          call h5dcreate_f(file_id, 'n2m', H5T_NATIVE_INTEGER, &
                          filespace, dset_var, hdf_error)
          call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n2m, &
                          dims,hdf_error)
          call h5dclose_f(dset_var, hdf_error)
          call h5sclose_f(filespace, hdf_error)

          call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error)
          call h5dcreate_f(file_id, 'n2mr', H5T_NATIVE_INTEGER, &
                          filespace, dset_var, hdf_error)
          call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n2mr, &
                          dims,hdf_error)
          call h5dclose_f(dset_var, hdf_error)
          call h5sclose_f(filespace, hdf_error)

          call h5gcreate_f(file_id, 'q1', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)

          call h5gcreate_f(file_id, 'q2', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)

          call h5gcreate_f(file_id, 'q3', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)

          call h5gcreate_f(file_id, 'te', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)

          call h5gcreate_f(file_id, 'sa', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)

          call h5fclose_f(file_id, hdf_error)
          endif
        enddo


        ndims = 1
        dims(1) = 1
        write(filename,'(a)')'data/slab_xmid.h5'
        inquire(file=filename,exist=tag_exist)
        if(.not.tag_exist)then
        call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, &
                        file_id, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'Ly', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, rext2, &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error)
        call h5dcreate_f(file_id, 'n2m', H5T_NATIVE_INTEGER, &
                          filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n2m, &
                          dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'n2mr', H5T_NATIVE_INTEGER, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n2mr, &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'n3m', H5T_NATIVE_INTEGER, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n3m, &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'n3mr', H5T_NATIVE_INTEGER, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n3mr, &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dims(1) = n2m
        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'ym', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, ym(1:n2m), &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dims(1) = n2mr
        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'ymr', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, ymr(1:n2mr), &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dims(1) = n3m
        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'zm', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, zm(1:n3m), &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dims(1) = n3mr
        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'zmr', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, zmr(1:n3mr), &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5gcreate_f(file_id, 'q1', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'q2', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'q3', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'te', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'sa', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5fclose_f(file_id, hdf_error)

        endif

        ndims = 1
        dims(1) = 1
        write(filename,'(a)')'data/slab_ymid.h5'
        inquire(file=filename,exist=tag_exist)

        if(.not.tag_exist)then
        call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, &
                        file_id, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'Lx', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, rext1, &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                          filespace, hdf_error)
        call h5dcreate_f(file_id, 'n1m', H5T_NATIVE_INTEGER, &
                          filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n1m, &
                          dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'n1mr', H5T_NATIVE_INTEGER, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n1mr, &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'n3m', H5T_NATIVE_INTEGER, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n3m, &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'n3mr', H5T_NATIVE_INTEGER, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_INTEGER, n3mr, &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dims(1) = n1m
        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'xm', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, xm(1:n1m), &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dims(1) = n1mr
        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'xmr', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, xmr(1:n1mr), &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dims(1) = n3m
        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'zm', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, zm(1:n3m), &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        dims(1) = n3mr
        call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
        call h5dcreate_f(file_id, 'zmr', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
        call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, zmr(1:n3mr), &
                        dims,hdf_error)
        call h5dclose_f(dset_var, hdf_error)
        call h5sclose_f(filespace, hdf_error)

        call h5gcreate_f(file_id, 'q1', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'q2', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'q3', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'te', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'sa', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5fclose_f(file_id, hdf_error)

        endif

      endif

      endif

      call h5close_f(hdf_error)

      return
      end subroutine slab_ini

