      subroutine profiles_rcd
      use param
      use local_arrays, only: dens, q1, q2, q3, dsal
      use mgrd_arrays, only: q1lr, q2lr, q3lr
#ifdef USE_BUDA
      use local_arrays, only: dens_d, q1_d, q2_d, q3_d, dsal_d
      use mgrd_arrays, only: q1lr_d, q2lr_d, q3lr_d
#endif
      use mpih
      use mpi_param
      USE hdf5

      implicit none
      integer i,j,k
      character filename*100
      character dataname*50

      integer :: hdf_error
      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: dset_var
      integer(HSIZE_T) :: dims(1)
      integer :: ndims
      integer :: comm, info

      real(DP),dimension(n3m)  :: uxavr,uyavr,uzavr
      real(DP),dimension(n3m)  :: uxrms,uyrms,uzrms
      real(DP),dimension(n3m)  :: terms,teavr,uzxte
      real(DP),dimension(n3mr) :: sarms,saavr,uzxsa
      real(DP),dimension(n3mr) :: uxrav,uyrav,uzrav
      real(DP),dimension(n3mr) :: uxrrm,uyrrm,uzrrm

      real(DP),dimension(n3m)  :: my_uxavr,my_uyavr,my_uzavr
      real(DP),dimension(n3m)  :: my_uxrms,my_uyrms,my_uzrms
      real(DP),dimension(n3m)  :: my_terms,my_teavr,my_uzxte
      real(DP),dimension(n3mr) :: my_sarms,my_saavr,my_uzxsa
      real(DP),dimension(n3mr) :: my_uxrav,my_uyrav,my_uzrav
      real(DP),dimension(n3mr) :: my_uxrrm,my_uyrrm,my_uzrrm

      real(DP) ucoefc,ucoefr
#ifdef USE_BUDA
      q1 = q1_d
      q2 = q2_d
      q3 = q3_d
      q1lr = q1lr_d
      q2lr = q2lr_d
      q3lr = q3lr_d
      dens = dens_d
      dsal = dsal_d
     
#endif

      call h5open_f(hdf_error)

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      ucoefc = 1.d0/dble(n1m*n2m)
      ucoefr = 1.d0/dble(n1mr*n2mr)

! calculate the profiles

      my_uxavr = 0.d0
      my_uyavr = 0.d0
      my_uzavr = 0.d0
      my_uxrms = 0.d0
      my_uyrms = 0.d0
      my_uzrms = 0.d0

      my_teavr = 0.d0
      my_terms = 0.d0
      my_uzxte = 0.d0
      do k=kstart,kend
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(i,j) &
!$OMP REDUCTION(+:my_uxavr,my_uyavr,my_uzavr) &
!$OMP REDUCTION(+:my_uxrms,my_uyrms,my_uzrms) &
!$OMP REDUCTION(+:my_terms,my_teavr,my_uzxte)
        do j=1,n2m
          do i=1,n1m
            my_uxavr(k) = my_uxavr(k) + q1(i,j,k)
            my_uyavr(k) = my_uyavr(k) + q2(i,j,k)
            my_uzavr(k) = my_uzavr(k) + q3(i,j,k)
            my_uxrms(k) = my_uxrms(k) + q1(i,j,k)**2
            my_uyrms(k) = my_uyrms(k) + q2(i,j,k)**2
            my_uzrms(k) = my_uzrms(k) + q3(i,j,k)**2
            my_teavr(k) = my_teavr(k) + dens(i,j,k)
            my_terms(k) = my_terms(k) + dens(i,j,k)**2 
            my_uzxte(k) = my_uzxte(k) + dens(i,j,k)* &
                         (q3(i,j,k)+q3(i,j,k+1))*0.5d0
          end do
        end do
!$OMP  END PARALLEL DO
      end do
      my_uxavr = my_uxavr*ucoefc
      my_uyavr = my_uyavr*ucoefc
      my_uzavr = my_uzavr*ucoefc
      my_uxrms = dsqrt(my_uxrms*ucoefc)
      my_uyrms = dsqrt(my_uyrms*ucoefc)
      my_uzrms = dsqrt(my_uzrms*ucoefc)

      my_teavr = my_teavr*ucoefc
      my_terms = dsqrt(my_terms*ucoefc)
      my_uzxte = my_uzxte*ucoefc

      my_uxrav = 0.d0
      my_uyrav = 0.d0
      my_uzrav = 0.d0
      my_uxrrm = 0.d0
      my_uyrrm = 0.d0
      my_uzrrm = 0.d0

      my_saavr = 0.d0
      my_sarms = 0.d0
      my_uzxsa = 0.d0
      do k=kstartr,kendr
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(i,j) &
!$OMP REDUCTION(+:my_uxrav,my_uyrav,my_uzrav) &
!$OMP REDUCTION(+:my_uxrrm,my_uyrrm,my_uzrrm) &
!$OMP REDUCTION(+:my_sarms,my_saavr,my_uzxsa)
        do j=1,n2mr
          do i=1,n1mr
            my_uxrav(k) = my_uxrav(k) + q1lr(i,j,k)
            my_uyrav(k) = my_uyrav(k) + q2lr(i,j,k)
            my_uzrav(k) = my_uzrav(k) + q3lr(i,j,k)
            my_uxrrm(k) = my_uxrrm(k) + q1lr(i,j,k)**2
            my_uyrrm(k) = my_uyrrm(k) + q2lr(i,j,k)**2
            my_uzrrm(k) = my_uzrrm(k) + q3lr(i,j,k)**2

            my_saavr(k) = my_saavr(k) + dsal(i,j,k)
            my_sarms(k) = my_sarms(k) + dsal(i,j,k)**2
            my_uzxsa(k) = my_uzxsa(k) + dsal(i,j,k)* &
                         (q3lr(i,j,k)+q3lr(i,j,k+1))*0.5d0
          end do
        end do
!$OMP  END PARALLEL DO
      end do

      my_uxrav = my_uxrav*ucoefr
      my_uyrav = my_uyrav*ucoefr
      my_uzrav = my_uzrav*ucoefr
      my_uxrrm = dsqrt(my_uxrrm*ucoefr)
      my_uyrrm = dsqrt(my_uyrrm*ucoefr)
      my_uzrrm = dsqrt(my_uzrrm*ucoefr)

      my_sarms = dsqrt(my_sarms*ucoefr)
      my_saavr = my_saavr*ucoefr
      my_uzxsa = my_uzxsa*ucoefr

      uxavr = 0.d0
      call MPI_REDUCE(my_uxavr,uxavr,n3m,MDP,MPI_SUM,0,comm,ierr)
      uyavr = 0.d0
      call MPI_REDUCE(my_uyavr,uyavr,n3m,MDP,MPI_SUM,0,comm,ierr)
      uzavr = 0.d0
      call MPI_REDUCE(my_uzavr,uzavr,n3m,MDP,MPI_SUM,0,comm,ierr)
      uxrms = 0.d0
      call MPI_REDUCE(my_uxrms,uxrms,n3m,MDP,MPI_SUM,0,comm,ierr)
      uyrms = 0.d0
      call MPI_REDUCE(my_uyrms,uyrms,n3m,MDP,MPI_SUM,0,comm,ierr)
      uzrms = 0.d0
      call MPI_REDUCE(my_uzrms,uzrms,n3m,MDP,MPI_SUM,0,comm,ierr)
      terms = 0.d0
      call MPI_REDUCE(my_terms,terms,n3m,MDP,MPI_SUM,0,comm,ierr)
      teavr = 0.d0
      call MPI_REDUCE(my_teavr,teavr,n3m,MDP,MPI_SUM,0,comm,ierr)
      uzxte = 0.d0
      call MPI_REDUCE(my_uzxte,uzxte,n3m,MDP,MPI_SUM,0,comm,ierr)
      sarms = 0.d0
      call MPI_REDUCE(my_sarms,sarms,n3mr,MDP,MPI_SUM,0,comm,ierr)
      saavr = 0.d0
      call MPI_REDUCE(my_saavr,saavr,n3mr,MDP,MPI_SUM,0,comm,ierr)
      uzxsa = 0.d0
      call MPI_REDUCE(my_uzxsa,uzxsa,n3mr,MDP,MPI_SUM,0,comm,ierr)

      uxrav = 0.d0
      call MPI_REDUCE(my_uxrav,uxrav,n3mr,MDP,MPI_SUM,0,comm,ierr)
      uyrav = 0.d0
      call MPI_REDUCE(my_uyrav,uyrav,n3mr,MDP,MPI_SUM,0,comm,ierr)
      uzrav = 0.d0
      call MPI_REDUCE(my_uzrav,uzrav,n3mr,MDP,MPI_SUM,0,comm,ierr)
      uxrrm = 0.d0
      call MPI_REDUCE(my_uxrrm,uxrrm,n3mr,MDP,MPI_SUM,0,comm,ierr)
      uyrrm = 0.d0
      call MPI_REDUCE(my_uyrrm,uyrrm,n3mr,MDP,MPI_SUM,0,comm,ierr)
      uzrrm = 0.d0
      call MPI_REDUCE(my_uzrrm,uzrrm,n3mr,MDP,MPI_SUM,0,comm,ierr)


! write to file from root processor
    
      IF(myid.eq.0)then

       write(filename,'(a)')'data/meanprofiles.h5'
       call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, &
                      file_id, hdf_error)

       ndims = 1
       dims(1) = n3m 

       write(dataname,'(a,i6.6)')'uxavr/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             uxavr(1:n3m), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'uyavr/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             uyavr(1:n3m), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'uzavr/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             uzavr(1:n3m), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'uxrms/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             uxrms(1:n3m), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'uyrms/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             uyrms(1:n3m), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'uzrms/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             uzrms(1:n3m), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'teavr/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             teavr(1:n3m), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'terms/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             terms(1:n3m), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'uzxte/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             uzxte(1:n3m), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)


       dims(1) = n3mr

       write(dataname,'(a,i6.6)')'saavr/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             saavr(1:n3mr), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'sarms/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             sarms(1:n3mr), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'uzxsa/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             uzxsa(1:n3mr), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'uxrav/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             uxrav(1:n3mr), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'uxrrm/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             uxrrm(1:n3mr), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'uyrav/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             uyrav(1:n3mr), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'uyrrm/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             uyrrm(1:n3mr), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'uzrav/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             uzrav(1:n3mr), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       write(dataname,'(a,i6.6)')'uzrrm/T',ntime
       call h5screate_simple_f(ndims, dims, filespace, hdf_error)
       call h5dcreate_f(file_id, trim(dataname), &
             H5T_NATIVE_DOUBLE, filespace, dset_var, hdf_error)
       call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, &
             uzrrm(1:n3mr), dims, hdf_error)
       call h5dclose_f(dset_var, hdf_error)
       call h5sclose_f(filespace, hdf_error)

       call h5fclose_f(file_id, hdf_error)
      
      ENDIF

      call h5close_f(hdf_error)

      return
      end subroutine profiles_rcd

!==============================================================
      subroutine profiles_ini
      use param
      use mpih
      use mpi_param
      use hdf5

      implicit none

      integer :: hdf_error
      integer(HID_T) :: file_id
      integer(HID_T) :: group_id
      integer(HID_T) :: filespace
      integer(HID_T) :: dset_var
      integer(HSIZE_T) :: dims(1)
      integer :: ndims
      character filename*100
      logical :: tag_exist, grp_exist

      call h5open_f(hdf_error)

!====================================================
!  create hdf5 file if reset or new simulation 

      if(myid.eq.0)then

!   reset file
      if(ireset.eq.1 .or. nread.eq.0) then

        ndims = 1
        dims(1) = 1
        write(filename,'(a)')'data/meanprofiles.h5'
        call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, &
                        file_id, hdf_error)

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

        call h5gcreate_f(file_id, 'uxavr', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uyavr', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uzavr', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uxrms', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uyrms', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uzrms', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'teavr', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'saavr', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'terms', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'sarms', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uzxte', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uzxsa', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uxrav', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uxrrm', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uyrav', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uyrrm', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uzrav', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uzrrm', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5fclose_f(file_id, hdf_error)

      endif

! continue and check if file exist
      if(ireset.eq.0 .and. nread.eq.1) then

        ndims = 1
        dims(1) = 1
        write(filename,'(a)')'data/meanprofiles.h5'
        
        inquire(file=filename,exist=tag_exist)

! if not exist, create file
        if(.not.tag_exist)then

        call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, &
                        file_id, hdf_error)

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

        call h5gcreate_f(file_id, 'uxavr', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uyavr', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uzavr', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uxrms', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uyrms', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uzrms', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'teavr', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'saavr', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'terms', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'sarms', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uzxte', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uzxsa', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uxrav', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uxrrm', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uyrav', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uyrrm', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uzrav', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5gcreate_f(file_id, 'uzrrm', group_id, hdf_error)
        call h5gclose_f(group_id, hdf_error)

        call h5fclose_f(file_id, hdf_error)

        else if(tag_exist)then

        call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, &
                      file_id, hdf_error)

        call h5lexists_f(file_id, 'zm', grp_exist, hdf_error)
        if(.not.grp_exist)then
          dims(1) = n3m
          call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
          call h5dcreate_f(file_id, 'zm', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
          call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, zm(1:n3m), &
                        dims,hdf_error)
          call h5dclose_f(dset_var, hdf_error)
          call h5sclose_f(filespace, hdf_error)
        endif

        call h5lexists_f(file_id, 'zmr', grp_exist, hdf_error)
        if(.not.grp_exist)then
          dims(1) = n3mr
          call h5screate_simple_f(ndims, dims, &
                        filespace, hdf_error)
          call h5dcreate_f(file_id, 'zmr', H5T_NATIVE_DOUBLE, &
                        filespace, dset_var, hdf_error)
          call h5dwrite_f(dset_var, H5T_NATIVE_DOUBLE, zmr(1:n3mr), &
                        dims,hdf_error)
          call h5dclose_f(dset_var, hdf_error)
          call h5sclose_f(filespace, hdf_error)
        endif

        call h5lexists_f(file_id, 'uxavr', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'uxavr', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5lexists_f(file_id, 'uyavr', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'uyavr', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5lexists_f(file_id, 'uzavr', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'uzavr', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5lexists_f(file_id, 'uxrms', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'uxrms', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5lexists_f(file_id, 'uyrms', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'uyrms', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5lexists_f(file_id, 'uzrms', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'uzrms', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5lexists_f(file_id, 'teavr', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'teavr', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5lexists_f(file_id, 'terms', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'terms', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5lexists_f(file_id, 'saavr', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'saavr', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5lexists_f(file_id, 'sarms', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'sarms', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5lexists_f(file_id, 'uzxte', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'uzxte', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5lexists_f(file_id, 'uzxsa', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'uzxsa', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5lexists_f(file_id, 'uxrav', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'uxrav', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5lexists_f(file_id, 'uyrav', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'uyrav', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5lexists_f(file_id, 'uzrav', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'uzrav', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5lexists_f(file_id, 'uxrrm', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'uxrrm', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5lexists_f(file_id, 'uyrrm', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'uyrrm', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5lexists_f(file_id, 'uzrrm', grp_exist, hdf_error)
        if(.not.grp_exist)then
          call h5gcreate_f(file_id, 'uzrrm', group_id, hdf_error)
          call h5gclose_f(group_id, hdf_error)
        endif

        call h5fclose_f(file_id, hdf_error)

        endif

      endif

      endif

      call h5close_f(hdf_error)

      return
      end subroutine profiles_ini

