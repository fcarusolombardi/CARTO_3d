!************************************************************
!     compute rms and average of velocity 
!

#ifdef USE_BUDA
      module gpu_stst
      contains
! we are  assuming that the block shape is 32 x BlocDim%y with BlocDim%y less than 32
      attributes(global) subroutine stst_kernel(kstart,kend,n1m,n2m, &
                         q1_me,q2_me,q3_me,dens_me,q3dens_me, & 
                         q1_rms,q2_rms,q3_rms,dens_rms)
       use local_arrays, only: q1=>q1_d,q2=>q2_d,q3=>q3_d,dens=>dens_d,dsal=>dsal_d
       implicit none
       integer,value:: kstart,kend,n1m,n2m
       real(DP):: q1_me(kstart:kend),q2_me(kstart:kend),q3_me(kstart:kend),dens_me(kstart:kend),q3dens_me(kstart:kend), & 
                         q1_rms(kstart:kend),q2_rms(kstart:kend),q3_rms(kstart:kend),dens_rms(kstart:kend)

       integer:: tx,ty,i,j,k
       real(DP):: usn12m,q3_loc
       real(DP):: loc_var(9),val
       real(DP),shared:: val_s(32)

       usn12m = 1.d0/(dble(n1m)*dble(n2m))
       tx = threadIdx%x
       ty = threadIdx%y
       k  = kstart+(blockIdx%x-1)


       do i=1,9
         loc_var(i)=0.d0
       end do

       do j=ty,n2m,blockDim%y
        do i=tx,n1m,blockDim%x
            !q1_me_loc = q1_me_loc + q1(i,j,k)*usn12m
            loc_var(1) = loc_var(1) + q1(i,j,k)*usn12m
            !q2_me_loc = q2_me_loc + q2(i,j,k)*usn12m
            loc_var(2) = loc_var(2) + q2(i,j,k)*usn12m
            q3_loc    = (q3(i,j,k)+q3(i,j,k+1))*0.5d0
            !q3_me_loc = q3_me_loc + q3_loc*usn12m
            loc_var(3) = loc_var(3) + q3_loc*usn12m
            !dens_me_loc = dens_me_loc + dens(i,j,k)*usn12m
            loc_var(4) = loc_var(4) + dens(i,j,k)*usn12m
            !q3dens_me_loc = q3dens_me_loc + dens(i,j,k)*q3_loc*usn12m
            loc_var(5) = loc_var(5) + dens(i,j,k)*q3_loc*usn12m
            !q1_rms_loc = q1_rms_loc + q1(i,j,k)**2*usn12m
            loc_var(6) = loc_var(6) + q1(i,j,k)*q1(i,j,k)*usn12m
            !q2_rms_loc = q2_rms_loc + q2(i,j,k)**2*usn12m
            loc_var(7) = loc_var(7) + q2(i,j,k)*q2(i,j,k)*usn12m
            !q3_rms_loc = q3_rms_loc + (q3(i,j,k)**2+q3(i,j,k+1)**2)*0.5d0*usn12m
            loc_var(8) = loc_var(8) + (q3(i,j,k)**2+q3(i,j,k+1)**2)*0.5d0*usn12m
            !dens_rms_loc = dens_rms_loc + dens(i,j,k)**2*usn12m
            loc_var(9) = loc_var(9) + dens(i,j,k)*dens(i,j,k)*usn12m
        end do
       end do
       ! Reduce inside each block
       do i=1,9
        val = __shfl_down(loc_var(i),16)
        loc_var(i) = loc_var(i) + val
        val = __shfl_down(loc_var(i),8)
        loc_var(i) = loc_var(i) + val
        val = __shfl_down(loc_var(i),4)
        loc_var(i) = loc_var(i) + val
        val = __shfl_down(loc_var(i),2)
        loc_var(i) = loc_var(i) + val
        val = __shfl_down(loc_var(i),1)
        loc_var(i) = loc_var(i) + val
        ! First thread in each warp writes to shared memory
         if (tx ==1) then
           val_s(ty)=loc_var(i)
         end if
         call syncthreads()
         ! first thread in first warp does final reduction
         if (tx==1 .and. ty ==1) then
            loc_var(i)=0.d0
            do j=1,blockDim%y
            loc_var(i)=loc_var(i)+val_s(j)
            end do
         endif
         call syncthreads()
       end do
       ! First thread in first warp writes back
       if ( tx == 1 .and. ty==1) then
        q1_me(k)     = q1_me(k)     +loc_var(1)
        q2_me(k)     = q2_me(k)     +loc_var(2)
        q3_me(k)     = q3_me(k)     +loc_var(3)
        dens_me(k)   = dens_me(k)   +loc_var(4) 
        q3dens_me(k) = q3dens_me(k) +loc_var(5)
        q1_rms(k)    = q1_rms(k)    +loc_var(6)
        q2_rms(k)    = q2_rms(k)    +loc_var(7)
        q3_rms(k)    = q3_rms(k)    +loc_var(8)
        dens_rms(k)   = dens_rms(k) +loc_var(9)
        end if


      end subroutine stst_kernel

       attributes(global) subroutine stst_kernel_r(kstartr,kendr,n1mr,n2mr, &
                         dsal_me,q3dsal_me, dsal_rms)
       use local_arrays, only: dsal=>dsal_d
       use mgrd_arrays, only: q3lr=>q3lr_d
       implicit none
       integer,value:: kstartr,kendr,n1mr,n2mr
       real(DP):: dsal_me(kstartr:kendr),q3dsal_me(kstartr:kendr),dsal_rms(kstartr:kendr)

       integer:: tx,ty,i,j,k, warpID,laneID
       real(DP):: usn12mr,q3_loc,dsal_loc
       real(DP):: loc_var(3),val
       real(DP),shared:: val_s(32)

       usn12mr = 1.d0/(dble(n1mr)*dble(n2mr))
       tx = threadIdx%x
       ty = threadIdx%y
       k  = kstartr+(blockIdx%x-1)
       warpID = ty
       laneId = tx

       do i=1,3
         loc_var(i)=0.d0
       end do

       do j=ty,n2mr,blockDim%y
        do i=tx,n1mr,blockDim%x
            dsal_loc=dsal(i,j,k)
            !dsal_me_loc = dsal_me_loc + dsal(i,j,k)*usn12mr
            !loc_var(1) = loc_var(1) + dsal(i,j,k)*usn12mr
            loc_var(1) = loc_var(1) + dsal_loc*usn12mr
            !q3dsal_me_loc = q3dsal_me_loc + dsal(i,j,k)*(q3lr(i,j,k)+q3lr(i,j,k+1))*0.5d0*usn12mr
            !loc_var(2) = loc_var(2) + dsal(i,j,k)*(q3lr(i,j,k)+q3lr(i,j,k+1))*0.5d0*usn12mr
            loc_var(2) = loc_var(2) + dsal_loc*(q3lr(i,j,k)+q3lr(i,j,k+1))*0.5d0*usn12mr
            !dsal_rms_loc = dsal_rms_loc + dsal(i,j,k)**2*usn12mr
            !loc_var(3) = loc_var(3) + dsal(i,j,k)**2*usn12mr
            loc_var(3) = loc_var(3) + dsal_loc*dsal_loc*usn12mr
        end do
       end do
       ! Reduce inside each block
       do i=1,3
        val = __shfl_down(loc_var(i),16)
        loc_var(i) = loc_var(i) + val
        val = __shfl_down(loc_var(i),8)
        loc_var(i) = loc_var(i) + val
        val = __shfl_down(loc_var(i),4)
        loc_var(i) = loc_var(i) + val
        val = __shfl_down(loc_var(i),2)
        loc_var(i) = loc_var(i) + val
        val = __shfl_down(loc_var(i),1)
        loc_var(i) = loc_var(i) + val
         if (laneID ==1) then
           val_s(warpID)=loc_var(i)
         end if
         call syncthreads()
         ! first warp does final reduction
         if (tx==1 .and. warpID ==1) then
            loc_var(i)=0.d0
            do j=1,blockDim%y
            loc_var(i)=loc_var(i)+val_s(j)
            end do
         endif
         call syncthreads()
       end do
       ! First thread in first warp write back
       if ( tx == 1 .and. warpID==1) then
        dsal_me(k)   = dsal_me(k)   +loc_var(1)
        q3dsal_me(k) = q3dsal_me(k) +loc_var(2)
        dsal_rms(k)   = dsal_rms(k) +loc_var(3)
        end if

      end subroutine stst_kernel_r

      subroutine stst_gpu
      use cudafor
      use param
      use mpi_param, only: kstart,kend,kstartr,kendr
      use stat_arrays
      implicit none

       integer    :: blocks
       type(dim3) :: threads
       threads = dim3(32, 8, 1)
       ! Main grid
       blocks =  kend - kstart +1 
       call stst_kernel<<<blocks,threads>>>( kstart,kend,n1m,n2m, &
                         q1_me_d,q2_me_d,q3_me_d,dens_me_d,q3dens_me_d, & 
                         q1_rms_d,q2_rms_d,q3_rms_d,dens_rms_d)

      ! Refined grid
       blocks =  kendr - kstartr +1 
       call stst_kernel_r<<<blocks,threads>>>( kstartr,kendr,n1mr,n2mr, &
                         dsal_me_d,q3dsal_me_d, dsal_rms_d)
       ! The stat arrays are copied back to CPU in stswr
       !q1_me     = q1_me_d
       !q2_me     = q2_me_d
       !q3_me     = q3_me_d
       !dens_me   = dens_me_d
       !q3dens_me = q3dens_me_d
       !!q1_rms    = q1_rms_d
       !q2_rms    = q2_rms_d
       !q3_rms    = q3_rms_d
       !dens_rms   = dens_rms_d
       !dsal_me   = dsal_me_d
       !q3dsal_me = q3dsal_me_d
       !dsal_rms  = dsal_rms_d
           
       
      end subroutine stst_gpu
      end module gpu_stst
#endif

      subroutine stst
      use param
      use mpi_param, only: kstart,kend,kstartr,kendr
#ifdef USE_BUDA
      use gpu_stst
      use local_arrays, only: q1=>q1_d,q2=>q2_d,q3=>q3_d,dens=>dens_d,dsal=>dsal_d
      use mgrd_arrays, only: q3lr=>q3lr_d
#else
      use local_arrays, only: q1,q2,q3,dens,dsal
      use mgrd_arrays, only: q3lr
#endif

      use stat_arrays
      use mpih
      implicit none
      real(DP) :: my_q1_rms_vol,my_q2_rms_vol,my_q3_rms_vol
      real(DP) :: my_q1q2_rms_vol, my_q1q2q3_rms_vol
      real(DP) :: q1_rms_vol,q2_rms_vol,q3_rms_vol
      real(DP) :: q1q2_rms_vol,q1q2q3_rms_vol
      real(DP) :: usn12m,usn12mr,lvol,volt
      real(DP) :: q1_me_loc, q2_me_loc, q3_me_loc, dens_me_loc, q3dens_me_loc
      real(DP) :: q1_rms_loc, q2_rms_loc, q3_rms_loc, dens_rms_loc
      real(DP) :: dsal_me_loc, dsal_rms_loc, q3dsal_me_loc
      integer :: i,j,k,istat


      usn12m = 1.d0/(dble(n1m)*dble(n2m))
      usn12mr = 1.d0/(dble(n1mr)*dble(n2mr))

      my_q1_rms_vol = 0.d0
      my_q2_rms_vol = 0.d0
      my_q3_rms_vol = 0.d0
      my_q1q2_rms_vol = 0.d0
      my_q1q2q3_rms_vol = 0.d0

#ifdef USE_BUDA
     ! Statistics on planes are computed with custom kernels
      call stst_gpu()
     ! Volume statistics are computed with cuf kernel
      !$cuf kernel do(3) <<<*,*>>>
      do k=kstart,kend
         do j=1,n2m
          do i=1,n1m
            lvol = g3rm_d(k)
            my_q1_rms_vol = my_q1_rms_vol + lvol*q1(i,j,k)**2
            my_q2_rms_vol = my_q2_rms_vol + lvol*q2(i,j,k)**2
            my_q3_rms_vol = my_q3_rms_vol  &
               + lvol*(q3(i,j,k)**2+q3(i,j,k+1)**2)*0.5d0
            my_q1q2_rms_vol = my_q1q2_rms_vol + lvol* &
                (q1(i,j,k)**2+q2(i,j,k)**2)
            my_q1q2q3_rms_vol = my_q1q2q3_rms_vol + lvol* &
                ( q1(i,j,k)**2+q2(i,j,k)**2 &
                 +(q3(i,j,k)**2+q3(i,j,k+1)**2)*0.5d0)
          end do
        end do
       end do

#else
      do k=kstart,kend
        lvol = g3rm(k)
!$OMP  PARALLEL DO &
!$OMP  DEFAULT(SHARED) &
!$OMP  PRIVATE(i,j) &
!$OMP  REDUCTION(+:my_q1_rms_vol,my_q2_rms_vol,my_q3_rms_vol) &
!$OMP  REDUCTION(+:my_q1q2_rms_vol,my_q1q2q3_rms_vol) &
!$OMP  REDUCTION(+:q1_me,q2_me,q3_me,dens_me,q3dens_me) &
!$OMP  REDUCTION(+:q1_rms,q2_rms,q3_rms,dens_rms)
        do j=1,n2m
          do i=1,n1m
            q1_me(k) = q1_me(k) + q1(i,j,k)*usn12m
            q2_me(k) = q2_me(k) + q2(i,j,k)*usn12m
            q3_me(k) = q3_me(k)  &
               + (q3(i,j,k)+q3(i,j,k+1))*0.5d0*usn12m
            dens_me(k) = dens_me(k) + dens(i,j,k)*usn12m

            q1_rms(k) = q1_rms(k) + q1(i,j,k)**2*usn12m
            q2_rms(k) = q2_rms(k) + q2(i,j,k)**2*usn12m
            q3_rms(k) = q3_rms(k)  &
               + (q3(i,j,k)**2+q3(i,j,k+1)**2)*0.5d0*usn12m
            dens_rms(k) = dens_rms(k) + dens(i,j,k)**2*usn12m
            q3dens_me(k) = q3dens_me(k) + dens(i,j,k)* &
                (q3(i,j,k)+q3(i,j,k+1))*0.5d0*usn12m

            my_q1_rms_vol = my_q1_rms_vol + lvol*q1(i,j,k)**2
            my_q2_rms_vol = my_q2_rms_vol + lvol*q2(i,j,k)**2
            my_q3_rms_vol = my_q3_rms_vol  &
               + lvol*(q3(i,j,k)**2+q3(i,j,k+1)**2)*0.5d0
            my_q1q2_rms_vol = my_q1q2_rms_vol + lvol* &
                (q1(i,j,k)**2+q2(i,j,k)**2)
            my_q1q2q3_rms_vol = my_q1q2q3_rms_vol + lvol* &
                ( q1(i,j,k)**2+q2(i,j,k)**2 &
                 +(q3(i,j,k)**2+q3(i,j,k+1)**2)*0.5d0)
          end do
        end do
!$OMP  END PARALLEL DO
      end do

      do k=kstartr,kendr
!$OMP  PARALLEL DO &
!$OMP  DEFAULT(SHARED) &
!$OMP  PRIVATE(i,j) &
!$OMP  REDUCTION(+:dsal_me,dsal_rms,q3dsal_me) 
        do j=1,n2mr
          do i=1,n1mr
            dsal_me(k) = dsal_me(k) + dsal(i,j,k)*usn12mr
            dsal_rms(k) = dsal_rms(k) + dsal(i,j,k)**2*usn12mr
            q3dsal_me(k) = q3dsal_me(k) + dsal(i,j,k)* &
               (q3lr(i,j,k)+q3lr(i,j,k+1))*0.5d0*usn12mr
          end do
        end do
!$OMP  END PARALLEL DO
      end do
#endif

      q1_rms_vol = 0.d0
      q2_rms_vol = 0.d0
      q3_rms_vol = 0.d0
      q1q2q3_rms_vol = 0.d0
      q1q2_rms_vol = 0.d0
      call MPI_REDUCE(my_q1_rms_vol,q1_rms_vol,1,MDP,MPI_SUM,0, MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_q2_rms_vol,q2_rms_vol,1,MDP,MPI_SUM,0, MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_q3_rms_vol,q3_rms_vol,1,MDP,MPI_SUM,0, MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_q1q2q3_rms_vol,q1q2q3_rms_vol,1,MDP, MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_q1q2_rms_vol,q1q2_rms_vol,1,MDP,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if(myid.eq.0) then
        volt = 1.d0/(dble(n1m)*dble(n2m)*dble(n3m))
        q1_rms_vol=dsqrt(q1_rms_vol*volt)/nu
        q2_rms_vol=dsqrt(q2_rms_vol*volt)/nu
        q3_rms_vol=dsqrt(q3_rms_vol*volt)/nu
        q1q2q3_rms_vol=dsqrt(q1q2q3_rms_vol*volt)/nu
        q1q2_rms_vol=dsqrt(q1q2_rms_vol*volt)/nu
        write(94,768) time,q1_rms_vol,q2_rms_vol,q3_rms_vol, q1q2q3_rms_vol, q1q2_rms_vol
      endif
768   format(1x,f10.4,5(1x,ES20.8))

      return
      end

!    
!***********************************************************************
      subroutine ststwr
      use mpih
      use param
      use mpi_param, only: kstart,kend,kstartr,kendr
      use stat_arrays
      use hdf5

      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_q1me
      integer(HID_T) :: dset_q2me
      integer(HID_T) :: dset_q3me
      integer(HID_T) :: dset_densme
      integer(HID_T) :: dset_dsalme

      integer(HID_T) :: dset_q1rms
      integer(HID_T) :: dset_q2rms
      integer(HID_T) :: dset_q3rms
      integer(HID_T) :: dset_densrms
      integer(HID_T) :: dset_dsalrms

      integer(HID_T) :: dset_q3dens
      integer(HID_T) :: dset_q3dsal

      integer(HID_T) :: dset_kindiss
      integer(HID_T) :: dset_thediss
      integer(HID_T) :: dset_saldiss
      integer(HID_T) :: dset_kedissr

      integer(HID_T) :: file_grid
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid

      integer(HID_T) :: plist_id

      integer(HSIZE_T) :: dims(1)
      integer(HSIZE_T) :: dims_grid(1)
      integer(HSIZE_T), dimension(1) :: data_count  
      integer(HSSIZE_T), dimension(1) :: data_offset 

      integer :: comm, info
      integer :: ndims

      character(30) flname

#ifdef USE_BUDA
       q1_me     = q1_me_d
       q2_me     = q2_me_d
       q3_me     = q3_me_d
       dens_me   = dens_me_d
       q3dens_me = q3dens_me_d
       q1_rms    = q1_rms_d
       q2_rms    = q2_rms_d
       q3_rms    = q3_rms_d
       dens_rms   = dens_rms_d

       dsal_me   = dsal_me_d
       q3dsal_me = q3dsal_me_d
       dsal_rms  = dsal_rms_d

       dissuc = dissuc_d
       disste = disste_d
       dissur = dissur_d
       disssa = disssa_d

#endif USE_BUDA
      call h5open_f(hdf_error)

!   Sort out MPI definitions
    comm = MPI_COMM_WORLD
    info = MPI_INFO_NULL

!   Form the name of the file
    flname = 'statistics_qte.h5'

!   Set offsets and element counts
 
    ndims = 1

    dims(1)=n3m

    data_count(1) = kend-kstart+1

    data_offset(1) = kstart-1

!   Open file and create dataspace

    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
    call h5fcreate_f(flname, H5F_ACC_TRUNC_F, file_id, hdf_error,  access_prp=plist_id)
    call h5pclose_f(plist_id, hdf_error)

!   Select hyperslab  and then write it

!   Q1me
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vx_mean', H5T_NATIVE_DOUBLE, &
                       filespace, dset_q1me, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q1me, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dwrite_f(dset_q1me, H5T_NATIVE_DOUBLE, &
              q1_me(kstart:kend), dims, hdf_error, &
              file_space_id = slabspace, mem_space_id = memspace,  &
              xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q1me, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)

!   Q2me
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vy_mean', H5T_NATIVE_DOUBLE, &
                       filespace, dset_q2me, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q2me, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dwrite_f(dset_q2me, H5T_NATIVE_DOUBLE, &
              q2_me(kstart:kend), dims, hdf_error, &
              file_space_id = slabspace, mem_space_id = memspace, &
              xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q2me, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)

!   Q3me
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vz_mean', H5T_NATIVE_DOUBLE, &
                       filespace, dset_q3me, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q3me, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dwrite_f(dset_q3me, H5T_NATIVE_DOUBLE, &
         q3_me(kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q3me, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)

!   densme
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'dens_mean', H5T_NATIVE_DOUBLE, &
                       filespace, dset_densme, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_densme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_densme, H5T_NATIVE_DOUBLE, &
         dens_me(kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_densme, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)

!   Q1rms
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vx_rms', H5T_NATIVE_DOUBLE, &
                       filespace, dset_q1rms, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q1rms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,&
                              hdf_error)
      call h5dwrite_f(dset_q1rms, H5T_NATIVE_DOUBLE,&
         q1_rms(kstart:kend), dims, &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q1rms, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)

!   Q2rms
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vy_rms', H5T_NATIVE_DOUBLE, &
                       filespace, dset_q2rms, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q2rms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dwrite_f(dset_q2rms, H5T_NATIVE_DOUBLE, &
         q2_rms(kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q2rms, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)

!   Q3rms
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vz_rms', H5T_NATIVE_DOUBLE, &
                       filespace, dset_q3rms, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q3rms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dwrite_f(dset_q3rms, H5T_NATIVE_DOUBLE, &
         q3_rms(kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q3rms, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)

!   densrms
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'dens_rms', H5T_NATIVE_DOUBLE, &
                       filespace, dset_densrms, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_densrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dwrite_f(dset_densrms, H5T_NATIVE_DOUBLE, &
         dens_rms(kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_densrms, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)

!   kindiss
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'kindiss_mean', H5T_NATIVE_DOUBLE, &
                       filespace, dset_kindiss, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_kindiss, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dwrite_f(dset_kindiss, H5T_NATIVE_DOUBLE, &
         dissuc(kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_kindiss, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)

!   thediss
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'thediss_mean', H5T_NATIVE_DOUBLE, &
                       filespace,dset_thediss, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_thediss, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dwrite_f(dset_thediss, H5T_NATIVE_DOUBLE, &
         disste(kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_thediss, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)

!     q3 x dens
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'q3dens_mean', H5T_NATIVE_DOUBLE, &
                       filespace, dset_q3dens, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q3dens, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dwrite_f(dset_q3dens, H5T_NATIVE_DOUBLE, &
         q3dens_me(kstart:kend), dims, &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q3dens, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)

!   Close properties and datasets
      call h5fclose_f(file_id, hdf_error)

!
!  Salinity statistics
!
      ndims = 1

      dims(1)=n3mr

      data_count(1) = kendr-kstartr+1

      data_offset(1) = kstartr-1

      flname = 'statistics_sal.h5'
!   Open file and create dataspace
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(flname, H5F_ACC_TRUNC_F, file_id, &
                       hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

!   dsalme
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'dsal_mean', H5T_NATIVE_DOUBLE, &
                       filespace, dset_dsalme, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_dsalme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dwrite_f(dset_dsalme, H5T_NATIVE_DOUBLE, &
              dsal_me(kstartr:kendr), dims, hdf_error, &
              file_space_id = slabspace, mem_space_id = memspace, &
              xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_dsalme, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)

!   dsalrms
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'dsal_rms', H5T_NATIVE_DOUBLE, &
                       filespace, dset_dsalrms, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_dsalrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dwrite_f(dset_dsalrms, H5T_NATIVE_DOUBLE, &
              dsal_rms(kstartr:kendr), dims, hdf_error, &
              file_space_id = slabspace, mem_space_id = memspace, &
              xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_dsalrms, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)

!   saldiss
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'saldiss_mean', H5T_NATIVE_DOUBLE, &
                       filespace, dset_saldiss, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_saldiss, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dwrite_f(dset_saldiss, H5T_NATIVE_DOUBLE, &
              disssa(kstartr:kendr), dims, hdf_error, &
              file_space_id = slabspace, mem_space_id = memspace, &
              xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_saldiss, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)

!   kedissr
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'kedissr_mean', H5T_NATIVE_DOUBLE, &
                       filespace, dset_kedissr, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_kedissr, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
     call h5dwrite_f(dset_kedissr, H5T_NATIVE_DOUBLE, &
              dissur(kstartr:kendr), dims, hdf_error, &
              file_space_id = slabspace, mem_space_id = memspace, &
              xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_kedissr, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)

!  q3 x dsal
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'q3dsal_mean', H5T_NATIVE_DOUBLE, &
                       filespace, dset_q3dsal, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q3dsal, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dwrite_f(dset_q3dsal, H5T_NATIVE_DOUBLE, &
              q3dsal_me(kstartr:kendr), dims, hdf_error,  &
              file_space_id = slabspace, mem_space_id = memspace, &
              xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q3dsal, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)

!   Close properties and datasets
      call h5fclose_f(file_id, hdf_error)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! WRITE MASTER FILE
      IF (myid.eq.0) then

      ndims=1
      flname='statistics_master.h5'
      call h5fcreate_f(flname,H5F_ACC_TRUNC_F, file_grid, hdf_error)

!   Write amount of averages 
      dims_grid(1)=1
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'averaging_time', H5T_NATIVE_INTEGER, &
                       dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_INTEGER, timeint_cdsp, &
                      dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!   Write Rayleigh number
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'RaT', H5T_NATIVE_DOUBLE, &
                       dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, rat, &
                      dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'RaS', H5T_NATIVE_DOUBLE, &
                       dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ras, &
                      dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!   Write Prandtl number
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'PrT', H5T_NATIVE_DOUBLE, &
                       dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, prt, &
                      dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!   Write Lewis number
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'Lew', H5T_NATIVE_DOUBLE, &
                       dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, lew, &
                      dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!   Write density ratio 
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'Rho_p', H5T_NATIVE_DOUBLE, &
                       dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, rhop, &
                      dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!   Write the grid information 
      dims_grid(1)=n2m
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'Y_cordin', H5T_NATIVE_DOUBLE, &
                       dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ym(1:n2m), &
                      dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n3m
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'Z_cordin', H5T_NATIVE_DOUBLE, &
                       dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zm(1:n3m), &
                      dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!   Close file
      call h5fclose_f(file_grid, hdf_error)

!  grid for dense field

      call h5fcreate_f(flname,H5F_ACC_TRUNC_F, file_grid, hdf_error)

!   Write amount of averages 
      dims_grid(1)=1
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'averaging_time', H5T_NATIVE_INTEGER, &
                       dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_INTEGER, timeint_cdsp, &
                      dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!   Write Rayleigh number
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'RaT', H5T_NATIVE_DOUBLE, &
                       dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, rat, &
                      dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(ndims, dims_grid,dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'RaS', H5T_NATIVE_DOUBLE, &
                       dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ras, &
                      dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!   Write Prandtl number
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'PrT', H5T_NATIVE_DOUBLE, &
                       dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, prt, &
                      dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!   Write Lewis number
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'Lew', H5T_NATIVE_DOUBLE, &
                       dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, lew, &
                      dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!   Write density ratio 
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'Rho_p', H5T_NATIVE_DOUBLE, &
                       dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, rhop, &
                      dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!   Write the grid information 
      dims_grid(1)=n2mr
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'Y_cordin', H5T_NATIVE_DOUBLE, &
                       dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ymr(1:n2mr), &
                      dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n3mr
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'Z_cordin', H5T_NATIVE_DOUBLE, &
                       dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zmr(1:n3mr), &
                      dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!   Close file
      call h5fclose_f(file_grid, hdf_error)

      ENDIF

      call h5close_f(hdf_error)
      return  
      end
! 
!***********************************************************************

      subroutine initstst
      use param
      use mpi_param, only: kstart,kend,kstartr,kendr
      use stat_arrays
      use mpih
      use hdf5

      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: memspace
      integer(HID_T) :: slabspace

      integer(HID_T) :: dset_q1me
      integer(HID_T) :: dset_q2me
      integer(HID_T) :: dset_q3me
      integer(HID_T) :: dset_densme
      integer(HID_T) :: dset_dsalme

      integer(HID_T) :: dset_q1rms
      integer(HID_T) :: dset_q2rms
      integer(HID_T) :: dset_q3rms
      integer(HID_T) :: dset_densrms
      integer(HID_T) :: dset_dsalrms

      integer(HID_T) :: dset_q3dens
      integer(HID_T) :: dset_q3dsal

      integer(HID_T) :: dset_kindiss
      integer(HID_T) :: dset_thediss
      integer(HID_T) :: dset_saldiss
      integer(HID_T) :: dset_kedissr

      integer(HSIZE_T) :: dims(1)

      integer(HSIZE_T) :: dims_grid(1)
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid

      integer(HID_T) :: plist_id
      integer(HID_T) :: plist_full
      integer(HSIZE_T), dimension(1) :: data_count  
      integer(HSSIZE_T), dimension(1) :: data_offset 

      integer :: comm, info
      integer :: ndims

      character flname*30

!   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      call h5open_f(hdf_error)

      IF(ireset.eq.0)then

!   Form the name of the file

      flname = 'statistics_qte.h5'

!   Set offsets and element counts
   
      ndims = 1

      dims(1)=n3m

      data_count(1)=kend-kstart+1

      data_offset(1) = kstart-1
         
!   Open file and create dataspace

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_full, hdf_error)
      call h5pset_fapl_mpio_f(plist_full, comm, info, hdf_error)
      call h5fopen_f(flname, H5F_ACC_RDONLY_F, file_id, hdf_error, &
                       access_prp=plist_full)
      call h5pclose_f(plist_full,hdf_error)

!   Q1me

      call h5dopen_f(file_id, 'Vx_mean', dset_q1me, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q1me, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dread_f(dset_q1me, H5T_NATIVE_DOUBLE, &
         q1_me(kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5dclose_f(dset_q1me, hdf_error)

!   Q2me

      call h5dopen_f(file_id, 'Vy_mean', dset_q2me, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q2me, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dread_f(dset_q2me, H5T_NATIVE_DOUBLE, &
         q2_me(kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5dclose_f(dset_q2me, hdf_error)

!   Q3me

      call h5dopen_f(file_id, 'Vz_mean', dset_q3me, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q3me, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dread_f(dset_q3me, H5T_NATIVE_DOUBLE, &
         q3_me(kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5dclose_f(dset_q3me, hdf_error)

!   densme

      call h5dopen_f(file_id, 'dens_mean', dset_densme, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_densme, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dread_f(dset_densme, H5T_NATIVE_DOUBLE, &
         dens_me(kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5dclose_f(dset_densme, hdf_error)

!   Q1rms

      call h5dopen_f(file_id, 'Vx_rms', dset_q1rms, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q1rms, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dread_f(dset_q1rms, H5T_NATIVE_DOUBLE, &
         q1_rms(kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5dclose_f(dset_q1rms, hdf_error)

!   Q2rms

      call h5dopen_f(file_id, 'Vy_rms', dset_q2rms, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q2rms, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dread_f(dset_q2rms, H5T_NATIVE_DOUBLE, &
         q2_rms(kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5dclose_f(dset_q2rms, hdf_error)

!   Q3rms

      call h5dopen_f(file_id, 'Vz_rms', dset_q3rms, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q3rms, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dread_f(dset_q3rms, H5T_NATIVE_DOUBLE, &
         q3_rms(kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5dclose_f(dset_q3rms, hdf_error)

!   densrms

      call h5dopen_f(file_id, 'dens_rms', dset_densrms, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_densrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dread_f(dset_densrms, H5T_NATIVE_DOUBLE, &
         dens_rms(kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5dclose_f(dset_densrms, hdf_error)

!   kindiss

      call h5dopen_f(file_id, 'kindiss_mean', dset_kindiss, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_kindiss, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dread_f(dset_kindiss, H5T_NATIVE_DOUBLE, &
         dissuc(kstart:kend), dims, &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5dclose_f(dset_kindiss, hdf_error)

!   thediss

      call h5dopen_f(file_id, 'thediss_mean', dset_thediss, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_thediss, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dread_f(dset_thediss, H5T_NATIVE_DOUBLE, &
         disste(kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5dclose_f(dset_thediss, hdf_error)

!   q3 x dens

      call h5dopen_f(file_id, 'q3dens_mean', dset_q3dens, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q3dens, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dread_f(dset_q3dens, H5T_NATIVE_DOUBLE, &
         q3dens_me(kstart:kend), dims, &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5dclose_f(dset_q3dens, hdf_error)

!   Close file

      call h5fclose_f(file_id, hdf_error)

!   Set offsets and element counts

      ndims = 1

      dims(1)=n3mr

      data_count(1)=kendr-kstartr+1

      data_offset(1) = kstartr-1

      flname = 'statistics_sal.h5'

!   Open file and create dataspace

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_full, hdf_error)
      call h5pset_fapl_mpio_f(plist_full, comm, info, hdf_error)

      call h5fopen_f(flname, H5F_ACC_RDONLY_F, file_id, hdf_error, &
                       access_prp=plist_full)

      call h5pclose_f(plist_full,hdf_error)

!   dsalme

      call h5dopen_f(file_id, 'dsal_mean', dset_dsalme, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_dsalme, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dread_f(dset_dsalme, H5T_NATIVE_DOUBLE, &
         dsal_me(kstartr:kendr), dims, &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5dclose_f(dset_dsalme, hdf_error)

!   dsalrms

      call h5dopen_f(file_id, 'dsal_rms', dset_dsalrms, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_dsalrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dread_f(dset_dsalrms, H5T_NATIVE_DOUBLE, &
         dsal_rms(kstartr:kendr), dims, &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5dclose_f(dset_dsalrms, hdf_error)

!   saldiss

      call h5dopen_f(file_id, 'saldiss_mean', dset_saldiss, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_saldiss, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dread_f(dset_saldiss, H5T_NATIVE_DOUBLE, &
         disssa(kstartr:kendr), dims, &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5dclose_f(dset_saldiss, hdf_error)

!   kedissr

      call h5dopen_f(file_id, 'kedissr_mean', dset_kedissr, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_kedissr, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dread_f(dset_kedissr, H5T_NATIVE_DOUBLE, &
         dissur(kstartr:kendr), dims, &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5dclose_f(dset_kedissr, hdf_error)

! q3 x dsal

      call h5dopen_f(file_id, 'q3dsal_mean', dset_q3dsal, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q3dsal, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dread_f(dset_q3dsal, H5T_NATIVE_DOUBLE, &
         q3dsal_me(kstartr:kendr), dims, &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5dclose_f(dset_q3dsal, hdf_error)

!  close file
      call h5fclose_f(file_id, hdf_error)

      flname = 'statistics_master.h5'
!   Read the grid & statistics information
      if(myid.eq.0)then
        ndims=1
        dims_grid(1)=1
        call h5fopen_f(flname, H5F_ACC_RDONLY_F, file_id, hdf_error)
        call h5screate_simple_f(ndims,dims_grid,dspace_grid,hdf_error)
        call h5dopen_f(file_id, 'averaging_time', dset_grid, hdf_error)
        call h5dread_f(dset_grid, H5T_NATIVE_INTEGER, timeint_cdsp, &
                       dims_grid,hdf_error)
        call h5dclose_f(dset_grid, hdf_error)
        call h5sclose_f(dspace_grid, hdf_error)
        call h5fclose_f(file_id, hdf_error)
      endif
      call MPI_BCAST(timeint_cdsp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      ENDIF

!m====================================================
!    if reset time, reset statistic averaging count

      IF(ireset.eq.1)then

        timeint_cdsp = 0

        q1_me = 0.0d0
        q2_me = 0.0d0
        q3_me = 0.0d0
        dens_me = 0.0d0
        q1_rms = 0.0d0
        q2_rms = 0.0d0
        q3_rms = 0.0d0
        dens_rms = 0.0d0
        disste = 0.0d0
        dissuc = 0.0d0
        q3dens_me = 0.d0
        dsal_me = 0.0d0
        dsal_rms = 0.0d0
        q3dsal_me = 0.d0
        disssa = 0.0d0
        dissur = 0.d0

      ENDIF

      call h5close_f(hdf_error)

      return
      end
