!===========================================
      subroutine PackZ_UnpackR(aa,bb)
      use mpih
      use mpi_param
      use param, only: n1m,n2m,n3m
#ifdef USE_CUDA
      use local_arrays, only: sbuf=>rwork1_d, rbuf=>rwork2_d
#ifdef CUDECOMPAVAIL
      use cudecomp
      use cudecomp_param
#endif
#else
      use local_arrays, only: sbuf=>rwork1, rbuf=>rwork2
#endif
      use nvtx
      implicit none
      real(DP), intent(in) :: aa(1:n1m,1:n2m,kstart:kend)
      real(DP), intent(out) :: bb(1:n3m,1:n1m,jstart:jend)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:)
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr, istat
#ifdef USE_CUDA
      attributes(managed) :: aa,bb
#endif

#ifndef CUDECOMPAVAIL
      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countj(i)*n1m
        aak(i)= dj* countk(i)*n1m
      end do
      
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 

      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do

      if(allocated(sbuf) .and. size(sbuf) < n1m*n2m*dk) deallocate(sbuf)
      if(.not. allocated(sbuf)) allocate(sbuf(n1m*n2m*dk))

      if(allocated(rbuf) .and. size(rbuf) < n1m*n3m*dj) deallocate(rbuf)
      if(.not. allocated(rbuf)) allocate(rbuf(n1m*n3m*dj))

      call nvtxStartRange("PACK", 2)
#ifndef USE_CUDA
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(nc,dr,offsetr,k,j,i)
#endif
      do kk = 0, numtasks-1
        nc = dispj(kk) + 1
        dr = countj(kk)
        offsetr = offsetj(kk)
#ifdef USE_CUDA
!$cuf kernel do(3)
#endif
        do k = kstart,kend
          do j=1,dr
            do i=1,n1m
              !sbuf(nc) = aa(i,j+offsetr,k)
              !nc=nc+1
              sbuf(nc+ (i-1) + (j-1)*n1m + (k-kstart)*n1m*dr) = aa(i,j+offsetr,k)
            enddo
          enddo
        enddo
      enddo
#ifndef USE_CUDA
!$OMP END PARALLEL DO
#endif
      call nvtxEndRange

      call nvtxStartRange("MPI", 1)

#ifdef NCCLAVAIL
      call alltoallv_isendrecv_nccl(sbuf, rbuf, dispj, aaj, dispk, aak)
#else
      call alltoallv_isendrecv(sbuf, rbuf, dispj, aaj, dispk, aak)
#endif
      !call MPI_ALLTOALLV(sbuf, aaj,dispj, MDP, &
      !                   rbuf, aak,dispk, MDP, &
      !                   MPI_COMM_WORLD, ierr)
      call nvtxEndRange


      call nvtxStartRange("UNPACK", 3)
#ifndef USE_CUDA
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(nc,dz,offsetz,k,j,i)
#endif
      do kk = 0, numtasks-1
        nc = dispk(kk) + 1
        dz = countk(kk)
        offsetz = offsetk(kk)
#ifdef USE_CUDA
       !$cuf kernel do(3)
#endif
        do k = 1,dz
          do j=jstart,jend
            do i=1,n1m
              !bb(k+offsetz,i,j) = rbuf(nc)
              !nc=nc+1
              bb(k+offsetz,i,j) = rbuf(nc+ (i-1) + (j-jstart)*n1m + (k-1)*n1m*(jend-jstart+1))
            enddo
          enddo
        enddo
      enddo
#ifndef USE_CUDA
!$OMP END PARALLEL DO
#endif
      call nvtxEndRange
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
#else
      istat = cudecompTransposeYToZ(cudecomp_handle, cudecomp_grid_desc,aa, bb, cudecomp_work_d, cudecomp_dtype)
#endif

      end subroutine PackZ_UnpackR

      subroutine PackZ_UnpackR_scaled(aa,bb,scal)
      use mpih
      use mpi_param
      use param, only: n1m,n2m,n3m
      use nvtx
      implicit none
      real(DP), intent(in) :: aa(1:n1m,1:n2m,kstart:kend)
      real(DP), intent(out) :: bb(1:n3m,1:n1m,jstart:jend)
      real(DP), intent(in) :: scal(1:n3m)
      real(DP), allocatable :: sbuf(:),rbuf(:)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:)
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countj(i)*n1m
        aak(i)= dj* countk(i)*n1m
      end do
      
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 

      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
     
      if(.not. allocated(sbuf)) allocate(sbuf(0:n1m*n2m*dk-1),stat=merr)
      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      
      if(.not. allocated(rbuf)) allocate(rbuf(0:n1m*n3m*dj-1),stat=merr)
      
      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(nc,dr,offsetr,k,j,i)
      do kk = 0, numtasks-1
        nc = dispj(kk)
        dr = countj(kk)
        offsetr = offsetj(kk)
        do k = kstart,kend
          do j=1,dr
            do i=1,n1m
              !sbuf(nc) = aa(i,j+offsetr,k)
              !nc=nc+1
              sbuf(nc+ (i-1) + (j-1)*n1m + (k-kstart)*n1m*dr) = aa(i,j+offsetr,k)
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      call nvtxStartRange("MPI", 1)
      call MPI_ALLTOALLV(sbuf, aaj,dispj, MDP, &
                         rbuf, aak,dispk, MDP, &
                         MPI_COMM_WORLD, ierr)
      call nvtxEndRange
     

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(nc,dz,offsetz,k,j,i)
      do kk = 0, numtasks-1
        nc = dispk(kk)
        dz = countk(kk)
        offsetz = offsetk(kk)
        do k = 1,dz
          do j=jstart,jend
            do i=1,n1m
              !bb(k+offsetz,i,j) = rbuf(nc)
              !nc=nc+1
              bb(k+offsetz,i,j) = rbuf(nc+ (i-1) + (j-jstart)*n1m + (k-1)*n1m*(jend-jstart+1)) * scal(k+offsetz)
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      if(allocated(sbuf)) deallocate(sbuf)
      if(allocated(rbuf)) deallocate(rbuf)
     
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackZ_UnpackR_scaled

#ifdef USE_BUDA

      subroutine PackZ_UnpackR_scaled_gpu(aa,bb,scal)
      use mpih
      use mpi_param
      use param, only: n1m,n2m,n3m
      use local_arrays, only: sbuf_h=>rwork1, rbuf_h=>rwork2, sbuf=>rwork1_d, rbuf=>rwork2_d
      use cudafor
      use nvtx
      implicit none
      real(DP), device, intent(in) :: aa(1:n1m,1:n2m,kstart:kend)
      real(DP), device, intent(out) :: bb(1:n3m,1:n1m,jstart:jend)
      real(DP), device, intent(in) :: scal(1:n3m)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:), srh(:)
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc,istat
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countj(i)*n1m
        aak(i)= dj* countk(i)*n1m
      end do
      
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 

      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do

      allocate(srh(2*(numtasks-1)))
     
      if(allocated(sbuf) .and. size(sbuf) < n1m*n2m*dk) deallocate(sbuf)
      if(.not. allocated(sbuf)) allocate(sbuf(n1m*n2m*dk))

      if(allocated(sbuf_h) .and. size(sbuf_h) < n1m*n2m*dk) deallocate(sbuf_h)
      if(.not. allocated(sbuf_h)) allocate(sbuf_h(n1m*n2m*dk))
      
      if(allocated(rbuf) .and. size(rbuf) < n1m*n3m*dj) deallocate(rbuf)
      if(.not. allocated(rbuf)) allocate(rbuf(n1m*n3m*dj))
      
      if(allocated(rbuf_h) .and. size(rbuf_h) < n1m*n3m*dj) deallocate(rbuf_h)
      if(.not. allocated(rbuf_h)) allocate(rbuf_h(n1m*n3m*dj))
      

      do kk = 0, numtasks-1
        nc = dispj(kk) + 1
        dr = countj(kk)
        offsetr = offsetj(kk)
        !$cuf kernel do(3)
        do k = kstart,kend
          do j=1,dr
            do i=1,n1m
              !sbuf(nc) = aa(i,j+offsetr,k)
              !nc=nc+1
              sbuf(nc+ (i-1) + (j-1)*n1m + (k-kstart)*n1m*dr) = aa(i,j+offsetr,k)
            enddo
          enddo
        enddo
      enddo


#ifdef NCCLAVAIL
      call alltoallv_isendrecv_nccl(sbuf, rbuf, dispj, aaj, dispk, aak)
#else
      call alltoallv_isendrecv(sbuf, rbuf, dispj, aaj, dispk, aak)
#endif


      do kk = 0, numtasks-1
        nc = dispk(kk) + 1
        dz = countk(kk)
        offsetz = offsetk(kk)
        !$cuf kernel do(3)
        do k = 1,dz
          do j=jstart,jend
            do i=1,n1m
              !bb(k+offsetz,i,j) = rbuf(nc)
              !nc=nc+1
              bb(k+offsetz,i,j) = rbuf(nc+ (i-1) + (j-jstart)*n1m + (k-1)*n1m*(jend-jstart+1)) * scal(k+offsetz)
            enddo
          enddo
        enddo
      enddo

      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)

      
      deallocate(srh)
      end subroutine PackZ_UnpackR_scaled_gpu
#endif
 
!============================================
      subroutine PackR_UnpackZ(aa,bb)
      use mpih
      use mpi_param
      use param, only: n1m,n2m,n3m
#ifdef USE_CUDA
      use local_arrays, only: sbuf=>rwork1_d, rbuf=>rwork2_d
#ifdef CUDECOMPAVAIL
      use cudecomp
      use cudecomp_param
#endif
#else
      use local_arrays, only: sbuf=>rwork1, rbuf=>rwork2
#endif
      use nvtx
      implicit none
      real(DP), intent(in) :: aa(1:n3m,1:n1m,jstart:jend)
      real(DP), intent(out) :: bb(1:n1m,1:n2m,kstart:kend)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:) 
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr, istat
#ifdef USE_CUDA
      attributes(managed) :: aa, bb
#endif

#ifndef CUDECOMPAVAIL
      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countj(i)*n1m
        aak(i)= dj* countk(i)*n1m
      end do
       
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 
      
      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
      
      if(allocated(rbuf) .and. size(rbuf) < n1m*n2m*dk) deallocate(rbuf)
      if(.not. allocated(rbuf)) allocate(rbuf(n1m*n2m*dk))

      if(allocated(sbuf) .and. size(sbuf) < n1m*n3m*dj) deallocate(sbuf)
      if(.not. allocated(sbuf)) allocate(sbuf(n1m*n3m*dj))

#ifndef USE_CUDA
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(nc,dz,offsetz,k,j,i)
#endif
      do kk = 0, numtasks-1
        nc = dispk(kk) + 1
        dz= countk(kk)
        offsetz = offsetk(kk)
#ifdef USE_CUDA
        !$cuf kernel do (3)
#endif
        do k = 1,dz
          do j=jstart,jend
            do i=1,n1m
              !sbuf(nc) = aa(k+offsetz,i,j) 
              !nc=nc+1
              sbuf(nc + (i-1) + (j-jstart)*n1m + (k-1)*n1m*(jend-jstart+1)) = aa(k+offsetz,i,j) 
            enddo
          enddo
        enddo
      enddo
#ifndef USE_CUDA
!$OMP END PARALLEL DO
#endif
 
      call nvtxStartRange("MPI", 1)
#ifdef NCCLAVAIL
call alltoallv_isendrecv_nccl(sbuf, rbuf, dispk, aak, dispj, aaj)
#else
      call alltoallv_isendrecv(sbuf, rbuf, dispk, aak, dispj, aaj)
#endif

      !call MPI_ALLTOALLV(sbuf, aak,dispk, MDP, &
      !                   rbuf, aaj,dispj, MDP, &
      !                   MPI_COMM_WORLD, ierr)
      call nvtxEndRange
#ifndef USE_CUDA
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(nc,dr,offsetr,k,j,i)
#endif
      do kk = 0, numtasks-1
        nc = dispj(kk) + 1
        dr= countj(kk)
        offsetr = offsetj(kk)
#ifdef USE_CUDA
        !$cuf kernel do (3)
#endif
        do k = kstart,kend
          do j=1,dr
            do i=1,n1m
              !bb(i,j+offsetr,k) = rbuf(nc)
              !nc=nc+1
              bb(i,j+offsetr,k) = rbuf(nc + (i-1) + (j-1)*n1m + (k-kstart)*n1m*dr)
            enddo
          enddo
        enddo
      enddo
#ifndef USE_CUDA
!$OMP END PARALLEL DO
#endif

      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
#else
      istat = cudecompTransposeZToY(cudecomp_handle, cudecomp_grid_desc,aa, bb, cudecomp_work_d, cudecomp_dtype)
#endif
      
      end subroutine PackR_UnpackZ

      subroutine PackR_UnpackZ_add(aa,cc)
      use mpih
      use mpi_param
      use param, only: n1m,n2m,n3m,n1,n2
      use nvtx
      implicit none
      real(DP), intent(in) :: aa(1:n3m,1:n1m,jstart:jend)
      real(DP), intent(inout) :: cc(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo)
      real(DP), allocatable :: sbuf(:),rbuf(:)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:) 
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countj(i)*n1m
        aak(i)= dj* countk(i)*n1m
      end do
       
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 
      
      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
      
      if(.not. allocated(rbuf)) allocate(rbuf(0:n1m*n2m*dk-1),stat=merr)

      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      
      if(.not. allocated(sbuf)) allocate(sbuf(0:n1m*n3m*dj-1),stat=merr)

      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for rbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(nc,dz,offsetz,k,j,i)
      do kk = 0, numtasks-1
        nc = dispk(kk)
        dz= countk(kk)
        offsetz = offsetk(kk)
        do k = 1,dz
          do j=jstart,jend
            do i=1,n1m
              !sbuf(nc) = aa(k+offsetz,i,j) 
              !nc=nc+1
              sbuf(nc + (i-1) + (j-jstart)*n1m + (k-1)*n1m*(jend-jstart+1)) = aa(k+offsetz,i,j) 
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
 
      call nvtxStartRange("MPI", 1)
      call MPI_ALLTOALLV(sbuf, aak,dispk, MDP, &
                         rbuf, aaj,dispj, MDP, &
                         MPI_COMM_WORLD, ierr)
      call nvtxEndRange
     
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(nc,dr,offsetr,k,j,i)
      do kk = 0, numtasks-1
        nc = dispj(kk)
        dr= countj(kk)
        offsetr = offsetj(kk)
        do k = kstart,kend
          do j=1,dr
            do i=1,n1m
              !bb(i,j+offsetr,k) = rbuf(nc)
              !nc=nc+1
              cc(i,j+offsetr,k) = cc(i,j+offsetr,k) + rbuf(nc + (i-1) + (j-1)*n1m + (k-kstart)*n1m*dr)
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      if(allocated(sbuf)) deallocate(sbuf)
      if(allocated(rbuf)) deallocate(rbuf)
     
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackR_UnpackZ_add

#ifdef USE_BUDA

      subroutine PackR_UnpackZ_gpu(aa,bb,cc)
      use mpih
      use mpi_param
      use param, only: n1m,n2m,n3m
      use local_arrays, only: sbuf_h=>rwork2, rbuf_h=>rwork1, sbuf=>rwork2_d, rbuf=>rwork1_d
      use nvtx
      implicit none
      real(DP), device, intent(in) :: aa(1:n3m,1:n1m,jstart:jend)
      real(DP), device, intent(out) :: bb(1:n1m,1:n2m,kstart:kend)
      real(DP), device :: cc(1:n1m,1:n2m,kstart:kend)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:) 
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countj(i)*n1m
        aak(i)= dj* countk(i)*n1m
      end do
       
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 
      
      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
      
      if(allocated(rbuf) .and. size(rbuf) < n1m*n2m*dk) deallocate(rbuf)
      if(.not. allocated(rbuf)) allocate(rbuf(n1m*n2m*dk))

      if(allocated(rbuf_h) .and. size(rbuf_h) < n1m*n2m*dk) deallocate(rbuf_h)
      if(.not. allocated(rbuf_h)) allocate(rbuf_h(n1m*n2m*dk))

      if(allocated(sbuf) .and. size(sbuf) < n1m*n3m*dj) deallocate(sbuf)
      if(.not. allocated(sbuf)) allocate(sbuf(n1m*n3m*dj))

      if(allocated(sbuf_h) .and. size(sbuf_h) < n1m*n3m*dj) deallocate(sbuf_h)
      if(.not. allocated(sbuf_h)) allocate(sbuf_h(n1m*n3m*dj))

      do kk = 0, numtasks-1
        nc = dispk(kk) + 1
        dz= countk(kk)
        offsetz = offsetk(kk)
        !$cuf kernel do (3)
        do k = 1,dz
          do j=jstart,jend
            do i=1,n1m
              !sbuf(nc) = aa(k+offsetz,i,j) 
              !nc=nc+1
              sbuf(nc + (i-1) + (j-jstart)*n1m + (k-1)*n1m*(jend-jstart+1)) = aa(k+offsetz,i,j) 
            enddo
          enddo
        enddo
      enddo
      
      !sbuf_h = sbuf
      !call MPI_ALLTOALLV(sbuf_h, aak,dispk, MDP, &
      !                   rbuf_h, aaj,dispj, MDP, &
      !                   MPI_COMM_WORLD, ierr)
      !rbuf = rbuf_h

      
      call nvtxStartRange("MPI", 1)
#ifdef NCCLAVAIL
      call alltoallv_isendrecv_nccl(sbuf, rbuf, dispk, aak, dispj, aaj)
#else
      call alltoallv_isendrecv(sbuf, rbuf, dispk, aak, dispj, aaj)
#endif
      call nvtxEndRange

      do kk = 0, numtasks-1
        nc = dispj(kk) + 1
        dr= countj(kk)
        offsetr = offsetj(kk)
        !$cuf kernel do (3)
        do k = kstart,kend
          do j=1,dr
            do i=1,n1m
              !bb(i,j+offsetr,k) = rbuf(nc)
              !nc=nc+1
              bb(i,j+offsetr,k) = rbuf(nc + (i-1) + (j-1)*n1m + (k-kstart)*n1m*dr)
              !cc(i,j+offsetr,k) = cc(i,j+offsetr,k) + bb(i,j+offsetr,k)
            enddo
          enddo
        enddo
      enddo

      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackR_UnpackZ_gpu

      subroutine PackR_UnpackZ_add_gpu(aa,cc)
      use mpih
      use mpi_param
      use param, only: n1m,n2m,n3m,n1,n2
      use local_arrays, only: sbuf_h=>rwork2, rbuf_h=>rwork1, sbuf=>rwork2_d, rbuf=>rwork1_d
      use nvtx
      implicit none
      real(DP), device, intent(in) :: aa(1:n3m,1:n1m,jstart:jend)
      real(DP), device, intent(inout) :: cc(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:) 
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countj(i)*n1m
        aak(i)= dj* countk(i)*n1m
      end do
       
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 
      
      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
      
      if(allocated(rbuf) .and. size(rbuf) < n1m*n2m*dk) deallocate(rbuf)
      if(.not. allocated(rbuf)) allocate(rbuf(n1m*n2m*dk))

      if(allocated(rbuf_h) .and. size(rbuf_h) < n1m*n2m*dk) deallocate(rbuf_h)
      if(.not. allocated(rbuf_h)) allocate(rbuf_h(n1m*n2m*dk))

      if(allocated(sbuf) .and. size(sbuf) < n1m*n3m*dj) deallocate(sbuf)
      if(.not. allocated(sbuf)) allocate(sbuf(n1m*n3m*dj))

      if(allocated(sbuf_h) .and. size(sbuf_h) < n1m*n3m*dj) deallocate(sbuf_h)
      if(.not. allocated(sbuf_h)) allocate(sbuf_h(n1m*n3m*dj))

      do kk = 0, numtasks-1
        nc = dispk(kk) + 1
        dz= countk(kk)
        offsetz = offsetk(kk)
        !$cuf kernel do (3)
        do k = 1,dz
          do j=jstart,jend
            do i=1,n1m
              !sbuf(nc) = aa(k+offsetz,i,j) 
              !nc=nc+1
              sbuf(nc + (i-1) + (j-jstart)*n1m + (k-1)*n1m*(jend-jstart+1)) = aa(k+offsetz,i,j) 
            enddo
          enddo
        enddo
      enddo
      
      !sbuf_h = sbuf
      !call MPI_ALLTOALLV(sbuf_h, aak,dispk, MDP, &
      !                   rbuf_h, aaj,dispj, MDP, &
      !                   MPI_COMM_WORLD, ierr)
      !rbuf = rbuf_h

      
      call nvtxStartRange("MPI", 1)
#ifdef NCCLAVAIL
      call alltoallv_isendrecv_nccl(sbuf, rbuf, dispk, aak, dispj, aaj) 
#else
      call alltoallv_isendrecv(sbuf, rbuf, dispk, aak, dispj, aaj)
#endif
      call nvtxEndRange

      do kk = 0, numtasks-1
        nc = dispj(kk) + 1
        dr= countj(kk)
        offsetr = offsetj(kk)
        !$cuf kernel do (3)
        do k = kstart,kend
          do j=1,dr
            do i=1,n1m
              !bb(i,j+offsetr,k) = rbuf(nc)
              !nc=nc+1
              cc(i,j+offsetr,k) = cc(i,j+offsetr,k) + rbuf(nc + (i-1) + (j-1)*n1m + (k-kstart)*n1m*dr)

            enddo
          enddo
        enddo
      enddo

      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackR_UnpackZ_add_gpu
#endif
!===========================================
      subroutine PackZ_UnpackR_refi(aa,bb)
      use mpih
      use mpi_param
      use param, only: n1mr,n2mr,n3mr
      use nvtx
      implicit none
      real(DP), intent(in) :: aa(1:n1mr,1:n2mr,kstartr:kendr)
      real(DP), intent(out) :: bb(1:n3mr,1:n1mr,jstartr:jendr)
      real(DP), allocatable :: sbuf(:),rbuf(:)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:)
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dkr* countjr(i)*n1mr
        aak(i)= djr* countkr(i)*n1mr
      end do
      
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 

      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
     
      if(.not. allocated(sbuf)) then
         allocate(sbuf(0:n1mr*n2mr*dkr-1),stat=merr)
      end if
 
      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      
      if(.not. allocated(rbuf)) then
        allocate(rbuf(0:n1mr*n3mr*djr-1),stat=merr)
      end if
      
      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(nc,dr,offsetr,k,j,i)
      do kk = 0, numtasks-1
        nc = dispj(kk)
        dr = countjr(kk)
        offsetr = offsetjr(kk)
        do k = kstartr,kendr
          do j=1,dr
            do i=1,n1mr
              sbuf(nc) = aa(i,j+offsetr,k)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      call nvtxStartRange("MPI", 1)
      call MPI_ALLTOALLV(sbuf, aaj, dispj, MDP, &
                         rbuf, aak, dispk, MDP, &
                         MPI_COMM_WORLD, ierr)
      call nvtxEndRange
     
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(nc,dz,offsetz,k,j,i)
      do kk = 0, numtasks-1
        nc = dispk(kk)
        dz = countkr(kk)
        offsetz = offsetkr(kk)
        do k = 1,dz
          do j=jstartr,jendr
            do i=1,n1mr
              bb(k+offsetz,i,j) = rbuf(nc)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      if(allocated(sbuf)) deallocate(sbuf)
      if(allocated(rbuf)) deallocate(rbuf)
     
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackZ_UnpackR_refi
 
!============================================
      subroutine PackR_UnpackZ_refi(aa,bb)
      use mpih
      use mpi_param
      use param, only: n1mr,n2mr,n3mr
      use nvtx
      implicit none
      real(DP), intent(in) ::  aa(1:n3mr,1:n1mr,jstartr:jendr)
      real(DP), intent(out) :: bb(1:n1mr,1:n2mr,kstartr:kendr)
      real(DP), allocatable :: sbuf(:),rbuf(:)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:) 
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dkr* countjr(i)*n1mr
        aak(i)= djr* countkr(i)*n1mr
      end do
       
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 
      
      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
      
      if(.not. allocated(rbuf)) then 
       allocate(rbuf(0:n1mr*n2mr*dkr-1),stat=merr)
      end if

      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      
      if(.not. allocated(sbuf)) then 
       allocate(sbuf(0:n1mr*n3mr*djr-1),stat=merr)
      end if

      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for rbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(nc,dz,offsetz,k,j,i)
      do kk = 0, numtasks-1
        nc = dispk(kk)
        dz = countkr(kk)
        offsetz = offsetkr(kk)
        do k = 1,dz
          do j=jstartr,jendr
            do i=1,n1mr
              sbuf(nc) = aa(k+offsetz,i,j)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
      
      call nvtxStartRange("MPI", 1)
      call MPI_ALLTOALLV(sbuf, aak, dispk, MDP, &
                         rbuf, aaj, dispj, MDP, &
                         MPI_COMM_WORLD, ierr)
      call nvtxEndRange
     
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(nc,dr,offsetr,k,j,i)
      do kk = 0, numtasks-1
        nc = dispj(kk)
        dr = countjr(kk)
        offsetr = offsetjr(kk)
        do k = kstartr,kendr
          do j=1,dr
            do i=1,n1mr
              bb(i,j+offsetr,k) = rbuf(nc)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      if(allocated(sbuf)) deallocate(sbuf)
      if(allocated(rbuf)) deallocate(rbuf)
     
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackR_UnpackZ_refi
!==========================================
      subroutine PackZ_UnpackRP(aa,bb)
      use mpih
      use mpi_param
      use param, only: n1m,n2m,n3m
      use nvtx
      implicit none
      real(DP),intent(in) :: aa(1:n1m,1:n2m+2,kstart:kend)
      real(DP),intent(out) :: bb(1:n3m,1:n1m,jstartp:jendp)
      real(DP), allocatable :: sbuf(:),rbuf(:)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:) 
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countjp(i)*n1m
        aak(i)= djp* countk(i)*n1m
      end do
      
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 

      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
 
     
      if(.not. allocated(sbuf)) allocate(sbuf(0:n1m*(n2m+2)*dk-1), stat=merr)

      if(.not. allocated(rbuf)) allocate(rbuf(0:n1m*n3m*djp-1), stat=merr)

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(nc,dr,offsetr,k,j,i)
      do kk = 0, numtasks-1
        nc = dispj(kk)
        dr = countjp(kk)
        offsetr = offsetjp(kk)
        do k = kstart,kend
          do j=1,dr
            do i=1,n1m
              sbuf(nc) = aa(i,j+offsetr,k)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
      
      call MPI_ALLTOALLV(sbuf, aaj,dispj, MDP, &
                        rbuf, aak,dispk, MDP, &
                        MPI_COMM_WORLD, ierr)
     
     
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(nc,dz,offsetz,k,j,i)
      do kk = 0, numtasks-1
        nc = dispk(kk)
        dz = countk(kk)
        offsetz = offsetk(kk)
        do k = 1,dz
          do j=jstartp,jendp
            do i=1,n1m
              bb(k+offsetz,i,j) = rbuf(nc)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      if(allocated(sbuf)) deallocate(sbuf)
      if(allocated(rbuf)) deallocate(rbuf)
     
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackZ_UnpackRP

      subroutine PackZ_complex_UnpackRP(aa,bb)
      use mpih
      use mpi_param
      use param, only: n1m,n2m,n3m,m2mh
#ifdef USE_CUDA
      use local_arrays, only: sbuf=>rwork1_d, rbuf=>rwork2_d
      use local_arrays, only: rwork3
#ifdef CUDECOMPAVAIL
      use cudecomp
      use cudecomp_param
#endif
#else
      use local_arrays, only: sbuf=>rwork1, rbuf=>rwork2
#endif
      use nvtx
      implicit none
      complex(DP),intent(in) :: aa(1:m2mh,1:n1m,kstart:kend)
      real(DP),intent(out) :: bb(1:n3m,1:n1m,jstartp:jendp)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:)
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr, n2mh, istat
      real(DP) :: coefnorm
#ifdef USE_CUDA
      real(DP), pointer, device, contiguous ::  aa_split(:,:,:)
      attributes(managed) :: aa,bb
#endif

      n2mh = n2m/2 + 1
      coefnorm = 1.d0/(dble(n1m)*dble(n2m))

#ifndef CUDECOMPAVAIL
      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countjp(i)*n1m
        aak(i)= djp* countk(i)*n1m
      end do

      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1))

      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do


      if(allocated(sbuf) .and. size(sbuf) < n1m*(n2m+2)*dk) deallocate(sbuf)
      if(.not. allocated(sbuf)) allocate(sbuf(n1m*(n2m+2)*dk))

      if(allocated(rbuf) .and. size(rbuf) < n1m*n3m*djp) deallocate(rbuf)
      if(.not. allocated(rbuf)) allocate(rbuf(n1m*n3m*djp))

      do kk = 0, numtasks-1
        nc = dispj(kk) + 1
        dr = countjp(kk)
        offsetr = offsetjp(kk)
#ifdef USE_CUDA
        !$cuf kernel do (3)
#endif
        do k = kstart,kend
          do j=1,dr
            do i=1,n1m
              !sbuf(nc) = aa(i,j+offsetr,k)
              !nc=nc+1
              if (j + offsetr <= n2mh) then
                sbuf(nc + (i-1) + (j-1)*n1m + (k-kstart)*n1m*dr) = dble(aa(j+offsetr,i,k))*coefnorm
              else if (j + offsetr <= 2*n2mh) then
                sbuf(nc + (i-1) + (j-1)*n1m + (k-kstart)*n1m*dr) = dimag(aa(j+offsetr-n2mh,i,k))*coefnorm
              endif
            enddo
          enddo
        enddo
      enddo

      !sbuf_h = sbuf
      !call MPI_ALLTOALLV(sbuf_h, aaj,dispj, MDP, &
      !                  rbuf_h, aak,dispk, MDP, &
      !                  MPI_COMM_WORLD, ierr)
      !rbuf = rbuf_h

      call nvtxStartRange("MPI", 1)
#ifdef NCCLAVAIL
      call alltoallv_isendrecv_nccl(sbuf, rbuf, dispj, aaj, dispk, aak)
#else
      call alltoallv_isendrecv(sbuf, rbuf, dispj, aaj, dispk, aak)
#endif

      call nvtxEndRange

      do kk = 0, numtasks-1
        nc = dispk(kk) + 1
        dz = countk(kk)
        offsetz = offsetk(kk)
#ifdef USE_CUDA
        !$cuf kernel do (3)
#endif
        do k = 1,dz
          do j=jstartp,jendp
            do i=1,n1m
              !bb(k+offsetz,i,j) = rbuf(nc)
              !nc=nc+1
              bb(k+offsetz,i,j) = rbuf(nc + (i-1) + (j-jstartp)*n1m + (k-1)*n1m*(jendp-jstartp+1))
            enddo
          enddo
        enddo
      enddo

      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)

      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
#else
      if (allocated(rwork3) .and. size(rwork3) < 2*n2mh*n1m*(kend-kstart+1)) deallocate(rwork3)
      if (.not. allocated(rwork3)) allocate(rwork3(2*n2mh*n1m*(kend-kstart+1)))
      aa_split(1:n1m, 1:2*n2mh, kstart:kend) => rwork3

      ! Deinterleave complex and transpose before cudecomp call
      ! JR Note: Can improve if input AA is already in (j,k,i)
      ! orientation, which is the expected y-pencil axis contiguous
      ! orientation from cudecomp
#ifdef USE_CUDA
      !$cuf kernel do (3)
#endif
      do k = kstart,kend
        do j = 1, n2mh
          do i = 1,n1m
            aa_split(i,j,k) = dble(aa(j,i,k)) * coefnorm
            aa_split(i,j+n2mh,k) = dimag(aa(j,i,k)) * coefnorm
          end do
        end do
      end do
      istat = cudecompTransposeYToZ(cudecomp_handle,cudecomp_grid_desc_ph, aa_split, bb, cudecomp_work_d, cudecomp_dtype)
#endif

      end subroutine PackZ_complex_UnpackRP

#ifdef USE_BUDA
      subroutine PackZ_UnpackRP_gpu(aa,bb)
      use mpih
      use mpi_param
      use param, only: n1m,n2m,n3m
      use local_arrays, only: sbuf_h=>rwork1, rbuf_h=>rwork2, sbuf=>rwork1_d, rbuf=>rwork2_d
      use nvtx
      implicit none
      real(DP),intent(in), device :: aa(1:n1m,1:n2m+2,kstart:kend)
      real(DP),intent(out), device :: bb(1:n3m,1:n1m,jstartp:jendp)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:) 
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countjp(i)*n1m
        aak(i)= djp* countk(i)*n1m
      end do
      
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 

      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
 
     
      if(allocated(sbuf) .and. size(sbuf) < n1m*(n2m+2)*dk) deallocate(sbuf)
      if(.not. allocated(sbuf)) allocate(sbuf(n1m*(n2m+2)*dk))
      if(allocated(sbuf_h) .and. size(sbuf_h) < n1m*(n2m+2)*dk) deallocate(sbuf_h)
      if(.not. allocated(sbuf_h)) allocate(sbuf_h(n1m*(n2m+2)*dk))

      if(allocated(rbuf) .and. size(rbuf) < n1m*n3m*djp) deallocate(rbuf)
      if(.not. allocated(rbuf)) allocate(rbuf(n1m*n3m*djp))
      if(allocated(rbuf_h) .and. size(rbuf_h) < n1m*n3m*djp) deallocate(rbuf_h)
      if(.not. allocated(rbuf_h)) allocate(rbuf_h(n1m*n3m*djp))

      do kk = 0, numtasks-1
        nc = dispj(kk) + 1
        dr = countjp(kk)
        offsetr = offsetjp(kk)
        !$cuf kernel do (3)
        do k = kstart,kend
          do j=1,dr
            do i=1,n1m
              !sbuf(nc) = aa(i,j+offsetr,k)
              !nc=nc+1
              sbuf(nc + (i-1) + (j-1)*n1m + (k-kstart)*n1m*dr) = aa(i,j+offsetr,k)
            enddo
          enddo
        enddo
      enddo
      
      !sbuf_h = sbuf
      !call MPI_ALLTOALLV(sbuf_h, aaj,dispj, MDP, &
      !                  rbuf_h, aak,dispk, MDP, &
      !                  MPI_COMM_WORLD, ierr)
      !rbuf = rbuf_h

      call nvtxStartRange("MPI", 1)
#ifdef NCCLAVAIL
      call alltoallv_isendrecv_nccl(sbuf, rbuf, dispj, aaj, dispk, aak)
#else
      call alltoallv_isendrecv(sbuf, rbuf, dispj, aaj, dispk, aak)
#endif
      call nvtxEndRange
     
      do kk = 0, numtasks-1
        nc = dispk(kk) + 1
        dz = countk(kk)
        offsetz = offsetk(kk)
        !$cuf kernel do (3)
        do k = 1,dz
          do j=jstartp,jendp
            do i=1,n1m
              !bb(k+offsetz,i,j) = rbuf(nc)
              !nc=nc+1
              bb(k+offsetz,i,j) = rbuf(nc + (i-1) + (j-jstartp)*n1m + (k-1)*n1m*(jendp-jstartp+1))
            enddo
          enddo
        enddo
      enddo

      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackZ_UnpackRP_gpu

#endif
 
!============================================
      subroutine PackR_UnpackZP(aa,bb)
      use mpih
      use mpi_param
      use param, only: n1m,n2m,n3m
      use nvtx
      implicit none
      real(DP), intent(in) :: aa(1:n3m,1:n1m,jstartp:jendp)
      real(DP), intent(out) :: bb(1:n1m,1:n2m+2,kstart:kend)
      real(DP), allocatable :: sbuf(:),rbuf(:)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:) 
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countjp(i)*n1m
        aak(i)= djp* countk(i)*n1m
      end do
       
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 
      
      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
      
      if(.not. allocated(rbuf)) allocate(rbuf(0:n1m*(n2m+2)*dk-1), stat=merr)

      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      
      if(.not. allocated(sbuf)) allocate(sbuf(0:n1m*n3m*djp-1), stat=merr)

      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for rbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(nc,dz,offsetz,k,j,i)
      do kk = 0, numtasks-1
        nc = dispk(kk)
        dz = countk(kk)
        offsetz = offsetk(kk)
        do k = 1,dz
          do j=jstartp,jendp
            do i=1,n1m
              sbuf(nc) = aa(k+offsetz,i,j) 
              nc=nc+1
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
      
      call nvtxStartRange("MPI", 1)
      call MPI_ALLTOALLV(sbuf, aak,dispk, MDP, &
                        rbuf, aaj,dispj, MDP, &
                        MPI_COMM_WORLD, ierr)
      call nvtxEndRange
     
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(nc,dr,offsetr,k,j,i)
      do kk = 0, numtasks-1
        nc = dispj(kk)
        dr = countjp(kk)
        offsetr = offsetjp(kk)
        do k = kstart,kend
          do j=1,dr
            do i=1,n1m
              bb(i,j+offsetr,k) = rbuf(nc)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      if(allocated(sbuf)) deallocate(sbuf)
      if(allocated(rbuf)) deallocate(rbuf)
     
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackR_UnpackZP

      subroutine PackR_UnpackZP_complex(aa,bb)
      use mpih
      use mpi_param
      use param, only: n1m,n2m,n3m,m2mh
#ifdef USE_CUDA
      use local_arrays, only: sbuf=>rwork1_d, rbuf=>rwork2_d
      use local_arrays, only: rwork3
#ifdef CUDECOMPAVAIL
      use cudecomp
      use cudecomp_param
#endif
#else
      use local_arrays, only: sbuf=>rwork1, rbuf=>rwork2
#endif
      use nvtx
      implicit none
      real(DP), intent(in) :: aa(1:n3m,1:n1m,jstartp:jendp)
      ! JR Note: Casting complex bb array to real with 2x length
      ! JR TODO: Is m2mh safe here? Should remove m vars.
      real(DP), intent(out) :: bb(1:2*m2mh,1:n1m,kstart:kend)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:)
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr, n2mh, istat
#ifdef USE_CUDA
      real(DP), pointer, device, contiguous ::  aa_split(:,:,:)
      attributes(managed) :: aa,bb
#endif

      n2mh = n2m/2 + 1

#ifndef CUDECOMPAVAIL
      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countjp(i)*n1m
        aak(i)= djp* countk(i)*n1m
      end do

      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1))

      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do

      if(allocated(rbuf) .and. size(rbuf) < n1m*(n2m+2)*dk) deallocate(rbuf)
      if(.not. allocated(rbuf)) allocate(rbuf(n1m*(n2m+2)*dk))

      if(allocated(sbuf) .and. size(sbuf) < n1m*n3m*djp) deallocate(sbuf)
      if(.not. allocated(sbuf)) allocate(sbuf(n1m*n3m*djp))

      do kk = 0, numtasks-1
        nc = dispk(kk) + 1
        dz = countk(kk)
        offsetz = offsetk(kk)
#ifdef USE_CUDA
        !$cuf kernel do (3)
#endif
        do k = 1,dz
          do j=jstartp,jendp
            do i=1,n1m
              !sbuf(nc) = aa(k+offsetz,i,j)
              !nc=nc+1
              sbuf(nc + (i-1) + (j-jstartp)*n1m + (k-1)*n1m*(jendp-jstartp+1)) = aa(k+offsetz,i,j)
            enddo
          enddo
        enddo
      enddo

      !sbuf_h = sbuf
      !call MPI_ALLTOALLV(sbuf_h, aak,dispk, MDP, &
      !                  rbuf_h, aaj,dispj, MDP, &
      !                  MPI_COMM_WORLD, ierr)
      !rbuf = rbuf_h


      call nvtxStartRange("MPI", 1)
#ifdef NCCLAVAIL
      call alltoallv_isendrecv_nccl(sbuf, rbuf, dispk, aak, dispj, aaj)
#else
      call alltoallv_isendrecv(sbuf, rbuf, dispk, aak, dispj, aaj)
#endif


      call nvtxEndRange

      do kk = 0, numtasks-1
        nc = dispj(kk) + 1
        dr = countjp(kk)
        offsetr = offsetjp(kk)
#ifdef USE_CUDA
        !$cuf kernel do (3)
#endif
        do k = kstart,kend
          do j=1,dr
            do i=1,n1m
              !bb(i,j+offsetr,k) = rbuf(nc)
              !nc=nc+1
              !JR Note: Write to odd/even indices for real and imaginary parts respectively
              if (j+offsetr <= n2mh) then
                bb(2*(j+offsetr) - 1,i,k) = rbuf(nc + (i-1) + (j-1)*n1m + (k-kstart)*n1m*dr)
              else if (j+offsetr <= 2*n2mh) then
                bb(2*(j+offsetr-n2mh),i,k) = rbuf(nc + (i-1) + (j-1)*n1m + (k-kstart)*n1m*dr)
              endif
            enddo
          enddo
        enddo
      enddo

      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)

      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
 #else
      if (allocated(rwork3) .and. size(rwork3) < 2*n2mh*n1m*(kend-kstart+1)) deallocate(rwork3)
      if (.not. allocated(rwork3)) allocate(rwork3(2*n2mh*n1m*(kend-kstart+1)))
      aa_split(1:n1m, 1:2*n2mh, kstart:kend) => rwork3

      istat = cudecompTransposeZToY(cudecomp_handle,cudecomp_grid_desc_ph, aa, aa_split, cudecomp_work_d, cudecomp_dtype)

      ! Reinterleave complex and transpose after cudecomp call
      ! JR Note: Can improve if output bb is already in (j,k,i)
      ! orientation, which is the expected y-pencil axis contiguous
      ! orientation from cudecomp
#ifdef USE_CUDA
      !$cuf kernel do (3)
#endif
      do k = kstart,kend
        do j = 1, n2mh
          do i = 1,n1m
            bb(2*j-1,i,k) = aa_split(i,j,k)
          end do
        end do
      end do

#ifdef USE_CUDA
      !$cuf kernel do (3)
#endif
      do k = kstart,kend
        do j = n2mh+1, 2*n2mh
          do i = 1,n1m
            bb(2*(j-n2mh),i,k) =  aa_split(i,j,k)
          end do
        end do
      end do
 #endif

      end subroutine PackR_UnpackZP_complex

#ifdef USE_BUDA
      subroutine PackR_UnpackZP_gpu(aa,bb)
      use mpih
      use mpi_param
      use param, only: n1m,n2m,n3m
      use local_arrays, only: sbuf_h=>rwork2, rbuf_h=>rwork1, sbuf=>rwork2_d, rbuf=>rwork1_d
      use nvtx
      implicit none
      real(DP), intent(in), device :: aa(1:n3m,1:n1m,jstartp:jendp)
      real(DP), intent(out), device :: bb(1:n1m,1:n2m+2,kstart:kend)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:) 
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countjp(i)*n1m
        aak(i)= djp* countk(i)*n1m
      end do
       
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 
      
      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
      
      if(allocated(rbuf) .and. size(rbuf) < n1m*(n2m+2)*dk) deallocate(rbuf)
      if(.not. allocated(rbuf)) allocate(rbuf(n1m*(n2m+2)*dk))
      if(allocated(rbuf_h) .and. size(rbuf_h) < n1m*(n2m+2)*dk) deallocate(rbuf_h)
      if(.not. allocated(rbuf_h)) allocate(rbuf_h(n1m*(n2m+2)*dk))

      if(allocated(sbuf) .and. size(sbuf) < n1m*n3m*djp) deallocate(sbuf)
      if(.not. allocated(sbuf)) allocate(sbuf(n1m*n3m*djp))
      if(allocated(sbuf_h) .and. size(sbuf_h) < n1m*n3m*djp) deallocate(sbuf_h)
      if(.not. allocated(sbuf_h)) allocate(sbuf_h(n1m*n3m*djp))

      do kk = 0, numtasks-1
        nc = dispk(kk) + 1
        dz = countk(kk)
        offsetz = offsetk(kk)
        !$cuf kernel do (3)
        do k = 1,dz
          do j=jstartp,jendp
            do i=1,n1m
              !sbuf(nc) = aa(k+offsetz,i,j) 
              !nc=nc+1
              sbuf(nc + (i-1) + (j-jstartp)*n1m + (k-1)*n1m*(jendp-jstartp+1)) = aa(k+offsetz,i,j) 
            enddo
          enddo
        enddo
      enddo
      
      !sbuf_h = sbuf
      !call MPI_ALLTOALLV(sbuf_h, aak,dispk, MDP, &
      !                  rbuf_h, aaj,dispj, MDP, &
      !                  MPI_COMM_WORLD, ierr)
      !rbuf = rbuf_h

      call nvtxStartRange("MPI", 1)

#ifdef NCCLAVAIL
      call alltoallv_isendrecv_nccl(sbuf, rbuf, dispk, aak, dispj, aaj)
#else
      call alltoallv_isendrecv(sbuf, rbuf, dispk, aak, dispj, aaj)
#endif

      call nvtxEndRange
     
      do kk = 0, numtasks-1
        nc = dispj(kk) + 1
        dr = countjp(kk)
        offsetr = offsetjp(kk)
        !$cuf kernel do (3)
        do k = kstart,kend
          do j=1,dr
            do i=1,n1m
              !bb(i,j+offsetr,k) = rbuf(nc)
              !nc=nc+1
              bb(i,j+offsetr,k) = rbuf(nc + (i-1) + (j-1)*n1m + (k-kstart)*n1m*dr)
            enddo
          enddo
        enddo
      enddo

      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackR_UnpackZP_gpu

#endif
!==========================================
      subroutine alltoallv_isendrecv(sbuf, rbuf, sdisp, ssize, rdisp, rsize)
      use constants
      use mpih
      use mpi_param
!@cuf use cudafor
      implicit none
      
      real(DP), intent(in) :: sbuf(*)
      real(DP), intent(out) :: rbuf(*)
      integer, intent(in) ::  sdisp(0:numtasks-1), ssize(0:numtasks-1), rdisp(0:numtasks-1), rsize(0:numtasks-1)
      integer, allocatable :: srh(:)
      integer :: i, istat, iter, sorc, dest, sdisp_me, ssize_me, rdisp_me, rsize_me, merr
#ifdef USE_CUDA
      attributes(device) :: sbuf, rbuf
#endif

      allocate(srh(2*(numtasks-1)))

!@cuf istat = cudaDeviceSynchronize

      ! custom alltoallv
      do iter = 1, numtasks-1
        if(iand(numtasks, numtasks-1) == 0) then
          sorc = ieor(myid, iter)
        else
          sorc = mod(myid + iter, numtasks)
        endif

        call MPI_IRECV(rbuf(rdisp(sorc)+1), rsize(sorc), MDP, sorc, 0, MPI_COMM_WORLD, srh(iter), merr)
      enddo

      do iter = 1, numtasks-1
        if(iand(numtasks, numtasks-1) == 0) then
          dest = ieor(myid, iter)
        else
          dest = mod(myid + iter, numtasks)
        endif

        call MPI_ISEND(sbuf(sdisp(dest)+1), ssize(dest), MDP, dest, 0, MPI_COMM_WORLD, srh(iter + numtasks - 1), merr)
      enddo

      rdisp_me = rdisp(myid)
      sdisp_me = sdisp(myid)
      ssize_me = ssize(myid)
      rsize_me = rsize(myid)

#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do i = 1, ssize_me
        rbuf(rdisp_me + i) = sbuf(sdisp_me + i)
      enddo

      call MPI_WAITALL(2*(numtasks-1), srh, MPI_STATUSES_IGNORE, merr)

!@cuf  istat = cudaDeviceSynchronize
      deallocate(srh)
      end subroutine alltoallv_isendrecv

#ifdef NCCLAVAIL
      subroutine alltoallv_isendrecv_nccl(sbuf, rbuf, sdisp, ssize, rdisp, rsize)
        use constants
        use mpih
        use nccl
        use mpi_param
        implicit none
      
        real(DP), intent(in) :: sbuf(*)
        real(DP), intent(out) :: rbuf(*)
        integer, intent(in) ::  sdisp(0:numtasks-1), ssize(0:numtasks-1), rdisp(0:numtasks-1), rsize(0:numtasks-1)
        integer :: i, istat, iter, sorc, dest, sdisp_me, ssize_me, rdisp_me, rsize_me, merr
      #ifdef USE_CUDA
        attributes(device) :: sbuf, rbuf
      #endif
      
      !@cuf istat = cudaDeviceSynchronize
      
        nccl_result = ncclGroupStart()
        ! custom alltoallv
        do iter = 1, numtasks-1
          if(iand(numtasks, numtasks-1) == 0) then
            sorc = ieor(myid, iter)
          else
            sorc = mod(myid + iter, numtasks)
          endif
      
          nccl_result = ncclRecv(rbuf(rdisp(sorc)+1), rsize(sorc), ncclDouble, sorc, nccl_comm, 0)
        enddo
      
        do iter = 1, numtasks-1
          if(iand(numtasks, numtasks-1) == 0) then
            dest = ieor(myid, iter)
          else
            dest = mod(myid + iter, numtasks)
          endif
      
          nccl_result = ncclSend(sbuf(sdisp(dest)+1), ssize(dest), ncclDouble, dest, nccl_comm, 0)
        enddo
      
        rdisp_me = rdisp(myid)
        sdisp_me = sdisp(myid)
        ssize_me = ssize(myid)
        rsize_me = rsize(myid)
      
        nccl_result = ncclGroupEnd()
      
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i = 1, ssize_me
          rbuf(rdisp_me + i) = sbuf(sdisp_me + i)
        enddo
      
        !call MPI_WAITALL(2*(numtasks-1), srh, MPI_STATUSES_IGNORE, merr)
      
      !@cuf  istat = cudaDeviceSynchronize
      end subroutine alltoallv_isendrecv_nccl

#endif
