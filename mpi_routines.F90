!===============================================
      subroutine block(n, p, irank, istart, iend, istartp, iendp, blcsz)
      use constants
      use param, only : n3m
      implicit none
      integer,intent(in) :: n,p,irank
      integer,intent(out) :: istart,iend
      integer,intent(out) :: istartp,iendp
      integer :: i
      integer,dimension(0:p-1),intent(out) :: blcsz
      
      do i=0,p-1
        blcsz(i) = floor(real((n+p-i-1)/p), kind=DP)
      enddo
      istart = sum(blcsz(0:irank))-blcsz(irank)+1
      iend = istart+blcsz(irank)-1

      if(istart.eq.1) then
            istartp = 2
      else
            istartp = istart
      end if

      if(iend.eq.n3m) then
            iendp = n3m-1
      else
            iendp = iend
      end if

      end subroutine block
!=================================================           
      subroutine mpi_workdistribution
      use param
      use mpih 
      use mpi_param
      implicit none
      integer :: i,jdums,jdume

      jdums= 1 ; jdume=1
      
      if(.not. allocated(countj)) allocate(countj(0:numtasks-1))
      if(.not. allocated(countjp)) allocate(countjp(0:numtasks-1))
      if(.not. allocated(countk)) allocate(countk(0:numtasks-1))

      if(.not. allocated(countjr)) allocate(countjr(0:numtasks-1))
      if(.not. allocated(countkr)) allocate(countkr(0:numtasks-1))

!   For PERIODIC pressure solver
      call block(n2m+2, numtasks, myid, jstartp, jendp,jdums,jdume,countjp)
      djp=jendp-jstartp+1

      call block(n2m, numtasks, myid, jstart, jend,jdums,jdume, countj)
      dj=jend-jstart+1
     
!     kstartp and kendp included for inflow-outflow routines 
      call block(n3m,numtasks,myid,kstart,kend,kstartp,kendp,countk)
      dk=kend-kstart+1

      countjr = countj*mref2
      countkr = countk*mref3
      jstartr = sum(countjr(0:myid))-countjr(myid)+1
      jendr = jstartr+countjr(myid)-1
      kstartr = sum(countkr(0:myid))-countkr(myid)+1
      kendr = kstartr+countkr(myid)-1

      djr=jendr-jstartr+1
      dkr=kendr-kstartr+1

#ifdef DEBUG
      if(myid.eq.0)then
        open(610,file='fact/count.out')
        do i=0,numtasks-1
          write(610,'(4i8)')i,countj(i),countk(i),countjp(i)
        enddo
        close(610)
        open(610,file='fact/countr.out')
        do i=0,numtasks-1
          write(610,'(3i8)')i,countjr(i),countkr(i)
        enddo
        close(610)
      endif

      write(*,'(i4,2(a,i5))') myid,' js :',jstart, ' je :', jend
      write(*,'(i4,2(a,i5))') myid,' jsp:',jstartp,' jep:', jendp
      write(*,'(i4,2(a,i5))') myid,' ks :',kstart, ' ke :', kend
      write(*,'(i4,2(a,i5))') myid,' jsr:',jstartr,' jer:', jendr
      write(*,'(i4,2(a,i5))') myid,' ksr:',kstartr,' ker:', kendr
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

      if( dj .lt. 1 ) then            
       write(*,*)'process ',myid,' has work load <1 cell in j direction'
       write(*,*)"Check grid dimensions and number of processes"
       
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif

      if( dk .lt. 1 ) then            
       write(*,*)'process ',myid,' has work load <1 cell in k direction'
       write(*,*)"Check grid dimensions and number of processes"
       
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif
  
      if(.not. allocated(offsetjp)) allocate(offsetjp(0:numtasks-1))
      if(.not. allocated(offsetj)) allocate(offsetj(0:numtasks-1))
      if(.not. allocated(offsetk)) allocate(offsetk(0:numtasks-1))
      
      offsetjp(:)=0
      offsetj(:)=0
      offsetk(:)=0
      do i=1,numtasks-1
        offsetjp(i)= offsetjp(i-1) + countjp(i-1)
        offsetj(i)= offsetj(i-1) + countj(i-1)
        offsetk(i)= offsetk(i-1) + countk(i-1)
      end do
      
      if(.not. allocated(offsetjr)) allocate(offsetjr(0:numtasks-1))
      if(.not. allocated(offsetkr)) allocate(offsetkr(0:numtasks-1))
      
      offsetjr(:)=0
      offsetkr(:)=0
      do i=1,numtasks-1
        offsetjr(i)= offsetjr(i-1) + countjr(i-1)
        offsetkr(i)= offsetkr(i-1) + countkr(i-1)
      end do

#ifdef DEBUG
      if(myid.eq.0)then
        open(610,file='fact/offset.out')
        do i=0,numtasks-1
          write(610,'(4i8)')i,offsetj(i),offsetk(i),offsetjp(i)
        enddo
        close(610)
        open(610,file='fact/offsetr.out')
        do i=0,numtasks-1
          write(610,'(3i8)')i,offsetjr(i),offsetkr(i)
        enddo
        close(610)
      endif
#endif

!-------For MPI-IO--------------------------------
      mydata= n2*dk*n1
      mydatam = n2m*dk*n1m

      if(myid .eq. numtasks-1) mydata = n2*(dk+1)*n1
      
      if(.not. allocated(countf)) allocate(countf(0:numtasks-1))
      if(.not. allocated(offsetf)) allocate(offsetf(0:numtasks-1))
       
      call MPI_ALLGATHER(mydata, 1, MPI_INTEGER, countf, 1, MPI_INTEGER,MPI_COMM_WORLD,ierr)
    
      offsetf(:)=0
      do i=1,numtasks-1
        offsetf(i)= offsetf(i-1) + countf(i-1)
      end do
      
      end subroutine mpi_workdistribution 
!===============================================
     subroutine update_both_ghosts(n1,n2,q1,ks,ke)
      use param, only: sendbuff,recvbuff
      use constants
      use mpih
#ifdef NCCLAVAIL
      use nccl
#endif
!@cuf use cudafor
      implicit none
      integer, intent(in) :: ks,ke
      real(DP),intent(inout) :: q1(n1,n2,ks-lvlhalo:ke+lvlhalo)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
      integer :: i,j,k
!@cuf integer :: istat
#ifdef USE_CUDA
      attributes(managed) :: q1
#endif

      if (numtasks > 1) then
!@cuf   istat = cudaDeviceSynchronize

        my_down=myid-1

        my_up=myid+1

        if(myid .eq. 0) my_down=MPI_PROC_NULL
        if(myid .eq. numtasks-1) my_up=MPI_PROC_NULL

        mydata= n1*n2*lvlhalo

#ifdef USE_CUDA
        if (allocated(sendbuff) .and. size(sendbuff) < 2*mydata) deallocate(sendbuff)
        if (.not. allocated(sendbuff)) allocate(sendbuff(2*mydata))
        if (allocated(recvbuff) .and. size(recvbuff) < 2*mydata) deallocate(recvbuff)
        if (.not. allocated(recvbuff)) allocate(recvbuff(2*mydata))

        ! Pack buffers
        if (my_up .ne. MPI_PROC_NULL) then
          !$cuf kernel do (3)
          do k = 1, lvlhalo
            do j = 1, n2
              do i = 1, n1
                sendbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 0 * n1*n2*lvlhalo) = q1(i,j,ke-lvlhalo+k)
              enddo
            enddo
          enddo
        endif

        if (my_down .ne. MPI_PROC_NULL) then
          !$cuf kernel do (3)
          do k = 1, lvlhalo
            do j = 1, n2
              do i = 1, n1
                sendbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 1 * n1*n2*lvlhalo) = q1(i,j,ks+k-1)
              enddo
            enddo
          enddo
        endif

!@cuf   istat = cudaDeviceSynchronize()

#ifndef NCCLAVAIL
        tag=1

        call MPI_ISEND(sendbuff(1), mydata, MDP, &
        my_up,tag,MPI_COMM_WORLD,req(1),ierr)

        call MPI_ISEND(sendbuff(n1*n2*lvlhalo+1), mydata,  MDP, &
        my_down,tag,MPI_COMM_WORLD,req(2), ierr)

        call MPI_IRECV(recvbuff(1), mydata,  MDP, &
        my_down,tag,MPI_COMM_WORLD,req(3),ierr)

        call MPI_IRECV(recvbuff(n1*n2*lvlhalo+1), mydata,  MDP, &
        my_up, tag,MPI_COMM_WORLD,req(4),ierr)

        call MPI_Waitall(4,req,status,ierr)
#else
        nccl_result = ncclGroupStart()
        if (my_up .ne. MPI_PROC_NULL) then
          nccl_result = ncclSend(sendbuff(1), mydata, ncclDouble, my_up, nccl_comm, 0)
          nccl_result = ncclRecv(recvbuff(n1*n2*lvlhalo+1), mydata, ncclDouble, my_up, nccl_comm, 0)
        endif
        if (my_down .ne. MPI_PROC_NULL) then
          nccl_result = ncclSend(sendbuff(n1*n2*lvlhalo+1), mydata, ncclDouble, my_down, nccl_comm, 0)
          nccl_result = ncclRecv(recvbuff(1), mydata, ncclDouble, my_down, nccl_comm, 0)
        endif
        nccl_result = ncclGroupEnd()
#endif

        ! Unpack buffers
        if (my_down .ne. MPI_PROC_NULL) then
          !$cuf kernel do (3)
          do k = 1, lvlhalo
            do j = 1, n2
              do i = 1, n1
                q1(i,j,ks-lvlhalo+k-1) = recvbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 0 * n1*n2*lvlhalo)
              enddo
            enddo
          enddo
        endif

        if (my_up .ne. MPI_PROC_NULL) then
          !$cuf kernel do (3)
          do k = 1, lvlhalo
            do j = 1, n2
              do i = 1, n1
                q1(i,j,ke+k) = recvbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 1 * n1*n2*lvlhalo)
              enddo
            enddo
          enddo
        endif
#else
        tag=1
        call MPI_ISEND(q1(1,1,ke-lvlhalo+1), mydata, MDP, &
        my_up,tag,MPI_COMM_WORLD,req(1),ierr)

        call MPI_ISEND(q1(1,1,ks), mydata,  MDP, &
        my_down,tag,MPI_COMM_WORLD,req(2), ierr)

        call MPI_IRECV(q1(1,1,ks-lvlhalo), mydata,  MDP, &
        my_down,tag,MPI_COMM_WORLD,req(3),ierr)

        call MPI_IRECV(q1(1,1,ke+1), mydata,  MDP, &
        my_up, tag,MPI_COMM_WORLD,req(4),ierr)

        call MPI_Waitall(4,req,status,ierr)
#endif

!@cuf   istat = cudaDeviceSynchronize
      endif

      end subroutine update_both_ghosts

      subroutine update_both_ghosts_3(n1,n2,q1,q2,q3,ks,ke)
      use param, only: sendbuff,recvbuff
      use constants
      use mpih
#ifdef NCCLAVAIL
      use nccl
#endif
!@cuf use cudafor
      implicit none
      integer, intent(in) :: ks,ke
      real(DP),intent(inout) :: q1(n1,n2,ks-lvlhalo:ke+lvlhalo)
      real(DP),intent(inout) :: q2(n1,n2,ks-lvlhalo:ke+lvlhalo)
      real(DP),intent(inout) :: q3(n1,n2,ks-lvlhalo:ke+lvlhalo)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
      integer :: i,j,k
!@cuf integer :: istat
#ifdef USE_CUDA
      attributes(managed) :: q1, q2, q3
#endif

      if (numtasks > 1) then
!@cuf   istat = cudaDeviceSynchronize


        my_down=myid-1

        my_up=myid+1

        if(myid .eq. 0) my_down=MPI_PROC_NULL
        if(myid .eq. numtasks-1) my_up=MPI_PROC_NULL

        mydata= n1*n2*lvlhalo
        if (allocated(sendbuff) .and. size(sendbuff) < 3*2*mydata) deallocate(sendbuff)
        if (.not. allocated(sendbuff)) allocate(sendbuff(3*2*mydata))
        if (allocated(recvbuff) .and. size(recvbuff) < 3*2*mydata) deallocate(recvbuff)
        if (.not. allocated(recvbuff)) allocate(recvbuff(3*2*mydata))

        ! Pack buffers
        if (my_up .ne. MPI_PROC_NULL) then
#ifdef USE_CUDA
          !$cuf kernel do (3)
#endif
          do k = 1, lvlhalo
            do j = 1, n2
              do i = 1, n1
                sendbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 0 * n1*n2*lvlhalo) = q1(i,j,ke-lvlhalo+k)
                sendbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 1 * n1*n2*lvlhalo) = q2(i,j,ke-lvlhalo+k)
                sendbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 2 * n1*n2*lvlhalo) = q3(i,j,ke-lvlhalo+k)
              enddo
            enddo
          enddo
        endif

        if (my_down .ne. MPI_PROC_NULL) then
#ifdef USE_CUDA
          !$cuf kernel do (3)
#endif
          do k = 1, lvlhalo
            do j = 1, n2
              do i = 1, n1
                sendbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 3 * n1*n2*lvlhalo) = q1(i,j,ks+k-1)
                sendbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 4 * n1*n2*lvlhalo) = q2(i,j,ks+k-1)
                sendbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 5 * n1*n2*lvlhalo) = q3(i,j,ks+k-1)
              enddo
            enddo
          enddo
        endif

!@cuf   istat = cudaDeviceSynchronize

#ifndef NCCLAVAIL
        tag=1

        call MPI_ISEND(sendbuff(1), 3*mydata, MDP, &
        my_up,tag,MPI_COMM_WORLD,req(1),ierr)

        call MPI_ISEND(sendbuff(3*n1*n2*lvlhalo+1), 3*mydata,  MDP, &
        my_down,tag,MPI_COMM_WORLD,req(2), ierr)

        call MPI_IRECV(recvbuff(1), 3*mydata,  MDP, &
        my_down,tag,MPI_COMM_WORLD,req(3),ierr)

        call MPI_IRECV(recvbuff(3*n1*n2*lvlhalo+1), 3*mydata,  MDP, &
        my_up, tag,MPI_COMM_WORLD,req(4),ierr)

        call MPI_Waitall(4,req,status,ierr)
#else
        nccl_result = ncclGroupStart()
        if (my_up .ne. MPI_PROC_NULL) then
          nccl_result = ncclSend(sendbuff(1), 3*mydata, ncclDouble, my_up, nccl_comm, 0)
          nccl_result = ncclRecv(recvbuff(3*n1*n2*lvlhalo+1), 3*mydata, ncclDouble, my_up, nccl_comm, 0)
        endif
        if (my_down .ne. MPI_PROC_NULL) then
          nccl_result = ncclSend(sendbuff(3*n1*n2*lvlhalo+1), 3*mydata, ncclDouble, my_down, nccl_comm, 0)
          nccl_result = ncclRecv(recvbuff(1), 3*mydata, ncclDouble, my_down, nccl_comm, 0)
        endif
        nccl_result = ncclGroupEnd()
#endif

        ! Unpack buffers
        if (my_down .ne. MPI_PROC_NULL) then
#ifdef USE_CUDA
          !$cuf kernel do (3)
#endif
          do k = 1, lvlhalo
            do j = 1, n2
              do i = 1, n1
                q1(i,j,ks-lvlhalo+k-1) = recvbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 0 * n1*n2*lvlhalo)
                q2(i,j,ks-lvlhalo+k-1) = recvbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 1 * n1*n2*lvlhalo)
                q3(i,j,ks-lvlhalo+k-1) = recvbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 2 * n1*n2*lvlhalo)
              enddo
            enddo
          enddo
        endif

        if (my_up .ne. MPI_PROC_NULL) then
#ifdef USE_CUDA
          !$cuf kernel do (3)
#endif
          do k = 1, lvlhalo
            do j = 1, n2
              do i = 1, n1
                q1(i,j,ke+k) = recvbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 3 * n1*n2*lvlhalo)
                q2(i,j,ke+k) = recvbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 4 * n1*n2*lvlhalo)
                q3(i,j,ke+k) = recvbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 5 * n1*n2*lvlhalo)
              enddo
            enddo
          enddo
        endif

!@cuf   istat = cudaDeviceSynchronize
      endif

      end subroutine update_both_ghosts_3
!=========================================
      subroutine update_upper_ghost(n1,n2,q1)
      use param, only: sendbuff,recvbuff
      use constants
      use mpih
      use mpi_param, only: kstart,kend,dk
#ifdef NCCLAVAIL
      use nccl
#endif
!@cuf use cudafor
      implicit none
      real(DP),intent(inout) :: q1(n1,n2,kstart-lvlhalo:kend+lvlhalo)
      integer,intent(in) :: n1,n2
      integer :: mydata, i, j, k
      integer :: my_down, my_up,tag
!@cuf integer :: istat
#ifdef USE_CUDA
      attributes(managed) :: q1
#endif

      if (numtasks > 1) then
!@cuf   istat = cudaDeviceSynchronize

        mydata= n1*n2*lvlhalo

        my_down= myid-1

        my_up= myid+1

        if(myid .eq. 0) my_down=MPI_PROC_NULL
        if(myid .eq. numtasks-1) my_up=MPI_PROC_NULL

#ifdef USE_CUDA
        if (allocated(sendbuff) .and. size(sendbuff) < mydata) deallocate(sendbuff)
        if (.not. allocated(sendbuff)) allocate(sendbuff(mydata))
        if (allocated(recvbuff) .and. size(recvbuff) < mydata) deallocate(recvbuff)
        if (.not. allocated(recvbuff)) allocate(recvbuff(mydata))

        ! Pack buffers
        if (my_down .ne. MPI_PROC_NULL) then
          !$cuf kernel do (3)
          do k = 1, lvlhalo
            do j = 1, n2
              do i = 1, n1
                sendbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2) = q1(i,j,kstart+k-1)
              enddo
            enddo
          enddo
        endif

!@cuf   istat = cudaDeviceSynchronize

#ifndef NCCLAVAIL
        tag=1

        call MPI_ISEND(sendbuff, mydata, MDP, &
        my_down, tag, MPI_COMM_WORLD, req(1), ierr)
        call MPI_IRECV(recvbuff, mydata, MDP, &
        my_up,tag, MPI_COMM_WORLD, req(2), ierr)

        call MPI_Waitall(2,req,status,ierr)
#else
        nccl_result = ncclGroupStart()
        if (my_down .ne. MPI_PROC_NULL) then
          nccl_result = ncclSend(sendbuff, mydata, ncclDouble, my_down, nccl_comm, 0)
        endif
        if (my_up .ne. MPI_PROC_NULL) then
          nccl_result = ncclRecv(recvbuff, mydata, ncclDouble, my_up, nccl_comm, 0)
        endif
        nccl_result = ncclGroupEnd()
#endif

        ! Unpack buffers
        if (my_up .ne. MPI_PROC_NULL) then
          !$cuf kernel do (3)
          do k = 1, lvlhalo
            do j = 1, n2
              do i = 1, n1
                q1(i,j,kend+k) = recvbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2)
              enddo
            enddo
          enddo
        endif
#else
        tag=1
        call MPI_ISEND(q1(1,1,kstart), mydata, MDP, &
        my_down, tag, MPI_COMM_WORLD, req(1), ierr)
        call MPI_IRECV(q1(1,1,kend+1), mydata, MDP, &
        my_up,tag, MPI_COMM_WORLD, req(2), ierr)

        call MPI_Waitall(2,req,status,ierr)
#endif

!@cuf   istat = cudaDeviceSynchronize
      endif

      end subroutine update_upper_ghost
!=========================================
      subroutine update_lower_ghost(n1,n2,q1)
      use constants
      use mpih
      use mpi_param, only: kstart,kend,dk
      implicit none
      real(DP),intent(inout) :: q1(n1,n2,kstart-lvlhalo:kend+lvlhalo)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
       
      mydata= n1*n2*lvlhalo
      
      my_down= myid-1
      
      my_up= myid+1

      if(myid .eq. 0) my_down=MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up=MPI_PROC_NULL
      
      
      tag=1
      
      call MPI_ISEND(q1(1,1,kend-lvlhalo+1), mydata,  MDP, &
      my_up, tag, MPI_COMM_WORLD, req(1), ierr)
      
      call MPI_IRECV(q1(1,1,kstart-lvlhalo), mydata,  MDP, &
      my_down,tag, MPI_COMM_WORLD, req(2), ierr)
       
      call MPI_Waitall(2,req,status,ierr)
    
      end subroutine update_lower_ghost
!=========================================
      subroutine sync_collision(n1,n2,q1)
      use param, only: sendbuff,recvbuff
      use mpih
      use mpi_param, only: kstart,kend,dk
#ifdef NCCLAVAIL
      use nccl
#endif
!@cuf use cudafor
      implicit none
      integer,intent(inout) :: q1(n1,n2,kstart-1:kend)
      !integer :: buf1(n1,n2), buf2(n1,n2)
      !integer, allocatable, dimension(:,:) :: buf1, buf2
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
      integer :: ic,jc
      integer :: i,j
!@cuf integer :: istat

#ifdef USE_CUDA
      attributes(managed) :: q1
#endif
       
      if (numtasks > 1) then
!@cuf   istat = cudaDeviceSynchronize
        mydata= n1*n2

        my_down= myid-1

        my_up= myid+1

        if(myid .eq. 0) my_down=MPI_PROC_NULL
        if(myid .eq. numtasks-1) my_up=MPI_PROC_NULL

        if (allocated(sendbuff) .and. size(sendbuff) < mydata) deallocate(sendbuff)
        if (.not. allocated(sendbuff)) allocate(sendbuff(mydata))
        if (allocated(recvbuff) .and. size(recvbuff) < mydata) deallocate(recvbuff)
        if (.not. allocated(recvbuff)) allocate(recvbuff(mydata))
     
        ! Pack buffers
        if (my_down .ne. MPI_PROC_NULL) then
#ifdef USE_CUDA
          !$cuf kernel do(2)
#endif
          do j = 1, n2
            do i = 1, n1
              sendbuff(i + (j-1) * n1) = q1(i,j,kstart-1)
            end do
          end do
        endif

!@cuf   istat = cudaDeviceSynchronize

#ifndef NCCLAVAIL
        tag=1

        call MPI_ISEND(sendbuff, mydata, MPI_INTEGER, &
        my_down, tag, MPI_COMM_WORLD, req(1), ierr)

        call MPI_IRECV(recvbuff, mydata, MPI_INTEGER, &
        my_up,tag, MPI_COMM_WORLD, req(2), ierr)

        call MPI_Waitall(2,req,status,ierr)
#else
        nccl_result = ncclGroupStart()
        if (my_down .ne. MPI_PROC_NULL) then
          nccl_result = ncclSend(sendbuff, mydata, ncclDouble, my_down, nccl_comm, 0)
        endif
        if (my_up .ne. MPI_PROC_NULL) then
          nccl_result = ncclRecv(recvbuff, mydata, ncclDouble, my_up, nccl_comm, 0)
        endif
        nccl_result = ncclGroupEnd()
#endif
    
        if (myid < numtasks-1) then
#ifdef USE_CUDA
          !$cuf kernel do (2)
#endif
          do jc=1,n2
           do ic=1,n1
             if(q1(ic,jc,kend).ne.recvbuff(ic + (jc-1) * n1).and. &
               q1(ic,jc,kend).ne.0.and.            &
               recvbuff(ic + (jc-1) * n1).ne.0) then
                q1(ic,jc,kend)=-1
             end if
           end do
          end do
        endif

        ! Pack buffers
        if (my_up .ne. MPI_PROC_NULL) then
#ifdef USE_CUDA
          !$cuf kernel do(2)
#endif
          do j = 1, n2
            do i = 1, n1
              sendbuff(i + (j-1) * n1) = q1(i,j,kend)
            end do
          end do
        endif


!@cuf   istat = cudaDeviceSynchronize

#ifndef NCCLAVAIL
        call MPI_ISEND(sendbuff, mydata, MPI_INTEGER, &
        my_up, tag, MPI_COMM_WORLD, req(1), ierr)

        call MPI_IRECV(recvbuff, mydata, MPI_INTEGER, &
        my_down,tag, MPI_COMM_WORLD, req(2), ierr)

        call MPI_Waitall(2,req,status,ierr)
#else
        nccl_result = ncclGroupStart()
        if (my_down .ne. MPI_PROC_NULL) then
          nccl_result = ncclRecv(recvbuff, mydata, ncclDouble, my_down, nccl_comm, 0)
        endif
        if (my_up .ne. MPI_PROC_NULL) then
          nccl_result = ncclSend(sendbuff, mydata, ncclDouble, my_up, nccl_comm, 0)
        endif
        nccl_result = ncclGroupEnd()
#endif

        if (myid > 0) then
#ifdef USE_CUDA
          !$cuf kernel do (2)
#endif
          do jc = 1,n2
            do ic = 1,n1
              q1(ic,jc,kstart-1) = recvbuff(ic + (jc-1) * n1)
            end do
          end do
        endif
!@cuf   istat = cudaDeviceSynchronize


      endif
 
      end subroutine sync_collision
!=========================================


!=========================================
      subroutine update_add_upper_ghost(n1,n2,q1)
      use param, only: sendbuff,recvbuff
      use constants
      use mpih
      use mpi_param, only: kstart,kend,dk
      use local_arrays, only: buf
!@cuf use cudafor
      implicit none
      real(DP),intent(inout) :: q1(n1,n2,kstart-lvlhalo:kend+lvlhalo-1)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag

      integer :: ic,jc,i,j,k
      real(DP) :: cksum
!@cuf integer :: istat

#ifdef USE_CUDA
      attributes(managed) :: q1
#endif

      if (numtasks > 1) then
!@cuf   istat = cudaDeviceSynchronize

        mydata= n1*n2*lvlhalo

        my_down= myid-1

        my_up= myid+1

        if(myid .eq. 0) my_down= MPI_PROC_NULL
        if(myid .eq. numtasks-1) my_up= MPI_PROC_NULL

        if (allocated(sendbuff) .and. size(sendbuff) < mydata) deallocate(sendbuff)
        if (.not. allocated(sendbuff)) allocate(sendbuff(mydata))
        if (allocated(recvbuff) .and. size(recvbuff) < mydata) deallocate(recvbuff)
        if (.not. allocated(recvbuff)) allocate(recvbuff(mydata))

        recvbuff = 0.d0

        ! Pack buffers
        if (my_down .ne. MPI_PROC_NULL) then
#ifdef USE_CUDA
          !$cuf kernel do (3)
#endif
          do k = 1, lvlhalo
            do j = 1, n2
              do i = 1, n1
                sendbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2) = q1(i,j,kstart-lvlhalo+k-1)
              enddo
            enddo
          enddo
        endif

!@cuf   istat = cudaDeviceSynchronize

        tag=1

        call MPI_ISEND(sendbuff,mydata,MDP, &
        my_down, tag, MPI_COMM_WORLD, req(1), ierr)

        call MPI_IRECV(recvbuff, mydata, MDP, &
        my_up,tag, MPI_COMM_WORLD, req(2), ierr)

        call MPI_Waitall(2,req,status,ierr)

        do ic=1,lvlhalo
         jc=kend-lvlhalo+ic
         if(jc.eq.kend)then
#ifdef USE_CUDA
           !$cuf kernel do (2)
#endif
           do j = 1, n2
             do i = 1, n1
               !q1(i,j,jc) = buf(i,j,ic)
               q1(i,j,jc) = recvbuff(i + (j-1) * n1 + (ic-1) * n1*n2)
             enddo
           enddo
         else
#ifdef USE_CUDA
           !$cuf kernel do (2)
#endif
           do j = 1, n2
             do i = 1, n1
               !q1(i,j,jc) = q1(i,j,jc) + buf(i,j,ic)
               q1(i,j,jc) = q1(i,j,jc) + recvbuff(i + (j-1) * n1 + (ic-1) * n1*n2)
             enddo
           enddo
         end if
        end do

!@cuf   istat = cudaDeviceSynchronize
      endif

      end subroutine update_add_upper_ghost

      subroutine update_add_upper_ghost_3(n1,n2,q1,q2,q3)
      use param, only: sendbuff,recvbuff
      use constants
      use mpih
      use mpi_param, only: kstart,kend,dk
      use local_arrays, only: buf
#ifdef NCCLAVAIL
      use nccl
#endif
!@cuf use cudafor
      implicit none
      real(DP),intent(inout) :: q1(n1,n2,kstart-lvlhalo:kend+lvlhalo-1)
      real(DP),intent(inout) :: q2(n1,n2,kstart-lvlhalo:kend+lvlhalo-1)
      real(DP),intent(inout) :: q3(n1,n2,kstart-lvlhalo:kend+lvlhalo-1)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag

      integer :: ic,jc,i,j,k
      real(DP) :: cksum
!@cuf integer :: istat

#ifdef USE_CUDA
      attributes(managed) :: q1,q2,q3
#endif

      if (numtasks > 1) then
!@cuf   istat = cudaDeviceSynchronize

        mydata= n1*n2*lvlhalo

        my_down= myid-1

        my_up= myid+1

        if(myid .eq. 0) my_down= MPI_PROC_NULL
        if(myid .eq. numtasks-1) my_up= MPI_PROC_NULL

        if (allocated(sendbuff) .and. size(sendbuff) < 3*mydata) deallocate(sendbuff)
        if (.not. allocated(sendbuff)) allocate(sendbuff(3*mydata))
        if (allocated(recvbuff) .and. size(recvbuff) < 3*mydata) deallocate(recvbuff)
        if (.not. allocated(recvbuff)) allocate(recvbuff(3*mydata))

        recvbuff = 0.d0

        ! Pack buffers
        if (my_down .ne. MPI_PROC_NULL) then
#ifdef USE_CUDA
          !$cuf kernel do (3)
#endif
          do k = 1, lvlhalo
            do j = 1, n2
              do i = 1, n1
                sendbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2                  ) = q1(i,j,kstart-lvlhalo+k-1)
                sendbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 1*n1*n2*lvlhalo) = q2(i,j,kstart-lvlhalo+k-1)
                sendbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 2*n1*n2*lvlhalo) = q3(i,j,kstart-lvlhalo+k-1)
              enddo
            enddo
          enddo
        endif

!@cuf   istat = cudaDeviceSynchronize

        tag=1

#ifndef NCCLAVAIL
        call MPI_ISEND(sendbuff,3*mydata,MDP, &
        my_down, tag, MPI_COMM_WORLD, req(1), ierr)

        call MPI_IRECV(recvbuff, 3*mydata, MDP, &
        my_up,tag, MPI_COMM_WORLD, req(2), ierr)

        call MPI_Waitall(2,req,status,ierr)
#else
        nccl_result = ncclGroupStart()
        if (my_down .ne. MPI_PROC_NULL) then
          nccl_result = ncclSend(sendbuff, 3*mydata, ncclDouble, my_down, nccl_comm, 0)
        endif
        if (my_up .ne. MPI_PROC_NULL) then
          nccl_result = ncclRecv(recvbuff, 3*mydata, ncclDouble, my_up, nccl_comm, 0)
        endif
        nccl_result = ncclGroupEnd()
#endif

        do ic=1,lvlhalo
         jc=kend-lvlhalo+ic
         if(jc.eq.kend)then
#ifdef USE_CUDA
           !$cuf kernel do (2)
#endif
           do j = 1, n2
             do i = 1, n1
               !q1(i,j,jc) = buf(i,j,ic)
               q1(i,j,jc) = recvbuff(i + (j-1) * n1 + (ic-1) * n1*n2                  )
               q2(i,j,jc) = recvbuff(i + (j-1) * n1 + (ic-1) * n1*n2 + 1*n1*n2*lvlhalo)
               q3(i,j,jc) = recvbuff(i + (j-1) * n1 + (ic-1) * n1*n2 + 2*n1*n2*lvlhalo)
             enddo
           enddo
         else
#ifdef USE_CUDA
           !$cuf kernel do (2)
#endif
           do j = 1, n2
             do i = 1, n1
               !q1(i,j,jc) = q1(i,j,jc) + buf(i,j,ic)
               q1(i,j,jc) = q1(i,j,jc) + recvbuff(i + (j-1) * n1 + (ic-1) * n1*n2                  )
               q2(i,j,jc) = q2(i,j,jc) + recvbuff(i + (j-1) * n1 + (ic-1) * n1*n2 + 1*n1*n2*lvlhalo)
               q3(i,j,jc) = q3(i,j,jc) + recvbuff(i + (j-1) * n1 + (ic-1) * n1*n2 + 2*n1*n2*lvlhalo)
             enddo
           enddo
         end if
        end do

!@cuf   istat = cudaDeviceSynchronize
      endif

      end subroutine update_add_upper_ghost_3
!=========================================
      subroutine update_add_lower_ghost(n1,n2,q1)
      use param, only: sendbuff,recvbuff
      use constants
      use mpih
      use mpi_param, only: kstart,kend,dk
      use local_arrays, only: buf
#ifdef NCCLAVAIL
      use nccl
#endif
!@cuf use cudafor
      implicit none
      real(DP),intent(inout) :: q1(n1,n2,kstart-lvlhalo:kend+lvlhalo-1)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
      integer :: ic,jc,i,j,k
!@cuf integer :: istat
#ifdef USE_CUDA
      attributes(managed) :: q1
#endif

      if (numtasks > 1) then
!@cuf   istat = cudaDeviceSynchronize

        mydata= n1*n2*(lvlhalo)

        my_down= myid-1

        my_up= myid+1

        if(myid .eq. 0) my_down= MPI_PROC_NULL
        if(myid .eq. numtasks-1) my_up= MPI_PROC_NULL

        if (allocated(sendbuff) .and. size(sendbuff) < mydata) deallocate(sendbuff)
        if (.not. allocated(sendbuff)) allocate(sendbuff(mydata))
        if (allocated(recvbuff) .and. size(recvbuff) < mydata) deallocate(recvbuff)
        if (.not. allocated(recvbuff)) allocate(recvbuff(mydata))

        recvbuff = 0.d0

        ! Pack buffers
        if (my_up .ne. MPI_PROC_NULL) then
#ifdef USE_CUDA
          !$cuf kernel do (3)
#endif
          do k = 1, lvlhalo
            do j = 1, n2
              do i = 1, n1
                sendbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2) = q1(i,j,kend+k-1)
              enddo
            enddo
          enddo
        endif

!@cuf   istat = cudaDeviceSynchronize

#ifndef NCCLAVAIL
        tag=1

        call MPI_ISEND(sendbuff,mydata,MDP, &
        my_up, tag, MPI_COMM_WORLD, req(1), ierr)

        call MPI_IRECV(recvbuff, mydata, MDP, &
        my_down,tag, MPI_COMM_WORLD, req(2), ierr)

        call MPI_Waitall(2,req,status,ierr)
#else
        nccl_result = ncclGroupStart()
        if (my_down .ne. MPI_PROC_NULL) then
          nccl_result = ncclRecv(recvbuff, 3*mydata, ncclDouble, my_down, nccl_comm, 0)
        endif
        if (my_up .ne. MPI_PROC_NULL) then
          nccl_result = ncclSend(sendbuff, 3*mydata, ncclDouble, my_up, nccl_comm, 0)
        endif
        nccl_result = ncclGroupEnd()
#endif


        do ic=1,lvlhalo
         jc=kstart+ic-2
#ifdef USE_CUDA
         !$cuf kernel do (2)
#endif
         do j = 1, n2
           do i = 1, n1
             !q1(i,j,jc) = q1(i,j,jc) + buf(i,j,ic)
             q1(i,j,jc) = q1(i,j,jc) + recvbuff(i + (j-1) * n1 +(ic-1) * n1*n2)
           enddo
         enddo
        end do

!@cuf   istat = cudaDeviceSynchronize
      endif

      end subroutine update_add_lower_ghost

      subroutine update_add_lower_ghost_3(n1,n2,q1,q2,q3)
      use param, only: sendbuff,recvbuff
      use constants
      use mpih
      use mpi_param, only: kstart,kend,dk
      use local_arrays, only: buf
!@cuf use cudafor
      implicit none
      real(DP),intent(inout) :: q1(n1,n2,kstart-lvlhalo:kend+lvlhalo-1)
      real(DP),intent(inout) :: q2(n1,n2,kstart-lvlhalo:kend+lvlhalo-1)
      real(DP),intent(inout) :: q3(n1,n2,kstart-lvlhalo:kend+lvlhalo-1)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
      integer :: ic,jc,i,j,k
!@cuf integer :: istat
#ifdef USE_CUDA
      attributes(managed) :: q1, q2, q3
#endif

      if (numtasks > 1) then
!@cuf   istat = cudaDeviceSynchronize

        mydata= n1*n2*(lvlhalo)

        my_down= myid-1

        my_up= myid+1

        if(myid .eq. 0) my_down= MPI_PROC_NULL
        if(myid .eq. numtasks-1) my_up= MPI_PROC_NULL

        if (allocated(sendbuff) .and. size(sendbuff) < 3*mydata) deallocate(sendbuff)
        if (.not. allocated(sendbuff)) allocate(sendbuff(3*mydata))
        if (allocated(recvbuff) .and. size(recvbuff) < 3*mydata) deallocate(recvbuff)
        if (.not. allocated(recvbuff)) allocate(recvbuff(3*mydata))

        recvbuff = 0.d0

        ! Pack buffers
        if (my_up .ne. MPI_PROC_NULL) then
#ifdef USE_CUDA
          !$cuf kernel do (3)
#endif
          do k = 1, lvlhalo
            do j = 1, n2
              do i = 1, n1
                sendbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2) = q1(i,j,kend+k-1)
                sendbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 1*n1*n2*lvlhalo) = q2(i,j,kend+k-1)
                sendbuff(i  + (j - 1) * n1 + (k - 1) * n1*n2 + 2*n1*n2*lvlhalo) = q3(i,j,kend+k-1)
              enddo
            enddo
          enddo
        endif

!@cuf   istat = cudaDeviceSynchronize

        tag=1

        call MPI_ISEND(sendbuff,3*mydata,MDP, &
        my_up, tag, MPI_COMM_WORLD, req(1), ierr)

        call MPI_IRECV(recvbuff, 3*mydata, MDP, &
        my_down,tag, MPI_COMM_WORLD, req(2), ierr)

        call MPI_Waitall(2,req,status,ierr)

        do ic=1,lvlhalo
         jc=kstart+ic-2
#ifdef USE_CUDA
         !$cuf kernel do (2)
#endif
         do j = 1, n2
           do i = 1, n1
             !q1(i,j,jc) = q1(i,j,jc) + buf(i,j,ic)
             q1(i,j,jc) = q1(i,j,jc) + recvbuff(i + (j-1) * n1 +(ic-1) * n1*n2)
             q2(i,j,jc) = q2(i,j,jc) + recvbuff(i + (j-1) * n1 +(ic-1) * n1*n2 + 1*n1*n2*lvlhalo)
             q3(i,j,jc) = q3(i,j,jc) + recvbuff(i + (j-1) * n1 +(ic-1) * n1*n2 + 2*n1*n2*lvlhalo)
           enddo
         enddo
        end do

!@cuf   istat = cudaDeviceSynchronize
      endif

      end subroutine update_add_lower_ghost_3

!==============================================
      subroutine mpi_globalsum_double_arr(var,nvar)
        use constants
        use mpih
!@cuf   use cudafor
        implicit none
        real(DP),intent(inout),dimension(nvar) :: var
        integer,intent(in) :: nvar
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: var
#endif

      if (numtasks > 1) then
!@cuf   istat = cudaDeviceSynchronize
        call MPI_ALLREDUCE(MPI_IN_PLACE,var,nvar,MPI_DOUBLE_PRECISION, &
             MPI_SUM,MPI_COMM_WORLD,ierr)
      endif

      end subroutine mpi_globalsum_double_arr

!==============================================
      subroutine mpi_globalsum_integer_arr(var,nvar)
        use constants
        use mpih
#if NCCLAVAIL
        use nccl
#endif
!@cuf   use cudafor
        implicit none
        integer,intent(inout),dimension(nvar) :: var
        integer,intent(in) :: nvar
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: var
#endif

      if (numtasks > 1) then
!@cuf   istat = cudaDeviceSynchronize
#ifndef NCCLAVAIL
        call MPI_ALLREDUCE(MPI_IN_PLACE,var,nvar,MPI_INTEGER, &
           MPI_SUM,MPI_COMM_WORLD,ierr)
#else
        nccl_result = ncclAllReduce(var, var, nvar, ncclInt, ncclSum, nccl_comm, 0)
#endif
      endif

      end subroutine mpi_globalsum_integer_arr

!==============================================

      subroutine mem_alloc
      use mpih
      use param
      use mpi_param
      use local_arrays
      use stat_arrays
      use rt_arrays

      implicit none
      integer :: merr, merr_all
     
      merr_all = 0 
!-------------------------------------------------
! Arrays with ghost cells
!-------------------------------------------------
      allocate(q1(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(q2(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(q3(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(pr(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(dens(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(dph(1:n1m,1:n2m+2,kstart-lvlhalo:kend+lvlhalo),stat=merr)
      merr_all = merr_all + merr
      allocate(dpht(1:n3m,1:n1m,jstartp:jendp),stat=merr)
      merr_all = merr_all + merr

#if defined(STRONG) || defined(STRONGRV)
!STRONG FV
      allocate(q1g(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(q2g(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(q3g(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(prg(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(qcapg(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(dphg(1:n1m,1:n2m+2,kstart-lvlhalo:kend+lvlhalo),stat=merr)
      merr_all = merr_all + merr
      allocate(dqg(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr                                
#endif

!PPRO FV
if (PPro.GE.1) then
      allocate(maskq1(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(maskq2(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(maskq3(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(maskpr(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
endif

!-------------------------------------------------
! Arrays for linear system
!-------------------------------------------------

      allocate(rhs(1:n1m,1:n2m,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(rhss(1:n1m,1:n2m,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(rhsss(1:n1m,1:n2m,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(rhs_t(1:n2m,1:n1m,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(dq(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr                    
      allocate(ru1(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(ru2(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(ru3(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr
#if defined(STRONG) || defined(STRONGRV)      
      allocate(ru1g(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(ru2g(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(ru3g(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr
#endif      
      allocate(qcap(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(hro(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(ruro(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(forclo(1:n1m,1:n2m,kstart:kend,3), stat=merr)
      merr_all = merr_all + merr
      allocate(forclo_t(1:n2m,1:n1m,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(rhst(1:n3m,1:n1m,jstart:jend), stat=merr)
      merr_all = merr_all + merr
      allocate(forclot(1:n3m,1:n1m,jstart:jend), stat=merr)
      merr_all = merr_all + merr

      !Scar arrays
!      allocate(scarpoint(1:nctot_3d,1:nctot_3d,1:nctot_3d), stat=merr)
!      merr_all = merr_all + merr
!      allocate(idscar(1:nctot_3d,1:nctot_3d,1:nctot_3d), stat=merr)
!      merr_all = merr_all + merr

!---------------------------------------------------
! Arrays for salinity
!---------------------------------------------------
      allocate(dsal(1:n1r,1:n2r,kstartr-lvlhalo:kendr+lvlhalo),stat=merr)
      merr_all = merr_all + merr
      allocate(hsal(1:n1r,1:n2r,kstartr:kendr), stat=merr)
      merr_all = merr_all + merr
      allocate(rusal(1:n1r,1:n2r,kstartr:kendr), stat=merr)
      merr_all = merr_all + merr
      allocate(rhsr(1:n1mr,1:n2mr,kstartr:kendr), stat=merr)
      merr_all = merr_all + merr

!========================================
      if(merr_all.ne.0)then
        write(*,*)myid, ' memory alloc error'
        write(*,*)merr_all
      endif
!-----------------------------------------------------
! statistical arrays
!-----------------------------------------------------

      allocate(q1_me(kstart:kend))
      allocate(q2_me(kstart:kend))
      allocate(q3_me(kstart:kend))
      allocate(q1_rms(kstart:kend))
      allocate(q2_rms(kstart:kend))
      allocate(q3_rms(kstart:kend))
      allocate(dens_me(kstart:kend))
      allocate(dens_rms(kstart:kend))
      allocate(dsal_me(kstartr:kendr))
      allocate(dsal_rms(kstartr:kendr))

      allocate(q3dens_me(kstart:kend))
      allocate(q3dsal_me(kstartr:kendr))

      allocate(dissuc(kstart:kend))
      allocate(disste(kstart:kend))
      allocate(disssa(kstartr:kendr))
      allocate(dissur(kstartr:kendr))

      allocate(buf(n1, n2, lvlhalo))

      return
      end subroutine mem_alloc

!==================================================      
      
      subroutine mem_dealloc
      use local_arrays
      use mpi_param
      use mpih
      use stat_arrays

      implicit none
      
      if(allocated(q1)) deallocate(q1)
      if(allocated(q2)) deallocate(q2)
      if(allocated(q3)) deallocate(q3)
      if(allocated(dens)) deallocate(dens)
      if(allocated(dsal)) deallocate(dsal)

      if(allocated(pr)) deallocate(pr)
      if(allocated(dph)) deallocate(dph)
      if(allocated(dq)) deallocate(dq)
      if(allocated(qcap)) deallocate(qcap)

      if(allocated(hro)) deallocate(hro)
      if(allocated(rhs)) deallocate(rhs)
      if(allocated(rhss)) deallocate(rhss)
      if(allocated(rhsss)) deallocate(rhsss)
      if(allocated(ru1)) deallocate(ru1)
      if(allocated(ru2)) deallocate(ru2)
      if(allocated(ru3)) deallocate(ru3)
      if(allocated(ruro)) deallocate(ruro)

      if(allocated(hsal)) deallocate(hsal)
      if(allocated(rusal)) deallocate(rusal)
      if(allocated(rhsr)) deallocate(rhsr)


#if defined(STRONG) || defined(STRONGRV)
!STRONG FV
      if(allocated(q1g)) deallocate(q1g)
      if(allocated(q2g)) deallocate(q2g)
      if(allocated(q3g)) deallocate(q3g)
      if(allocated(prg)) deallocate(prg)
#endif

!PPRO FV
if (PPro.GE.1) then
      if(allocated(maskq1)) deallocate(maskq1)
      if(allocated(maskq2)) deallocate(maskq2)
      if(allocated(maskq3)) deallocate(maskq3)
      if(allocated(maskpr)) deallocate(maskpr)
endif

!---------------------------------------
      if(allocated(countj)) deallocate(countj)
      if(allocated(countk)) deallocate(countk)
      if(allocated(countjp)) deallocate(countjp)

      if(allocated(offsetj)) deallocate(offsetj)
      if(allocated(offsetk)) deallocate(offsetk)
      if(allocated(offsetjp)) deallocate(offsetjp)

      if(allocated(countjr)) deallocate(countjr)
      if(allocated(countkr)) deallocate(countkr)

      if(allocated(offsetjr)) deallocate(offsetjr)
      if(allocated(offsetkr)) deallocate(offsetkr)
      
      if(allocated(countf)) deallocate(countf)
      if(allocated(offsetf)) deallocate(offsetf)

      if(allocated(q1_me)) deallocate(q1_me)
      if(allocated(q2_me)) deallocate(q2_me)
      if(allocated(q3_me)) deallocate(q3_me)
      if(allocated(q1_rms)) deallocate(q1_rms)
      if(allocated(q2_rms)) deallocate(q2_rms)
      if(allocated(q3_rms)) deallocate(q3_rms)
      if(allocated(dens_me)) deallocate(dens_me)
      if(allocated(dens_rms)) deallocate(dens_rms)
      if(allocated(dsal_me)) deallocate(dsal_me)
      if(allocated(dsal_rms)) deallocate(dsal_rms)

      if(allocated(q3dens_me)) deallocate(q3dens_me)
      if(allocated(q3dsal_me)) deallocate(q3dsal_me)

      if(allocated(disste)) deallocate(disste)
      if(allocated(dissuc)) deallocate(dissuc)
      if(allocated(disssa)) deallocate(disssa)
      if(allocated(dissur)) deallocate(dissur)

      return    
      end subroutine mem_dealloc
!================================================
!================================================
      
      subroutine mpi_clean_ibm
      use mpih
      use param
      use hdf5
      use mpi_param, only : kstart, kend
      use ibm_param
      implicit none
      integer i, n, ntemp,m
 
       do i=1,4

! !     Solid body 
!         ntemp=0
!         do n=1,npunfx(i)
!          if(kstart.le.indgeo(i,n,3).and.kend.ge.indgeo(i,n,3)) then
!           ntemp = ntemp + 1
!           do m=1,3
!            indgeo(i,ntemp,m)=indgeo(i,n,m)
!            indgeoe(i,ntemp,m)=indgeoe(i,n,m)
!           end do
!           distb(i,ntemp)=distb(i,n)
!          end if
!         enddo
!         npunfx(i) = ntemp
        
!         ntemp=0 
!         do n=1,npunifx(i)
!          if(kstart.le.indgeoee(i,n,3).and.kend.ge.indgeoee(i,n,3)) then
!          ntemp = ntemp + 1
!            do m=1,3
!             indgeoee(i,ntemp,m)=indgeoee(i,n,m)
!            enddo
!          end if
!         end do
!         npunifx(i) = ntemp

! !     Solid body MI
!         ntemp=0
!         do n=1,npunfxMI(i)
!          if(kstart.le.indgeoMI(i,n,3).and.kend.ge.indgeoMI(i,n,3)) then
!           ntemp = ntemp + 1
!           do m=1,3
!            indgeoMI(i,ntemp,m)=indgeoMI(i,n,m)
!            indgeoeMI(i,ntemp,m)=indgeoeMI(i,n,m)
!           end do
!           distbMI(i,ntemp)=distbMI(i,n)
!          end if
!         enddo
!         npunfxMI(i) = ntemp
        
!         ntemp=0 
!         do n=1,npunifxMI(i)
!          if(kstart.le.indgeoeeMI(i,n,3).and.kend.ge.indgeoeeMI(i,n,3)) then
!          ntemp = ntemp + 1
!            do m=1,3
!             indgeoeeMI(i,ntemp,m)=indgeoeeMI(i,n,m)
!            enddo
!          end if
!         end do
!         npunifxMI(i) = ntemp

! !     Solid body AO
!         ntemp=0
!         do n=1,npunfxAO(i)
!          if(kstart.le.indgeoAO(i,n,3).and.kend.ge.indgeoAO(i,n,3)) then
!           ntemp = ntemp + 1
!           do m=1,3
!            indgeoAO(i,ntemp,m)=indgeoAO(i,n,m)
!            indgeoeAO(i,ntemp,m)=indgeoeAO(i,n,m)
!           end do
!           distbAO(i,ntemp)=distbAO(i,n)
!          end if
!         enddo
!         npunfxAO(i) = ntemp
        
!         ntemp=0 
!         do n=1,npunifxAO(i)
!          if(kstart.le.indgeoeeAO(i,n,3).and.kend.ge.indgeoeeAO(i,n,3)) then
!          ntemp = ntemp + 1
!            do m=1,3
!             indgeoeeAO(i,ntemp,m)=indgeoeeAO(i,n,m)
!            enddo
!          end if
!         end do
!         npunifxAO(i) = ntemp

!     Solid body AS
        ntemp=0
        do n=1,npunfxAS(i)
         if(kstart.le.indgeoAS(i,n,3).and.kend.ge.indgeoAS(i,n,3)) then
          ntemp = ntemp + 1
          do m=1,3
           indgeoAS(i,ntemp,m)=indgeoAS(i,n,m)
           indgeoeAS(i,ntemp,m)=indgeoeAS(i,n,m)
          end do
          distbAS(i,ntemp)=distbAS(i,n)
         end if
        enddo
        npunfxAS(i) = ntemp
        
        ntemp=0 
        do n=1,npunifxAS(i)
         if(kstart.le.indgeoeeAS(i,n,3).and.kend.ge.indgeoeeAS(i,n,3)) then
         ntemp = ntemp + 1
           do m=1,3
            indgeoeeAS(i,ntemp,m)=indgeoeeAS(i,n,m)
           enddo
         end if
        end do
        npunifxAS(i) = ntemp


       end do

      end subroutine mpi_clean_ibm

#ifdef NCCLAVAIL
      subroutine setup_nccl
      use mpih
!@cuf use cudafor
!@cuf use nccl
      implicit none

#ifdef USE_CUDA
      if (myid == 0) then
        ! Rank 0 generates unique id
        nccl_result = ncclGetUniqueId(nccl_id)
      endif

      call MPI_Bcast(nccl_id%internal, sizeof(nccl_id%internal), &
                     MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

      nccl_result = ncclCommInitRank(nccl_comm, numtasks, nccl_id, myid)
#endif
      end subroutine
#endif

#ifdef CUDECOMPAVAIL
      subroutine setup_cudecomp
      use param
      use mpih
      use mpi_param
      use cudecomp_param
      implicit none
      integer :: istat
      integer(8) :: work_num_elements, work_num_elements_ph
      type(cudecompGridDescConfig) :: config
      type(cudecompGridDescAutotuneOptions) :: at_options

      istat = cudecompInit(cudecomp_handle, MPI_COMM_WORLD)

      ! Autotuning options
      istat = cudecompGridDescAutotuneOptionsSetDefaults(at_options)
      at_options%dtype = CUDECOMP_DOUBLE

      ! JR Note: autotuning can consume additional memory due to
      ! leak issues with OpenMPI. If you run out of memory, disable
      ! autotuning and use a hardcoded setting.
      at_options%autotune_transpose_backend = .true.
      !at_options%disable_nccl_backends = .true.
      !at_options%disable_nvshmem_backends = .true.

      ! Setup desciptor for sol routine transposes
      istat = cudecompGridDescConfigSetDefaults(config)
      config%pdims = [1, numtasks] ! slabs
      config%gdims = [n1m, n2m, n3m] ! slabs

      ! JR Note: hardcoding autotuned setting for 1 Marconi node due to
      ! memory leaking from MPI
      if (numtasks <= 4) then
        config%transpose_comm_backend = CUDECOMP_TRANSPOSE_COMM_NVSHMEM_PL
        at_options%autotune_transpose_backend = .false.
      else
        config%transpose_comm_backend = CUDECOMP_TRANSPOSE_COMM_NCCL
      endif
      config%transpose_axis_contiguous = [.false., .false., .true.]
      istat = cudecompGridDescCreate(cudecomp_handle, cudecomp_grid_desc, config, at_options)

      ! Setup desciptor for phcalc routine transposes
      !istat = cudecompGridDescConfigSetDefaults(config)
      config%pdims = [1, numtasks] ! slabs
      config%gdims = [n1m, n2m+2, n3m] ! slabs
      !config%transpose_comm_backend = CUDECOMP_TRANSPOSE_COMM_NCCL
      config%transpose_axis_contiguous = [.false., .false., .true.]
      istat = cudecompGridDescCreate(cudecomp_handle, cudecomp_grid_desc_ph, config)

      istat = cudecompGetTransposeWorkspaceSize(cudecomp_handle, cudecomp_grid_desc, work_num_elements)
      istat = cudecompGetTransposeWorkspaceSize(cudecomp_handle, cudecomp_grid_desc_ph, work_num_elements_ph)
      istat = cudecompMalloc(cudecomp_handle, cudecomp_grid_desc, cudecomp_work_d, max(work_num_elements, work_num_elements_ph))

      end subroutine setup_cudecomp
#endif
