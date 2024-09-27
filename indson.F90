!***********************************************************************
      subroutine indson
!
!     This routnine finds the indices in the computational grid
!     close to the physical location of the probe.
!     This is necessary since the location of the probe is given
!     in terms of physical coordinates (sonda.in)
!
       use param
       use probes
       implicit none

       integer:: ierr, ll, j,jp,k,kp,i,ip
!
      do ll = 1,nson
         ierr=0
         if((coson(ll,1).lt.xc(1)).or.(coson(ll,1).gt.xc(n1))) ierr=1
         if((coson(ll,2).lt.yc(1)).or.(coson(ll,2).gt.yc(n2))) ierr=1
         if((coson(ll,3).lt.zc(1)).or.(coson(ll,3).gt.zc(n3))) ierr=1
      if(ierr.eq.1) then
         write( 6,*) 'WARNING: probe ',ll,' is outside the domain'
         write(32,*) 'WARNING: probe ',ll,' is outside the domain'
         mason(ll,1)=1
         mason(ll,2)=1
         mason(ll,3)=1
      else
!
!     SEARCHING FOR THE INDICES IN THE GRID CLOSE TO THE PROBE
!
         do i=1,n1m
            ip = i+1
            if((xc(i).le.coson(ll,1)).and. (xc(ip).gt.coson(ll,1))) then
               mason(ll,1)=i
            end if
         end do
         do j=1,n2m
            jp = j+1
            if((yc(j).le.coson(ll,2)).and. (yc(jp).gt.coson(ll,2))) then
               mason(ll,2)=j
            end if
         end do
         do k=1,n1m
            kp = k+1
            if((zc(k).le.coson(ll,3)).and. (zc(kp).gt.coson(ll,3))) then
               mason(ll,3)=k
            end if
         end do
      end if
      end do

      return
      end     

!***********************************************************************
      subroutine indsonAorta 

       use param
       use probes
       use mpih
       use mpi_param
       implicit none 

       integer:: j,jp,k,kp,i,ip
       integer :: tmp(numtasks)


!LEFT HEART
!     Searching for the indices in the grid close to the probe
!
!     ----
         do i=1,n1m
            ip = i+1
            if((xc(i).le.cosonAO(1)).and.  (xc(ip).gt.cosonAO(1))) then
               masonAO(1)=i
            end if
         end do

         do j=1,n2m
            jp = j+1
            if((yc(j).le.cosonAO(2)).and. (yc(jp).gt.cosonAO(2))) then
               masonAO(2)=j
            end if
         end do

         rankAO = -1
         do k=kstart,kend
            kp = k+1
            if((zc(k).le.cosonAO(3)).and. (zc(kp).gt.cosonAO(3))) then
               masonAO(3)=k
               rankAO = myid
            end if
         end do

         call MPI_Allgather(rankAO, 1, MPI_INT, tmp, 1, MPI_INT, MPI_COMM_WORLD, ierr)
         do i = 1, numtasks
           if (tmp(i) .eq. i-1) rankAO = i-1
         end do
!     ----
         do i=1,n1m
            ip = i+1
            if((xc(i).le.cosonLV(1)).and.  (xc(ip).gt.cosonLV(1))) then
               masonLV(1)=i
            end if
         end do

         do j=1,n2m
            jp = j+1
            if((yc(j).le.cosonLV(2)).and.  (yc(jp).gt.cosonLV(2))) then
               masonLV(2)=j
            end if
         end do

         rankLV = -1
         do k=kstart,kend
            kp = k+1
            if((zc(k).le.cosonLV(3)).and.  (zc(kp).gt.cosonLV(3))) then
               masonLV(3)=k
               rankLV=myid
            end if
         end do
      
         call MPI_Allgather(rankLV, 1, MPI_INT, tmp, 1, MPI_INT, MPI_COMM_WORLD, ierr)
         do i = 1, numtasks
           if (tmp(i) .eq. i-1) rankLV = i-1
         end do
!     ----
         do i=1,n1m
            ip = i+1
            if((xc(i).le.cosonLA(1)).and.  (xc(ip).gt.cosonLA(1))) then
               masonLA(1)=i
            end if
         end do

         do j=1,n2m
            jp = j+1
            if((yc(j).le.cosonLA(2)).and.  (yc(jp).gt.cosonLA(2))) then
               masonLA(2)=j
            end if
         end do

         rankLA = -1
         do k=kstart,kend
            kp = k+1
            if((zc(k).le.cosonLA(3)).and.  (zc(kp).gt.cosonLA(3))) then
               masonLA(3)=k
               rankLA=myid
            end if
         end do

         call MPI_Allgather(rankLA, 1, MPI_INT, tmp, 1, MPI_INT, MPI_COMM_WORLD, ierr)
         do i = 1, numtasks
           if (tmp(i) .eq. i-1) rankLA = i-1
         end do




!RIGHT HEART
!     Searching for the indices in the grid close to the probe
!
!     ----
         do i=1,n1m
            ip = i+1
            if((xc(i).le.cosonAP(1)).and.  (xc(ip).gt.cosonAP(1))) then
               masonAP(1)=i
            end if
         end do

         do j=1,n2m
            jp = j+1
            if((yc(j).le.cosonAP(2)).and. (yc(jp).gt.cosonAP(2))) then
               masonAP(2)=j
            end if
         end do

         rankAP = -1
         do k=kstart,kend
            kp = k+1
            if((zc(k).le.cosonAP(3)).and. (zc(kp).gt.cosonAP(3))) then
               masonAP(3)=k
               rankAP = myid
            end if
         end do

         call MPI_Allgather(rankAP, 1, MPI_INT, tmp, 1, MPI_INT, MPI_COMM_WORLD, ierr)
         do i = 1, numtasks
           if (tmp(i) .eq. i-1) rankAP = i-1
         end do
!     ----
         do i=1,n1m
            ip = i+1
            if((xc(i).le.cosonRV(1)).and.  (xc(ip).gt.cosonRV(1))) then
               masonRV(1)=i
            end if
         end do

         do j=1,n2m
            jp = j+1
            if((yc(j).le.cosonRV(2)).and.  (yc(jp).gt.cosonRV(2))) then
               masonRV(2)=j
            end if
         end do

         rankRV = -1
         do k=kstart,kend
            kp = k+1
            if((zc(k).le.cosonRV(3)).and.  (zc(kp).gt.cosonRV(3))) then
               masonRV(3)=k
               rankRV=myid
            end if
         end do
      
         call MPI_Allgather(rankRV, 1, MPI_INT, tmp, 1, MPI_INT, MPI_COMM_WORLD, ierr)
         do i = 1, numtasks
           if (tmp(i) .eq. i-1) rankRV = i-1
         end do
!     ----
         do i=1,n1m
            ip = i+1
            if((xc(i).le.cosonRA(1)).and.  (xc(ip).gt.cosonRA(1))) then
               masonRA(1)=i
            end if
         end do

         do j=1,n2m
            jp = j+1
            if((yc(j).le.cosonRA(2)).and.  (yc(jp).gt.cosonRA(2))) then
               masonRA(2)=j
            end if
         end do

         rankRA = -1
         do k=kstart,kend
            kp = k+1
            if((zc(k).le.cosonRA(3)).and.  (zc(kp).gt.cosonRA(3))) then
               masonRA(3)=k
               rankRA=myid
            end if
         end do

         call MPI_Allgather(rankRA, 1, MPI_INT, tmp, 1, MPI_INT, MPI_COMM_WORLD, ierr)
         do i = 1, numtasks
           if (tmp(i) .eq. i-1) rankRA = i-1
         end do




         
      return
      end     


 !***********************************************************************
      subroutine search_nodicheckEF()
      use param
      use mls_param
      implicit none 

      real(DP), dimension (3,n_EFcheck) :: coords
      real(DP) :: D, Dmin
      integer :: i,n

      !VENTRICOLO
      !Nodo1
      coords(1,1) =  0.02714
      coords(2,1) =  0.005196
      coords(3,1) =  3.06138
      !Nodo2
      coords(1,2) =  -0.3828
      coords(2,2) =  0.0
      coords(3,2) =  2.859
      !Nodo3
      coords(1,3) =  -0.8537
      coords(2,3) =  0.0
      coords(3,3) =  2.207
      !Nodo4
      coords(1,4) =  -0.9817
      coords(2,4) =  0.0
      coords(3,4) =  1.464
      !Nodo5
      coords(1,5) =  0.5632
      coords(2,5) =  0.0
      coords(3,5) =  2.77
      !Nodo6
      coords(1,6) =  1.103
      coords(2,6) =  0.0
      coords(3,6) =  2.062
      !Nodo7
      coords(1,7) =  1.403
      coords(2,7) =  0.0
      coords(3,7) =  1.174
      !Nodo8
      coords(1,8) =  1.222
      coords(2,8) =  0.0
      coords(3,8) =  0.4365
      !ATRIO
      !Nodo9
      coords(1,9) =  0.766
      coords(2,9) =  0.0
      coords(3,9) = -1.794
      !Nodo10
      coords(1,10) =  0.065
      coords(2,10) =  0.0
      coords(3,10) = -1.455
      !Nodo11
      coords(1,11) =  0.065
      coords(2,11) =  0.0
      coords(3,11) = -0.838
      !Nodo12
      coords(1,12) = -0.033
      coords(2,12) =  0.0
      coords(3,12) = -0.263
      !Nodo13
      coords(1,13) =  1.094
      coords(2,13) =  0.0
      coords(3,13) = -1.589
      !Nodo14
      coords(1,14) =  1.409
      coords(2,14) =  0.0
      coords(3,14) = -0.823
      !Nodo15
      coords(1,15) =  1.312
      coords(2,15) =  0.0
      coords(3,15) = -0.298

!    Apply translation
      
      do n=1,n_EFcheck
      Dmin=100._DP
         do i=1,nvtot
            D= sqrt((xyz0(1,i)-coords(1,n))**2 + &
                    (xyz0(2,i)-coords(2,n))**2 + &
                    (xyz0(3,i)-coords(3,n))**2)
            if (D.lt.Dmin) then
               Dmin = D
               nodicheckEF(n)  = i
            endif
         enddo
      enddo

!      do n=1,n_EFcheck
!         write(*,'(a,3(2x,f12.8))') "Looking for:",coords(1:3,n)
!         write(*,'(a,3(2x,f12.8))') "Found      :",xyz0(1:3, nodicheckEF(n))
!             !   ,"body:",   LabelV(nodicheckEF(n))
!      enddo


      return
      end

