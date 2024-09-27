! **********************  set boundary conditions  *********
!     Upper and lower wall boundary conditions
!
      subroutine boucdq
      use constants
      use param
      use local_arrays, only: q1,q2,q3,pr
      use mpi_param, only: kstart,kend
      use mpih
      use outflow_vars
#ifdef NCCLAVAIL
      use nccl
#endif      
!@cuf use cudafor
      implicit none
      integer i,j
      real(DP) dq1x2, dq2x2, dq3x2
      real(DP) area, darea
      real(DP) qout2, qoutf, cor
      real(DP) g3rcn3

      real(DP) :: velin(n1,n2)
!@cuf integer :: istat

      g3rcn3 = g3rc(n3)


!   Radiative B.C. at outflow (k=1)
!   A wave equation is solved to try to make the solution as
!   clean as possible. Courant number gives speed of waves.

      if(kend.eq.n3m) then
#ifdef USE_CUDA
      !$cuf kernel do (2)
#endif
      do j=1,n2m
        do i=1,n1m

!    radiative b.c. for x1 velocity at the outflow    
          dq1x2=(qb1n(i,j)-q1(i,j,n3m))*2.*dx3/g3rcn3
          dqb1n(i,j)=-dt*(ga*dq1x2+ro*dq1x2o(i,j))*cou

!    radiative b.c. for x2 velocity at the outflow    
          dq2x2=(qb2n(i,j)-q2(i,j,n3m))*2.*udx3c(n3)
          dqb2n(i,j)=-dt*(ga*dq2x2+ro*dq2x2o(i,j))*cou

!    radiative b.c. for x3 velocity at the outflow
          dq3x2=(qb3n(i,j)-q3(i,j,n3m))*udx3c(n3)
          dqb3n(i,j)=-dt*(ga*dq3x2+ro*dq3x2o(i,j))*cou

          dqb3s(i,j) = -qb3s(i,j)
          dqb2s(i,j) = 0.
          dqb1s(i,j) = 0.

!    save data for next time step
          dq1x2o(i,j)=dq1x2
          dq2x2o(i,j)=dq2x2
          dq3x2o(i,j)=dq3x2



        enddo   
      enddo   
!@cuf istat = cudaDeviceSynchronize !JDR TMP
      end if

      if (numtasks .gt. 1) then
#ifdef USE_CUDA
        if (allocated(sendbuff) .and. size(sendbuff) < n1*n2*3) deallocate(sendbuff)
        if (.not. allocated(sendbuff)) allocate(sendbuff(n1*n2*3))

        if (myid .eq. numtasks - 1) then
          !$cuf kernel do(2)
          do j = 1, n2
            do i = 1, n1
              sendbuff(i + (j-1) * n1 + 0 * n1*n2) = dqb1n(i,j)
              sendbuff(i + (j-1) * n1 + 1 * n1*n2) = dqb2n(i,j)
              sendbuff(i + (j-1) * n1 + 2 * n1*n2) = dqb3n(i,j)
            enddo
          enddo
        endif

!@cuf   istat = cudaDeviceSynchronize
#ifndef NCCLAVAIL        
        call MPI_BCAST(sendbuff,n1*n2*3,MDP,numtasks-1,MPI_COMM_WORLD,ierr)
#else
        nccl_result = ncclBroadcast(sendbuff,sendbuff,n1*n2*3,ncclDouble,numtasks-1,nccl_comm,0)
#endif
        
        if (myid .ne. numtasks - 1) then
          !$cuf kernel do(2)
          do j = 1, n2
            do i = 1, n1
              dqb1n(i,j) =  sendbuff(i + (j-1) * n1 + 0 * n1*n2)
              dqb2n(i,j) =  sendbuff(i + (j-1) * n1 + 1 * n1*n2)
              dqb3n(i,j) =  sendbuff(i + (j-1) * n1 + 2 * n1*n2)
            enddo
          enddo
        endif
#else
        call MPI_BCAST(dqb1n,n1*n2,MDP,numtasks-1,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dqb2n,n1*n2,MDP,numtasks-1,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dqb3n,n1*n2,MDP,numtasks-1,MPI_COMM_WORLD,ierr)
#endif
      end if


!    Flow is adjusted so that it is globally divergence free
!    First mass imbalance is calculated

      area = 0.
      qout=0.
      qinf=0.
      darea=1.0/dx2/dx1

#ifdef USE_CUDA
      !$cuf kernel do (2)
#endif
      do j=1,n2m
        do i=1,n1m
          area = area+darea
          qout=qout+dqb3n(i,j)*darea
          qinf=qinf+dqb3s(i,j)*darea
        end do
      end do

!@cuf istat = cudaDeviceSynchronize !JDR TMP
!   Then the correction is calculated and implemented at the outflow

      cor=(qinf-qout)/area
!      qout2=0.
!      qoutf=0.
!      darea=1.0/dx2/dx1
#ifdef USE_CUDA
      !$cuf kernel do (2)
#endif
      do j=1,n2m
        do i=1,n1m
         dqb3n(i,j)=dqb3n(i,j)+cor
        end do
      end do

!@cuf istat = cudaDeviceSynchronize !JDR TMP
      return
      end
!     =======================================================================


!     =======================================================================
      subroutine boucqt
      use constants
      use param
      use outflow_vars
!@cuf use cudafor
      IMPLICIT NONE
      real(DP) :: darea
      integer :: i,j
!@cuf integer :: istat
 

!   Inflow/outflow velocities are updated, and mass imbalance is checked
#ifdef USE_CUDA
      !$cuf kernel do(2)
#endif
      do j=1,n2m
        do i=1,n1m
          qb1s(i,j) = qb1s(i,j)+dqb1s(i,j)
          qb2s(i,j) = qb2s(i,j)+dqb2s(i,j)
          qb3s(i,j) = qb3s(i,j)+dqb3s(i,j)
          qb1n(i,j) = qb1n(i,j)+dqb1n(i,j)
          qb2n(i,j) = qb2n(i,j)+dqb2n(i,j)
          qb3n(i,j) = qb3n(i,j)+dqb3n(i,j)
        enddo
      enddo

!@cuf istat = cudaDeviceSynchronize !JDR TMP

      return
      end

!     =======================================================================
      subroutine FixedBodiesPRib
      use param
      use outflow_vars
      use ibm_param
      use local_arrays, only: pr
!@cuf use cudafor

      integer :: n, i, j, k, ie, je, ke
!@cuf integer :: istat

! #ifdef USE_CUDA
!       !$cuf kernel do(1)
! #endif
!       do n=1,npunfx(4)
!            i=indgeo(4,n,1)
!            j=indgeo(4,n,2)
!            k=indgeo(4,n,3)
!            ie=indgeoe(4,n,1)
!            je=indgeoe(4,n,2)
!            ke=indgeoe(4,n,3)
!            pr(i,j,k)=pr(ie,je,ke)
!       end do

! #ifdef USE_CUDA
!       !$cuf kernel do(1)
! #endif
!       do n=1,npunifx(4)
!          i=indgeoee(4,n,1)
!          j=indgeoee(4,n,2)
!          k=indgeoee(4,n,3)
!          pr(i,j,k)=0.
!       end do

! !francesco--opening/closing valves:Mitral
! #ifdef USE_CUDA
!       !$cuf kernel do(1)
! #endif
!       do n=1,npunfxMI(4)
!            i=indgeoMI(4,n,1)
!            j=indgeoMI(4,n,2)
!            k=indgeoMI(4,n,3)
!            ie=indgeoeMI(4,n,1)
!            je=indgeoeMI(4,n,2)
!            ke=indgeoeMI(4,n,3)
!            pr(i,j,k)=(1.-ftMI)*pr(i,j,k) + &
!                           ftMI*pr(ie,je,ke)
!       end do

! #ifdef USE_CUDA
!       !$cuf kernel do(1)
! #endif
!       do n=1,npunifxMI(4)
!          i=indgeoeeMI(4,n,1)
!          j=indgeoeeMI(4,n,2)
!          k=indgeoeeMI(4,n,3)
!          pr(i,j,k)=(1.-ftMI)*pr(i,j,k)
!       end do

! !francesco--opening/closing valves:Aorta
! #ifdef USE_CUDA
!       !$cuf kernel do(1)
! #endif
!       do n=1,npunfxAO(4)
!            i=indgeoAO(4,n,1)
!            j=indgeoAO(4,n,2)
!            k=indgeoAO(4,n,3)
!            ie=indgeoeAO(4,n,1)
!            je=indgeoeAO(4,n,2)
!            ke=indgeoeAO(4,n,3)
!            pr(i,j,k)=(1.-ftAO)*pr(i,j,k) + &
!                           ftAO*pr(ie,je,ke)
!       end do

! #ifdef USE_CUDA
!       !$cuf kernel do(1)
! #endif
!       do n=1,npunifxAO(4)
!          i=indgeoeeAO(4,n,1)
!          j=indgeoeeAO(4,n,2)
!          k=indgeoeeAO(4,n,3)
!          pr(i,j,k)=(1.-ftAO)*pr(i,j,k)
!       end do

!francesco--opening/closing valves:AortaSistema
#ifdef USE_CUDA
      !$cuf kernel do(1)
#endif
      do n=1,npunfxAS(4)
           i=indgeoAS(4,n,1)
           j=indgeoAS(4,n,2)
           k=indgeoAS(4,n,3)
           ie=indgeoeAS(4,n,1)
           je=indgeoeAS(4,n,2)
           ke=indgeoeAS(4,n,3)
           pr(i,j,k)=(1.-ftAS)*pr(i,j,k) + &
                          ftAS*pr(ie,je,ke)
      end do

#ifdef USE_CUDA
      !$cuf kernel do(1)
#endif
      do n=1,npunifxAS(4)
         i=indgeoeeAS(4,n,1)
         j=indgeoeeAS(4,n,2)
         k=indgeoeeAS(4,n,3)
         pr(i,j,k)=(1.-ftAS)*pr(i,j,k)
      end do

!@cuf istat = cudaDeviceSynchronize !JDR TMP

      return
      end

!     -----------------------------------------------------------------
!     -----------------------------------------------------------------

      subroutine ContactModel(snv,env,xyzL,xyzvL)
      use mls_param
      use mls_local, only: coll
!@cuf use cudafor
      implicit none
      integer :: snv,env
!      real(DP),dimension(3,snv:env),intent(in) :: xyzL,xyzvL
      real(DP),dimension(3,snv:env),intent(inout) :: xyzL,xyzvL
!      real(DP),dimension(3,snv:env):: xyzL,xyzvL
      integer :: ntr,ps,pe
      integer :: ci,cj,ck
      integer :: v1,v2,v3
      integer :: tr, inp
!@cuf integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: xyzL, xyzvL
#endif

      fcol = 0
      
!     needs to be improved
      if(ContactA.eq.1.and.ContactM.eq.1) then
        ps = 3
        pe = 7
      elseif(ContactA.eq.1.and.ContactM.eq.0) then
        ps = 3
        pe = 5
      elseif(ContactA.eq.0.and.ContactM.eq.1)then
        ps = 6
        pe = 7
      end if

      if(ContactA.gt.0.or.ContactM.gt.0) then

        !JDR TODO: move setting of fcol into mlsForceTiled
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do tr = 1, trcnt_h
          ntr = mytr(tr)

          inp = face_to_part(ntr)

          if (inp .ge. ps .and. inp .le. pe) then
            ci = pind(4,ntr)
            cj = pind(5,ntr)
            ck = pind(6,ntr)
            if(coll(ci,cj,ck).eq.-1) fcol(ntr) = 1
          endif
        end do

!@cuf   istat = cudaDeviceSynchronize

        call mpi_globalsum_integer_arr(fcol,nftot)

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=fstart(ps),fend(pe)
          if(fcol(ntr).eq.1) then

             v1 = vert_of_face(1,ntr)
             v2 = vert_of_face(2,ntr)
             v3 = vert_of_face(3,ntr)

             !move back all vertices of face to old posn
             xyzL(1:3,v1) = xyzold(1:3,v1)
             xyzL(1:3,v2) = xyzold(1:3,v2)
             xyzL(1:3,v3) = xyzold(1:3,v3)

             xyzvL(1:3,v1) = 0.d0
             xyzvL(1:3,v2) = 0.d0
             xyzvL(1:3,v3) = 0.d0
           end if
        end do

!@cuf   istat = cudaDeviceSynchronize

      end if

      return
      end

!     -----------------------------------------------------------------
      subroutine ContactModelV(snv,env,xyzL,xyzvL)
      use mls_param
      use mls_local, only: coll
!@cuf use cudafor
      implicit none
      integer :: snv,env
!      real(DP),dimension(3,snv:env),intent(in) :: xyzL,xyzvL
      real(DP),dimension(3,snv:env) :: xyzL,xyzvL
      integer :: ntr,ps,pe,inp,leaflet
      integer :: ci,cj,ck
      integer :: v1,v2,v3,vertice
!@cuf integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: xyzL, xyzvL
#endif


        fcolv = 0 

!AORTIC VALVE      
      if (ContactA.EQ.1) then      
!PRIMA
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(2),vend(2) 
          if(dum_forv(ntr).eq.0)then

!            inp = vert_to_part(ntr)
            ci = pindv(4,ntr)
            cj = pindv(5,ntr)
            ck = pindv(6,ntr)

            leaflet = LabelVaortic(ntr)
            if (leaflet.EQ.1) then
               if(coll(ci,cj,ck).ne.leaflet.and.coll(ci,cj,ck).gt.0) then
                  fcolv(ntr) = 1
                  coll(ci,cj,ck) = leaflet
               endif
            endif
            
          end if
        end do
!@cuf   istat = cudaDeviceSynchronize

!SECONDA
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(2),vend(2) 
          if(dum_forv(ntr).eq.0)then

!            inp = vert_to_part(ntr)
            ci = pindv(4,ntr)
            cj = pindv(5,ntr)
            ck = pindv(6,ntr)

            leaflet = LabelVaortic(ntr)
            if (leaflet.EQ.2) then
               if(coll(ci,cj,ck).ne.leaflet.and.coll(ci,cj,ck).gt.0) then
                  fcolv(ntr) = 1
                  coll(ci,cj,ck) = leaflet
               endif
            endif
            
          end if
        end do
!@cuf   istat = cudaDeviceSynchronize

!TERZA
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(2),vend(2) 
          if(dum_forv(ntr).eq.0)then

!            inp = vert_to_part(ntr)
            ci = pindv(4,ntr)
            cj = pindv(5,ntr)
            ck = pindv(6,ntr)

            leaflet = LabelVaortic(ntr)
            if (leaflet.EQ.3) then
               if(coll(ci,cj,ck).ne.leaflet.and.coll(ci,cj,ck).gt.0) then
                  fcolv(ntr) = 1
                  coll(ci,cj,ck) = leaflet
               endif
            endif
            
          end if
        end do
!@cuf   istat = cudaDeviceSynchronize

        call mpi_globalsum_integer_arr(fcolv,nvtot)      

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(2),vend(2)
          if(fcolv(ntr).eq.1) then
             !move back all vertices of face to old posn
             xyzL(1:3,ntr) = xyzold(1:3,ntr)
             xyzvL(1:3,ntr) = 0.d0
          end if
        end do
!@cuf   istat = cudaDeviceSynchronize

     endif !ContactA



!MITRAL VALVE      
      if (ContactM.EQ.1) then      
!PRIMA
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(3),vend(3) 
          if(dum_forv(ntr).eq.0)then

!            inp = vert_to_part(ntr)
            ci = pindv(4,ntr)
            cj = pindv(5,ntr)
            ck = pindv(6,ntr)

            leaflet = LabelVmitral(ntr)
            if (leaflet.EQ.1) then
               if(coll(ci,cj,ck).ne.leaflet.and.coll(ci,cj,ck).gt.0) then
                  fcolv(ntr) = 1
                  coll(ci,cj,ck) = leaflet
               endif
            endif
            
          end if
        end do
!@cuf   istat = cudaDeviceSynchronize

!SECONDA
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(3),vend(3) 
          if(dum_forv(ntr).eq.0)then

!            inp = vert_to_part(ntr)
            ci = pindv(4,ntr)
            cj = pindv(5,ntr)
            ck = pindv(6,ntr)

            leaflet = LabelVmitral(ntr)
            if (leaflet.EQ.2) then
               if(coll(ci,cj,ck).ne.leaflet.and.coll(ci,cj,ck).gt.0) then
                  fcolv(ntr) = 1
                  coll(ci,cj,ck) = leaflet
               endif
            endif
            
          end if
        end do
!@cuf   istat = cudaDeviceSynchronize

        call mpi_globalsum_integer_arr(fcolv,nvtot)      

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(3),vend(3)
          if(fcolv(ntr).eq.1) then
             !move back all vertices of face to old posn
             xyzL(1:3,ntr) = xyzold(1:3,ntr)
             xyzvL(1:3,ntr) = 0.d0
          end if
        end do
!@cuf   istat = cudaDeviceSynchronize

        
     endif !ContactM
     


!PULMONARY VALVE      
      if (ContactP.EQ.1) then      
!PRIMA
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(4),vend(4) 
          if(dum_forv(ntr).eq.0)then

!            inp = vert_to_part(ntr)
            ci = pindv(4,ntr)
            cj = pindv(5,ntr)
            ck = pindv(6,ntr)

            leaflet = LabelVpulmo(ntr)
            if (leaflet.EQ.1) then
               if(coll(ci,cj,ck).ne.leaflet.and.coll(ci,cj,ck).gt.0) then
                  fcolv(ntr) = 1
                  coll(ci,cj,ck) = leaflet
               endif
            endif
            
          end if
        end do
!@cuf   istat = cudaDeviceSynchronize

!SECONDA
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(4),vend(4) 
          if(dum_forv(ntr).eq.0)then

!            inp = vert_to_part(ntr)
            ci = pindv(4,ntr)
            cj = pindv(5,ntr)
            ck = pindv(6,ntr)

            leaflet = LabelVpulmo(ntr)
            if (leaflet.EQ.2) then
               if(coll(ci,cj,ck).ne.leaflet.and.coll(ci,cj,ck).gt.0) then
                  fcolv(ntr) = 1
                  coll(ci,cj,ck) = leaflet
               endif
            endif
            
          end if
        end do
!@cuf   istat = cudaDeviceSynchronize

!TERZA
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(4),vend(4) 
          if(dum_forv(ntr).eq.0)then

!            inp = vert_to_part(ntr)
            ci = pindv(4,ntr)
            cj = pindv(5,ntr)
            ck = pindv(6,ntr)

            leaflet = LabelVpulmo(ntr)
            if (leaflet.EQ.3) then
               if(coll(ci,cj,ck).ne.leaflet.and.coll(ci,cj,ck).gt.0) then
                  fcolv(ntr) = 1
                  coll(ci,cj,ck) = leaflet
               endif
            endif
            
          end if
        end do
!@cuf   istat = cudaDeviceSynchronize

        call mpi_globalsum_integer_arr(fcolv,nvtot)      

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(4),vend(4)
          if(fcolv(ntr).eq.1) then
             !move back all vertices of face to old posn
             xyzL(1:3,ntr) = xyzold(1:3,ntr)
             xyzvL(1:3,ntr) = 0.d0
          end if
        end do
!@cuf   istat = cudaDeviceSynchronize
        
     endif !ContactP


!TRICUSPID VALVE      
      if (ContactT.EQ.1) then      
!PRIMA
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(5),vend(5) 
          if(dum_forv(ntr).eq.0)then

!            inp = vert_to_part(ntr)
            ci = pindv(4,ntr)
            cj = pindv(5,ntr)
            ck = pindv(6,ntr)

            leaflet = LabelVtricu(ntr)
            if (leaflet.EQ.1) then
               if(coll(ci,cj,ck).ne.leaflet.and.coll(ci,cj,ck).gt.0) then
                  fcolv(ntr) = 1
                  coll(ci,cj,ck) = leaflet
               endif
            endif
            
          end if
        end do
!@cuf   istat = cudaDeviceSynchronize

!SECONDA
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(5),vend(5) 
          if(dum_forv(ntr).eq.0)then

!            inp = vert_to_part(ntr)
            ci = pindv(4,ntr)
            cj = pindv(5,ntr)
            ck = pindv(6,ntr)

            leaflet = LabelVtricu(ntr)
            if (leaflet.EQ.2) then
               if(coll(ci,cj,ck).ne.leaflet.and.coll(ci,cj,ck).gt.0) then
                  fcolv(ntr) = 1
                  coll(ci,cj,ck) = leaflet
               endif
            endif
            
          end if
        end do
!@cuf   istat = cudaDeviceSynchronize

!TERZA
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(5),vend(5) 
          if(dum_forv(ntr).eq.0)then

!            inp = vert_to_part(ntr)
            ci = pindv(4,ntr)
            cj = pindv(5,ntr)
            ck = pindv(6,ntr)

            leaflet = LabelVtricu(ntr)
            if (leaflet.EQ.3) then
               if(coll(ci,cj,ck).ne.leaflet.and.coll(ci,cj,ck).gt.0) then
                  fcolv(ntr) = 1
                  coll(ci,cj,ck) = leaflet
               endif
            endif
            
          end if
        end do
!@cuf   istat = cudaDeviceSynchronize

        call mpi_globalsum_integer_arr(fcolv,nvtot)      

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(5),vend(5)
          if(fcolv(ntr).eq.1) then
             !move back all vertices of face to old posn
             xyzL(1:3,ntr) = xyzold(1:3,ntr)
             xyzvL(1:3,ntr) = 0.d0
          end if
        end do
!@cuf   istat = cudaDeviceSynchronize
        
     endif !ContactT



!AORTICA -> AORTA
      if (ContactA.EQ.0) then      
!PRIMA
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(2),vend(2) 
          if(dum_forv(ntr).eq.0)then

!            inp = vert_to_part(ntr)
            ci = pindv(4,ntr)
            cj = pindv(5,ntr)
            ck = pindv(6,ntr)

            if(coll(ci,cj,ck).EQ.-5) then
               fcolv(ntr) = 1
!               coll(ci,cj,ck) = leaflet
            endif

          end if
        end do
!@cuf   istat = cudaDeviceSynchronize


        call mpi_globalsum_integer_arr(fcolv,nvtot)      

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(2),vend(2)
          if(fcolv(ntr).eq.1) then
             !move back all vertices of face to old posn
             xyzL(1:3,ntr) = xyzold(1:3,ntr)
             xyzvL(1:3,ntr) = 0.d0
          end if
        end do
!@cuf   istat = cudaDeviceSynchronize
        
     endif

        
     !POLMONARE -> ARTERIA POLMONARE
!PULMONARY VALVE      
      if (ContactP.EQ.0) then      
!PRIMA
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(4),vend(4) 
          if(dum_forv(ntr).eq.0)then

!            inp = vert_to_part(ntr)
            ci = pindv(4,ntr)
            cj = pindv(5,ntr)
            ck = pindv(6,ntr)

            if(coll(ci,cj,ck).EQ.-7) then
               fcolv(ntr) = 1
!               coll(ci,cj,ck) = leaflet
            endif

          end if
        end do
!@cuf   istat = cudaDeviceSynchronize

        call mpi_globalsum_integer_arr(fcolv,nvtot)      

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do ntr=vstart(4),vend(4)
          if(fcolv(ntr).eq.1) then
             !move back all vertices of face to old posn
             xyzL(1:3,ntr) = xyzold(1:3,ntr)
             xyzvL(1:3,ntr) = 0.d0
          end if
        end do
!@cuf   istat = cudaDeviceSynchronize
        
     endif !ContactP pulmo art po 
        
      return
      end


!     -----------------------------------------------------------------

