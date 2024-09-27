!
!***********************************************************************
!
!           SUBROUTINE  TSCHEM
!
!   This subroutine manages the whole integration scheme.
!   The following equations are solved:          
!   
!    ~~     n
!   Q  -  Q                n         n       n-1   alp       2  ~~   n 
!  --------- = -alp*grad (P ) + gam*H + rho*H   + ----- nabla ( Q + Q )
!    d t                                          2 Re
!
!          i                           i               i
!   where H  are the nonlinear terms, P  the pressure Q  the velocities
!       ~~
!   and Q  the provisional non solenoidal velocity field.
!   The superscripts (~~, n, n-1) indicate the time step level. 
!                        n
!   The nonlinear terms H  are computed in the routines HDNL*, while
!   in the routines INVTR* are computed the remaining terms, updated
!   the non linear terms and inverted the equation to find the provisional
!   field at the new time step.
!       ~~
!   The Q  velocity field is projected onto a solenoidal field by a 
!   scalar Phi computed through the equation
!
!                         2            1          ~~
!                    nabla (Phi ) =  ------ div ( Q  )
!                                    alp dt
!
!   The right hand side of this equation is computed in the routine
!   DIVG, while the equation is solved in PHCALC.
!
!   In the routine UPDVP the solenoidal velocity field at the new time
!   step is then computed through
!
!                n+1  ~~
!               Q   = Q  - alt*dt grad (Phi)
!
!   Finally in the routine PRCALC is updated the pressure field
!
!                n+1   n        alp dt      2
!               P   = P + Phi - ------ nabla (Phi)
!                                2 Re
!
!   When the scalar field is computed (density, concentration,
!   temperature) the routines HDNLRO and INVTRRO are used. The same
!   strategy at the velocity field is used, except that the scalar
!   field does not need any correction.
!
!   All variables are located on a staggered grid with the velocities
!   on the faces of the computational cell and all the scalars at the
!   centre. This is important when terms belonging to different equations
!   are avaluated.
!
!   Further details of the scheme can be found in the paper
!   "A finite-difference scheme for three-dimensional incompressible
!    flows in cylindrical coordinates" by R. Verzicco and P. Orlandi
!    J. of Comp. Phys. 1996.
!
!
      subroutine tschemStrong
!@cuf use cudafor
      use param
      use local_arrays
      use mgrd_arrays
      use mpih
      use mpi_param, only: kstart,kend,kstartr,kendr
      use mls_param     
      use mls_local, only: for_xc,for_yc,for_zc, for_sc, coll
      use probes
      use nvtx
#ifdef NCCLAVAIL
      use nccl
#endif      
      implicit none
!@cuf integer :: istat
      integer :: ns,inp,ntr,nsub
      integer :: j,k,i,v1,ic,jc,kc
      integer :: vsi, vei, esi, eei, fsi, fei
      real(DP) :: dinvm
      integer :: jjj,jjjdiv
      character*50 filename,strucfilename
#ifdef DEBUG
      integer :: kstartp
      real(DP) :: cksum2,cksum3,mck2,mck3,cksum1,mck1
#endif


      call preambolotsch

      dts = dt/dble(lmax)
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!!!!!!
!STRONG Store n-th configuration before iterating
#ifdef USE_CUDA
      !$cuf kernel do (3)
      do kc=kstart-lvlhalo,kend+lvlhalo      
        do jc=1,n2
          do ic=1,n1
#else
      do kc=kstart-lvlhalo,kend+lvlhalo
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,ic)
        do jc=1,n2m
          do ic=1,n1m
#endif
          q1g(ic,jc,kc)=q1(ic,jc,kc)
          q2g(ic,jc,kc)=q2(ic,jc,kc)
          q3g(ic,jc,kc)=q3(ic,jc,kc)
          prg(ic,jc,kc)=pr(ic,jc,kc)

          dphg(ic,jc,kc)=dph(ic,jc,kc)                              
        enddo 
       enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      enddo

!@cuf istat = cudaDeviceSynchronize() !JDR TMP


#ifdef USE_CUDA
      !$cuf kernel do (3)
      do kc=kstart,kend
        do jc=1,n2
          do ic=1,n1
#else
      do kc=kstart,kend
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,ic)
        do jc=1,n2
          do ic=1,n1
#endif
          ru1g(ic,jc,kc)=ru1(ic,jc,kc)
          ru2g(ic,jc,kc)=ru2(ic,jc,kc)
          ru3g(ic,jc,kc)=ru3(ic,jc,kc)

          qcapg(ic,jc,kc)=qcap(ic,jc,kc)
          dqg(ic,jc,kc)=dq(ic,jc,kc)                    
        enddo 
       enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      enddo

!@cuf istat = cudaDeviceSynchronize() !JDR TMP

      

#ifdef USE_CUDA
      !$cuf kernel do (2)
        do jc=1,n2
          do ic=1,n1
#else
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,ic)
        do jc=1,n2
          do ic=1,n1
#endif
          qb1sg(ic,jc)=qb1s(ic,jc)
          qb2sg(ic,jc)=qb2s(ic,jc)
          qb3sg(ic,jc)=qb3s(ic,jc)
          qb1ng(ic,jc)=qb1n(ic,jc)
          qb2ng(ic,jc)=qb2n(ic,jc)
          qb3ng(ic,jc)=qb3n(ic,jc)
          dq1x2og(ic,jc)=dq1x2o(ic,jc)
          dq2x2og(ic,jc)=dq2x2o(ic,jc)
          dq3x2og(ic,jc)=dq3x2o(ic,jc)
        enddo 
       enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif

!@cuf istat = cudaDeviceSynchronize() !JDR TMP
      
!!!!!!!!!!!!!!!!!!!!



!====== Predictor ======================================================== 
!     Compute pressure and viscous loads
        call nvtxStartRange("mlsStruc", 29)
        call mlsStrucPredictorNewton
        call nvtxEndRange


      if(imlsfor.eq.1)then
!QUI DEVO USARE xyzp xyzvp xyza
!       Computing normals and triangle properties
        call nvtxStartRange("convert_geo", 32)
        call convert_geo(1, nvtot, 1, netot, 1, nftot, xyzp, xyzvp, xyza, &
                vert_of_face, tri_ver, tri_vel, tri_bar, vel_tri, acc_tri)
        call nvtxEndRange

        do inp = 1,Nparticle

          call calculate_area(Surface(inp), vstart(inp), vend(inp), &
                   fstart(inp), fend(inp), xyzp(:,vstart(inp):vend(inp)), &
                   vert_of_face(:,fstart(inp):fend(inp)), &
                   sur(fstart(inp):fend(inp)))

        end do
!       Find indices of all centroids of particles and bounding box
        call nvtxStartRange("findindices", 31)
        call findindicesFaces
        call findindicesNodes
        call findindicesTiled
#ifdef OFFSETBODY        
        call indindicesTiledOff
#endif
        call findindicesFacesProbeP
        call findindicesFacesProbeN        
        call nvtxEndRange
        endif



!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!      REAL timef !etime_,t(2)
!
!   TIME INTEGRATION : implicit viscous, 3rd order RK (Adams Bashfort)  
!                                                                       

      DO ns=1,nsst

        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)

        if(ntime.le.1)then
            aldto = alm(1)*dt
        else
            aldto = al*dt
        end if

        call nvtxStartRange("boucdq", 30)
        call boucdq
        call nvtxEndRange

        !print*, "q1:", sum(q1(:,:,kstart:kend)), minval(q1(:,:,kstart:kend)), maxval(q1(:,:,kstart:kend))
        !print*, "q2:", sum(q2(:,:,kstart:kend)), minval(q2(:,:,kstart:kend)), maxval(q2(:,:,kstart:kend))
        !print*, "q3:", sum(q3(:,:,kstart:kend)), minval(q3(:,:,kstart:kend)), maxval(q3(:,:,kstart:kend))
!=========================================
!  explicit terms
        call nvtxStartRange("hdnlq123", 9)
        call hdnlq1


        call hdnlq2


        call hdnlq3
        call nvtxEndRange

!@cuf   istat = cudaDeviceSynchronize() !JDR TMP

#ifdef CALTEMP
        call hdnlte
#endif

#ifdef CALSCAL
        call hdnlsa
#endif


!====================================
!  implicit terms
        call nvtxStartRange("invtrq123", 12)
        call invtrq1


        call invtrq2                     !! DQ2HAT=Q2HAT-Q2(N)


        call invtrq3     !! DQ3HAT=Q3HAT-Q3(N)
        call nvtxEndRange
!        call update_upper_ghost(n1,n2,q3)


#ifdef CALTEMP
        call invtrte                     !! DENSITY  ==>>>
        call update_both_ghosts(n1,n2,dens,kstart,kend)
#endif


#ifdef CALSCAL
        call invtrsa                  !! SALINITY  ==>>>
        call update_both_ghosts(n1r,n2r,dsal,kstartr,kendr)
#endif


      call nvtxStartRange("update_both_ghosts", 20)
#ifdef USE_CUDA
      call update_both_ghosts_3(n1,n2,q1,q2,q3,kstart,kend)
#else
      call update_both_ghosts(n1,n2,q1,kstart,kend)
      call update_both_ghosts(n1,n2,q2,kstart,kend)
      call update_both_ghosts(n1,n2,q3,kstart,kend)
#endif
      call nvtxEndRange
!     ============================================
!     Forcing for fluid-structure interaction
!     ============================================

      do jjjdiv=1,nsstepdiv
       if(imlsfor.eq.1) then

!     simple mlsforce without tiling
        !call mlsForce

!     call this instead of mlsforce for tiled forcing
      do jjj=1,nsstep
        call nvtxStartRange("zero forces", 20)
        for_xc = 0.0 ; for_yc = 0.0 ; for_zc = 0.0
        call nvtxEndRange


#ifdef OFFSETBODY
!     ========OffBody Face-Centered-Forcing=========
        call mlsForceTiledOff3Comp
#endif
!     ========Face-Centered-Forcing=========
        call nvtxStartRange("mlsForceTiled3Comp", 16)
        call mlsForceTiled3Comp
        call nvtxEndRange

        call nvtxStartRange("velforce3Comp", 17)
        call velforce3Comp
        call nvtxEndRange

      end do

        !for_sc = 0.0
        ! call mlsForce_scalar
        ! call velforce_scalar

        call nvtxStartRange("sync_collision", 18)
        call sync_collision(n1,n2,coll)
        call nvtxEndRange
       end if
!   ======================================================
!   End of forcing for Fluid-Structure interaction
!   ======================================================


!====================================
! pressure and velocity correction

        call update_upper_ghost(n1,n2,q3)

        call nvtxStartRange("divg", 18)
        call divg        !here only q3 ghost(up) cell are needed
        call nvtxEndRange

        
        call nvtxStartRange("phcalc", 19)
        call phcalc
        call nvtxEndRange

        call update_both_ghosts(n1m,n2m+2,dph,kstart,kend)
        
        call nvtxStartRange("updvp", 21)
        call updvp                 !! SOLENOIDAL VEL FIELD
        call nvtxEndRange


        call nvtxStartRange("velbc", 22)
        call velbc                 !! velocity boundary condition 
        call nvtxEndRange

#ifdef USE_CUDA
        call update_both_ghosts_3(n1,n2,q1,q2,q3,kstart,kend)
#else
        call update_both_ghosts(n1,n2,q1,kstart,kend)
        call update_both_ghosts(n1,n2,q2,kstart,kend)
        call update_both_ghosts(n1,n2,q3,kstart,kend)
#endif
        
        call nvtxStartRange("prcalc", 23)
        call prcalc                         !! PRESSURE FIELD
        call nvtxEndRange

        call update_both_ghosts(n1,n2,pr,kstart,kend)

  

!==================================================
!     multi resolution information exchange
!   expanding vel field from q* to q*lr
        !call mgrd_velitp

#ifdef CALSCAL
!      salinity from refined mesh to base mesh
        call mgrd_dsalc
#endif

      if (infig .eq. 1) then
        call nvtxStartRange("FixedBodies", 25)
        call FixedBodiesPRib
        call nvtxEndRange
        call update_both_ghosts(n1,n2,pr,kstart,kend)
      endif
 
   enddo!jjjdiv nsstepdiv

      call nvtxStartRange("boucqt", 24)
      call boucqt
      call nvtxEndRange

      ENDDO         ! end of RK loop

!================================
!     CALLS FOR STRUCTURAL SOLVER
     
      if(imlsstr.eq.1.and.ntime.ge.3)then

        ! Note: not quite sure what this does?
!         if (time .lt. 2.5 .and. pressAO .lt. 28.2) then
!           prete = prete + 10.*1.7* 5d-5/3.d0 * (dt/5d-5)
! #ifdef USE_CUDA
!           !$cuf kernel do (1)
! #endif
!           do i = estart(2), eend(2)
!              v1 = vert_of_edge(1,i)
!              if ((xyz(2,v1).LT.-1.d0).AND.(xyz(3,v1).GT.1.d0)) then
!                dist0(i) = (1.d0 - prete) * dist00(i)
!              endif
!           end do
! !@cuf     istat = cudaDeviceSynchronize
!         endif



        ! NOTE: calculating normalized normal here for fpxyz in mlsStruc
        call nvtxStartRange("calculate_normal", 33)
        call calculate_normal(1, nvtot, 1, nftot, xyzp, &
                vert_of_face, tri_nor)
        call nvtxEndRange

        call mlsStruc3CompStrong

!     Compute internal force
        call internalforce(1,nvtot,xyzp,xyzvp)
        call internalforce_3d(1,nvtot_3d,xyzp_3d,xyzvp_3d)



!     Reduce the forces from all processors over each particle
!       do ntr=1,maxnv
!       call mpi_globalsum_double_arr(fxyz(1,ntr,:),Nparticle)
!       call mpi_globalsum_double_arr(fxyz(2,ntr,:),Nparticle)
!       call mpi_globalsum_double_arr(fxyz(3,ntr,:),Nparticle)
!       end do

        ! Reduce the forces from all processors over each particle
        if (numtasks > 1) then
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#ifndef NCCLAVAIL
          call MPI_ALLREDUCE(MPI_IN_PLACE,fxyz,3*nvtot,MPI_DOUBLE_PRECISION, &
                  MPI_SUM,MPI_COMM_WORLD,ierr)
#else
          nccl_result = ncclAllReduce(fxyz,fxyz,3*nvtot,ncclDouble,ncclSum,nccl_comm,0)
#endif
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)

#ifndef NCCLAVAIL
          call MPI_ALLREDUCE(MPI_IN_PLACE,fxyz_3d,3*nvtot_3d,MPI_DOUBLE_PRECISION, &
                  MPI_SUM,MPI_COMM_WORLD,ierr)
#else
          nccl_result = ncclAllReduce(fxyz_3d,fxyz_3d,3*nvtot_3d,ncclDouble,ncclSum,nccl_comm,0)
#endif
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        endif


!store stuff
#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = 1, nvtot_3d
         dinvm = 1.d0/mass_of_vert_3d(ntr)
         xyzak_3d(1,ntr)=(fpxyz_3d(1,ntr)+fxyz_3d(1,ntr))*dinvm
         xyzak_3d(2,ntr)=(fpxyz_3d(2,ntr)+fxyz_3d(2,ntr))*dinvm
         xyzak_3d(3,ntr)=(fpxyz_3d(3,ntr)+fxyz_3d(3,ntr))*dinvm

         xyzakm1_3d(1,ntr)=xyzak_3d(1,ntr)
         xyzakm1_3d(2,ntr)=xyzak_3d(2,ntr)
         xyzakm1_3d(3,ntr)=xyzak_3d(3,ntr)

         xyzvk_3d(1,ntr)=xyzvp_3d(1,ntr)
         xyzvk_3d(2,ntr)=xyzvp_3d(2,ntr)
         xyzvk_3d(3,ntr)=xyzvp_3d(3,ntr)

         xyzvkm1_3d(1,ntr)=xyzvk_3d(1,ntr)
         xyzvkm1_3d(2,ntr)=xyzvk_3d(2,ntr)
         xyzvkm1_3d(3,ntr)=xyzvk_3d(3,ntr)

         xyzk_3d(1,ntr)=xyzp_3d(1,ntr)
         xyzk_3d(2,ntr)=xyzp_3d(2,ntr)
         xyzk_3d(3,ntr)=xyzp_3d(3,ntr)
      enddo

#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = nvstart_2dwet,nvend_2dwet
         v1 = tag_2dwet(ntr)
         xyzak(1,ntr)=xyzak_3d(1,v1)
         xyzak(2,ntr)=xyzak_3d(1,v1)
         xyzak(3,ntr)=xyzak_3d(1,v1)

         xyzakm1(1,ntr)=xyzakm1_3d(1,v1)
         xyzakm1(2,ntr)=xyzakm1_3d(2,v1)
         xyzakm1(3,ntr)=xyzakm1_3d(3,v1)

         xyzvk(1,ntr)=xyzvk_3d(1,v1)
         xyzvk(2,ntr)=xyzvk_3d(2,v1)
         xyzvk(3,ntr)=xyzvk_3d(3,v1)

         xyzvkm1(1,ntr)=xyzvkm1_3d(1,v1)
         xyzvkm1(2,ntr)=xyzvkm1_3d(2,v1)
         xyzvkm1(3,ntr)=xyzvkm1_3d(3,v1)

         xyzk(1,ntr)=xyzk_3d(1,v1)
         xyzk(2,ntr)=xyzk_3d(2,v1)
         xyzk(3,ntr)=xyzk_3d(3,v1)
      enddo

#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = nvstart_2dstr,nvend_2dstr
         dinvm = 1.d0/mass_of_vert(ntr)
         xyzak(1,ntr)=(fpxyz(1,ntr)+fxyz(1,ntr))*dinvm
         xyzak(2,ntr)=(fpxyz(2,ntr)+fxyz(2,ntr))*dinvm
         xyzak(3,ntr)=(fpxyz(3,ntr)+fxyz(3,ntr))*dinvm

         xyzakm1(1,ntr)=xyzak(1,ntr)
         xyzakm1(2,ntr)=xyzak(2,ntr)
         xyzakm1(3,ntr)=xyzak(3,ntr)

         xyzvk(1,ntr)=xyzvp(1,ntr)
         xyzvk(2,ntr)=xyzvp(2,ntr)
         xyzvk(3,ntr)=xyzvp(3,ntr)

         xyzvkm1(1,ntr)=xyzvk(1,ntr)
         xyzvkm1(2,ntr)=xyzvk(2,ntr)
         xyzvkm1(3,ntr)=xyzvk(3,ntr)

         xyzk(1,ntr)=xyzp(1,ntr)
         xyzk(2,ntr)=xyzp(2,ntr)
         xyzk(3,ntr)=xyzp(3,ntr)
      enddo
!end store stuff

!====== Corrector ========================================================                  
       do iter=1,n_max_iter
           !------ Recall n-th configuration -------------------------------   
#ifdef USE_CUDA
      !$cuf kernel do (3)
      do kc=kstart-lvlhalo,kend+lvlhalo 
        do jc=1,n2
          do ic=1,n1
#else
      do kc=kstart-lvlhalo,kend+lvlhalo 
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,ic)
        do jc=1,n2
          do ic=1,n1
#endif
          q1(ic,jc,kc)=q1g(ic,jc,kc)
          q2(ic,jc,kc)=q2g(ic,jc,kc)
          q3(ic,jc,kc)=q3g(ic,jc,kc)
          pr(ic,jc,kc)=prg(ic,jc,kc)

          dph(ic,jc,kc)=dphg(ic,jc,kc)                              
        enddo 
       enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      enddo

!@cuf istat = cudaDeviceSynchronize() !JDR TMP


#ifdef USE_CUDA
      !$cuf kernel do (3)
      do kc=kstart,kend 
        do jc=1,n2
          do ic=1,n1
#else
      do kc=kstart,kend
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,ic)
        do jc=1,n2
          do ic=1,n1
#endif
          ru1(ic,jc,kc)=ru1g(ic,jc,kc)
          ru2(ic,jc,kc)=ru2g(ic,jc,kc)
          ru3(ic,jc,kc)=ru3g(ic,jc,kc)          

          qcap(ic,jc,kc)=qcapg(ic,jc,kc)
          dq(ic,jc,kc)=dqg(ic,jc,kc)          
       enddo
       enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      enddo

!@cuf istat = cudaDeviceSynchronize() !JDR TMP
      

#ifdef USE_CUDA
      !$cuf kernel do (2)
        do jc=1,n2
          do ic=1,n1
#else
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,ic)
        do jc=1,n2
          do ic=1,n1
#endif
          qb1s(ic,jc)=qb1sg(ic,jc)
          qb2s(ic,jc)=qb2sg(ic,jc)
          qb3s(ic,jc)=qb3sg(ic,jc)
          qb1n(ic,jc)=qb1ng(ic,jc)
          qb2n(ic,jc)=qb2ng(ic,jc)
          qb3n(ic,jc)=qb3ng(ic,jc)
          dq1x2o(ic,jc)=dq1x2og(ic,jc)
          dq2x2o(ic,jc)=dq2x2og(ic,jc)
          dq3x2o(ic,jc)=dq3x2og(ic,jc)
       enddo
       enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif

!@cuf istat = cudaDeviceSynchronize() !JDR TMP


!mi serve per il modello di contatto mi sembra con xyz
!       Computing normals and triangle properties
        call nvtxStartRange("convert_geo", 32)
        call convert_geo(1, nvtot, 1, netot, 1, nftot, xyz, xyzv, xyza, &
                vert_of_face, tri_ver, tri_vel, tri_bar, vel_tri, acc_tri)
        call nvtxEndRange

        do inp = 1,Nparticle

          call calculate_area(Surface(inp), vstart(inp), vend(inp), &
                   fstart(inp), fend(inp), xyz(:,vstart(inp):vend(inp)), &
                   vert_of_face(:,fstart(inp):fend(inp)), &
                   sur(fstart(inp):fend(inp)))

        end do
!       Find indices of all centroids of particles and bounding box
        call nvtxStartRange("findindices", 31)
        call findindicesFaces
        call findindicesNodes
        call findindicesTiled
#ifdef OFFSETBODY        
        call indindicesTiledOff
#endif
        call findindicesFacesProbeP
        call findindicesFacesProbeN
        call nvtxEndRange



        call mlsStrucCorrectorNewton

        call nvtxEndRange


!       Computing normals and triangle properties
        call nvtxStartRange("convert_geo", 32)
        call convert_geo(1, nvtot, 1, netot, 1, nftot, xyzkp1, xyzvkp1, xyzauk, &
                vert_of_face, tri_ver, tri_vel, tri_bar, vel_tri, acc_tri)
        call nvtxEndRange

        do inp = 1,Nparticle

          call calculate_area(Surface(inp), vstart(inp), vend(inp), &
                   fstart(inp), fend(inp), xyzkp1(:,vstart(inp):vend(inp)), &
                   vert_of_face(:,fstart(inp):fend(inp)), &
                   sur(fstart(inp):fend(inp)))

        end do
!       Find indices of all centroids of particles and bounding box
        call nvtxStartRange("findindices", 31)
        call findindicesFaces
        call findindicesNodes
        call findindicesTiled        
#ifdef OFFSETBODY        
        call indindicesTiledOff
#endif
        call findindicesFacesProbeP
        call findindicesFacesProbeN        
        call nvtxEndRange


!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!      REAL timef !etime_,t(2)
!
!   TIME INTEGRATION : implicit viscous, 3rd order RK (Adams Bashfort)  
!                                                                       

      DO ns=1,nsst

        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)

        if(ntime.le.1)then
            aldto = alm(1)*dt
        else
            aldto = al*dt
        end if

        call nvtxStartRange("boucdq", 30)
        call boucdq
        call nvtxEndRange

        !print*, "q1:", sum(q1(:,:,kstart:kend)), minval(q1(:,:,kstart:kend)), maxval(q1(:,:,kstart:kend))
        !print*, "q2:", sum(q2(:,:,kstart:kend)), minval(q2(:,:,kstart:kend)), maxval(q2(:,:,kstart:kend))
        !print*, "q3:", sum(q3(:,:,kstart:kend)), minval(q3(:,:,kstart:kend)), maxval(q3(:,:,kstart:kend))
!=========================================
!  explicit terms
        call nvtxStartRange("hdnlq123", 9)
        call hdnlq1


        call hdnlq2


        call hdnlq3
        call nvtxEndRange

!@cuf   istat = cudaDeviceSynchronize() !JDR TMP

#ifdef CALTEMP
        call hdnlte
#endif

#ifdef CALSCAL
        call hdnlsa
#endif


!====================================
!  implicit terms
        call nvtxStartRange("invtrq123", 12)
        call invtrq1


        call invtrq2                     !! DQ2HAT=Q2HAT-Q2(N)


        call invtrq3     !! DQ3HAT=Q3HAT-Q3(N)
        call nvtxEndRange
!        call update_upper_ghost(n1,n2,q3)


#ifdef CALTEMP
        call invtrte                     !! DENSITY  ==>>>
        call update_both_ghosts(n1,n2,dens,kstart,kend)
#endif


#ifdef CALSCAL
        call invtrsa                  !! SALINITY  ==>>>
        call update_both_ghosts(n1r,n2r,dsal,kstartr,kendr)
#endif


      call nvtxStartRange("update_both_ghosts", 20)
#ifdef USE_CUDA
      call update_both_ghosts_3(n1,n2,q1,q2,q3,kstart,kend)
#else
      call update_both_ghosts(n1,n2,q1,kstart,kend)
      call update_both_ghosts(n1,n2,q2,kstart,kend)
      call update_both_ghosts(n1,n2,q3,kstart,kend)
#endif
      call nvtxEndRange
!     ============================================
!     Forcing for fluid-structure interaction
!     ============================================

      do jjjdiv=1,nsstepdiv
       if(imlsfor.eq.1) then

!     simple mlsforce without tiling
        !call mlsForce

!     call this instead of mlsforce for tiled forcing
      do jjj=1,nsstep
        call nvtxStartRange("zero forces", 20)
        for_xc = 0.0 ; for_yc = 0.0 ; for_zc = 0.0
        call nvtxEndRange


#ifdef OFFSETBODY
!     ========OffBody Face-Centered-Forcing=========
        call mlsForceTiledOff3Comp
#endif
!     ========Face-Centered-Forcing=========
        call nvtxStartRange("mlsForceTiled3Comp", 16)
        call mlsForceTiled3Comp
        call nvtxEndRange

        call nvtxStartRange("velforce3Comp", 17)
        call velforce3Comp
        call nvtxEndRange

      end do

        !for_sc = 0.0
        ! call mlsForce_scalar
        ! call velforce_scalar

        call nvtxStartRange("sync_collision", 18)
        call sync_collision(n1,n2,coll)
        call nvtxEndRange
       end if
!   ======================================================
!   End of forcing for Fluid-Structure interaction
!   ======================================================


!====================================
! pressure and velocity correction

        call update_upper_ghost(n1,n2,q3)

        call nvtxStartRange("divg", 18)
        call divg        !here only q3 ghost(up) cell are needed
        call nvtxEndRange

        
        call nvtxStartRange("phcalc", 19)
        call phcalc
        call nvtxEndRange

        call update_both_ghosts(n1m,n2m+2,dph,kstart,kend)
        
        call nvtxStartRange("updvp", 21)
        call updvp                 !! SOLENOIDAL VEL FIELD
        call nvtxEndRange


        call nvtxStartRange("velbc", 22)
        call velbc                 !! velocity boundary condition 
        call nvtxEndRange

#ifdef USE_CUDA
        call update_both_ghosts_3(n1,n2,q1,q2,q3,kstart,kend)
#else
        call update_both_ghosts(n1,n2,q1,kstart,kend)
        call update_both_ghosts(n1,n2,q2,kstart,kend)
        call update_both_ghosts(n1,n2,q3,kstart,kend)
#endif
        
        call nvtxStartRange("prcalc", 23)
        call prcalc                         !! PRESSURE FIELD
        call nvtxEndRange

        call update_both_ghosts(n1,n2,pr,kstart,kend)

  

!==================================================
!     multi resolution information exchange
!   expanding vel field from q* to q*lr
        !call mgrd_velitp

#ifdef CALSCAL
!      salinity from refined mesh to base mesh
        call mgrd_dsalc
#endif

      if (infig .eq. 1) then
        call nvtxStartRange("FixedBodies", 25)
        call FixedBodiesPRib
        call nvtxEndRange
        call update_both_ghosts(n1,n2,pr,kstart,kend)
      endif
   enddo !jjjdiv

      call nvtxStartRange("boucqt", 24)
      call boucqt
      call nvtxEndRange
   
      ENDDO         ! end of RK loop

!================================

        ! NOTE: calculating normalized normal here for fpxyz in mlsStruc
        call nvtxStartRange("calculate_normal", 33)
        call calculate_normal(1, nvtot, 1, nftot, xyzkp1, &
                vert_of_face, tri_nor)
        call nvtxEndRange

        call mlsStruc3CompStrong


!     Compute internal force
        call internalforce(1,nvtot,xyzkp1,xyzvkp1)
        call internalforce_3d(1,nvtot_3d,xyzkp1_3d,xyzvkp1_3d)

        ! Reduce the forces from all processors over each particle
        if (numtasks > 1) then
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#ifndef NCCLAVAIL
          call MPI_ALLREDUCE(MPI_IN_PLACE,fxyz,3*nvtot,MPI_DOUBLE_PRECISION, &
                  MPI_SUM,MPI_COMM_WORLD,ierr)
#else
          nccl_result = ncclAllReduce(fxyz,fxyz,3*nvtot,ncclDouble,ncclSum,nccl_comm,0)
#endif

          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#ifndef NCCLAVAIL
          call MPI_ALLREDUCE(MPI_IN_PLACE,fxyz_3d,3*nvtot_3d,MPI_DOUBLE_PRECISION, &
                  MPI_SUM,MPI_COMM_WORLD,ierr)
#else
          nccl_result = ncclAllReduce(fxyz_3d,fxyz_3d,3*nvtot_3d,ncclDouble,ncclSum,nccl_comm,0)
#endif
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        endif


!store stuff
#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = 1, nvtot_3d
         dinvm = 1.d0/mass_of_vert_3d(ntr)
         xyzakp1_3d(1,ntr)=(fpxyz_3d(1,ntr)+fxyz_3d(1,ntr))*dinvm
         xyzakp1_3d(2,ntr)=(fpxyz_3d(2,ntr)+fxyz_3d(2,ntr))*dinvm
         xyzakp1_3d(3,ntr)=(fpxyz_3d(3,ntr)+fxyz_3d(3,ntr))*dinvm
      enddo

#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = nvstart_2dwet,nvend_2dwet
         v1 = tag_2dwet(ntr)
         xyzakp1(1,ntr)=xyzakp1_3d(1,v1)
         xyzakp1(2,ntr)=xyzakp1_3d(2,v1)
         xyzakp1(3,ntr)=xyzakp1_3d(3,v1)
      enddo

#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = nvstart_2dstr,nvend_2dstr
         dinvm = 1.d0/mass_of_vert(ntr)
         xyzakp1(1,ntr)=(fpxyz(1,ntr)+fxyz(1,ntr))*dinvm
         xyzakp1(2,ntr)=(fpxyz(2,ntr)+fxyz(2,ntr))*dinvm
         xyzakp1(3,ntr)=(fpxyz(3,ntr)+fxyz(3,ntr))*dinvm
      enddo
!end store stuff


           !--------- Check errors --------------------------------------------                                                                                                             
            max_error=-1.0d0
#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = 1, nvtot
         max_error=max(max_error, abs(xyzkp1(1,ntr)-xyzk(1,ntr)))
         max_error=max(max_error, abs(xyzkp1(2,ntr)-xyzk(2,ntr)))
         max_error=max(max_error, abs(xyzkp1(3,ntr)-xyzk(3,ntr)))
      enddo


!      write(*,*) 'FSI iter and err: ',iter,max_error
      if (max_error.LE.err_max) then
!store stuff
#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = 1, nvtot_3d
         xyza0_3d(1,ntr)=xyza_3d(1,ntr)
         xyza0_3d(2,ntr)=xyza_3d(2,ntr)
         xyza0_3d(3,ntr)=xyza_3d(3,ntr)

         xyzv0_3d(1,ntr)=xyzv_3d(1,ntr)
         xyzv0_3d(2,ntr)=xyzv_3d(2,ntr)
         xyzv0_3d(3,ntr)=xyzv_3d(3,ntr)

         xyza_3d(1,ntr)=xyzakp1_3d(1,ntr)
         xyza_3d(2,ntr)=xyzakp1_3d(2,ntr)
         xyza_3d(3,ntr)=xyzakp1_3d(3,ntr)

         xyzv_3d(1,ntr)=xyzvkp1_3d(1,ntr)
         xyzv_3d(2,ntr)=xyzvkp1_3d(2,ntr)
         xyzv_3d(3,ntr)=xyzvkp1_3d(3,ntr)

         xyz_3d(1,ntr)=xyzkp1_3d(1,ntr)
         xyz_3d(2,ntr)=xyzkp1_3d(2,ntr)
         xyz_3d(3,ntr)=xyzkp1_3d(3,ntr)
      enddo

#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = nvstart_2dwet,nvend_2dwet
         v1 = tag_2dwet(ntr)
         xyza0(1,ntr)=xyza0_3d(1,v1)
         xyza0(2,ntr)=xyza0_3d(2,v1)
         xyza0(3,ntr)=xyza0_3d(3,v1)

         xyzv0(1,ntr)=xyzv0_3d(1,v1)
         xyzv0(2,ntr)=xyzv0_3d(2,v1)
         xyzv0(3,ntr)=xyzv0_3d(3,v1)

         xyza(1,ntr)=xyza_3d(1,v1)
         xyza(2,ntr)=xyza_3d(2,v1)
         xyza(3,ntr)=xyza_3d(3,v1)

         xyzv(1,ntr)=xyzv_3d(1,v1)
         xyzv(2,ntr)=xyzv_3d(2,v1)
         xyzv(3,ntr)=xyzv_3d(3,v1)

         xyz(1,ntr)=xyz_3d(1,v1)
         xyz(2,ntr)=xyz_3d(2,v1)
         xyz(3,ntr)=xyz_3d(3,v1)
      enddo

#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = nvstart_2dstr,nvend_2dstr
         xyza0(1,ntr)=xyza(1,ntr)
         xyza0(2,ntr)=xyza(2,ntr)
         xyza0(3,ntr)=xyza(3,ntr)

         xyzv0(1,ntr)=xyzv(1,ntr)
         xyzv0(2,ntr)=xyzv(2,ntr)
         xyzv0(3,ntr)=xyzv(3,ntr)

         xyza(1,ntr)=xyzakp1(1,ntr)
         xyza(2,ntr)=xyzakp1(2,ntr)
         xyza(3,ntr)=xyzakp1(3,ntr)

         xyzv(1,ntr)=xyzvkp1(1,ntr)
         xyzv(2,ntr)=xyzvkp1(2,ntr)
         xyzv(3,ntr)=xyzvkp1(3,ntr)

         xyz(1,ntr)=xyzkp1(1,ntr)
         xyz(2,ntr)=xyzkp1(2,ntr)
         xyz(3,ntr)=xyzkp1(3,ntr)
      enddo
!end store stuff
         return !break iterations

      else
!store stuff
#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = 1, nvtot_3d
         xyzakm1_3d(1,ntr)=xyzak_3d(1,ntr)
         xyzakm1_3d(2,ntr)=xyzak_3d(2,ntr)
         xyzakm1_3d(3,ntr)=xyzak_3d(3,ntr)

         xyzvkm1_3d(1,ntr)=xyzvk_3d(1,ntr)
         xyzvkm1_3d(2,ntr)=xyzvk_3d(2,ntr)
         xyzvkm1_3d(3,ntr)=xyzvk_3d(3,ntr)

         xyzak_3d(1,ntr)=xyzakp1_3d(1,ntr)
         xyzak_3d(2,ntr)=xyzakp1_3d(2,ntr)
         xyzak_3d(3,ntr)=xyzakp1_3d(3,ntr)

         xyzvk_3d(1,ntr)=xyzvkp1_3d(1,ntr)
         xyzvk_3d(2,ntr)=xyzvkp1_3d(2,ntr)
         xyzvk_3d(3,ntr)=xyzvkp1_3d(3,ntr)

         xyzk_3d(1,ntr)=xyzkp1_3d(1,ntr)
         xyzk_3d(2,ntr)=xyzkp1_3d(2,ntr)
         xyzk_3d(3,ntr)=xyzkp1_3d(3,ntr)
      enddo

#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = nvstart_2dwet,nvend_2dwet
         v1 = tag_2dwet(ntr)
         xyzakm1(1,ntr)=xyzakm1_3d(1,v1)
         xyzakm1(2,ntr)=xyzakm1_3d(2,v1)
         xyzakm1(3,ntr)=xyzakm1_3d(3,v1)

         xyzvkm1(1,ntr)=xyzvkm1_3d(1,v1)
         xyzvkm1(2,ntr)=xyzvkm1_3d(2,v1)
         xyzvkm1(3,ntr)=xyzvkm1_3d(3,v1)

         xyzak(1,ntr)=xyzak_3d(1,v1)
         xyzak(2,ntr)=xyzak_3d(2,v1)
         xyzak(3,ntr)=xyzak_3d(3,v1)

         xyzvk(1,ntr)=xyzvk_3d(1,v1)
         xyzvk(2,ntr)=xyzvk_3d(2,v1)
         xyzvk(3,ntr)=xyzvk_3d(3,v1)

         xyzk(1,ntr)=xyzk_3d(1,v1)
         xyzk(2,ntr)=xyzk_3d(2,v1)
         xyzk(3,ntr)=xyzk_3d(3,v1)
      enddo


#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = nvstart_2dstr,nvend_2dstr
         xyzakm1(1,ntr)=xyzak(1,ntr)
         xyzakm1(2,ntr)=xyzak(2,ntr)
         xyzakm1(3,ntr)=xyzak(3,ntr)

         xyzvkm1(1,ntr)=xyzvk(1,ntr)
         xyzvkm1(2,ntr)=xyzvk(2,ntr)
         xyzvkm1(3,ntr)=xyzvk(3,ntr)

         xyzak(1,ntr)=xyzakp1(1,ntr)
         xyzak(2,ntr)=xyzakp1(2,ntr)
         xyzak(3,ntr)=xyzakp1(3,ntr)

         xyzvk(1,ntr)=xyzvkp1(1,ntr)
         xyzvk(2,ntr)=xyzvkp1(2,ntr)
         xyzvk(3,ntr)=xyzvkp1(3,ntr)

         xyzk(1,ntr)=xyzkp1(1,ntr)
         xyzk(2,ntr)=xyzkp1(2,ntr)
         xyzk(3,ntr)=xyzkp1(3,ntr)
      enddo

      endif


      enddo !enddo iter corrector





!do this in any case
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do ntr = 1, nvtot_3d
         xyza0_3d(1,ntr)=xyza_3d(1,ntr)
         xyza0_3d(2,ntr)=xyza_3d(2,ntr)
         xyza0_3d(3,ntr)=xyza_3d(3,ntr)

         xyzv0_3d(1,ntr)=xyzv_3d(1,ntr)
         xyzv0_3d(2,ntr)=xyzv_3d(2,ntr)
         xyzv0_3d(3,ntr)=xyzv_3d(3,ntr)

         xyza_3d(1,ntr)=xyzakp1_3d(1,ntr)
         xyza_3d(2,ntr)=xyzakp1_3d(2,ntr)
         xyza_3d(3,ntr)=xyzakp1_3d(3,ntr)

         xyzv_3d(1,ntr)=xyzvkp1_3d(1,ntr)
         xyzv_3d(2,ntr)=xyzvkp1_3d(2,ntr)
         xyzv_3d(3,ntr)=xyzvkp1_3d(3,ntr)

         xyz_3d(1,ntr)=xyzkp1_3d(1,ntr)
         xyz_3d(2,ntr)=xyzkp1_3d(2,ntr)
         xyz_3d(3,ntr)=xyzkp1_3d(3,ntr)
      enddo

#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = nvstart_2dwet,nvend_2dwet
         v1 = tag_2dwet(ntr)
         xyza0(1,ntr)=xyza0_3d(1,v1)
         xyza0(2,ntr)=xyza0_3d(2,v1)
         xyza0(3,ntr)=xyza0_3d(3,v1)

         xyzv0(1,ntr)=xyzv0_3d(1,v1)
         xyzv0(2,ntr)=xyzv0_3d(2,v1)
         xyzv0(3,ntr)=xyzv0_3d(3,v1)

         xyza(1,ntr)=xyza_3d(1,v1)
         xyza(2,ntr)=xyza_3d(2,v1)
         xyza(3,ntr)=xyza_3d(3,v1)

         xyzv(1,ntr)=xyzv_3d(1,v1)
         xyzv(2,ntr)=xyzv_3d(2,v1)
         xyzv(3,ntr)=xyzv_3d(3,v1)

         xyz(1,ntr)=xyz_3d(1,v1)
         xyz(2,ntr)=xyz_3d(2,v1)
         xyz(3,ntr)=xyz_3d(3,v1)
      enddo


#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = nvstart_2dstr,nvend_2dstr
         xyza0(1,ntr)=xyza(1,ntr)
         xyza0(2,ntr)=xyza(2,ntr)
         xyza0(3,ntr)=xyza(3,ntr)

         xyzv0(1,ntr)=xyzv(1,ntr)
         xyzv0(2,ntr)=xyzv(2,ntr)
         xyzv0(3,ntr)=xyzv(3,ntr)

         xyza(1,ntr)=xyzakp1(1,ntr)
         xyza(2,ntr)=xyzakp1(2,ntr)
         xyza(3,ntr)=xyzakp1(3,ntr)

         xyzv(1,ntr)=xyzvkp1(1,ntr)
         xyzv(2,ntr)=xyzvkp1(2,ntr)
         xyzv(3,ntr)=xyzvkp1(3,ntr)

         xyz(1,ntr)=xyzkp1(1,ntr)
         xyz(2,ntr)=xyzkp1(2,ntr)
         xyz(3,ntr)=xyzkp1(3,ntr)
      enddo
!end store stuff


      end if ! if(imlsstr.eq.1.and.ntime.ge.3)


!     END OF  STRUCTURAL SOLVER
!     ==================================================================       


      return                                                            
      end                                                               

