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
      subroutine tschem
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
      integer :: j,k,i,v1
      integer :: vsi, vei, esi, eei, fsi, fei
      integer :: jjj, jjjdiv
      character*50 filename,strucfilename
#ifdef DEBUG
      integer :: kstartp
      real(DP) :: cksum2,cksum3,mck2,mck3,cksum1,mck1
#endif
      real(DP)    :: ti(2)
      real(DP)    :: xpr,ypr,zpr !for testing

      call preambolotsch
      dts = dt/dble(lmax)
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef FLUID_STRUCTURE_SOLVER 
      if(imlsfor.eq.1)then
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
        ! call findindices
        call findindicesFaces
        call findindicesNodes
        call findindicesTiled
#ifdef OFFSETBODY        
        call findindicesTiledOff
#endif
        call findindicesFacesProbeP
        call findindicesFacesProbeN
        call nvtxEndRange
      end if
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
!=========================================
!  explicit terms

#ifndef EXPLICIT
        call nvtxStartRange("hdnlq123", 9)
        call hdnlq1


        call hdnlq2


        call hdnlq3
        call nvtxEndRange
#endif
!@cuf   istat = cudaDeviceSynchronize() !JDR TMP

#ifdef CALTEMP
        call hdnlte
#endif
#ifdef CALSCAL
        call hdnlsa
#endif

!====================================
!  implicit terms
#ifndef EXPLICIT
        call nvtxStartRange("invtrq123", 12)
        call invtrq1


        call invtrq2                     !! DQ2HAT=Q2HAT-Q2(N)


        call invtrq3     !! DQ3HAT=Q3HAT-Q3(N)
        call nvtxEndRange
!        call update_upper_ghost(n1,n2,q3)
#else
#ifdef SCALINGTEST
        ti(1) = MPI_WTIME()    ! CPU time        
#endif
        call explicitNS
#ifdef SCALINGTEST
        ti(2) = MPI_WTIME()    ! CPU time        
        timesforscaling(1)=(ti(2)-ti(1))
#endif
#endif


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
#ifdef SCALINGTEST
      ti(1) = MPI_WTIME()    ! CPU time        
#endif
      call update_both_ghosts_3(n1,n2,q1,q2,q3,kstart,kend)
#ifdef SCALINGTEST      
      ti(2) = MPI_WTIME()    ! CPU time        
      timesforscaling(2)=(ti(2)-ti(1))
#endif
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
#ifdef SCALINGTEST
       ti(1) = MPI_WTIME()    ! CPU time        
#endif
       call mlsForceTiledOff3Comp
#ifdef SCALINGTEST
       ti(2) = MPI_WTIME()    ! CPU time 
       timesforscaling(3)=(ti(2)-ti(1))       
#endif
#endif
!     ========Face-Centered-Forcing=========
        call nvtxStartRange("mlsForceTiled3Comp", 16)
#ifdef SCALINGTEST
        ti(1) = MPI_WTIME()    ! CPU time        
#endif
        call mlsForceTiled3Comp
#ifdef SCALINGTEST
        ti(2) = MPI_WTIME()    ! CPU time     
        timesforscaling(4)=(ti(2)-ti(1))          
#endif
        call nvtxEndRange
        call nvtxStartRange("velforce3Comp", 17)
#ifdef SCALINGTEST
        ti(1) = MPI_WTIME()    ! CPU time        
#endif
        call velforce3Comp
#ifdef SCALINGTEST
        ti(2) = MPI_WTIME()    ! CPU time        
        timesforscaling(5)=(ti(2)-ti(1))          
#endif
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
#ifdef SCALINGTEST
        ti(1) = MPI_WTIME()    ! CPU time        
#endif
        call divg        !here only q3 ghost(up) cell are needed
#ifdef SCALINGTEST
        ti(2) = MPI_WTIME()    ! CPU time        
        timesforscaling(6)=(ti(2)-ti(1))          
#endif
        call nvtxEndRange

        call nvtxStartRange("phcalc", 19)
#ifdef SCALINGTEST
        ti(1) = MPI_WTIME()    ! CPU time        
#endif
        call phcalc
#ifdef SCALINGTEST
        ti(2) = MPI_WTIME()    ! CPU time        
        timesforscaling(7)=(ti(2)-ti(1))          
#endif
        call nvtxEndRange
        call update_both_ghosts(n1m,n2m+2,dph,kstart,kend)
        
        call nvtxStartRange("updvp", 21)
#ifdef SCALINGTEST
        ti(1) = MPI_WTIME()    ! CPU time        
#endif
        call updvp                 !! SOLENOIDAL VEL FIELD
#ifdef SCALINGTEST
        ti(2) = MPI_WTIME()    ! CPU time        
        timesforscaling(8)=(ti(2)-ti(1))          
#endif
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
#ifdef SCALINGTEST
        ti(1) = MPI_WTIME()    ! CPU time        
#endif
        call prcalc                         !! PRESSURE FIELD
#ifdef SCALINGTEST
        ti(2) = MPI_WTIME()    ! CPU time        
        timesforscaling(9)=(ti(2)-ti(1))          
#endif
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
      !     CALLS FOR STRUCTURAL SOLVER
      
      if ((imlsstr.eq.1.and.ntime.ge.3).OR.(Ppro.GT.0))then
         ! NOTE: calculating normalized normal here for fpxyz in mlsStruc
         
        call nvtxStartRange("calculate_normal", 33)
        call calculate_normal(1, nvtot, 1, nftot, xyz, &
                vert_of_face, tri_nor)
        call nvtxEndRange
!     Compute pressure and viscous loads
        call nvtxStartRange("mlsStruc", 29)
#ifdef SCALINGTEST
        ti(1) = MPI_WTIME()    ! CPU time        
#endif
        call mlsStruc3Comp
#ifdef SCALINGTEST
        ti(2) = MPI_WTIME()    ! CPU time        
        timesforscaling(10)=(ti(2)-ti(1))          
#endif
        call nvtxEndRange

!     Compute internal force
#ifdef SCALINGTEST
        ti(1) = MPI_WTIME()    ! CPU time        
#endif

#ifndef SOLOUNO
        call internalforce(1,nvtot,xyz,xyzv)
#endif

        call internalforce_3d(1,nvtot_3d,xyz_3d,xyzv_3d)
#ifdef SCALINGTEST      
        ti(2) = MPI_WTIME()    ! CPU time        
        timesforscaling(11)=(ti(2)-ti(1))          
#endif
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

      end if

!     END OF  STRUCTURAL SOLVER
!     ==================================================================       

!================================
      !     CALLS FOR MLS PROBES
     !  do k=kstart-lvlhalo,kend+lvlhalo
     !    do j=1,n2
     !       do i=1,n1
     !        q1(i,j,k)=xc(i)
     !        q2(i,j,k)=yc(j)
     !        q3(i,j,k)=zc(k)
     !        pr(i,j,k)=xm(i)+ym(j)+zm(k)
              
     !        ! q1(i,j,k)=cos(30.D0*xc(i))
     !        ! q2(i,j,k)=cos(30.D0*yc(j))
     !        ! q3(i,j,k)=cos(30.D0*zc(k))
     !        ! pr(i,j,k)=cos(30.D0*xm(i))*cos(30.D0*ym(j))*cos(30.D0*zm(k))
     !      enddo
     !    enddo
     ! enddo
     
     !  do i=1,Nmlsprobe !for testing, probes distributed along the main diagonal
     !     xyzMLSprobe(1,i)= -xyz_tr(1) + rext1/real(Nmlsprobe)*real(i-1) + 0.5D0*rext1/real(Nmlsprobe)
     !     xyzMLSprobe(2,i)= -xyz_tr(2) + rext2/real(Nmlsprobe)*real(i-1) + 0.5D0*rext2/real(Nmlsprobe)
     !     xyzMLSprobe(3,i)= -xyz_tr(3) + alx3/real(Nmlsprobe)*real(i-1) + 0.5D0*alx3/real(Nmlsprobe)
     !  enddo
     !  call mlsInterp3Comp
      
     !  do i=1,Nmlsprobe
     !     write(*,*) i,xyzMLSprobe(1:3,i)
     !     !         write(*,*) outMLSprobe(1:4,i)
     !     xpr = xyzMLSprobe(1,i)
     !     ypr = xyzMLSprobe(2,i)
     !     zpr = xyzMLSprobe(3,i)
     !     write(*,*) outMLSprobe(1,i) , xpr, abs(outMLSprobe(1,i)  -xpr )    
     !     write(*,*) outMLSprobe(2,i) , ypr,  abs(outMLSprobe(2,i) -ypr )
     !     write(*,*) outMLSprobe(3,i) , zpr, abs(outMLSprobe(3,i)  -zpr )
     !     write(*,*) outMLSprobe(4,i) , xpr+ypr+zpr, abs(outMLSprobe(4,i) - (xpr+ypr+zpr) )

     !     ! write(*,*) outMLSprobe(1,i) , cos(30.D0*xpr), abs(outMLSprobe(1,i) -cos(30.D0*xpr) )
     !     ! write(*,*) outMLSprobe(2,i) , cos(30.D0*ypr), abs(outMLSprobe(2,i) -cos(30.D0*ypr) )
     !     ! write(*,*) outMLSprobe(3,i) , cos(30.D0*zpr), abs(outMLSprobe(3,i) -cos(30.D0*zpr) )
     !     ! write(*,*) outMLSprobe(4,i) , cos(30.D0*xpr)*cos(30.D0*ypr)*cos(30.D0*zpr), abs(outMLSprobe(4,i) -cos(30.D0*xpr)*cos(30.D0*ypr)*cos(30.D0*zpr) )          
     !  enddo      
     !  stop
      
!     END OF  STRUCTURAL SOLVER
!     ==================================================================       
#endif
      
      

      return                                                            
      end                                                               

