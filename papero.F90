!==========================================================================
!     Cartesian version - MLS - IBM deforming objects (MPI)
!     Vamsi Spandan, Phys. Fluids, Univ.Twente
!     Marco D de Tullio, Politecnico di Bari
!     Single phase code: E.P.van der Poel, R.O.Monico
!     Multiple resolution code: Yantao Yang
!     Cardiac Simulation : Francesco Viola 
!     GPU implementation : Josh Romero, Massimiliano Fatica (NVIDIA)
!     Group Lead: Roberto Verzicco
!==========================================================================
!     This code is made for simulating three-dimensional flows in
!     cartesian coordinates.
!
!      navier-stokes equations are solved by a fractional step method
!     ( Kim and Moin ) with the pressure in the first step.
!
!     The time advancement of the solution is obtained by a
!     Runge-Kutta 3rd order low storage scheme (Wray) or a 2nd
!     order Adams-Bashfort scheme.
!
!     The Poisson  equation for the pressure is solved directly
!     introducing FFT in the azimutal and vertical direction.
!     Because of the presence of the walls in the vertical direction
!     the cosFFT in that direction is required
!
!                 Roberto Verzicco and Paolo Orlandi
!                 dipartimento di meccanica ed aeronautica
!                 universita' la sapienza di roma
!
!     All variables are calculated in a staggered grid:
!
!        q2= vr, q3= vz
!
!        dq2,dq3 : velocity correction
!
!        qcap :divergence of the  non free divergent velocity field
!
!        non linear terms:
!
!        ru2, ru3, ruro : old step
!        h2,h3, hro : new step
!        dens : density
!        pr : pressure
!        dph : pressure correction
!       pressure solver is in the package press.f
!       non linear terms are calculated in hdnl routines
!       the invertions of momentum equations is performed in invtr
!       routines
!
!======================================================================
      program papero
      use mpih
      use mpi_param
      use param
      use mls_param
#ifdef USE_CUDA
      use cudafor
#endif
      use probes
      implicit none
      character(len=4)   :: dummy
      integer :: tfini,tfin,n,ns
      real(DP) :: ts,te,npproR,qtrest,qtframe
      integer :: omp_get_num_threads
      integer :: inp, i , j,nppro
#ifdef USE_CUDA
      integer :: istat, local_comm, mydev
#endif


      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)

#ifdef USE_CUDA
#ifdef NCCLAVAIL
      ! Set devices to local rank
      call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, &
                               myid, MPI_INFO_NULL, local_comm, ierr)
      call MPI_Comm_rank(local_comm, mydev, ierr)

      istat = cudaSetDevice(mydev)
      call setup_nccl
#else
      istat = cudaSetDevice(myid) !c'era questo ma probabilmente ci va l'altro
#endif
#endif

      numthreads = 1
!$OMP PARALLEL
      numthreads = omp_get_num_threads()
!$OMP END PARALLEL

      if (myid .eq. 0) then
        write(*,'(3x,a,i8)') 'MPI tasks=', numtasks
        write(*,'(3x,a,i8)') 'OMP threads per task=', numthreads
        write(*,'(3x,a,i8)') 'No. of processors=', numthreads*numtasks
      endif

!     ------------------------------------------
!     Initial allocations for geometries 
      allocate(geofile(Nparticle),iopen_pos(Nparticle),iopen_neg(Nparticle),wet_3d(Nparticle))
      allocate(body_nsl(Nparticle),body_flux(Nparticle))
      allocate(geofile_3d(Nparticle_3d))
      allocate(geofileEF_2d(NparticleEF_2d))
      allocate(geofileEF_1d(NparticleEF_1d))
      allocate(geofileCV_1d(NparticleCV_1d))

      allocate(rke(Nparticle),rkc(Nparticle),rkb(Nparticle),kv(Nparticle))
      allocate(kat(Nparticle),kal(Nparticle),kae(Nparticle))
      allocate(rhos(Nparticle),thck(Nparticle))
      allocate(cv1(Nparticle),cv2(Nparticle))
      ! allocate(rke_3d(Nparticle),rkc_3d(Nparticle),kv_3d(Nparticle))
      allocate(rke_3d(Nparticle),kv_3d(Nparticle))
      allocate(rkc_3d(9)) !nine chambers for hyperelastic
      allocate(kae_3d(Nparticle))
      allocate(rhos_3d(Nparticle))
      allocate(cv1_3d(Nparticle),cv2_3d(Nparticle))
      allocate(rpm(Nparticle))
      allocate(xyz_tr(3))
!     ------------------------------------------

!*******************************************************
!******* Read input file bou.in by all processes********
!*******************************************************
301     format(a4)
      open(unit=15,file='bou.in',status='old')
        read(15,301) dummy
        read(15,*) n1,n2,n3
        read(15,301) dummy
        read(15,*) ireset,nread,nwrit,folderload
        read(15,301) dummy
        read(15,*) Ppro,nreadstart,nreadend,ECG_Ppro
        read(15,301) dummy
        read(15,*) nsst,nsstep,nsstepdiv
        read(15,301) dummy
        read(15,*) ntst,tmax,qtrest,qtframe,tpin,tecg
        read(15,301) dummy
        read(15,*) alx3,rext1,rext2,istr3,str3,lmax
        read(15,301) dummy
        read(15,*) xyz_tr(3),xyz_tr(1),xyz_tr(2)
        read(15,301) dummy
        read(15,*) rey,ras,prs,rat,prt,sbl
        read(15,301) dummy
        read(15,*) idtv,dtmax,dt,cflmax,cfllim,resid
        read(15,301) dummy
        read(15,*) undr, err_max, n_max_iter
        read(15,301) dummy
        read(15,*) ubctop,ubcbot
        read(15,301) dummy
        read(15,*) tbctop, tbcbot, sbctop, sbcbot
        read(15,301) dummy
        read(15,*) infig,bcs,cou,tilemax
        read(15,301) dummy
        read(15,*) imlsfor, imlsstr, imlsref
        read(15,301) dummy
        read(15,*) wcheck, fcheck, wcub, wexp, nel
        read(15,301) dummy
        read(15,*) wcon,wscl
        read(15,301) dummy
        read(15,*) LSTAR,USTAR,densitySTAR,SMZ_S,HR,Nperiod,dyssy
        read(15,301) dummy
        read(15,*) stepEF_3d,stepEF_2d,stepEF_1d,tolEFbido, maxitEFbido, maxrestEFbido
        read(15,301) dummy
        read(15,*) ECG_D,ECG_nson,EscSettoPurk,EscCellConn,EscZonSet
        read(15,301) dummy
        read(15,*) nat_pac,bsn,leads,minf

      close(15)

#ifdef SCALINGTEST
      infig=0
#endif
      
! ============================================
!     Non dimensional quantities for heart simulation
      TSTAR = LSTAR/USTAR !sec
      periodDim = (60.d0/HR) !sec
      period = (60.d0/HR)/TSTAR !nondim
      tmax = min(Nperiod*period,tmax)
      tmax = tmax + 0.01d0*period

      !tpin tecg in bou.in are in ms, qtframe frame and qtrest rest
      trest=period/qtrest
      tframe=period/qtframe
      tpin=tpin/(1000d0*TSTAR)
      tecg=tecg/(1000d0*TSTAR)

      !Added for correction of HR induced by restarting
      HR_S1S2 = HR
      write(*,*) "HR = ", HR, &
                 " periodDim = ",periodDim, &
                 " period = ",period
      if (PPro.GT.0) then
         nread = 1
         ireset = 0
         if(myid.eq.0) call system('mkdir cfield') 
      endif
!     Clean stuff eventually
      if (((nread.EQ.1).AND.(folderload.EQ.'./restart/')).OR.(PPro.GT.0)) then
      else
         if(myid.eq.0) call system('rm *Mask.gts')
         if(myid.eq.0) call system('rm *_3d.out')
         if(myid.eq.0) call system('rm -r restart')
         if(myid.eq.0) call system('rm -r vtkfiles')
         if(myid.eq.0) call system('rm pressure.out')
         if(myid.eq.0) call system('rm electroEF.out')
         if(myid.eq.0) call system('rm ECG*.out')
         if(myid.eq.0) call system('rm EGM*.out')
         if(myid.eq.0) call system('mkdir data fact vtkfiles restart')
      endif

! ============================================

!      open(unit=15,file='ir3files.in',status='old')
      open(unit=15,file='vtkfiles.in',status='old')
      do inp = 1,Nparticle_3d
      read(15,301)dummy
      read(15,*)geofile_3d(inp)
      read(15,301)dummy
      read(15,*) rke_3d(inp), &
                 kv_3d(inp),kae_3d(inp), &
                 rhos_3d(inp), &
                 cv1_3d(inp),cv2_3d(inp)
      end do

      open(unit=15,file='gtsfiles.in',status='old')
      do inp = 1,Nparticle
      read(15,301)dummy
      read(15,*)geofile(inp),iopen_pos(inp),iopen_neg(inp),wet_3d(inp), &
                body_nsl(inp),body_flux(inp)
      read(15,301)dummy
      read(15,*) rke(inp), rkb(inp), &
                 kv(inp),kat(inp),kal(inp),kae(inp), &
                 rhos(inp),thck(inp),cv1(inp),cv2(inp)
      end do
      close(15)

#ifdef ELECTRO
      open(unit=15,file='gtsfilesEF_1d.in',status='old')
      do inp = 1,NparticleEF_1d
         read(15,301)dummy
         read(15,*)geofileEF_1d(inp)
      enddo
      close(15)

      open(unit=15,file='gtsfilesEF_2d.in',status='old')
      do inp = 1,NparticleEF_2d
         read(15,301)dummy
         read(15,*)geofileEF_2d(inp)
      enddo
      close(15)
#endif

!CORONARY VEINS
      open(unit=15,file='gtsfilesCV_1d.in',status='old')
      do inp = 1,NparticleCV_1d
         read(15,301)dummy
         read(15,*)geofileCV_1d(inp)
      enddo
      close(15)

!*******************************************************
!******* Read input file leads.in by all processes********
!*******************************************************
      allocate(ECG_PHI(ECG_nson))
      allocate(ECG_E(3,ECG_nson))
      open(unit=15,file='leads.in',status='old')
          do i=1,ECG_NSON
             read(15,*) ECG_E(1,i),ECG_E(2,i),ECG_E(3,i)
             ECG_E(1,i)=ECG_E(1,i)/(LSTAR*1000.d0)
             ECG_E(2,i)=ECG_E(2,i)/(LSTAR*1000.d0)
             ECG_E(3,i)=ECG_E(3,i)/(LSTAR*1000.d0)
             write(*,*) 'leads pos Xs',i,ECG_E(1,i),ECG_E(2,i),ECG_E(3,i)
          enddo
      close(15)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)


!     probe positions
      open(unit=15,file='sonda.in',status='old')
      read(15,*) nson
      do i=1,nson
       read(15,*) (coson(i,j),j=1,3)
      end do
! ! Apply translation
!       coson(:,1)=coson(:,1)
!       coson(:,2)=coson(:,2)
!       coson(:,3)=coson(:,3)
      close(15)


      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!============================================
!************ End of input file**************
!============================================
      if( n1 .ne. m1 ) then
      if(myid.eq.0) then
        write(*,*) 'Error: n1 must be = m1'
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif
      endif

      if( n2 .ne. m2 ) then
      if(myid.eq.0) then
        write(*,*) 'Error: n2 must be = m2'
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif
      endif

      if( n3 .ne. m3 ) then
      if(myid.eq.0) then
        write(*,*) 'Error: n3 must be = m3'
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif
      endif
!============================================
      rint = 0.d0
!
!     DEFINITIONS FOR THE NATURAL CONVECTION
!
      lew = prs/prt
      rhop = rat*lew/ras

!m ============================================
!    coefficients in equations

!    nondimensionalized by S free-fall vel
!     viscosity
      !nu = dsqrt(prs/ras)
      nu = 1.d0/rey ! hardcoded to be consistent with Francesco's code
!     diffusivity
      kps = 1.d0 / dsqrt(Ras*Prs)
      kpt = Lew / dsqrt(Ras*PrS)
!     buoyancy coefs
      byct = rhop
      bycs = 1.d0

!============================================
!    boundary values for scalars
      denstop = 0.d0
      densbot = 0.d0

      dsaltop = 0.d0
      dsalbot = 0.d0

!================================================
      call openfi
!================================================

      tfini=dt*ntst
      tfin=tfini
      n1m=n1-1
      n2m=n2-1
      n3m=n3-1

!======================
      n1mr=n1m*mref1
      n1r=n1mr+1
      n2mr=n2m*mref2
      n2r=n2mr+1
      n3mr=n3m*mref3
      n3r=n3mr+1

      usref1=1.d0/dble(mref1)
      usref2=1.d0/dble(mref2)
      usref3=1.d0/dble(mref3)

!===================================================
!    refinement of time step

!      lmax = max(dble(mref1),dble(mref2),dble(mref3))

!======================
!
!====================================================
      if(myid.eq.0) then
        write(*,112)rext1/alx3, rext2/alx3
        write(*,202)rat,prt
        write(*,203)ras,prs
        write(*,213)rhop
        if(idtv.eq.1)then
          write(*,204)cflmax
        else
          write(*,205)dtmax,cfllim
        endif
        write(*,206)tmax
      endif
!====================================================
112   format(//,6x,'CARDIOVASCULAR FLOW SIMULATOR',//, &
       3x, '3D Cell with aspect-ratios:  Lx/H=',f6.3, ' Ly/H=', f6.3)
142   format(3x,'Periodic lateral wall boundary condition')
202   format(3x,'Parameters temperature: ', &
                 ' RaT=',es10.3,' PrT= ',es10.3)
203   format(3x,'Parameters salinity:    ' &
                ,' RaS=',es10.3,' PrS= ',es10.3)
213   format(3x,'Density ratio:    Rhop=',es10.3,/)
204   format(3x,'Variable dt and fixed cfl= ',es11.4)
205   format(3x,'Fixed dt= ',es11.4,' and maximum cfl=',es11.4)
206   format(3x,'Total integration time: ',f8.2)
!====================================================
!
!     assign coefficients for time marching schemes
!
      if(nsst.gt.1) then
        gam(1)=8.d0/15.d0
        gam(2)=5.d0/12.d0
        gam(3)=3.d0/4.d0
        rom(1)=0.d0
        rom(2)=-17.d0/60.d0
        rom(3)=-5.d0/12.d0
!======================================================
        if(myid.eq.0) then
          write(*,100) (gam(n),n=1,nsst),(rom(n),n=1,nsst)
        endif
100     format(3x,'The time scheme is a III order Runge-Kutta', &
              /,3x,'gam= ',3f8.3,4x,'ro= ',3f8.3)
!======================================================
      else
        gam(1)=1.5d0
        gam(2)=0.d0
        gam(3)=0.d0
        rom(1)=-0.5d0
        rom(2)=0.d0
        rom(3)=0.d0
!======================================================
        if(myid.eq.0) then
          write(*,110) gam(1),rom(1)
        endif
110     format(3x,'The time scheme is the Adams-Bashfort', &
              /,3x,'gam= ',f8.3,4x,'ro= ',f8.3)

!======================================================
      endif

      do ns=1,nsst
        alm(ns)=(gam(ns)+rom(ns))
      enddo

!======================================================

#ifdef DEBUG
      if(myid.eq.0) write(*,*)myid, 'start work distribution'
#endif
      call mpi_workdistribution

#ifdef USE_CUDA
#ifdef CUDECOMPAVAIL
      call setup_cudecomp
#endif
#endif
      
#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0) write(*,*)myid, 'start mem alloc'
#endif
      call mem_alloc
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!======================================================
!
!     the solution of the problem starts
!
!======================================================

      ts=MPI_WTIME()

#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0) write(*,*)myid, 'start gcurv'
#endif




      if (PPro.LE.0) then
         iread=0
         !running mode
         call gcurv
      else
         !postprocessing mode
!         ntst = 0 !save directly after blanking
         ntst  = 1 !make 1step before saving (it computes fpxyz)
         iread = 1
         nPProR = (nreadend-nreadstart+PPro)/PPro
         nPPro = nint( nPProR )
         call allocate_stuffPpro
         do nwrit=nreadstart,nreadend,PPro            
            write(*,*) "ANALYZING", iread,"/",nPPro
            write(*,*) "NWRITE", nwrit


            call gcurv
            call deallocate_stuff
            call deallocate_mls_local 
            call mgrd_mem_dealloc
            iread = iread + 1
 
        enddo
         call deallocate_stuffPpro

      endif




#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0) write(*,*)myid, ' gcurv exited'
#endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      te=MPI_WTIME()

      if(myid.eq.0) then
        open(27,file="fact/Total_time.out")
        write(27,*)"Total simulation time in sec.: ", te-ts
        close(27)
      endif

!================================================

      call closefi

!==================================================

      call mem_dealloc
      call mgrd_mem_dealloc
#ifdef DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0)write(*,*)'memory released'
#endif

      call MPI_FINALIZE(ierr)
!==================================================

      stop
      end

