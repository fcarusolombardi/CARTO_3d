!===========================================================
! Declaration of global variables
!***********************************************************      
      module constants
        integer, parameter :: DP = kind(1.d0)
        real(DP), parameter :: PI = acos(-1.d0)
      end module

      module param
        use constants
        implicit none
!m=============================================================
!m      maximal wall time in second
        real(DP),parameter :: walltimemax = 8500000.d0
!m===========================================================
!m      grid size
        integer,parameter :: m1=3,  mref1=1
        integer,parameter :: m2=3,  mref2=mref1
        integer,parameter :: m3=3,  mref3=1
!m
        integer :: m2m,m3m,m2mh,m1m
        integer :: m1r,m2r,m3r
        integer :: m1mr,m2mr,m3mr,m2mhr
        parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
        parameter (m1mr=m1m*mref1,m1r=m1mr+1)
        parameter (m2mr=m2m*mref2,m2r=m2mr+1)
        parameter (m3mr=m3m*mref3,m3r=m3mr+1)
        parameter (m2mh=m2m/2+1,m2mhr=m2mr/2+1)
!m==============================================================
!m      inital condition
!m  0:  both uniform S=T=0.5, uni
!m  1:  both linear, lin
!m  2:  S=0.5, T linear, mix
!m  3:  interior interface(s), itf*
!m      total homogeneous layers: Nitf+1
!m      total sharp inital interface: Nitf+2 (two at bounaries)
        integer,parameter :: tag_ini = 1
        integer,parameter :: Nitf = 1
        ! for sharp interface, number of interior interface
        real(DP), dimension(1:m3)  :: Tiniprof
        real(DP), dimension(1:m3r) :: Siniprof
!==========================================================			
!       read from input file bou.in
!==========================================================
        integer   :: n1, n2, n3, nsst, nwrit, nread
        integer   :: ECG_Ppro,Ppro,nreadstart,nreadend,iread
        integer   :: ntst, ireset,nsstep,nsstepdiv
        real(DP)      :: trest, tpin, tmax
        real(DP)      :: alx3, str3, rext1, rext2
        integer   :: istr3, lmax
        real(DP)      :: rey, rat, prt, ras, prs, lew, rhop, sbl
        integer   :: idtv
        real(DP)      :: dtmax, resid, cflmax, dt,dt_o, cfllim, dts
        real(DP)   :: undr, err_max, Lgrid
        integer   :: stepEF_3d,stepEF_2d,stepEF_1d,n_max_iter,maxitEFbido,maxrestEFbido,EscSettoPurk,EscCellConn
        integer   :: ubctop, ubcbot, tbctop, tbcbot, sbctop, sbcbot
        integer   :: nat_pac,bsn,leads,minf
        real(DP)      :: tframe, tecg,cou,tolEFbido,EscZonSet
!==========================================================			
!       read from input file part.in
!==========================================================
        integer   :: infig,bcs,tilemax
        integer   :: imlsfor,imlsstr,imlsref
        integer   :: wcheck,fcheck,wcub,wexp,nel
        real(DP)      :: wcon, wscl
        real(DP):: LSTAR,USTAR,TSTAR,densitySTAR,SMZ,SMZ_S,HR,HR_S1S2,Nperiod,dyssy
        real(DP),dimension(:),allocatable ::xyz_tr
        character*50  gtsfx,datfx,stlfx,stlfxAS,stlfxPpro
        character*100 folderload
!=================================================
!       end of input file
!=================================================
        real(DP)      :: time,timeBeat,periodDim,period
        real(DP)      :: denstop, densbot, dsaltop, dsalbot
        integer :: nBeat,switchEXPLV,switchEXPLA,switchSMZ

!******* Grid parameters**************************
        real(DP) :: dx2,dx3,dx1
        real(DP) :: dx2q,dx3q,dx1q
!
        real(DP) :: dx2r,dx3r,dx1r
        real(DP) :: dx2qr,dx3qr,dx1qr
        
        real(DP), dimension(1:m3) :: g3rc,g3rm
        real(DP), dimension(:), allocatable :: xc, yc, zc
        real(DP), dimension(:), allocatable :: xm, ym, zm
        
        real(DP), dimension(1:m1r) :: xcr,xmr
        real(DP), dimension(1:m2r) :: ycr,ymr
        real(DP), dimension(1:m3r) :: zcr,zmr,g3rcr,g3rmr

!====================================================
!******* QUANTITIES FOR DERIVATIVES******************
        real(DP), dimension(1:m3r) :: udx3cr,udx3mr
        real(DP), dimension(:), allocatable :: udx3c,udx3m

#ifdef USE_CUDA
        attributes(managed) :: xyz_tr
        attributes(managed) :: xc, yc, zc
        attributes(managed) :: xm, ym, zm
        attributes(managed) :: udx3c, udx3m
#endif

!==========================================================
!******* Grid indices**************************************
        integer, dimension(:), allocatable :: jmv,jpv
        integer, dimension(:), allocatable :: imv,ipv
        integer, dimension(:), allocatable :: kmv,kpv
        integer, dimension(:), allocatable :: jmhv
        integer, dimension(1:m3) :: kmc,kpc,kup,kum

        integer, dimension(1:m1r) :: imvr,ipvr
        integer, dimension(1:m2r) :: jmvr,jpvr
        integer, dimension(1:m3r) :: kmvr,kpvr

#ifdef USE_CUDA
        attributes(managed) :: jmv, jpv
        attributes(managed) :: imv, ipv
        attributes(managed) :: kmv, kpv
        attributes(managed) :: jmhv
#endif
  
!===========================================================
!******* Metric coefficients *******************************
        real(DP), dimension(1:m2) :: ap3j,ac3j,am3j
        real(DP), dimension(:), allocatable :: ap3ck,ac3ck,am3ck
        real(DP), dimension(1:m3) :: ap3ssk,ac3ssk,am3ssk   
        real(DP), dimension(1:m3r) :: ap3sskr,ac3sskr,am3sskr
        real(DP), dimension(:), allocatable :: ap3sk,ac3sk,am3sk

#ifdef USE_CUDA
        attributes(managed) :: ap3ck,ac3ck,am3ck
        attributes(managed) :: ap3sk,ac3sk,am3sk
#endif

!============================================================
!******* Variables for FFTW and Poisson solver****************
        real(DP), dimension(13) :: ifx1
        real(DP), dimension(3*m2/2+1) :: trigx1
        real(DP), dimension(1:m2) :: ap
        real(DP), dimension(1:m1) :: ao
        real(DP), dimension(:), allocatable :: ak1, ak2
        real(DP), dimension(:), allocatable :: amphk,acphk,apphk
        integer*8 :: fwd_plan,bck_plan

#ifdef USE_CUDA
        integer :: cufft_fwd_plan,cufft_bck_plan
        attributes(managed) :: ak1, ak2
        attributes(managed) :: amphk,acphk,apphk
#endif
        
!===========================================================
!******* Other variables ***********************************
        integer  :: n2m, n3m, n1m
        integer :: n1r,n2r,n3r,n1mr,n2mr,n3mr
        integer  :: iaxsy
        real(DP) :: rint
        real(DP) :: nu, kps, kpt
        real(DP) :: byct, bycs    ! coefs for buoyancy force term
        real(DP) :: al,ga,ro,aldto
        real(DP) :: beta
        real(DP) :: qqmax,qqtot
        real(DP) :: re
        real(DP) :: anutlow, anutupp
        real(DP) :: anuslow, anusupp
        real(DP) :: denmax, denmin, densm
        real(DP) :: max_error,scarto1gmres,scarto2gmres
        integer :: ntime, iter,countjgmres
        integer, parameter:: ndv=3
        real(DP), dimension(1:ndv) :: vmax
        real(DP), dimension(1:3) :: gam,rom,alm
        real(DP), dimension(1:m1,1:m2) :: densbs,densbn
        real(DP), dimension(1:m1r,1:m2r) :: dsalbs,dsalbn
        real(DP) :: usref1, usref2, usref3              

        real(DP), allocatable, dimension(:) :: aml, acl, apl
        real(DP), allocatable, dimension(:,:,:) :: buf
        real(DP),dimension(:),allocatable ::timesforscaling
#ifdef USE_CUDA
        attributes(managed) :: aml, acl, apl
        attributes(managed) :: buf
#endif

       ! Send/recv buffers for use with MPI
       real(DP), allocatable, dimension(:) :: sendbuff
       real(DP), allocatable, dimension(:) :: recvbuff
#ifdef USE_CUDA
       attributes(device) :: sendbuff, recvbuff
#endif

       real(DP) :: pressOUT(6)

#ifdef USE_CUDA
  integer(int_ptr_kind()):: freeMem,totalMem
  real(DP):: rFreeMem,rUsedMem,rTotalMem,minUsedMem,maxUsedMem,minFreeMem,maxFreeMem
  ! integer:: istat
#endif
      end module param

      
!************* End of param module******************************


!===============================================================
!******* 3D arrays, dynamically allocated by each process*******
      module local_arrays
        use constants
        use param
        implicit none
        real(DP),allocatable,dimension(:,:,:) :: q1,q2,q3,dens,dsal
        real(DP),allocatable,dimension(:,:,:) :: q1g,q2g,q3g,prg
        real(DP),allocatable,dimension(:,:,:) :: hro,hsal,rhs,rhss,rhsss,rhsr,rhs_t
        real(DP),allocatable,dimension(:,:,:) :: ru1,ru2,ru3,ruro,rusal
        real(DP),allocatable,dimension(:,:,:) :: ru1g,ru2g,ru3g        
        real(DP),allocatable,dimension(:,:,:) :: pr,qcap,dph,dq,dpht
        real(DP),allocatable,dimension(:,:,:) :: qcapg,dphg,dqg               
        real(DP),allocatable,dimension(:,:,:) :: forclo_t
        real(DP),allocatable,dimension(:,:,:,:) :: forclo
        real(DP),allocatable,dimension(:),target :: rwork1,rwork2,rwork3
!        complex(DP),allocatable,dimension(:,:,:) :: cwork1_3D
        real(DP),allocatable,dimension(:,:,:) :: rhst,forclot
        real(DP),allocatable,dimension(:,:,:) :: maskq1,maskq2,maskq3,maskpr

#ifdef USE_CUDA
        attributes(managed) :: q1,q2,q3,dens,dsal
        attributes(managed) :: rhs,rhss,rhsss,rhs_t
        attributes(managed) :: ru1,ru2,ru3
        attributes(managed) :: ru1g,ru2g,ru3g        
        attributes(managed) :: pr,qcap,dph,dq,dpht
        attributes(managed) :: qcapg,dphg,dqg        
        attributes(managed) :: q1g,q2g,q3g,prg
        attributes(managed) :: forclo, forclo_t
        attributes(managed) :: rwork1,rwork2,rwork3
!        attributes(managed) :: cwork1_3D
        attributes(managed) :: rhst,forclot

        real(DP),allocatable,dimension(:),device :: rwork1_d,rwork2_d
#endif
      end module local_arrays

!===============================================================
      module stat_arrays
        use constants
        implicit none
        real(DP),allocatable, dimension(:) :: q1_me,q1_rms 
        real(DP),allocatable, dimension(:) :: q2_me,q2_rms
        real(DP),allocatable, dimension(:) :: q3_me,q3_rms 
        real(DP),allocatable, dimension(:) :: dens_me,dens_rms 
        real(DP),allocatable, dimension(:) :: dsal_me,dsal_rms
        real(DP),allocatable, dimension(:) :: q3dens_me,q3dsal_me 
        real(DP),allocatable, dimension(:) :: disste,disssa
        real(DP),allocatable, dimension(:) :: dissuc,dissur
        integer :: timeint_cdsp
      end module stat_arrays

!=====================================================       
      module mpih
#ifdef NCCLAVAIL 
        !@cuf use nccl  
#endif
        implicit none
        include 'mpif.h'
        integer :: myid, numtasks, numthreads, ierr
        integer, parameter :: master=0
        integer, parameter :: lvlhalo=3
        integer :: MDP = MPI_DOUBLE_PRECISION
        integer :: STATUS(MPI_STATUS_SIZE,4)
        integer :: req(1:4)
        integer(kind=MPI_OFFSET_KIND) :: disp, offset
#ifdef NCCLAVAIL 
        type(ncclUniqueId) :: nccl_id
        type(ncclResult) :: nccl_result
        type(ncclComm) :: nccl_comm
#endif
      end module mpih

      module cudecomp_param
#ifdef CUDECOMPAVAIL
        use constants
        use cudafor
        use cudecomp
        implicit none
        type(cudecompHandle) :: cudecomp_handle
        type(cudecompGridDesc) :: cudecomp_grid_desc
        type(cudecompGridDesc) :: cudecomp_grid_desc_ph

        integer :: cudecomp_dtype=CUDECOMP_DOUBLE
        real(DP), pointer, device, contiguous :: cudecomp_work_d(:)
#endif
      end module cudecomp_param
      
      module mpi_param
        implicit none
        integer :: istart,iend, jstart,jend, kstart,kend
        integer :: kstartp,kendp
        integer :: jstartr,jendr, kstartr,kendr
        integer :: jstartp,jendp
        integer :: dj,dk,mydata,mydatam
        integer :: djp,djr,dkr
        integer, allocatable, dimension(:) :: offsetj,offsetk
        integer, allocatable, dimension(:) :: offsetjr,offsetkr
        integer, allocatable, dimension(:) :: offsetjp
        integer, allocatable, dimension(:) :: countj,countk
        integer, allocatable, dimension(:) :: countjr,countkr
        integer, allocatable, dimension(:) :: countjp
        integer, allocatable, dimension(:) :: countf
        integer(8), allocatable, dimension(:) :: offsetf 
      end module mpi_param

!==========================================================
      module slab_param
        use constants
        implicit none
        integer, parameter :: nslab = 15
        integer :: idslab(1:nslab), idslabr(1:nslab)
        integer :: kslab(1:nslab), kslabr(1:nslab)
        real(DP) :: zslab(1:nslab)
      end module slab_param
!
!===============================================================
!******* for multi grid arrays*******
      module mgrd_arrays
        use constants
        use param
        implicit none
        integer, parameter :: mrefa=mref1*mref2*mref3
        integer irangs(0:m1),jrangs(0:m2),krangs(0:m3)
        integer indc1sal(m1),indc2sal(m2),indc3sal(m3)
        real(DP), dimension(4,0:m1r) :: cxq1, cxq2, cxq3, cxrs
        real(DP), dimension(4,0:m2r) :: cyq1, cyq2, cyq3, cyrs
        real(DP), dimension(4,0:m3r) :: czq1, czq2, czq3, czrs
        real(DP),allocatable,dimension(:,:,:) :: q1lr,q2lr,q3lr,dsalc
      end module mgrd_arrays
!
!===============================================================
!******* initial interpoaltion *******
      module init_itp
        use constants
        implicit none
        integer n1o,n2o,n3o,n1om,n2om,n3om
        integer n1or,n2or,n3or,n1omr,n2omr,n3omr
        integer mref1o,mref2o,mref3o
        integer istr3o
        real(DP) str3o

        real(DP),allocatable,dimension(:) :: xcold,xmold
        real(DP),allocatable,dimension(:) :: ycold,ymold
        real(DP),allocatable,dimension(:) :: zcold,zmold
        real(DP),allocatable,dimension(:) :: xcrold,xmrold
        real(DP),allocatable,dimension(:) :: ycrold,ymrold
        real(DP),allocatable,dimension(:) :: zcrold,zmrold

        integer,allocatable,dimension(:) :: iiq1c1,iiq1c2,iiq1c3
        integer,allocatable,dimension(:) :: iiq2c1,iiq2c2,iiq2c3
        integer,allocatable,dimension(:) :: iiq3c1,iiq3c2,iiq3c3
        integer,allocatable,dimension(:) :: iidec1,iidec2,iidec3
        integer,allocatable,dimension(:) :: iisac1,iisac2,iisac3

        real(DP),allocatable,dimension(:,:,:) :: varold 

      end module init_itp 
!====================================================

!====================================================
      module outflow_vars
        use constants
        use param
        implicit none
        real(DP) :: qinf, qout
        real(DP),allocatable,  dimension(:,:) :: qb1s, qb1n
        real(DP),allocatable,  dimension(:,:) :: qb1sg, qb1ng                
        real(DP),allocatable,  dimension(:,:) :: dqb1s, dqb1n
        real(DP),allocatable,  dimension(:,:) :: qb2s, qb2n
        real(DP),allocatable,  dimension(:,:) :: qb2sg, qb2ng        
        real(DP),allocatable,  dimension(:,:) :: dqb2s, dqb2n
        real(DP),allocatable,  dimension(:,:) :: qb3s, qb3n
        real(DP),allocatable,  dimension(:,:) :: qb3sg, qb3ng        
        real(DP),allocatable,  dimension(:,:) :: dqb3s, dqb3n
        real(DP),allocatable,  dimension(:,:) :: dq1x2o, dq2x2o, dq3x2o
        real(DP),allocatable,  dimension(:,:) :: dq1x2og, dq2x2og, dq3x2og        
#ifdef USE_CUDA
        attributes(managed) :: qb1s, qb1n
        attributes(managed) :: qb1sg, qb1ng        
        attributes(managed) :: dqb1s, dqb1n
        attributes(managed) :: qb2s, qb2n
        attributes(managed) :: qb2sg, qb2ng        
        attributes(managed) :: dqb2s, dqb2n
        attributes(managed) :: qb3s, qb3n
        attributes(managed) :: qb3sg, qb3ng 
        attributes(managed) :: dqb3s, dqb3n
        attributes(managed) :: dq1x2o, dq2x2o, dq3x2o
        attributes(managed) :: dq1x2og, dq2x2og, dq3x2og        
#endif

      end module outflow_vars
!====================================================

!====================================================
      module ibm_param
        use constants
        use param
        implicit none
        integer :: npunt,npunr,npunz,mpun,mpunPpro
        integer :: nbfx,mpugeo,mpugeoPpro
        integer :: nbfxAS,nbfxPpro
        integer :: niix, nifx, njix, njfx, nkix, nkfx
        parameter (mpun=1000000)
        parameter (mpugeo=1000000)

        real(DP),dimension(mpugeo,9) :: xyzbfx,xyzbfxAS
        real(DP),dimension(mpugeo,3) :: barfx, nxyzfx
        real(DP),dimension(mpugeo,3) :: barfxAS, nxyzfxAS
        real(DP),dimension(mpugeo) :: areafx,areafxAS

        real(DP), dimension(4,mpun) :: vbfx, vbifx
        real(DP), dimension(4,mpun) :: vbfxAS, vbifxAS

        real(DP),dimension(m1,m2) :: plth1, plth2, kpllo,kplup


        real(DP),dimension(:,:),allocatable :: xyzbfxPpro,barfxPpro, nxyzfxPpro,vbfxPpro,vbifxPpro
        real(DP),dimension(:),allocatable :: areafxPpro
        integer,dimension(:,:,:),allocatable :: indgeo, indgeoe, indgeoee
        integer,dimension(:,:,:),allocatable :: indgeoAS, indgeoeAS, indgeoeeAS
        integer,dimension(:,:,:),allocatable :: indgeoPpro, indgeoePpro,indgeoeePpro
        integer,dimension(:,:),allocatable :: ntribfxPpro
#ifdef USE_CUDA
        attributes(managed) :: indgeo, indgeoe, indgeoee
        attributes(managed) :: indgeoAS, indgeoeAS, indgeoeeAS
        attributes(managed) :: indgeoPpro, indgeoePpro, indgeoeePpro
#endif

        integer,dimension(4) :: npunifx,npunfx
        integer,dimension(4) :: npunifxAS,npunfxAS
        integer,dimension(4) :: npunifxPpro,npunfxPpro

        real(DP),dimension(:,:),allocatable :: distb,distbAS,distbPpro
#ifdef USE_CUDA
        attributes(managed) :: distb,distbAS,distbPpro
#endif

        real(DP) :: ftAS, ftAS2, ftAS3

        integer,dimension(4,mpun)::ntribfx,ntribfxAS

        real(DP),dimension(mpun) :: q1bo,q2bo,q3bo
      end module ibm_param
!====================================================

!====================================================
      module mls_param
        use constants
        use param
        implicit none

        integer :: Nparticle, Nparticle_3d, Npx, Npy, Npz, Nmlsprobe
        integer :: NparticleEF_2d, NparticleEF_1d,NparticleCV_1d
        parameter (Nparticle=1,Nparticle_3d=1,Npx=1,Npy=1,Npz=1,Nmlsprobe=10)
        parameter (NparticleEF_2d=1,NparticleEF_1d=3,NparticleCV_1d=1)
        integer, parameter:: n_EFcheck=15   
        real(DP) :: h_probe
        parameter (h_probe=1.50D0)
        real(DP) ECG_D

        integer, dimension (n_EFcheck) :: nodicheckEF
!     --------variables for structural solver------------------------
      integer, parameter :: max_n_edge_of_vert=50
      integer, parameter :: max_n_edge_of_vert_1d=5
      integer, parameter :: max_n_edge_of_vert_3d=100

      integer, dimension(:),allocatable :: n_edge_of_vert
      integer, dimension(:,:), allocatable :: vert_of_edge
      integer, dimension(:,:), allocatable :: vert_of_face
      integer, dimension(:,:), allocatable :: edge_of_face
      integer, dimension(:,:), allocatable :: vert_of_vert
      integer, dimension(:,:), allocatable :: edge_of_vert
      integer, dimension(:,:), allocatable :: face_of_edge
      integer, dimension(:,:), allocatable :: v1234
      integer, dimension(:), allocatable:: vert_to_part
      integer, dimension(:), allocatable:: face_to_part
      integer, dimension(:), allocatable:: edge_to_part
      integer, dimension(:), allocatable:: vert_to_chamb
      integer, dimension(:), allocatable:: vert_to_chamb4V
      integer, dimension(:), allocatable:: edge_to_chamb
      integer, dimension(:), allocatable:: face_to_chamb
      integer, dimension(:), allocatable:: face_to_chamb4V

      integer, dimension(:),allocatable :: n_edge_of_vert_3d
      integer, dimension(:),allocatable :: n_cell_of_vert_3d
      integer, dimension(:),allocatable :: n_edge_of_cell_3d !2
      integer, dimension(:),allocatable :: n_cell_of_edge_3d !2
      integer, dimension(:,:), allocatable :: vert_of_cell_3d
      integer, dimension(:,:), allocatable :: vert_of_vert_3d
      integer, dimension(:,:), allocatable :: cell_of_vert_3d
      integer, dimension(:,:), allocatable :: vert_of_edge_3d !2
      integer, dimension(:,:), allocatable :: edge_of_cell_3d !2
      integer, dimension(:,:), allocatable :: cell_of_edge_3d !2
      integer, dimension(:,:), allocatable :: vert_of_face_3d !3
      integer, dimension(:,:), allocatable :: edge_of_vert_3d !3
      integer, dimension(:,:), allocatable :: face_of_cell_3d
      integer, dimension(:,:), allocatable :: cell_of_face_3d
      integer, dimension(:), allocatable :: count_edge_of_vert_3d !3

      integer, dimension(:), allocatable:: vert_to_part_3d
      integer, dimension(:), allocatable:: face_to_part_3d
      integer, dimension(:), allocatable:: edge_to_part_3d
      integer, dimension(:), allocatable:: cell_to_part_3d
      integer, dimension(:), allocatable:: vert_to_chamb_3d
      integer, dimension(:), allocatable:: vert_to_chamb_3dBou
      integer, dimension(:), allocatable:: vert_to_chamb_3dBouOld
      integer, dimension(:), allocatable:: vert_to_chamb_3d4V
      integer, dimension(:), allocatable:: face_to_chamb_3d
      integer, dimension(:), allocatable:: edge_to_chamb_3d
      integer, dimension(:), allocatable:: cell_to_chamb_3d
      integer, dimension(:), allocatable:: cell_to_chamb_3dBou
      integer, dimension(:), allocatable:: tag_2dwet


!      integer, dimension(:), allocatable:: LabelH
      ! integer, dimension(:), allocatable:: LabelBound
      integer, dimension(:), allocatable:: Segments_node_3d
      real(DP), dimension(3,17)  :: xyz_seg
      real(DP), dimension(4)  :: xyz_inf
      integer, dimension(:), allocatable:: LabelStenosi
      real(DP), dimension(:), allocatable:: LabelSmooth
      real(DP), dimension(:), allocatable:: LabelSmooth_cell
      ! real(DP), dimension(:), allocatable:: LabelDelay
      integer, dimension(:), allocatable:: LabelVaortic
      integer, dimension(:), allocatable:: LabelVmitral
      integer, dimension(:), allocatable:: LabelVpulmo
      integer, dimension(:), allocatable:: LabelVtricu      
      integer, dimension(:), allocatable:: LabelPurkSetto
      integer, dimension(:), allocatable:: LabelSkipConn

      integer, dimension(:),allocatable::dum_for, dum_forv
      integer, dimension(:,:),allocatable::pind,pindr,pindt,pindv,pindtOff
      integer, dimension(:,:),allocatable::pind_probeP,pind_probeN

      integer,dimension(:,:,:),allocatable::bboxind
      
      real(DP),dimension(:),allocatable::theta0,theta
      real(DP),dimension(:),allocatable::dist0,dist,dist00
      real(DP),dimension(:),allocatable::sur0,sur
      real(DP),dimension(:,:),allocatable::aalpha0,aalpha

      real(DP),dimension(:),allocatable::dist0_3d,dist_3d
      real(DP),dimension(:),allocatable::vol0_3d,vol_3d
      real(DP),dimension(:),allocatable::sur0_3d,sur_3d
      real(DP),dimension(:,:),allocatable::aalpha0_3d,aalpha_3d
      real(DP),dimension(:,:),allocatable::versCFface_3d
      real(DP),dimension(:),allocatable::distCFface_3d,g1interpface_3d
      real(DP),dimension(:,:),allocatable::cell_bar
      real(DP),dimension(:,:,:),allocatable::normalfaceofcells_3d
      real(DP),dimension(:,:,:),allocatable::AmatrFibers_cell_3d
      real(DP),dimension(:,:,:),allocatable::AmatrFibers_node_3d
      real(DP),dimension(:,:,:),allocatable::AmatrFibers_edge_3d
      real(DP),dimension(:),allocatable::edgefiber_cosangle_3d      
      real(DP),dimension(:,:,:),allocatable::Mintcells_3d
      real(DP),dimension(:,:,:),allocatable::Mextcells_3d
      real(DP),dimension(:,:,:),allocatable::Mintfaces_3d
      real(DP),dimension(:,:,:),allocatable::Mextfaces_3d
      ! Fibrillation and Scar variables
      real(DP),dimension(:),allocatable::IstimEF_3d,IstimEF_3dS1,IstimEF_3dS2
      real(DP),dimension(:,:),allocatable::xyz_gz,xyz_sc
      integer,dimension(:,:),allocatable::vert_of_face_gz,vert_of_face_sc
      integer,dimension(:,:),allocatable::vert_of_edge_gz,vert_of_edge_sc
      integer,dimension(:,:),allocatable::edge_of_face_gz,edge_of_face_sc

      real(DP),dimension(:),allocatable::mvol
      real(DP),dimension(:,:),allocatable::dismax
      real(DP),dimension(:,:),allocatable::meanxyz,meanvel

      real(DP),dimension(:),allocatable::Volume0,Volume
      real(DP),dimension(:),allocatable::Surface0,Surface
      real(DP),dimension(:),allocatable::Surface0_3d,Surface_3d
      real(DP),dimension(:),allocatable::Volume_chamb0,Volume_chamb,Volume_chambShift
      real(DP),dimension(:),allocatable::Surface_chamb0,Surface_chamb,Surface_chambShift
      real(DP),dimension(:),allocatable::celvol,Hboxx
      real(DP),dimension(:),allocatable::celvolr,Hboxxr
      real(DP),dimension(:),allocatable::shwtx,shwty,shwtz1,shwtz2

      real(DP), dimension(:,:), allocatable :: xyz,xyz0,xyzold
      real(DP), dimension(:,:), allocatable :: xyzv,xyzv0,xyza
      real(DP), dimension(:,:), allocatable :: xyza0

      real(DP), dimension(:,:), allocatable :: xyzMLSprobe,outMLSprobe
      integer, dimension(:,:),allocatable:: pindMLSprobe
      
      real(DP), dimension(:,:), allocatable :: xyz_3d,xyz0_3d,xyzold_3d
      real(DP), dimension(:,:), allocatable :: xyzv_3d,xyzv0_3d,xyza_3d
      real(DP), dimension(:,:), allocatable :: xyza0_3d

      real(DP), dimension(:,:), allocatable :: XEF_3d
      real(DP), dimension(:), allocatable :: potEFcell_3d,XEF_3ddt
      real(DP), dimension(:), allocatable :: potEFnode_3d
      real(DP), dimension(:), allocatable :: potEFface_3d
      !--------------------------------------------------- 
      !APD
      real(DP), dimension(:,:), allocatable :: t_apd_3d,t_apd_2d
      real(DP), dimension(:), allocatable :: old_1,old_1_2d
      integer, dimension(:), allocatable :: f_apd,f_apd_2d,scar_tag_2d
      real(DP), dimension(:), allocatable :: LATface_3d,LATnode_3d
      !--------------------------------------------------- 
      !CV calc
      real(DP), dimension(:), allocatable :: CVdiv
      real(DP), dimension(:,:), allocatable :: CVcell,CVnode,CVface_3d,CVrot
      real(DP), dimension(:,:,:), allocatable :: CVgrad_cell,CVgrad_node
      !Laplace interp
      !---------------------------------------------------
      integer,dimension(:),allocatable :: myotag
      real(DP),dimension(:),allocatable ::meshquality_3d
      !---------------------------------------------------
      ! CARTO diffusivity
      !--------------------------------------------------- 
      real(DP),dimension(:),allocatable ::CARTO_Dnode,CARTO_Dnode3d,CARTO_Dface,CARTO_Dcell3d,CARTO_Dface3d
      real(DP)::cartodelta
      !---------------------------------------------------
      !EGM
      !---------------------------------------------------
#ifdef EGM
      integer:: Negm
      real(DP),dimension(:,:),allocatable ::xyz_egm
      real(DP),dimension(:),allocatable ::EGMvec
#endif
      !---------------------------------------------------
      !Scar 3D
      integer,dimension(:),allocatable::scar_cell
      !--------------------------------------------------- 
      real(DP), dimension(:,:), allocatable :: gradcell_3d
      real(DP), dimension(:,:), allocatable :: gradface_3d

      real(DP), dimension(:), allocatable :: ECG_PHI
      real(DP), dimension(:,:), allocatable :: ECG_E

#ifdef BIDOMAIN      
      real(DP), dimension(:), allocatable :: potextEFcell_3d
      real(DP), dimension(:), allocatable :: potextEFnode_3d
      real(DP), dimension(:), allocatable :: potextEFface_3d
      real(DP), dimension(:,:), allocatable :: gradextcell_3d
      real(DP), dimension(:,:), allocatable :: gradextface_3d
      real(DP), dimension(:,:), allocatable :: QEFbido
      real(DP), dimension(:,:), allocatable :: HEFbido
      real(DP), dimension(:,:), allocatable :: checkEFbido
      integer, dimension(:), allocatable :: IPIVEFbido
#endif
      
      integer, dimension(:),allocatable :: n_edge_of_vertEF_2d
      integer, dimension(:,:), allocatable :: vert_of_edgeEF_2d
      integer, dimension(:,:), allocatable :: vert_of_faceEF_2d
      integer, dimension(:,:), allocatable :: edge_of_faceEF_2d
      integer, dimension(:,:), allocatable :: vert_of_vertEF_2d
      integer, dimension(:,:), allocatable :: edge_of_vertEF_2d
      integer, dimension(:,:), allocatable :: face_of_edgeEF_2d
      integer, dimension(:,:), allocatable :: face_of_vertEF_2d
      integer, dimension(:), allocatable :: n_face_of_vertEF_2d
      integer, dimension(:),allocatable :: n_edge_of_vertEF_1d
      integer, dimension(:,:), allocatable :: vert_of_edgeEF_1d
      integer, dimension(:,:), allocatable :: vert_of_vertEF_1d
      integer, dimension(:,:), allocatable :: edge_of_vertEF_1d
      integer, dimension(:),allocatable :: n_edge_of_vertCV_1d
      integer, dimension(:,:), allocatable :: vert_of_edgeCV_1d
      integer, dimension(:,:), allocatable :: vert_of_vertCV_1d
      integer, dimension(:,:), allocatable :: edge_of_vertCV_1d

      real(DP), dimension(:,:), allocatable :: xyzEF_2d,xyz0EF_2d,xyzSEF_2d,xyzS1EF_2d
      !CV
      real(DP), dimension(:), allocatable :: CVface_EF_2d
      
      real(DP),dimension(:),allocatable::dist0EF_2d,distEF_2d
      real(DP),dimension(:),allocatable::sur0EF_2d,surEF_2d
      real(DP),dimension(:,:),allocatable::tri_barEF_2d
      real(DP),dimension(:,:),allocatable::tri_norEF_2d
      real(DP),dimension(:),allocatable::Surface0EF_2d,SurfaceEF_2d
      real(DP),dimension(:,:),allocatable::versCFedgeEF_2d
      real(DP),dimension(:),allocatable::distCFedgeEF_2d,g1interpedgeEF_2d
      real(DP),dimension(:,:,:),allocatable::normaledgeoffacesEF_2d
      real(DP),dimension(:,:,:),allocatable::AmatrFibersEF_2d
      real(DP), dimension(:,:), allocatable :: XEF_2d
      real(DP),dimension(:,:,:),allocatable::MintfacesEF_2d
      real(DP),dimension(:,:,:),allocatable::MextfacesEF_2d
      real(DP),dimension(:,:,:),allocatable::MintedgesEF_2d
      real(DP),dimension(:,:,:),allocatable::MextedgesEF_2d
      real(DP), dimension(:,:), allocatable :: gradfaceEF_2d
      real(DP), dimension(:,:), allocatable :: gradedgeEF_2d
      real(DP), dimension(:), allocatable :: potEFface_2d
      real(DP), dimension(:), allocatable :: potEFedge_2d
      real(DP), dimension(:), allocatable :: potEFnode_2d
      real(DP),dimension(:),allocatable::IstimEF_2dS1,IstimEF_2dS2

      real(DP), dimension(:,:), allocatable :: xyzEF_1d,xyz0EF_1d,xyzSEF_1d      
      real(DP), dimension(:,:), allocatable :: xyzCV_1d,xyz0CV_1d,xyzSCV_1d      
      real(DP), dimension(:,:), allocatable :: XEF_1d
      real(DP),dimension(:),allocatable::dist0EF_1d,distEF_1d
      real(DP),dimension(:,:),allocatable::edg_barEF_1d
      real(DP),dimension(:),allocatable::IstimEF_1d
      real(DP),dimension(:),allocatable::MintedgesEF_1d
!      real(DP),dimension(:,:,:),allocatable::MextedgesEF_1d                                                                                         
      real(DP),dimension(:),allocatable::MintnodesEF_1d
 !     real(DP),dimension(:,:,:),allocatable::MextnodesEF_1d                                                                                         
      real(DP), dimension(:), allocatable :: potEFnode_1d
      integer :: nBundAV,nBundPurk,nBundAtri,nPurkVentr
      integer, dimension(:), allocatable:: vert_to_partEF_1d
      integer, dimension(:), allocatable:: edge_to_partEF_1d
      integer, dimension(:),allocatable :: vAVmaster
      integer, dimension(:),allocatable :: vAVslave
      integer, dimension(:),allocatable :: vPurkmaster
      integer, dimension(:),allocatable :: fPurkslave
      integer, dimension(:),allocatable :: vAtrimaster
      integer, dimension(:),allocatable :: cAtrislave
      integer, dimension(:),allocatable :: fVentrmaster
      integer, dimension(:),allocatable :: cVentrslave
      integer, dimension(:),allocatable :: fcheckEF_2d



!STRONG
#ifdef STRONG
      real(DP), dimension(:,:), allocatable :: xyzp,xyzvp
      real(DP), dimension(:,:), allocatable :: xyzak,xyzakm1,xyzauk,xyzakp1
      real(DP), dimension(:,:), allocatable :: xyzvk,xyzvkm1,xyzvuk,xyzvkp1
      real(DP), dimension(:,:), allocatable :: xyzk,xyzkp1
      real(DP), dimension(:,:), allocatable :: xyzp_3d,xyzvp_3d
      real(DP), dimension(:,:), allocatable :: xyzak_3d,xyzakm1_3d,xyzauk_3d,xyzakp1_3d
      real(DP), dimension(:,:), allocatable :: xyzvk_3d,xyzvkm1_3d,xyzvuk_3d,xyzvkp1_3d
      real(DP), dimension(:,:), allocatable :: xyzk_3d,xyzkp1_3d
#endif

      real(DP),dimension(:,:),allocatable::fpxyz
      real(DP),dimension(:,:),allocatable::fxyz
      real(DP),dimension(:,:),allocatable::fxyzm
      real(DP),dimension(:),allocatable::tauface,pressface,taufaceAV,pressfaceAV
      real(DP),dimension(:),allocatable:: mass_of_vert
      real(DP),dimension(:),allocatable:: AFung
      real(DP),dimension(:),allocatable:: StiffV

      real(DP),dimension(:,:),allocatable::fpxyz_3d,fxyz_3d,fxyzm_3d
      real(DP),dimension(:),allocatable:: mass_of_vert_3d
      real(DP),dimension(:),allocatable:: AFung_3d
      real(DP),dimension(:),allocatable:: StiffV_3d
 
      real(DP),dimension(:,:),allocatable::tri_ver,tri_vel,vel_tri,acc_tri
      real(DP),dimension(:,:),allocatable::tri_bar,tri_nor

      integer, dimension(:), allocatable :: mytr, mytl, mypr
      integer, dimension(:), allocatable :: trcnt, tlcnt, tpcnt
      integer :: trcnt_h,vertLV,vertLA

#ifdef USE_CUDA
      attributes(managed) :: n_edge_of_vert
      attributes(managed) :: vert_of_edge
      attributes(managed) :: vert_of_face
      attributes(managed) :: vert_of_vert
      attributes(managed) :: edge_of_vert
      attributes(managed) :: face_of_edge
      attributes(managed) :: v1234
      attributes(managed) :: vert_to_part
      attributes(managed) :: face_to_part
      attributes(managed) :: edge_to_part
      attributes(managed) :: vert_to_chamb
      attributes(managed) :: vert_to_chamb4V
      attributes(managed) :: edge_to_chamb
      attributes(managed) :: face_to_chamb
      attributes(managed) :: face_to_chamb4V


      attributes(managed) :: n_edge_of_vert_3d
      attributes(managed) :: n_cell_of_vert_3d
      attributes(managed) :: n_edge_of_cell_3d
      attributes(managed) :: n_cell_of_edge_3d
      attributes(managed) :: vert_of_cell_3d
      attributes(managed) :: vert_of_vert_3d
      attributes(managed) :: cell_of_vert_3d
      attributes(managed) :: vert_of_edge_3d
      attributes(managed) :: edge_of_cell_3d
      attributes(managed) :: cell_of_edge_3d
      attributes(managed) :: vert_of_face_3d
      attributes(managed) :: edge_of_vert_3d
      attributes(managed) :: vert_to_part_3d
      attributes(managed) :: face_to_part_3d
      attributes(managed) :: edge_to_part_3d
      attributes(managed) :: cell_to_part_3d
      attributes(managed) :: vert_to_chamb_3d
      attributes(managed) :: vert_to_chamb_3dBou
      attributes(managed) :: vert_to_chamb_3d4V
      attributes(managed) :: face_to_chamb_3d
      attributes(managed) :: edge_to_chamb_3d
      attributes(managed) :: cell_to_chamb_3d
      attributes(managed) :: cell_to_chamb_3dBou,vert_to_chamb_3dBouold
      attributes(managed) :: tag_2dwet
      attributes(managed) :: face_of_cell_3d
      attributes(managed) :: cell_of_face_3d
      attributes(managed) :: versCFface_3d
      attributes(managed) :: distCFface_3d
      attributes(managed) :: g1interpface_3d


!      attributes(managed) :: LabelH
      ! attributes(managed) :: LabelBound
      attributes(managed) :: LabelStenosi
      attributes(managed) :: LabelSmooth
      attributes(managed) :: LabelSmooth_cell
      ! attributes(managed) :: LabelDelay
      attributes(managed) :: LabelVaortic
      attributes(managed) :: LabelVmitral
      attributes(managed) :: LabelVpulmo
      attributes(managed) :: LabelVtricu
      attributes(managed) :: LabelPurkSetto
      attributes(managed) :: LabelSkipConn

      attributes(managed) :: dum_for, dum_forv
      attributes(managed) :: pind, pindt, pindv, pindtOff
      attributes(managed) :: pind_probeP,pind_probeN
      attributes(managed) :: AFung
      attributes(managed) :: StiffV
      attributes(managed) :: AFung_3d
      attributes(managed) :: StiffV_3d

      attributes(managed) :: theta0,theta
      attributes(managed) :: dist00,dist0,dist
      attributes(managed) :: sur0, sur
      attributes(managed) :: aalpha0,aalpha
      attributes(managed) :: celvol, Hboxx

      attributes(managed) :: dist0_3d,dist_3d
      attributes(managed) :: vol0_3d, vol_3d
      attributes(managed) :: sur0_3d, sur_3d
      attributes(managed) :: aalpha0_3d,aalpha_3d
      attributes(managed) :: cell_bar
      attributes(managed) :: normalfaceofcells_3d
      attributes(managed) :: AmatrFibers_cell_3d
      attributes(managed) :: AmatrFibers_node_3d
      attributes(managed) :: AmatrFibers_edge_3d
      attributes(managed) :: edgefiber_cosangle_3d
      attributes(managed) :: Mintcells_3d
      attributes(managed) :: Mextcells_3d
      attributes(managed) :: Mintfaces_3d
      attributes(managed) :: Mextfaces_3d
      attributes(managed) :: IstimEF_3d,IstimEF_3dS1,IstimEF_3dS2
      attributes(managed) :: xyz_gz,xyz_sc
      attributes(managed) :: vert_of_edge_gz,vert_of_edge_sc
      attributes(managed) :: edge_of_face_gz,edge_of_face_sc
      attributes(managed) :: vert_of_face_gz,vert_of_face_sc
      attributes(managed) :: edge_of_face,count_edge_of_vert_3d !parallel_ini
      
      attributes(managed) :: Volume0, Volume
      attributes(managed) :: Surface0, Surface
      attributes(managed) :: Surface0_3d, Surface_3d
      attributes(managed) :: Volume_chamb0, Volume_chamb, Volume_chambShift
      attributes(managed) :: Surface_chamb0, Surface_chamb, Surface_chambShift

      attributes(managed) :: xyz,xyz0,xyzold
      attributes(managed) :: xyzv,xyzv0,xyza
      attributes(managed) :: xyza0

      attributes(managed) :: xyzMLSprobe,outMLSprobe
      attributes(managed) :: pindMLSprobe

      attributes(managed) :: xyz_3d,xyz0_3d,xyzold_3d
      attributes(managed) :: xyzv_3d,xyzv0_3d,xyza_3d
      attributes(managed) :: xyza0_3d
      attributes(managed) :: XEF_3d
      attributes(managed) :: potEFcell_3d,XEF_3ddt
      attributes(managed) :: potEFnode_3d
      attributes(managed) :: potEFface_3d
      attributes(managed) :: gradcell_3d
      attributes(managed) :: gradface_3d
      !APD
      attributes(managed) :: t_apd_3d,old_1,f_apd
      attributes(managed) :: t_apd_2d,old_1_2d,f_apd_2d
      attributes(managed) :: LATface_3d,LATnode_3d
      !CV valc
      attributes(managed) :: CVcell,CVnode,CVgrad_cell,CVgrad_node,CVface_3d
      attributes(managed) :: CVdiv,CVrot
      !Laplace interp
      attributes(managed) :: myotag,meshquality_3d
      !CARTO variables
      attributes(managed) :: CARTO_Dnode,CARTO_Dnode3d,CARTO_Dface,CARTO_Dcell3d,scar_tag_2d
      attributes(managed) :: CARTO_Dface3d
      !EGM
#ifdef EGM
      attributes(managed) :: EGMvec,xyz_egm
#endif
      !Scar 3D
      attributes(managed) :: scar_cell
      
      attributes(managed) :: ECG_PHI 
      attributes(managed) :: ECG_E 

#ifdef BIDOMAIN
      attributes(managed) :: potextEFcell_3d
      attributes(managed) :: potextEFnode_3d
      attributes(managed) :: potextEFface_3d
      attributes(managed) :: gradextcell_3d
      attributes(managed) :: gradextface_3d
      attributes(managed) :: QEFbido
      attributes(managed) :: HEFbido
      attributes(managed) :: checkEFbido
#endif
! #endif                                                                                                                                             
! #ifdef BUNDLE                                                                                                                                      
      attributes(managed) :: xyzEF_2d,xyz0EF_2d,xyzSEF_2d,xyzS1EF_2d
      attributes(managed) :: n_edge_of_vertEF_2d
      attributes(managed) :: vert_of_edgeEF_2d
      attributes(managed) :: vert_of_faceEF_2d
      attributes(managed) :: edge_of_faceEF_2d
      attributes(managed) :: vert_of_vertEF_2d
      attributes(managed) :: edge_of_vertEF_2d
      attributes(managed) :: face_of_edgeEF_2d
      attributes(managed) :: face_of_vertEF_2d
      attributes(managed) :: n_face_of_vertEF_2d
      attributes(managed) :: dist0EF_2d,distEF_2d
      attributes(managed) :: sur0EF_2d,surEF_2d
      attributes(managed) :: tri_barEF_2d,tri_norEF_2d
      attributes(managed) :: SurfaceEF_2d,Surface0EF_2d
      attributes(managed) :: versCFedgeEF_2d
      attributes(managed) :: distCFedgeEF_2d,g1interpedgeEF_2d
      attributes(managed) :: normaledgeoffacesEF_2d
      attributes(managed) :: AmatrFibersEF_2d
      attributes(managed) :: XEF_2d
      
      attributes(managed) :: CVface_EF_2d
      
      attributes(managed) :: MintfacesEF_2d
      attributes(managed) :: MextfacesEF_2d
      attributes(managed) :: MintedgesEF_2d
      attributes(managed) :: MextedgesEF_2d
      attributes(managed) :: gradfaceEF_2d
      attributes(managed) :: gradedgeEF_2d
      attributes(managed) :: potEFface_2d
      attributes(managed) :: potEFedge_2d
      attributes(managed) :: potEFnode_2d
      attributes(managed) :: IstimEF_2dS1,IstimEF_2dS2
      
      attributes(managed) :: xyzEF_1d,xyz0EF_1d,xyzSEF_1d      
      attributes(managed) :: xyzCV_1d,xyz0CV_1d,xyzSCV_1d      
      attributes(managed) :: n_edge_of_vertEF_1d
      attributes(managed) :: vert_of_edgeEF_1d
      attributes(managed) :: vert_of_vertEF_1d
      attributes(managed) :: edge_of_vertEF_1d
      attributes(managed) :: n_edge_of_vertCV_1d
      attributes(managed) :: vert_of_edgeCV_1d
      attributes(managed) :: vert_of_vertCV_1d
      ! attributes(managed) :: edge_of_vertCV_1d
      attributes(managed) :: edg_barEF_1d
      attributes(managed) :: dist0EF_1d,distEF_1d
      attributes(managed) :: XEF_1d
      attributes(managed) :: IstimEF_1d
      attributes(managed) :: MintedgesEF_1d
!      attributes(managed) :: MextedgesEF_1d                                                                                                        
      attributes(managed) :: MintnodesEF_1d
!      attributes(managed) :: MextnodesEF_1d                                                                                                        
      attributes(managed) :: potEFnode_1d
      attributes(managed) :: vert_to_partEF_1d
      attributes(managed) :: edge_to_partEF_1d

!STRONG
#ifdef STRONG
      attributes(managed) :: xyzp,xyzvp
      attributes(managed) :: xyzak,xyzakm1,xyzauk,xyzakp1
      attributes(managed) :: xyzvk,xyzvkm1,xyzvuk,xyzvkp1
      attributes(managed) :: xyzk,xyzkp1
      attributes(managed) :: xyzp_3d,xyzvp_3d
      attributes(managed) :: xyzak_3d,xyzakm1_3d,xyzauk_3d,xyzakp1_3d
      attributes(managed) :: xyzvk_3d,xyzvkm1_3d,xyzvuk_3d,xyzvkp1_3d
      attributes(managed) :: xyzk_3d,xyzkp1_3d
#endif

      attributes(managed) :: fpxyz, fxyz,fxyzm
      attributes(managed) :: mass_of_vert
      attributes(managed) :: tri_ver,tri_vel, vel_tri, acc_tri
      attributes(managed) :: tri_bar, tri_nor
      attributes(managed) :: tauface,pressface,taufaceAV,pressfaceAV

      attributes(managed) :: fpxyz_3d,fxyz_3d,fxyzm_3d
      attributes(managed) :: mass_of_vert_3d

      attributes(managed) :: mytr, mytl, mypr
      attributes(managed) :: trcnt, tlcnt, tpcnt
#endif
      
     
      integer :: nvtot, netot, nftot, nttot, count2,count2EF_2d
      integer :: nvtot_3d, netot_3d, nftot_3d, nctot_3d,count2_3d
      integer :: inpstart_2dstr,inpend_2dstr,inpstart_2dwet,inpend_2dwet !NON SO SE PASSARLE A CUDA
      integer :: nvstart_2dstr,nvend_2dstr,nestart_2dstr,neend_2dstr,nfstart_2dstr,nfend_2dstr
      integer :: nvstart_2dwet,nvend_2dwet,nestart_2dwet,neend_2dwet,nfstart_2dwet,nfend_2dwet
      integer :: maxnv,maxne,maxnf
      integer :: ntOffstart,ntOffend,ntOfftot
      integer :: inpstartOff,inpendOff
! #ifdef BUNDLE                                                                                                                                     
      integer :: nvtotEF_2d, netotEF_2d, nftotEF_2d
      integer :: nvtotEF_1d, netotEF_1d
      integer :: nvtotCV_1d, netotCV_1d
! #endif 
      integer   :: lead_LV,lead_RV,lead_RA

      real(DP)    :: usVolbub,usVolele,dOff
      real(DP)    :: tstartEF_LV,tstartEF_LA
!     inclusions for objects with different properties
      character *50,dimension(:),allocatable :: geofile,geofile_3d,geofileEF_2d,geofileEF_1d,geofileCV_1d

      integer,dimension(:),allocatable:: nvi,nei,nfi,ntilei
      integer,dimension(:),allocatable:: nvi_3d,nei_3d,nfi_3d,nci_3d
      integer,dimension(:),allocatable:: vstart,vend    !vertex start and end of particles  
      integer,dimension(:),allocatable:: estart,eend    !edge start and end of particles 
      integer,dimension(:),allocatable:: fstart,fend    !face start and end of particles 
      integer,dimension(:),allocatable:: tilestart,tileend    !tile start and end of particles 
      integer,dimension(:),allocatable:: vstart_3d,vend_3d    !vertex start and end of particles 3d  
      integer,dimension(:),allocatable:: estart_3d,eend_3d    !edge start and end of particles 3d
      integer,dimension(:),allocatable:: fstart_3d,fend_3d    !face start and end of particles 3d
      integer,dimension(:),allocatable:: cstart_3d,cend_3d    !face start and end of particles 3d
      integer,dimension(:),allocatable:: tstart,tend    !tile start and end

! #ifdef BUNDLE                                                                                                                                     
      integer,dimension(:),allocatable:: nviEF_2d,neiEF_2d,nfiEF_2d
      integer,dimension(:),allocatable:: nviEF_1d,neiEF_1d
      integer,dimension(:),allocatable:: nviCV_1d,neiCV_1d
      integer,dimension(:),allocatable:: vstartEF_2d,vendEF_2d    !vertex start and end of particles                                                
      integer,dimension(:),allocatable:: estartEF_2d,eendEF_2d    !edge start and end of particles                                                  
      integer,dimension(:),allocatable:: fstartEF_2d,fendEF_2d   !face start and end of particles                                                   
      integer,dimension(:),allocatable:: vstartEF_1d,vendEF_1d    !vertex start and end of particles                                                
      integer,dimension(:),allocatable:: estartEF_1d,eendEF_1d    !edge start and end of particless                                                 
      integer,dimension(:),allocatable:: vstartCV_1d,vendCV_1d    !vertex start and end of particles                                                
      integer,dimension(:),allocatable:: estartCV_1d,eendCV_1d    !edge start and end of particless                                                 
! #endif

      integer,dimension(:,:),allocatable:: boundary1   
      integer,dimension(:,:),allocatable:: boundary2   
      integer,dimension(:,:),allocatable:: boundary1_3d   
      integer,dimension(:,:),allocatable:: boundary2_3d  
      integer,dimension(:),allocatable:: boundary2EF_2d
      integer,dimension(:),allocatable:: vert_boundEF_2d  
      integer,dimension(:,:),allocatable:: vert_of_vert_boundEF_2d  

      integer,dimension(:),allocatable:: boundary2EF_1d  
      integer,dimension(:),allocatable:: boundary2CV_1d  
      integer,dimension(:),allocatable::iopen_pos,iopen_neg
      integer,dimension(:),allocatable::body_nsl
      integer,dimension(:),allocatable::wet_3d
      integer,dimension(:),allocatable:: faceid_t
      integer,dimension(:),allocatable:: fcol,fcolv
      integer :: ContactA, ContactM, ContactP, ContactT, Npresstarget

      real(DP),dimension(:),allocatable :: rke,rkc,rkb,kv,kat,kal,kae,kb,ke
      real(DP),dimension(:),allocatable :: rhos,thck,cv1,cv2
      real(DP),dimension(:),allocatable :: rke_3d,rkc_3d,kv_3d,kae_3d,ke_3d,f_act_3d
      real(DP),dimension(:),allocatable :: rhos_3d,cv1_3d,cv2_3d
      real(DP),dimension(:),allocatable :: rpm,body_flux
      
      real(DP),dimension(:),allocatable :: EFtstart_3d,EFtstart_2d,EFtstart_1d
      real(DP),dimension(:),allocatable :: astressEFcell_3d,astressEFnode_3d,astressEFedge_3d
      real(DP),dimension(:),allocatable :: Purk_tstim
      real(DP),dimension(:),allocatable :: Bund_tstim
      real(DP),dimension(:,:),allocatable:: BCoffset   
      real(DP),dimension(:,:),allocatable:: BCoffset_3d
      real(DP),dimension(:,:),allocatable:: BCoffsetEF_2d   
      real(DP),dimension(:,:),allocatable:: BCoffsetEF_1d     
      real(DP),dimension(:,:),allocatable:: BCoffsetCV_1d     
      real(DP):: timestartEF, timestartEFLA
      real(DP):: minpotLV,minpotLA,minpotRV,minpotRA

#ifdef USE_CUDA
      attributes(managed) :: nfi
      attributes(managed) :: tstart
!      attributes(managed) :: boundary1
      attributes(managed) :: boundary2
      attributes(managed) :: BCoffset
      attributes(managed) :: iopen_pos, iopen_neg, body_nsl
      attributes(managed) :: faceid_t
      attributes(managed) :: fcol, fcolv
      attributes(managed) :: rke, rkc, rkb, kv, kat, kal, kae, kb, ke
      attributes(managed) :: rhos,thck, cv1, cv2, body_flux
      attributes(managed) :: rke_3d, rkc_3d, kv_3d, kae_3d, ke_3d,f_act_3d
      attributes(managed) :: rhos_3d,cv1_3d, cv2_3d
      attributes(managed) :: EFtstart_3d,EFtstart_2d,EFtstart_1d,Purk_tstim,Bund_tstim
      attributes(managed) :: astressEFcell_3d,astressEFnode_3d,astressEFedge_3d
!      attributes(managed) :: boundary1_3d
      attributes(managed) :: boundary2_3d
      attributes(managed) :: BCoffset_3d

      attributes(managed) :: boundary2EF_2d
      attributes(managed) :: boundary2EF_1d
      attributes(managed) :: boundary2CV_1d
      attributes(managed) :: vert_boundEF_2d  
      attributes(managed) :: vert_of_vert_boundEF_2d        
      attributes(managed) :: BCoffsetEF_2d   
      attributes(managed) :: BCoffsetEF_1d   
      attributes(managed) :: BCoffsetCV_1d   
#endif

      end module mls_param
!====================================================
      module mls_local
        use constants
        use param
        implicit none
        real(DP), dimension(:,:,:), allocatable :: for_xc, for_yc, for_zc
        real(DP), dimension(:,:,:), allocatable :: for_te, for_sc
        integer, dimension(:,:,:), allocatable :: coll, phfield
#ifdef USE_CUDA
        attributes(managed) :: for_xc, for_yc, for_zc
        attributes(managed) :: coll
#endif
      end module mls_local
!====================================================
      module tile_arrays
        implicit none
      
        real,dimension(:,:,:), allocatable :: albegaBar
        integer, dimension(:,:), allocatable :: tri_tiling
        integer :: Navamax
#ifdef USE_CUDA
        attributes(managed) :: albegaBar, tri_tiling
#endif
      end module tile_arrays
!====================================================
      module probes
        use constants
       integer :: nson,mason(65,3),masonAO(3), masonLV(3), masonLA(3)
       integer :: ECG_nson
       integer :: masonAP(3), masonRV(3), masonRA(3)
       real(DP) :: coson(65,3), cosonAO(3), cosonLV(3),cosonLA(3)
       real(DP) :: cosonAP(3), cosonRV(3),cosonRA(3)
       integer :: rankAO, rankLV, rankLA, rankAP, rankRV, rankRA
       integer :: nodesMovingProbes(3,2)
      end module probes
!====================================================                                                                               
      module rt_arrays
        use param
        ! ray tracing arrays                                                                                                           

        ! color maps                                                                                                                   
        integer,allocatable,dimension(:,:,:,:) :: color
        integer, allocatable, dimension(:,:,:) :: ibpointX, ibpointY, ibpointZ, ibpoint
        integer, allocatable,dimension(:,:) :: intersect

        #ifdef USE_CUDA
        ! color maps                                                                                                                   
        attributes(managed) :: color
        attributes(managed) :: ibpointX, ibpointY, ibpointZ, ibpoint
        attributes(managed) :: intersect
#endif
      end module rt_arrays
