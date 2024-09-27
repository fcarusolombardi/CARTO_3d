!***********************************************************************
      subroutine inirea
      use mpih
      use mpi_param, only: kstart,kend,kstartr,kendr
      use local_arrays, only: dens,q2,q3,q1,dsal,pr
      use param
      use init_itp
      IMPLICIT NONE
      real(DP) :: dummyr 
      character(100)  :: filcnw2
      character(100) :: ipfi
      integer :: i,j
      integer, parameter :: ghosts = 4
      integer  :: kstartog,kendog,kstartogr,kendogr

      write(ipfi,199) nwrit
 199   format(i8.8)
      
!  Reading old grid information by rank0
      if (myid .eq. 0) then
        filcnw2 = trim(folderload)//'continua_grid'//trim(ipfi)//'.dat'
        open(13,file=filcnw2,status='old')
        rewind(13)                                                      
        read(13,*) n1o,n2o,n3o
        read(13,*) dummyr,dummyr,time
        read(13,*) istr3o,str3o
        read(13,*) mref1o, mref2o, mref3o
        close(13)
        write(*,'(5x,a)') 'old grid info read'
      endif
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!   Bcasting old grid information and time
      call MPI_BCAST(n1o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(n2o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(n3o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(istr3o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(str3o,1,MDP,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(time,1,MDP,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mref1o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mref2o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mref3o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
      
!   Check whether grid specifications have been updated
      IF( n2o.ne.n2 .or. n3o.ne.n3 .or. n1o.ne.n1  &
         .or. istr3o.ne.istr3 .or. str3o.ne.str3   &
         .or. mref1.ne.mref1o .or. mref2.ne.mref2o &
         .or. mref3.ne.mref3o) then
        if(myid.eq.0) write(*,'(5x,a)') "Interpolating new grid"
        if(n1.gt.n1o*2.or.n2.gt.n2o*2.or.n3.gt.n3o*2) then
          if(myid.eq.0) write(*,*) "New grid resolution can not be more &
                       & than twice the old resolution"
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
        endif

        n1om = n1o - 1
        n2om = n2o - 1
        n3om = n3o - 1       
        n1omr = n1om*mref1o
        n2omr = n2om*mref2o
        n3omr = n3om*mref3o
        n1or = n1omr + 1
        n2or = n2omr + 1
        n3or = n3omr + 1

        allocate(xcold(1:n1o),xmold(1:n1o))
        allocate(ycold(1:n2o),ymold(1:n2o))
        allocate(zcold(1:n3o),zmold(1:n3o))
        allocate(xcrold(1:n1or),xmrold(1:n1or))
        allocate(ycrold(1:n2or),ymrold(1:n2or))
        allocate(zcrold(1:n3or),zmrold(1:n3or))

        call iniitp_grid

        allocate(iiq1c1(1:n1),iiq1c2(1:n2),iiq1c3(1:n3))
        allocate(iiq2c1(1:n1),iiq2c2(1:n2),iiq2c3(1:n3))
        allocate(iiq3c1(1:n1),iiq3c2(1:n2),iiq3c3(1:n3))
        allocate(iidec1(1:n1),iidec2(1:n2),iidec3(1:n3))
        allocate(iisac1(1:n1r),iisac2(1:n2r),iisac3(1:n3r)) 

        call iniitp_indc
      
!   dens
#ifdef CALTEMP
        kstartog = max(iidec3(kstart),1)
        kendog   = min(iidec3(kend)+1,n3om)
        allocate(varold(1:n1o,1:n2o,kstartog-1:kendog+1))

        call mpi_read_continua(n1o,n2o,n3o,kstartog,kendog,4, &
               varold(1:n1o,1:n2o,kstartog-1:kendog+1))
        call interp_dens(kstartog,kendog)

        deallocate(varold)
        if(myid.eq.0)write(*,'(5x,a)')'Density interpolated!'
#endif

!   q1
        kstartog = max(iiq1c3(kstart),1)
        kendog   = min(iiq1c3(kend)+1,n3om)
        allocate(varold(1:n1o,1:n2o,kstartog-1:kendog+1))

        call mpi_read_continua(n1o,n2o,n3o,kstartog,kendog,1, &
               varold(1:n1o,1:n2o,kstartog-1:kendog+1))
        call interp_velx(kstartog,kendog)

        deallocate(varold)

!   q2
        kstartog = max(iiq2c3(kstart),1)
        kendog   = min(iiq2c3(kend)+1,n3om)
        allocate(varold(1:n1o,1:n2o,kstartog-1:kendog+1))

        call mpi_read_continua(n1o,n2o,n3o,kstartog,kendog,2, &
               varold(1:n1o,1:n2o,kstartog-1:kendog+1))
        call interp_vely(kstartog,kendog)

        deallocate(varold)

!   q3
        kstartog = max(iiq3c3(kstart),1)
        kendog   = min(iiq3c3(kend)+1,n3om)
        allocate(varold(1:n1o,1:n2o,kstartog-1:kendog+1))

        call mpi_read_continua(n1o,n2o,n3o,kstartog,kendog,3, &
               varold(1:n1o,1:n2o,kstartog-1:kendog+1))
        call interp_velz(kstartog,kendog)

        deallocate(varold)
        if(myid.eq.0)write(*,'(5x,a)')'Velocity interpolated!'

!   dsal
#ifdef CALSCAL
        kstartogr = max(iisac3(kstartr),1)
        kendogr   = min(iisac3(kendr)+1,n3omr)
        allocate(varold(1:n1or,1:n2or,kstartogr-1:kendogr+1))

        call mpi_read_continua(n1or,n2or,n3or,kstartogr,kendogr,5, &
               varold(1:n1or,1:n2or,kstartogr-1:kendogr+1))
        call interp_dsal(kstartogr,kendogr)

        deallocate(varold)
        if(myid.eq.0)write(*,'(5x,a)')'Salinity interpolated!'
#endif

!  deallocate old grid

        deallocate(xcold,xmold)
        deallocate(ycold,ymold)
        deallocate(zcold,zmold)
        deallocate(xcrold,xmrold)
        deallocate(ycrold,ymrold)
        deallocate(zcrold,zmrold)

        deallocate(iiq1c1,iiq1c2,iiq1c3)
        deallocate(iiq2c1,iiq2c2,iiq2c3)
        deallocate(iiq3c1,iiq3c2,iiq3c3)
        deallocate(iidec1,iidec2,iidec3)
        deallocate(iisac1,iisac2,iisac3)

      ELSE

        if(myid.eq.0)write(*,'(5x,a)')'Same mesh, reading data'

!     One to one HDF read

        call mpi_read_continua(n1,n2,n3,kstart,kend,1,q1)
        call mpi_read_continua(n1,n2,n3,kstart,kend,2,q2)
        call mpi_read_continua(n1,n2,n3,kstart,kend,3,q3)
#ifdef CALTEMP
        call mpi_read_continua(n1,n2,n3,kstart,kend,4,dens)
#endif
#ifdef CALSCAL
        call mpi_read_continua(n1r,n2r,n3r,kstartr,kendr,5,dsal)
#endif
        call mpi_read_continua(n1r,n2r,n3r,kstartr,kendr,6,pr)!FV

      ENDIF

      if(kstart.eq.1) then
        do j=1,n2
          do i=1,n1
            dens(i,j,0) = densbot 
          enddo
        enddo
        do j=1,n2r
          do i=1,n1r
            dsal(i,j,0) = dsalbot
          enddo
        enddo
      endif

      if(kend.eq.n3m) then
        do j=1,n2
          do i=1,n1
            dens(i,j,n3) = denstop
          enddo
        enddo
        do j=1,n2r
          do i=1,n1r
            dsal(i,j,n3r) = dsaltop
          enddo
        enddo
      endif

      if (ireset.eq.1) then                                             
        time=0.d0
      endif                                                             

      return                                                            
      end                                                               


!===================================================================                                                                                                                          
      subroutine inistruc
      use constants
      use param
      use mpih
      use mpi_param, only: kstart, kend
      use mls_param
!@cuf   use cudafor
      implicit none                                                                       


      character*100 filename,strucfilename
      character*100 ipfi,ipfip

      real(DP) :: OH,HP,OP,A,B,C,sogliola
      real(DP) :: x0,y0,z0,xH,yH,zH,xP,yP,zP

      integer :: i,contaste
!@cuf   integer :: istat

      write(ipfi,199) nwrit !FV                                                                                                                                                              
199   format(i8.8)                                                                                                                                                                          
      strucfilename = trim(folderload)//'continua_str'//trim(ipfi)//'.dat'                                                                                                                    
      write(*,*) "LEGGO", strucfilename                                                                                                                                                      
           open(919,file=strucfilename,form='unformatted')                                                                                                                                   
                 read(919)xyz(1:3,1:nvtot)                                                                                                                                                   
                 read(919)xyzv(1:3,1:nvtot)                                                                                                                                                  
                 read(919)xyzv0(1:3,1:nvtot)                                                                                                                                                 
                 read(919)xyza(1:3,1:nvtot)                                                                                                                                                  
                 read(919)xyza0(1:3,1:nvtot)                                                                                                                                                 
                 read(919)fxyz(1:3,1:nvtot)                                                 
                 read(919)HR
                 read(919)xyz_3d(1:3,1:nvtot_3d)                                             
                 read(919)xyzv_3d(1:3,1:nvtot_3d)                                            
                 read(919)xyzv0_3d(1:3,1:nvtot_3d)                                           
                 read(919)xyza_3d(1:3,1:nvtot_3d)                                            
                 read(919)xyza0_3d(1:3,1:nvtot_3d)                                           
                 read(919)fxyz_3d(1:3,1:nvtot_3d)                                            
                 read(919)astressEFcell_3d(1:nctot_3d)                                                 
                 read(919)EFtstart_3d(1:nvtot_3d)
#ifdef ELECTRO
                 read(919)XEF_1d(1:21,1:nvtotEF_1d)
                 read(919)XEF_2d(1:20,1:nftotEF_2d)
#ifdef TP06
                 read(919)XEF_3d(1:22,1:nctot_3d)
#endif
#ifdef MINIMAL_MODEL
                 read(919)XEF_3d(1:4,1:nctot_3d)
#endif
                 read(919)IstimEF_3dS1(1:nctot_3d)
                 read(919)astressEFcell_3d(1:nctot_3d)
                 read(919)Purk_tstim(1:nPurkVentr)
                 read(919)Bund_tstim(1:nBundAtri)
#endif

           close(919)                                                                        

#ifdef ELECTRO
      ! open(unit=15,file='meshes/myotag.txt',status='old')
      ! do i=1,nctot_3d
      !    read(15,*) myotag(i)
      ! enddo
      ! close(15) 
      
#endif
           
      HR=HR_S1S2     
      periodDim = (60.d0/HR)
      period = (60.d0/HR)/TSTAR
      write(*,*) "restart: HR = ", HR, &
                 " periodDim = ",periodDim, &
                 " period = ",period
      timeBeat= modulo(time,period)
      nBeat= floor(time/period)

#ifdef ELECTRO 

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=1,nvtotEF_1d 
         potEFnode_1d(i) = XEF_1d(1,i)  
      enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=1,nftotEF_2d 
         potEFface_2d(i) = XEF_2d(1,i)  
      enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
      do i=1,nctot_3d
#ifdef MIMIMAL_MODEL
         potEFcell_3d(i) = 2.7D0*XEF_3d(1,i)-83.0D0
!#else
!         potEFcell_3d(i) = XEF_3d(1,i)
#endif
      enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

   
#ifdef BIDOMAIN

      call boxEFbido(checkEFbido(1:nctot_3d,1),potEFcell_3d(1:nctot_3d),0) !check vero sistemalin
      call boxEFbido(checkEFbido(1:nctot_3d,2),potextEFcell_3d(1:nctot_3d),1)  !check vero sistemalin

      countjgmres=0
      do j=1,maxrestEFbido
          call SolveVextGmres(errorEFbido)
          countjgmres=countjgmres+1
!          write(*,*) 'GMRES', ntime, j, errorEFbido 
          if (errorEFbido.LT.tolEFbido) go to 667
       enddo
667 continue 

      call boxEFbido(checkEFbido(1:nctot_3d,3),potextEFcell_3d(1:nctot_3d),1) !check vero sistemalin

       scarto1gmres=0.0!check vero sistemalin
       scarto2gmres=0.0!check vero sistemalin
#ifdef USE_CUDA 
       !$cuf kernel do (1)
#endif
       do i=1,nctot_3d
          scarto1gmres = scarto1gmres + ( checkEFbido(i,1) - checkEFbido(i,2) )**2
          scarto2gmres = scarto2gmres + ( checkEFbido(i,1) - checkEFbido(i,3) )**2
       enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

#endif


#endif


      return
      end

!===================================================================                                                                                                                          
      subroutine iniecg
      use constants
      use param
      use mpih
      use mpi_param, only: kstart, kend
      use mls_param
!@cuf   use cudafor
      implicit none                                                                       


      character*100 filename,strucfilename
      character*100 ipfi,ipfip

      real(DP) :: OH,HP,OP,A,B,C,sogliola
      real(DP) :: x0,y0,z0,xH,yH,zH,xP,yP,zP

      integer :: i,contaste,inp
!@cuf   integer :: istat

      write(ipfi,199) nwrit !FV                                                                                                                                                  
199   format(i8.8)                                                                                                                                                               
      strucfilename = trim(folderload)//'continua_ecg'//trim(ipfi)//'.dat'                                                                                                       
      write(*,*) "LEGGO", strucfilename                                                                                                                                          
           open(919,file=strucfilename,form='unformatted')                                                                                                                       
                 read(919)time
                 read(919)xyz(1:3,1:nvtot)                                                                                                                                       
                 read(919)xyz_3d(1:3,1:nvtot_3d)                                                        
                 read(919)potEFnode_1d(1:nvtotEF_1d)
                 read(919)potEFface_2d(1:nftotEF_2d)
                 read(919)potEFcell_3d(1:nctot_3d)
           close(919)                                                                        
           nBeat= floor(time/period)

#ifdef ELECTRO 

#ifndef ELEGEO0
      !ALTRE UTILIZZATE PER EF
      call calculate_distance(dist_3d, 1, nvtot_3d, 1, netot_3d, &
               xyz_3d, vert_of_edge_3d)

      call calculate_volume_cells(1,nvtot_3d,1,nctot_3d,xyz_3d(:,1:nvtot_3d),  &
               vert_of_cell_3d(:,1:nctot_3d),vol_3d(1:nctot_3d))

      
      call calculate_bar_cells(1,nvtot_3d,1,nctot_3d,xyz_3d(:,1:nvtot_3d),  &
           vert_of_cell_3d(:,1:nctot_3d),cell_bar(:,1:nctot_3d))

      do inp=1,1
         call calculate_area(Surface_3d(inp),vstart_3d(inp),vend_3d(inp),fstart_3d(inp),fend_3d(inp),         &
                    xyz_3d(:,vstart_3d(inp):vend_3d(inp)),vert_of_face_3d(:,fstart_3d(inp):fend_3d(inp)),    &
                           sur_3d(fstart_3d(inp):fend_3d(inp)))
      end do

      call calculate_ginterp_face_cells(1,nvtot_3d,1,nftot_3d,1,nctot_3d,         &
           xyz_3d(:,1:nvtot_3d),vert_of_face_3d(:,1:nftot_3d),    &
           cell_of_face_3d(:,1:nftot_3d),cell_bar(:,1:nctot_3d), &
           versCFface_3d(:,1:nftot_3d),distCFface_3d(1:nftot_3d),g1interpface_3d(1:nftot_3d))

      call calculate_normals_face_cells(1,nvtot_3d,1,nftot_3d,1,nctot_3d,xyz_3d(:,1:nftot_3d),  &
           vert_of_face_3d(:,1:nftot_3d),face_of_cell_3d(:,1:nctot_3d),   &
           cell_bar(:,1:nctot_3d),normalfaceofcells_3d(:,:,1:nctot_3d))
!END ALTRE UTILIZZATE PER EF
#endif

      


#endif





      return
      end

