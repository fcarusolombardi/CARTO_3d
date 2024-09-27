! Trilinear interpolation to new grid for MPI
! continua.dat files.
! assumes alx3 = 1
!=====================================================
!  salinity 
      subroutine interp_dsal(kstartor,kendor)
      use param
      use mpih
      use local_arrays, only: dsal
      use mpi_param, only: kstartr,kendr
      use init_itp
      implicit none
      integer,intent(in) :: kstartor,kendor

      real(DP),allocatable,dimension(:) :: xold,yold,zold
      real(DP),allocatable,dimension(:) :: xnew,ynew,znew
      integer,allocatable,dimension(:) :: idcx,idcy,idcz

      integer i,j,k,ic,jc,kc
      real(DP) :: zla,xlc,ylc,zlc
      real(DP) :: v111,v112,v121,v122
      real(DP) :: v211,v212,v221,v222
      real(DP) :: v11,v12,v21,v22,v1,v2

      allocate(xold(0:n1or),yold(0:n2or),zold(0:n3or))
      allocate(xnew(1:n1r),ynew(1:n2r),znew(1:n3r))
      allocate(idcx(1:n1r),idcy(1:n2r),idcz(1:n3r))

!============================================================
!    setup old new grids

      xold(1:n1or) = xmrold(1:n1or)
      xold(0) = 2.d0*xold(1) - xold(2)
      yold(1:n2or) = ymrold(1:n2or)
      yold(0) = 2.d0*yold(1) - yold(2)
      zold(1:n3or) = zmrold(1:n3or)
      zold(0) = 0.d0

      xnew = xmr
      ynew = ymr
      znew = zmr
      idcx=iisac1
      idcy=iisac2
      idcz=iisac3
      if(kstartor.eq.1) then
        do i=1,n1or
          do j=1,n2or
            varold(i,j,0)=dsalbot
          enddo
        enddo
      endif
      if(kendor.eq.n3omr) then
        do i=1,n1or
          do j=1,n2or
            varold(i,j,n3or)=dsaltop
          enddo
        enddo
      endif

!==========================================================
!    INTERP
      do k=kstartr,kendr
       kc = idcz(k)
       zla = 1.d0/(zold(kc+1)-zold(kc))
       zlc = (znew(k)-zold(kc))*zla
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(i,j,ic,jc) &
!$OMP PRIVATE(v111,v112,v121,v122,v211,v212,v221,v222) &
!$OMP PRIVATE(xlc,ylc,v11,v12,v21,v22,v1,v2)
       do j=1,n2mr
        do i=1,n1mr
!  setup local coefficients
          ic = idcx(i)
          xlc = (xnew(i)-xold(ic))/(xold(ic+1)-xold(ic)) 

          jc = idcy(j) 
          ylc = (ynew(j)-yold(jc))/(yold(jc+1)-yold(jc))

!  setup data at eight corners
          if(ic.eq.0) ic = n1omr
          if(jc.eq.0) jc = n2omr

          v121 = varold(ic,jc+1,kc)
          v111 = varold(ic,jc,kc)
          v221 = varold(ic+1,jc+1,kc)
          v211 = varold(ic+1,jc,kc)

          v122 = varold(ic,jc+1,kc+1)
          v112 = varold(ic,jc,kc+1)
          v222 = varold(ic+1,jc+1,kc+1)
          v212 = varold(ic+1,jc,kc+1)

! trilinear interpolation
          v11 = v111*(1.d0-zlc)+v112*zlc
          v12 = v121*(1.d0-zlc)+v122*zlc
          v21 = v211*(1.d0-zlc)+v212*zlc
          v22 = v221*(1.d0-zlc)+v222*zlc 

          v1 = v11*(1.d0-ylc)+v12*ylc
          v2 = v21*(1.d0-ylc)+v22*ylc

          dsal(i,j,k) = v1*(1.d0-xlc)+v2*xlc

        enddo
       enddo
!$OMP END PARALLEL DO
      enddo

!================================================================
!  periodic boundary condition
      do kc=kstartr-1,kendr+1
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc)
        do jc=1,n2mr
          dsal(n1r,jc,kc) = dsal(1,jc,kc)
        enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ic)
        do ic=1,n1r
          dsal(ic,n2r,kc) = dsal(ic,1,kc)
        enddo
!$OMP END PARALLEL DO
      enddo

!================================================================
!  cleanup memory

      deallocate(xold,yold,zold)
      deallocate(xnew,ynew,znew)
      deallocate(idcx,idcy,idcz)

      return
      end subroutine interp_dsal
