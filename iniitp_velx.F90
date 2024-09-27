! Trilinear interpolation to new grid for MPI
! continua.dat files.
! assumes alx3 = 1
!=====================================================
!  x velocity 
      subroutine interp_velx(kstarto,kendo)
      use param
      use mpih
      use local_arrays, only: q1
      use mpi_param, only: kstart,kend
      use init_itp
      implicit none
      integer,intent(in) :: kstarto,kendo

      real(DP),allocatable,dimension(:) :: xold,yold,zold
      real(DP),allocatable,dimension(:) :: xnew,ynew,znew
      integer,allocatable,dimension(:) :: idcx,idcy,idcz

      integer i,j,k,ic,jc,kc
      real(DP) :: zla,xlc,ylc,zlc
      real(DP) :: v111,v112,v121,v122
      real(DP) :: v211,v212,v221,v222
      real(DP) :: v11,v12,v21,v22,v1,v2

      allocate(xold(0:n1o),yold(0:n2o),zold(0:n3o))
      allocate(xnew(1:n1),ynew(1:n2),znew(1:n3))
      allocate(idcx(1:n1),idcy(1:n2),idcz(1:n3))

!============================================================
!    setup old new grids

      xold(1:n1o) = xcold(1:n1o)
      xold(0) = 2.d0*xold(1) - xold(2)
      yold(1:n2o) = ymold(1:n2o)
      yold(0) = 2.d0*yold(1) - yold(2)
      zold(1:n3o) = zmold(1:n3o)
      zold(0) = 0.d0

      xnew = xc
      ynew = ym
      znew = zm
      idcx=iiq1c1
      idcy=iiq1c2
      idcz=iiq1c3
      if(kstarto.le.1) then
        do i=1,n1o
          do j=1,n2o
            varold(i,j,0)=dble(ubcbot)*varold(i,j,1)
          enddo
        enddo
      endif
      if(kendo.ge.n3o-1) then
        do i=1,n1o
          do j=1,n2o
            varold(i,j,n3o)=dble(ubctop)*varold(i,j,n3om)
          enddo
        enddo
      endif

!==========================================================
!    INTERP
      do k=kstart,kend
       kc = idcz(k)
       zla = 1.d0/(zold(kc+1)-zold(kc))
       zlc = (znew(k)-zold(kc))*zla
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(i,j,ic,jc) &
!$OMP PRIVATE(v111,v112,v121,v122,v211,v212,v221,v222) &
!$OMP PRIVATE(xlc,ylc,v11,v12,v21,v22,v1,v2)
       do j=1,n2m
        do i=1,n1m
!  setup local coefficients
          ic = idcx(i)
          xlc = (xnew(i)-xold(ic))/(xold(ic+1)-xold(ic)) 

          jc = idcy(j) 
          ylc = (ynew(j)-yold(jc))/(yold(jc+1)-yold(jc))

!  setup data at eight corners
          if(ic.eq.0) ic = n1o-1
          if(jc.eq.0) jc = n2o-1

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

          q1(i,j,k) = v1*(1.d0-xlc)+v2*xlc

        enddo
       enddo
!$OMP END PARALLEL DO
      enddo

!================================================================
!  cleanup memory

      deallocate(xold,yold,zold)
      deallocate(xnew,ynew,znew)
      deallocate(idcx,idcy,idcz)

      return
      end subroutine interp_velx
