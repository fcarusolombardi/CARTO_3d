!***********************************************************************
!                                                                      *
!                       INITIAL CONDITION                              *
!                                                                      *
!***********************************************************************
      subroutine inqpr
      use param
      use local_arrays, only: q2,q3,dens,q1,dsal
      use outflow_vars
      use mpi_param, only: kstart,kend,kstartr,kendr
      use mpih
      implicit none
      integer :: j,k,i
      real(DP) :: eps, varptb

! prescribed initial distributions of T and S

      call initprofs

      do k=kstart,kend
       do j=1,n2m
        do i=1,n1m
          dens(i,j,k) = Tiniprof(k)
        enddo
       enddo
      end do

      do k=kstartr,kendr
       do j=1,n2mr
        do i=1,n1mr
          dsal(i,j,k) = Siniprof(k)
        enddo
       enddo
      end do

! random noise

      eps = 0.001d0              ! amplitude of pertubation
      eps = 0.0
      call random_seed()

      do k=kstartr,kendr
        do j=1,n2mr
          do i=1,n1mr
            call random_number(varptb)
            dens(i,j,k) = dens(i,j,k) + eps*(2.d0*varptb-1.d0)
            dsal(i,j,k) = dsal(i,j,k) + eps*(2.d0*varptb-1.d0)
          enddo
        enddo
      end do

!  velocity field

      do k=kstart-lvlhalo,kend+lvlhalo
        do j=1,n2
          do i=1,n1
            q1(i,j,k)=0.d0
            q2(i,j,k)=0.0d0
            q3(i,j,k)=0.d0
          enddo
        enddo
      enddo

      qb1s = 0.0 ; qb1n = 0.0
      qb2s = 0.0 ; qb2n = 0.0
      qb3s = 0.0 ; qb3n = 0.0
      
      dqb1s = 0.0 ; dqb1n = 0.0
      dqb2s = 0.0 ; dqb2n = 0.0
      dqb3s = 0.0 ; dqb3n = 0.0

      dq1x2o=0.0 ; dq2x2o=0.0 ; dq3x2o=0.0 

      return
      end


!==========================================================
!    initial distributions of T and S
!==========================================================
      subroutine initprofs
      use param
      implicit none
      integer :: k
      real(DP) :: densmean, dsalmean
      real(DP) :: dvar, dhl

! prescribed initial distributions of T and S

      if(tag_ini.eq.0)then
        densmean = (densbot+denstop)/2.d0
        do k=1,n3
          Tiniprof(k) = densmean
        end do
        dsalmean = (dsalbot+dsaltop)/2.d0
        do k=1,n3r
          Siniprof(k) = dsalmean
        end do
      endif

      if(tag_ini.eq.1)then
        do k=1,n3
          Tiniprof(k) = densbot + zm(k) * (denstop-densbot)/alx3
        end do
        do k=1,n3r
          Siniprof(k) = dsalbot + zmr(k) * (dsaltop-dsalbot)/alx3
        end do
      endif

      if(tag_ini.eq.2)then
        do k=1,n3
          Tiniprof(k) = densbot + zm(k) * (denstop-densbot)/alx3
        end do
        dsalmean = (dsalbot+dsaltop)/2.d0
        do k=1,n3r
          Siniprof(k) = dsalmean
        end do
      endif

      if(tag_ini.eq.3)then
        dhl = alx3/dble(Nitf+1)
        dvar = (denstop-densbot)/dble(Nitf+2)
        do k=1,n3
          Tiniprof(k) = densbot + dble(floor(zm(k)/dhl)+1)*dvar
        end do
        dvar = (dsaltop-dsalbot)/dble(Nitf+2)
        do k=1,n3r
          Siniprof(k) = dsalbot + dble(floor(zmr(k)/dhl)+1)*dvar
        end do
      endif

      return
      end

