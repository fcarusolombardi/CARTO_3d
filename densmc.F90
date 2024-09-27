!***********************************************************************
!***********************************************************************
      subroutine densmc
      use param
      use local_arrays, only: dens, dsal
      use mpih
      implicit none
      integer :: j,i
      real(DP) :: fcder,fcdern
      real(DP) :: surf,fac2
      real(DP) :: my_anutupp, my_anusupp
      real(DP) :: adens0,adens1,adens2, adsal0,adsal1,adsal2
      real(DP) :: del1,del2,del1n,del2n,udel1q,udel2q,udel1qn,udel2qn

!
!     COMPUTATION OF THE NUSSELT NUMBER AT THE 
!     LOWER AND UPPER WALLS
!

      anutlow = 0.d0
      anutupp = 0.d0
      anuslow = 0.d0
      anusupp = 0.d0

      my_anutupp = 0.d0
      my_anusupp = 0.d0

      if(myid.eq.0) then

!  heat flux at lower plate
        del1 = zm(1)-zc(1)
        del2  = zm(2)-zc(1)
        udel1q = 1.d0/del1**2
        udel2q = 1.d0/del2**2
        fcder = 1.d0/(1.d0/del1 - 1.d0/del2)
        surf = rext1*rext2
        fac2 = 1.d0/(dx1*dx2)
        adens0 = 0.d0
        adens1 = 0.d0
        adens2 = 0.d0
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(j,i) &
!$OMP REDUCTION(+:adens0,adens1,adens2)
        do j=1,n2m
          do i=1,n1m
            adens0 = adens0 + densbs(i,j)*fac2 
            adens1 = adens1 + dens(i,j,1)*fac2
            adens2 = adens2 + dens(i,j,2)*fac2
          enddo
        enddo
!$OMP  END PARALLEL DO
        adens0 = adens0 / surf
        adens1 = adens1 / surf
        adens2 = adens2 / surf

        anutlow = (adens1-adens0)/del1

! salinity flux at lower plate

        del1 = zmr(1)-zcr(1)
        del2  = zmr(2)-zcr(1)
        udel1q = 1.d0/del1**2
        udel2q = 1.d0/del2**2
        fcder = 1.d0/(1.d0/del1 - 1.d0/del2)
        surf = rext1*rext2
        fac2 = 1.d0/(dx1r*dx2r)
        adsal0 = 0.d0
        adsal1 = 0.d0
        adsal2 = 0.d0 
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(j,i) &
!$OMP REDUCTION(+:adsal0,adsal1,adsal2)
        do j=1,n2mr
          do i=1,n1mr
            adsal0 = adsal0 + dsalbs(i,j)*fac2
            adsal1 = adsal1 + dsal(i,j,1)*fac2
            adsal2 = adsal2 + dsal(i,j,2)*fac2
          enddo
        end do
!$OMP  END PARALLEL DO
        adsal0 = adsal0 / surf
        adsal1 = adsal1 / surf
        adsal2 = adsal2 / surf
        anuslow = (adsal1-adsal0)/del1
      endif

      if(myid.eq.numtasks-1) then

! heat flux at upper plate
        del1n = zc(n3)-zm(n3m)
        del2n  = zc(n3)-zm(n3m-1)
        udel1qn = 1.d0/del1n**2
        udel2qn = 1.d0/del2n**2
        fcdern = 1.d0/(1.d0/del2n - 1.d0/del1n)
        surf = rext1*rext2
        fac2 = 1.d0/(dx1*dx2)
        adens0 = 0.d0
        adens1 = 0.d0
        adens2 = 0.d0
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(j,i) &
!$OMP REDUCTION(+:adens0,adens1,adens2)
        do j=1,n2m
          do i=1,n1m
            adens0 = adens0 + densbn(i,j)*fac2
            adens1 = adens1 + dens(i,j,n3m)*fac2
            adens2 = adens2 + dens(i,j,n3m-1)*fac2
          enddo
        enddo
!$OMP  END PARALLEL DO
        adens0 = adens0 / surf
        adens1 = adens1 / surf
        adens2 = adens2 / surf
        my_anutupp = (adens0-adens1)/del1n

! salinity flux at upper plate

        del1n = zcr(n3r) - zmr(n3mr)
        del2n = zcr(n3r) - zmr(n3mr-1)
        udel1qn = 1.d0/del1n**2
        udel2qn = 1.d0/del2n**2
        fcdern = 1.d0/(1.d0/del2n - 1.d0/del1n)
        surf = rext1*rext2
        fac2 = 1.d0/(dx1r*dx2r)
        adsal0 = 0.d0
        adsal1 = 0.d0
        adsal2 = 0.d0
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(j,i) &
!$OMP REDUCTION(+:adsal0,adsal1,adsal2)
        do j=1,n2mr
          do i=1,n1mr
            adsal0 = adsal0 + dsalbn(i,j)*fac2
            adsal1 = adsal1 + dsal(i,j,n3mr)*fac2
            adsal2 = adsal2 + dsal(i,j,n3mr-1)*fac2
          enddo
        end do
!$OMP  END PARALLEL DO
        adsal0 = adsal0 / surf
        adsal1 = adsal1 / surf
        adsal2 = adsal2 / surf
        my_anusupp = (adsal0-adsal1)/del1n
      endif

      call MPI_REDUCE(my_anutupp,anutupp,1,MDP, &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_anusupp,anusupp,1,MDP, &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if(myid.eq.0)then
        write(97,546) time, dabs(anutlow), dabs(anutupp),  &
                           dabs(anuslow), dabs(anusupp)
      endif
546   format(1x,f10.4,4(1x,ES20.8))

      return         
      end               
