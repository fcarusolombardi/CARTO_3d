
!***********************************************************************
!                                                                      *
!                       CONDITIONS                                     *
!                                                                      *
!***********************************************************************
      subroutine dsalbo
      use param
      implicit none
      integer :: j,i
      
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(i,j)
      do j=1,n2r
        do i=1,n1r
          dsalbn(i,j) = dsaltop
          dsalbs(i,j) = dsalbot
        enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end
