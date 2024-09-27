
!***********************************************************************
!                                                                      *
!                       CONDIZIONI AL CONTORNO                         *
!                                                                      *
!***********************************************************************
      subroutine densbo
      use param
      implicit none
      integer :: j,i

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(i,j)
      do j=1,n2
        do i=1,n1
          densbn(i,j)=denstop
          densbs(i,j)=densbot
        enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end
