      subroutine ntrvb_line(am,ac,ap,f,n)
      use mpih
      implicit none
      integer :: n
      real(DP), dimension(n) :: am, ac, ap, f
      integer :: i,nm,ii
!
!  ******** reduction of trid. matrix to an upper rigth matrix
!
      do i=2,n
        ac(i)=ac(i)*ac(i-1)-ap(i-1)*am(i)
        ap(i)=ap(i)*ac(i-1)
        f(i)=f(i)*ac(i-1)-f(i-1)*am(i)
      enddo

!  ******** calculation of the unknown by backward elimination
!

      f(n)=f(n)/ac(n)
      nm=n-1
      
      do ii=1,nm
        i=n-ii
        f(i)=(f(i)-ap(i)*f(i+1))/ac(i)
      enddo
   
      return
      end 
