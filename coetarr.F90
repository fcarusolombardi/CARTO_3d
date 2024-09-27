!************************************************************************
!                                                                       *
! ****************************** subrout coetar  ********************** *
!                                                                       *
!    this subroutine calculates the coefficients for the                *
!    integration in the radial direction with non-uniform coor. trasf.  *
!                                                                       *
!************************************************************************
      subroutine coetarr
      use param
      use mpih
      implicit none
      integer :: km,kc,kp,k
      real(DP):: a33,a33m,a33p

!  **coefficients for q1, q2 ap3sk,am3sk,ac3sk
!   ap3ssk,am3ssk,ac3ssk, psc   for x3 differentation
!  s means staggered that is at k+1/2 location
!
      do kc=2,n3mr-1
        kp=kc+1
        km=kc-1
        a33 = dx3qr/g3rmr(kc)
        a33p= +a33/g3rcr(kp)
        a33m= +a33/g3rcr(kc)
        ap3sskr(kc)=a33p
        am3sskr(kc)=a33m
        ac3sskr(kc)=-(a33p+a33m)
      enddo
!    
!    lower wall  bound. conditions 
!    differemtiation of sec. derivative at 1+1/2
!    
      kc=1
      kp=kc+1
      a33 = dx3qr/g3rmr(kc)
      a33p= +a33/g3rcr(kp)
      a33m= +a33/g3rcr(kc)
      ap3sskr(kc)=a33p
      am3sskr(kc)=0.d0
      ac3sskr(kc)=-(a33p+a33m*2.d0)
!    
!    upper wall  bound. conditions
!    differentiation of sec. derivative at n3-1/2
!    
      kc=n3mr
      kp=kc+1
      a33=dx3qr/g3rmr(kc)
      a33p= +a33/g3rcr(kp)
      a33m= +a33/g3rcr(kc)
      ap3sskr(kc)=0.d0
      am3sskr(kc)=a33m
      ac3sskr(kc)=-(a33m+a33p*2.d0)

#ifdef DEBUG 
      if(myid.eq.0)then
        open(401,file='fact/csskr.out')
        do k=1, n3mr
          write(401,'(i5,3f25.10)')k, ap3sskr(k), am3sskr(k), ac3sskr(k)
        enddo
        close(401)
      endif
#endif

      return
      end

