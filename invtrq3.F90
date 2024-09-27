
!************************************************************************
!                       SUBROUTINE INVTR3
!   This subroutine performs the computation of Q~~ for the q3 momentum 
!   equation (axial direction) by a factored implicit scheme.
!   Viscous terms are treated implicitly, nonlinear terms explicitly.
!   For details see the introduction of INVTR1 
!
      subroutine invtrq3
      use param
      use ibm_param
      use local_arrays, only: q3,qcap,pr,ru3,rhs,forclo
      use mpi_param, only: kstart,kend,kstartp,kendp
!@cuf use cudafor
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,ip,im
      real(DP)    :: udx3
      real(DP)    :: dq32,dq33,dcq3,dpx33,dq31
      real(DP)    :: app,acc,amm
      real(DP)    :: alre,udx1q,udx2q

      integer :: i,j,k,ie,je,ke,n
      real(DP) :: q3e
!@cuf integer :: istat
      
      alre=al*nu

      udx1q=dx1q
      udx2q=dx2q
!
!  compute the rhs of the factored equation
!  everything at i+1/2,j+1/2,k
!


#ifdef USE_CUDA
      !$cuf kernel do (3)
      do kc=kstartp,kend
        do jc=1,n2m
            do ic=1,n1m
              km=kmv(kc)
              kp=kc+1
              udx3 = al*udx3c(kc)
              amm=am3ck(kc)
              acc=ac3ck(kc)
              app=ap3ck(kc)
              jm=jmv(jc)
              jp=jpv(jc)
              im=imv(ic)
              ip=ipv(ic)
#else
      do kc=kstartp,kend
        km=kmv(kc)
        kp=kc+1
        udx3 = al*udx3c(kc)
        amm=am3ck(kc)
        acc=ac3ck(kc)
        app=ap3ck(kc)
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,jm,jp,ic,ip,im) &
!$OMP PRIVATE(dq31,dq32,dq33,dcq3,dpx33)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
            do ic=1,n1m
              im=imv(ic)
              ip=ipv(ic)
#endif
!
!   11 second derivatives of q3
!
            dq31=(q3(im,jc,kc)      &
                -2.d0*q3(ic,jc,kc)  &
                +q3(ip,jc,kc))*udx1q
!
!   22 second derivatives of q3
!
            dq32=(q3(ic,jm,kc)      &
                -2.d0*q3(ic,jc,kc)  &
                +q3(ic,jp,kc))*udx2q
!
!   33 second derivatives of q3
!
            dq33=q3(ic,jc,kp)*app   &
               +q3(ic,jc,kc)*acc    &
               +q3(ic,jc,km)*amm
!
!   viscous terms
!
            dcq3=dq32+dq33+dq31
!
!  component of grad(pr) along x3 direction
!
            dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*udx3
!=======================================================     
            rhs(ic,jc,kc)=(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc) &
                         +alre*dcq3-dpx33)*dt 
!=======================================================
!
!  updating of the non-linear terms
!
            ru3(ic,jc,kc)=qcap(ic,jc,kc)
         enddo
       enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      enddo

!@cuf istat = cudaDeviceSynchronize() ! JDR TMP
!     --------------------------------------------
!     Direct forcing
!     TODO: Need to handle parallel case access to k variables
 
      ! if(infig.eq.1) then

! !Solid Body
! #ifdef USE_CUDA
!         !$cuf kernel do(1)
! #endif
!         do n=1,npunfx(3)
!           i=indgeo(3,n,1)
!           j=indgeo(3,n,2)
!           k=indgeo(3,n,3)
!           ie=indgeoe(3,n,1)
!           je=indgeoe(3,n,2)
!           ke=indgeoe(3,n,3)
! !         q3e=((al*dt+aldto)*q3(ie,je,ke)-al*dt*q3bo(n))*usaldto
!           q3e = q3(ie,je,ke)

!           rhs(i,j,k) = -q3(i,j,k) + q3e*distb(3,n)
!           !q3bo(n)= q3(ie,je,ke) !JR: not used anywhere
!           forclo(i,j,k,3)= 0.d0
!         end do

! #ifdef USE_CUDA
!         !$cuf kernel do(1)
! #endif
!         do n=1,npunifx(3)
!           i=indgeoee(3,n,1)
!           j=indgeoee(3,n,2)
!           k=indgeoee(3,n,3)
!           rhs(i,j,k) = -q3(i,j,k) !+ vbifx(1,n)
!           forclo(i,j,k,3)= 0.d0
!         end do

! !Solid Body MI
! !Mitral open/close
! #ifdef USE_CUDA
!         !$cuf kernel do(1)
! #endif
!         do n=1,npunfxMI(3)
!           i=indgeoMI(3,n,1)
!           j=indgeoMI(3,n,2)
!           k=indgeoMI(3,n,3)
!           ie=indgeoeMI(3,n,1)
!           je=indgeoeMI(3,n,2)
!           ke=indgeoeMI(3,n,3)
!           !q3e=((al*dt+aldto)*q3(ie,je,ke)-al*dt*q3o(ie,je,ke))*usaldto
!           q3e=q3(ie,je,ke)
!           rhs(i,j,k) = (1.-ftMI)*rhs(i,j,k) + &
!                 ftMI*( -q3(i,j,k) + q3e*distbMI(3,n) )
!           forclo(i,j,k,3)= (1.-ftMI)
!         end do

! #ifdef USE_CUDA
!         !$cuf kernel do(1)
! #endif
!         do n=1,npunifxMI(3)
!           i=indgeoeeMI(3,n,1)
!           j=indgeoeeMI(3,n,2)
!           k=indgeoeeMI(3,n,3)
!           rhs(i,j,k) = (1.-ftMI)*rhs(i,j,k) + ftMI*(-q3(i,j,k))
!           forclo(i,j,k,3)= (1.-ftMI)
!         enddo

! !Aorta open/close
! !Solid Body AO
! #ifdef USE_CUDA
!         !$cuf kernel do(1)
! #endif
!         do n=1,npunfxAO(3)
!           i=indgeoAO(3,n,1)
!           j=indgeoAO(3,n,2)
!           k=indgeoAO(3,n,3)
!           ie=indgeoeAO(3,n,1)
!           je=indgeoeAO(3,n,2)
!           ke=indgeoeAO(3,n,3)
!           !q3e=((al*dt+aldto)*q3(ie,je,ke)-al*dt*q3o(ie,je,ke))*usaldto
!           q3e=q3(ie,je,ke)
!           rhs(i,j,k) = (1.-ftAO)*rhs(i,j,k) + &
!                 ftAO*( -q3(i,j,k) + q3e*distbAO(3,n) )
!           forclo(i,j,k,3)= (1.-ftAO)
!         end do

! #ifdef USE_CUDA
!         !$cuf kernel do(1)
! #endif
!         do n=1,npunifxAO(3)
!           i=indgeoeeAO(3,n,1)
!           j=indgeoeeAO(3,n,2)
!           k=indgeoeeAO(3,n,3)
!           rhs(i,j,k) = (1.-ftAO)*rhs(i,j,k) + ftAO*(-q3(i,j,k))
!           forclo(i,j,k,3)= (1.-ftAO)
!         enddo

!Aorta system open/close
!Solid Body AS
#ifdef USE_CUDA
        !$cuf kernel do(1)
#endif
        do n=1,npunfxAS(3)
          i=indgeoAS(3,n,1)
          j=indgeoAS(3,n,2)
          k=indgeoAS(3,n,3)
          ie=indgeoeAS(3,n,1)
          je=indgeoeAS(3,n,2)
          ke=indgeoeAS(3,n,3)
!         q3e=((al*dt+aldto)*q3(ie,je,ke)-al*dt*q3bo(n))*usaldto
          q3e = q3(ie,je,ke)


          rhs(i,j,k) = (1.-ftAS)*rhs(i,j,k) +             &
                  ftAS*(-q3(i,j,k) + q3e*distbAS(3,n))+   &
                  (1.-ftAS)*dt*ftAS2*(-q3(i,j,k)) + dt*ftAS3
          !q1bo(n)= q1(ie,je,ke)
          forclo(i,j,k,3)= (1.-ftAS)
        end do

#ifdef USE_CUDA
        !$cuf kernel do(1)
#endif
        do n=1,npunifxAS(3)
          i=indgeoeeAS(3,n,1)
          j=indgeoeeAS(3,n,2)
          k=indgeoeeAS(3,n,3)
          rhs(i,j,k) = (1.-ftAS)*rhs(i,j,k)+ftAS*(-q3(i,j,k)) +  &
                  (1.-ftAS)*dt*ftAS2*(-q3(i,j,k)) + dt*ftAS3
          forclo(i,j,k,3)= (1.-ftAS)
        end do

!@cuf istat = cudaDeviceSynchronize() ! JDR TMP

      ! end if


!     --------------------------------------------


      call solxi(beta*al*dx1q,3)
      call solxj(beta*al*dx2q,3)
      call solq3k(3)

      return
      end


