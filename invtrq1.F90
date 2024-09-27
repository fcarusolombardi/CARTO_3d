
!************************************************************************
!                       SUBROUTINE INVTR1
!   This subroutine performs the computation of Q~~ for the q1 momentum 
!   equation (radial direction) by a factored implicit scheme.
!   For details see the introduction of INVTR1
!   
!
      subroutine invtrq1
      use param
      use outflow_vars
      use ibm_param
      use local_arrays, only: pr,rhs,ru1,q1,dq,forclo
      use mpi_param, only: kstart,kend,kstartp,kendp
!@cuf use cudafor

      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,im,ip
      real(DP)    :: udx1,amm,acc,app
      real(DP)    :: dcq1,dpx11
      real(DP)    :: d22q1,d33q1,d11q1
      real(DP)    :: alre,udx1q,udx2q
      real(DP)    :: codeas, codebs, codeae, codebe
      real(DP)    :: usaldto,q1e
      
      integer :: i,j,k,ie,je,ke,n
!@cuf integer :: istat
!
      alre=al*nu

      udx1=dx1*al
      udx1q=dx1q
      udx2q=dx2q

      codeas = dx3q/(g3rc(2)*g3rm(1))
      codebs = dx3q/(g3rc(1)*g3rm(1))
      codeae = dx3q/(g3rc(n3m+1)*g3rm(n3m))
      codebe = dx3q/(g3rc(n3m)*g3rm(n3m))
!
!  compute the rhs of the factored equation
!  everything at i,j+1/2,k+1/2
!
!    points inside the flowfield
!
#ifdef USE_CUDA
      !$cuf kernel do(3)
      do kc=kstart,kend
        do jc=1,n2m
          do ic=1,n1m
            km=kmv(kc)
            kp=kpv(kc)
            amm=am3sk(kc)
            acc=ac3sk(kc)
            app=ap3sk(kc)
            jm=jmv(jc)
            jp=jpv(jc)
            im=imv(ic)
            ip=ipv(ic)
#else
      do kc=kstart,kend
        km=kmv(kc)
        kp=kpv(kc)
        amm=am3sk(kc)
        acc=ac3sk(kc)
        app=ap3sk(kc)
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,jm,jp,ic,ip,im) &
!$OMP PRIVATE(d11q1,d22q1,d33q1,dcq1,dpx11)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
#endif

!
!   11 second derivative of q1
!
            d11q1=(q1(ip,jc,kc)       &
                 -2.d0*q1(ic,jc,kc)   &
                 +q1(im,jc,kc))*udx1q
!
!   22 second derivative of q1
!
            d22q1=(q1(ic,jp,kc)       &
                  -2.d0*q1(ic,jc,kc)  &
                 +q1(ic,jm,kc))*udx2q
!
!   33 second derivative of q1
!
            if (kc .eq. n3m) then
              d33q1=((qb1n(ic,jc)-q1(ic,jc,kc))*codeae*2._DP -  &
                    (q1(ic,jc,kc)-q1(ic,jc,km))*codebe )
            else if (kc .eq. 1) then
              d33q1=((q1(ic,jc,kp)-q1(ic,jc,kc))*codeas   -  &
                    (q1(ic,jc,kc)-qb1s(ic,jc) )*codebs*2._DP )
            else
              d33q1=q1(ic,jc,kp)*app &
                   +q1(ic,jc,kc)*acc &
                   +q1(ic,jc,km)*amm
            endif

!
!    viscid terms
!
            dcq1=d11q1+d22q1+d33q1
!
!
!   component of grad(pr) along 2 direction
!
            dpx11=(pr(ic,jc,kc)-pr(im,jc,kc))*udx1

!
            rhs(ic,jc,kc)=(ga*dq(ic,jc,kc)+ro*ru1(ic,jc,kc) &
                         +alre*dcq1-dpx11)*dt

!===========================================================
!
            ru1(ic,jc,kc)=dq(ic,jc,kc)
          enddo
        enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      enddo

!@cuf istat = cudaDeviceSynchronize() ! JDR TMP

      usaldto = 1./aldto
!     --------------------------------------------
!     Direct forcing
!     TODO: Need to handle parallel case access to k variables
 
!       if(infig.eq.1) then

! !Solid Body
! #ifdef USE_CUDA
!         !$cuf kernel do(1)
! #endif
!         do n=1,npunfx(1)
!           i=indgeo(1,n,1)
!           j=indgeo(1,n,2)
!           k=indgeo(1,n,3)
!           ie=indgeoe(1,n,1)
!           je=indgeoe(1,n,2)
!           ke=indgeoe(1,n,3)
! !          q1e=((al*dt+aldto)*q1(ie,je,ke)-al*dt*q1bo(n))*usaldto
!           q1e = q1(ie,je,ke)

!           rhs(i,j,k) = -q1(i,j,k) + q1e*distb(1,n)
!           ! q1bo(n)= q1(ie,je,ke) !JR: not used anywhere
!           forclo(i,j,k,1) = 0.d0
!         end do

! #ifdef USE_CUDA
!         !$cuf kernel do(1)
! #endif
!         do n=1,npunifx(1)
!           i=indgeoee(1,n,1)
!           j=indgeoee(1,n,2)
!           k=indgeoee(1,n,3)
!           rhs(i,j,k) = -q1(i,j,k) !+ vbifx(1,n)
!           forclo(i,j,k,1) = 0.d0
!         end do

! !Mitral open/close
! !Solid Body MI
! #ifdef USE_CUDA
!         !$cuf kernel do (1)
! #endif
!         do n=1,npunfxMI(1)
!           i=indgeoMI(1,n,1)
!           j=indgeoMI(1,n,2)
!           k=indgeoMI(1,n,3)
!           ie=indgeoeMI(1,n,1)
!           je=indgeoeMI(1,n,2)
!           ke=indgeoeMI(1,n,3)
! !         q1e=((al*dt+aldto)*q1(ie,je,ke)-al*dt*q1o(ie,je,ke))*usaldto
!           q1e=q1(ie,je,ke)

!           rhs(i,j,k) = (1.-ftMI)*rhs(i,j,k) + &
!                 ftMI*( -q1(i,j,k) + q1e*distbMI(1,n) )
!           forclo(i,j,k,1)= (1.-ftMI)
!          end do

! #ifdef USE_CUDA
!         !$cuf kernel do (1)
! #endif
!         do n=1,npunifxMI(1)
!           i=indgeoeeMI(1,n,1)
!           j=indgeoeeMI(1,n,2)
!           k=indgeoeeMI(1,n,3)
!           rhs(i,j,k) = (1.-ftMI)*rhs(i,j,k) + ftMI*(-q1(i,j,k))
!           forclo(i,j,k,1)= (1.-ftMI)
!         enddo

! !Aortra open/close
! !Solid Body AO
! #ifdef USE_CUDA
!         !$cuf kernel do (1)
! #endif
!         do n=1,npunfxAO(1)
!           i=indgeoAO(1,n,1)
!           j=indgeoAO(1,n,2)
!           k=indgeoAO(1,n,3)
!           ie=indgeoeAO(1,n,1)
!           je=indgeoeAO(1,n,2)
!           ke=indgeoeAO(1,n,3)
! !         q1e=((al*dt+aldto)*q1(ie,je,ke)-al*dt*q1o(ie,je,ke))*usaldto
!           q1e=q1(ie,je,ke)
!           rhs(i,j,k) = (1.-ftAO)*rhs(i,j,k) + &
!                 ftAO*( -q1(i,j,k) + q1e*distbAO(1,n) )
!           forclo(i,j,k,1)= (1.-ftAO)
!         end do

! #ifdef USE_CUDA
!         !$cuf kernel do (1)
! #endif
!         do n=1,npunifxAO(1)
!           i=indgeoeeAO(1,n,1)
!           j=indgeoeeAO(1,n,2)
!           k=indgeoeeAO(1,n,3)
!           rhs(i,j,k) = (1.-ftAO)*rhs(i,j,k) + ftAO*(-q1(i,j,k))
!           forclo(i,j,k,1)= (1.-ftAO)
!         enddo

!Aorta system open/close
!Solid Body AS
#ifdef USE_CUDA
        !$cuf kernel do(1)
#endif
        do n=1,npunfxAS(1)
          i=indgeoAS(1,n,1)
          j=indgeoAS(1,n,2)
          k=indgeoAS(1,n,3)
          ie=indgeoeAS(1,n,1)
          je=indgeoeAS(1,n,2)
          ke=indgeoeAS(1,n,3)
!         q1e=((al*dt+aldto)*q1(ie,je,ke)-al*dt*q1bo(n))*usaldto
          q1e = q1(ie,je,ke)


          rhs(i,j,k) = (1.-ftAS)*rhs(i,j,k) +              &
                  ftAS*(-q1(i,j,k) + q1e*distbAS(1,n))+    &
                  (1.-ftAS)*dt*ftAS2*(-q1(i,j,k))
          !q1bo(n)= q1(ie,je,ke)
          forclo(i,j,k,1) = (1.-ftAS)
        end do

#ifdef USE_CUDA
        !$cuf kernel do(1)
#endif
        do n=1,npunifxAS(1)
          i=indgeoeeAS(1,n,1)
          j=indgeoeeAS(1,n,2)
          k=indgeoeeAS(1,n,3)
          rhs(i,j,k) = (1.-ftAS)*rhs(i,j,k)+ftAS*(-q1(i,j,k)) +  &
                  (1.-ftAS)*dt*ftAS2*(-q1(i,j,k))
          forclo(i,j,k,1) = (1.-ftAS)
        end do

!@cuf istat = cudaDeviceSynchronize() ! JDR TMP

      ! end if


!     --------------------------------------------


      call solxi(beta*al*dx1q,1)

      call solxj(beta*al*dx2q,1)
      
      call solq1k(1)
     
      return
      end

