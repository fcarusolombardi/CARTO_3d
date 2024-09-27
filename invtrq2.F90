
!***********************************************************************
!                       SUBROUTINE INVTR2
!   This subroutine performs the computation of Q~~ for the q2 momentum 
!   equation (radial direction) by a factored implicit scheme.
!   For details see the introduction of INVTR1
!   
!
      subroutine invtrq2
      use param
      use ibm_param
      use outflow_vars
      use local_arrays, only: q2,pr,rhs,dph,ru2,forclo
      use mpi_param, only: kstart,kend,kstartp,kendp
!@cuf use cudafor
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,im,ip
      real(DP)    :: udx2,amm,app,acc
      real(DP)    :: dcq2,dpx22
      real(DP)    :: d22q2,d33q2,d11q2
      real(DP)    :: alre,udx1q,udx2q

      real(DP)    :: codeas,codebs,codeae,codebe
      real(DP)    :: usaldto,q2e


      integer :: i,j,k,ie,je,ke,n
!@cuf integer :: istat

      

      alre=al*nu

      udx2=dx2*al
      udx1q=dx1q
      udx2q=dx2q

      codeas = dx3q/(g3rc(2)*g3rm(1))
      codebs = dx3q/(g3rc(1)*g3rm(1))
      codeae = dx3q/(g3rc(n3m+1)*g3rm(n3m))
      codebe = dx3q/(g3rc(n3m)*g3rm(n3m))
!
!  compute the rhs of the factored equation
!  everything at i+1/2,j,k+1/2
!
!    points inside the flowfield
!
#ifdef USE_CUDA
      !$cuf kernel do (3)
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
!$OMP PRIVATE(d11q2,d22q2,d33q2,dcq2,dpx22)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
#endif
!
!
!   11 second derivative of q2
!
            d11q2=(q2(ip,jc,kc)        &
                 -2.d0*q2(ic,jc,kc)    &
                 +q2(im,jc,kc))*udx1q

!
!   22 second derivative of q2
!
            d22q2=(q2(ic,jp,kc)       &
                 -2.d0*q2(ic,jc,kc)   &
                 +q2(ic,jm,kc))*udx2q
!
!   33 second derivative of q2
!
            if (kc .eq. n3m) then
              d33q2=((qb2n(ic,jc)-q2(ic,jc,kc))*codeae*2. - &
                    (q2(ic,jc,kc)-q2(ic,jc,km))*codebe )
            elseif (kc .eq. 1) then
              d33q2=((q2(ic,jc,kp)-q2(ic,jc,kc))*codeas -  &
                    (q2(ic,jc,kc)-qb2s(ic,jc) )*codebs*2. )
            else
              d33q2=q2(ic,jc,kp)*app  &
                  +q2(ic,jc,kc)*acc   &
                  +q2(ic,jc,km)*amm
            endif

!
!    viscid terms
!
            dcq2=d22q2+d33q2+d11q2
!
!
!   component of grad(pr) along 2 direction
!
            dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*udx2


            rhs(ic,jc,kc)=(ga*dph(ic,jc,kc)+ro*ru2(ic,jc,kc) &
                         +alre*dcq2-dpx22)*dt

!===========================================================
!
            ru2(ic,jc,kc)=dph(ic,jc,kc)
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
 
      ! if(infig.eq.1) then

! !Solid Body
! #ifdef USE_CUDA
!         !$cuf kernel do(1)
! #endif
!         do n=1,npunfx(2)
!           i=indgeo(2,n,1)
!           j=indgeo(2,n,2)
!           k=indgeo(2,n,3)
!           ie=indgeoe(2,n,1)
!           je=indgeoe(2,n,2)
!           ke=indgeoe(2,n,3)
! !         q2e=((al*dt+aldto)*q2(ie,je,ke)-al*dt*q2bo(n))*usaldto
!           q2e = q2(ie,je,ke)

!           rhs(i,j,k) = -q2(i,j,k) + q2e*distb(2,n)
!           !q2bo(n)= q2(ie,je,ke) !JR: not used anywhere
!           forclo(i,j,k,2)= 0.d0
!         end do

! #ifdef USE_CUDA
!         !$cuf kernel do(1)
! #endif
!         do n=1,npunifx(2)
!           i=indgeoee(2,n,1)
!           j=indgeoee(2,n,2)
!           k=indgeoee(2,n,3)
!           rhs(i,j,k) = -q2(i,j,k) !+ vbifx(1,n)
!           forclo(i,j,k,2)= 0.d0
!         end do

! !Mitral open/close
! !Solid Body MI
! #ifdef USE_CUDA
!         !$cuf kernel do(1)
! #endif
!          do n=1,npunfxMI(2)
!            i=indgeoMI(2,n,1)
!            j=indgeoMI(2,n,2)
!            k=indgeoMI(2,n,3)
!            ie=indgeoeMI(2,n,1)
!            je=indgeoeMI(2,n,2)
!            ke=indgeoeMI(2,n,3)
!            !q2e=((al*dt+aldto)*q2(ie,je,ke)-al*dt*q2o(ie,je,ke))*usaldto
!            q2e=q2(ie,je,ke)
!            rhs(i,j,k) = (1.-ftMI)*rhs(i,j,k) + &
!                  ftMI*( -q2(i,j,k) + q2e*distbMI(2,n) )
!            forclo(i,j,k,2)= (1.-ftMI)
!          end do

! #ifdef USE_CUDA
!         !$cuf kernel do(1)
! #endif
!          do n=1,npunifxMI(2)
!             i=indgeoeeMI(2,n,1)
!             j=indgeoeeMI(2,n,2)
!             k=indgeoeeMI(2,n,3)
!             rhs(i,j,k) = (1.-ftMI)*rhs(i,j,k) + ftMI*(-q2(i,j,k))
!             forclo(i,j,k,2)= (1.-ftMI)
!          enddo

! !Aorta open/close
! !Solid Body AO
! #ifdef USE_CUDA
!         !$cuf kernel do(1)
! #endif
!          do n=1,npunfxAO(2)
!            i=indgeoAO(2,n,1)
!            j=indgeoAO(2,n,2)
!            k=indgeoAO(2,n,3)
!            ie=indgeoeAO(2,n,1)
!            je=indgeoeAO(2,n,2)
!            ke=indgeoeAO(2,n,3)
!            !q2e=((al*dt+aldto)*q2(ie,je,ke)-al*dt*q2o(ie,je,ke))*usaldto
!            q2e=q2(ie,je,ke)
!            rhs(i,j,k) = (1.-ftAO)*rhs(i,j,k) + &
!                  ftAO*( -q2(i,j,k) + q2e*distbAO(2,n) )
!            forclo(i,j,k,2)= (1.-ftAO)
!          end do

! #ifdef USE_CUDA
!         !$cuf kernel do(1)
! #endif
!          do n=1,npunifxAO(2)
!             i=indgeoeeAO(2,n,1)
!             j=indgeoeeAO(2,n,2)
!             k=indgeoeeAO(2,n,3)
!             rhs(i,j,k) = (1.-ftAO)*rhs(i,j,k) + ftAO*(-q2(i,j,k))
!             forclo(i,j,k,2)= (1.-ftAO)
!          enddo

!Aorta system open/close
!Solid Body AS
#ifdef USE_CUDA
        !$cuf kernel do(1)
#endif
        do n=1,npunfxAS(2)
          i=indgeoAS(2,n,1)
          j=indgeoAS(2,n,2)
          k=indgeoAS(2,n,3)
          ie=indgeoeAS(2,n,1)
          je=indgeoeAS(2,n,2)
          ke=indgeoeAS(2,n,3)
!         q2e=((al*dt+aldto)*q2(ie,je,ke)-al*dt*q2bo(n))*usaldto
          q2e = q2(ie,je,ke)


          rhs(i,j,k) = (1.-ftAS)*rhs(i,j,k) +    &
                  ftAS*(-q2(i,j,k) + q2e*distbAS(2,n))+  &
                  (1.-ftAS)*dt*ftAS2*(-q2(i,j,k))
          !q1bo(n)= q1(ie,je,ke)
          forclo(i,j,k,2)= (1.-ftAS)
        end do

#ifdef USE_CUDA
        !$cuf kernel do(1)
#endif
        do n=1,npunifxAS(2)
          i=indgeoeeAS(2,n,1)
          j=indgeoeeAS(2,n,2)
          k=indgeoeeAS(2,n,3)
          rhs(i,j,k) = (1.-ftAS)*rhs(i,j,k)+ftAS*(-q2(i,j,k)) +  &
                  (1.-ftAS)*dt*ftAS2*(-q2(i,j,k))
          forclo(i,j,k,2)= (1.-ftAS)
        end do

!@cuf istat = cudaDeviceSynchronize() ! JDR TMP

      ! end if


!     --------------------------------------------



      call solxi(beta*al*dx1q,2)

      call solxj(beta*al*dx2q,2)
      
      call solq2k(2)
     
      return
      end

