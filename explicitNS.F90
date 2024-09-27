

!************************************************************************
!                       SUBROUTINE EXPLICITNS
!   This subroutine performs the computation of Q~~ for the q1,2,3 momentum 
!   equation (radial direction) by an explicit AB scheme.
!   
!
      subroutine explicitNS
      use param
      use outflow_vars
      use ibm_param
      use local_arrays, only: rhs,rhss,rhsss,q1,q2,q3,pr,ru1,ru2,ru3
      use mpi_param, only: kstart,kend,kstartp,kendp
!@cuf use cudafor

      implicit none
      integer :: jc,kc,km,kp,kp3,jp,jm,ic,im,ip
      real(DP)    :: udx1,udx2,udx3
      real(DP)    :: amm,acc,app,amm3,acc3,app3
      real(DP)    :: dcq1,dpx11
      real(DP)    :: dcq2,dpx22
      real(DP)    :: dcq3,dpx33
      real(DP)    :: d22q1,d33q1,d11q1
      real(DP)    :: d22q2,d33q2,d11q2
      real(DP)    :: dq32,dq33,dq31
      real(DP)    :: alre,udx1q,udx2q
      real(DP)    :: codeas, codebs, codeae, codebe
      real(DP)    :: usaldto,q1e,q2e,q3e
      real(DP)    :: h11,h12,h13,h21,h22,h23,h31,h32,h33
      real(DP)    :: densit,dsalit

      integer :: i,j,k,ie,je,ke,n
!@cuf integer :: istat
!
      alre=al*nu

      udx1=dx1*al
      udx1q=dx1q
      udx2=dx2*al
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
            kp3=kc+1
            udx3 = al*udx3c(kc)              
            amm3=am3ck(kc)
            acc3=ac3ck(kc)
            app3=ap3ck(kc)
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
        kp3=kc+1
        udx3 = al*udx3c(kc)              
        amm3=am3ck(kc)
        acc3=ac3ck(kc)
        app3=ap3ck(kc)
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

!from hdnlq1
!     q1 q1 term
!
!
!                 d  q_t q_t 
!                ------------
!                 d   t      
!
            h11=( (q1(ip,jc,kc)+q1(ic,jc,kc)) &
                *(q1(ip,jc,kc)+q1(ic,jc,kc))  &
                -(q1(im,jc,kc)+q1(ic,jc,kc))  &
                *(q1(im,jc,kc)+q1(ic,jc,kc))  &
               )*udx1*0.25d0

!     q1 q2 term
!
!
!                 d  q_t q_r 
!                ------------
!                 d   r      
!
            h12=( (q2(ic,jp,kc)+q2(im,jp,kc)) &
                *(q1(ic,jp,kc)+q1(ic,jc,kc))  &
                -(q2(ic,jc,kc)+q2(im,jc,kc))  &
                *(q1(ic,jc,kc)+q1(ic,jm,kc))  &
               )*udx2*0.25d0
!
!     q1 q3 term
!
!
!                 d  q_t q_x 
!                -----------
!                 d   x      
!
!     inflow/outflow boundary condition
            if (kc .eq. n3m) then
              h13=(qb1n(ic,jc)*(qb3n(ic,jc)+qb3n(im,jc))                &
              -(q3(ic,jc,kc)+q3(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jc,km)) &
              *0.5)*udx3m(kc)*0.5
            elseif (kc .eq. 1) then
              h13=-(qb1s(ic,jc)*(qb3s(ic,jc)+qb3s(im,jc))  &
                 -(q3(ic,jc,kp3)+q3(im,jc,kp3))*(q1(ic,jc,kc)+q1(ic,jc,kp3)) &
                 *0.5)*udx3m(kc)*0.5
            else
              h13=( (q3(ic,jc,kp3)+q3(im,jc,kp3)) &
                    *(q1(ic,jc,kp3)+q1(ic,jc,kc)) &
                   -(q3(ic,jc,kc)+q3(im,jc,kc)) &
                    *(q1(ic,jc,kc)+q1(ic,jc,km)) &
                  )*udx3m(kc)*0.25d0
            endif
!
            dcq1=-(h11+h12+h13)


!from hdnlq2
!     q1 q2 term
!
!
!                 d  q_t q_r 
!                ------------
!                 d   t      
!
            h21=( (q2(ip,jc,kc)+q2(ic,jc,kc)) &
                *(q1(ip,jc,kc)+q1(ip,jm,kc))  &
                -(q2(ic,jc,kc)+q2(im,jc,kc))  &
                *(q1(ic,jc,kc)+q1(ic,jm,kc))  &
               )*udx1*0.25d0
      
!     q2 q2 term
!
!
!                 d  q_r q_r 
!                ------------
!                 d   r      
!
            h22=( (q2(ic,jp,kc)+q2(ic,jc,kc)) &
                *(q2(ic,jp,kc)+q2(ic,jc,kc))  &
                -(q2(ic,jm,kc)+q2(ic,jc,kc))  &
                *(q2(ic,jm,kc)+q2(ic,jc,kc))  &
               )*udx2*0.25d0
!
!     q2 q3 term
!
!
!                 d  q_x q_r
!                -----------
!                 d   x
!
!     inflow/outflow boundary condition
          if(kc .eq. n3m) then
            h23=( qb2n(ic,jc)*(qb3n(ic,jc)+qb3n(ic,jm)) &
               -(q3(ic,jc,kc)+q3(ic,jm,kc))*(q2(ic,jc,kc)+q2(ic,jc,km)) &
               *0.5 )*udx3m(kc)*0.5
          elseif(kc .eq. 1) then
            h23= ( (q3(ic,jc,kp3)+q3(ic,jm,kp3))*(q2(ic,jc,kc)+q2(ic,jc,kp3))*0.5 &
                -qb2s(ic,jc)*(qb3s(ic,jc)+qb3s(ic,jm)))*udx3m(kc)*0.5
          else
            h23=( (q3(ic,jc,kp3)+q3(ic,jm,kp3)) &
                  *(q2(ic,jc,kp3)+q2(ic,jc,kc)) &
                 -(q3(ic,jc,kc)+q3(ic,jm,kc)) &
                  *(q2(ic,jc,kc)+q2(ic,jc,km)) &
                )*udx3m(kc)*0.25d0
          endif
!
!
!      the buoyancy term
!
!      temperature           
#ifdef CALTEMP_VC
            densit=(dens(ic,jc,kc)+dens(ic,jm,kc))*0.5d0
#else
            densit=0.0d0
#endif

!      salinity 
#ifdef CALSCAL_VC
            dsalit=(dsalc(ic,jc,kc)+dsalc(ic,jm,kc))*0.5d0
#else
            dsalit=0.0d0
#endif


            dcq2=-(h21+h22+h23) + byct * densit - bycs * dsalit

!from hdnlq3
!
!    q3 q1 term
!
!                d  q_x q_t 
!              -------------
!                d   t      
!
            h31=( (q1(ip,jc,kc)+q1(ip,jc,km))   &
                *(q3(ip,jc,kc)+q3(ic,jc,kc))    &
                -(q1(ic,jc,kc)+q1(ic,jc,km))    &
                *(q3(ic,jc,kc)+q3(im,jc,kc)))*udx1*0.25d0
!
!    q3 q2 term
!
!                d  q_x q_r 
!              -------------
!                d   r      
!
            h32=( (q2(ic,jp,kc)+q2(ic,jp,km)) &
                *(q3(ic,jp,kc)+q3(ic,jc,kc))  &
                -(q2(ic,jc,kc)+q2(ic,jc,km))  &
                *(q3(ic,jc,kc)+q3(ic,jm,kc)))*udx2*0.25d0
!
!    q3 q3 term
!
!               d  q_x q_x 
!              ------------
!               d   x      
!
            h33=( (q3(ic,jc,kp3)+q3(ic,jc,kc)) &
                *(q3(ic,jc,kp3)+q3(ic,jc,kc))  &
                -(q3(ic,jc,kc)+q3(ic,jc,km))  &
                *(q3(ic,jc,kc)+q3(ic,jc,km))  &
               )*udx3c(kc)*0.25d0
!
!      the buoyancy term
!
!      temperature           
#ifdef CALTEMP_HC
            densit=(dens(ic,jc,kc)+dens(ic,jc,km))*0.5d0
#else
            densit=0.0d0
#endif

!      salinity 
#ifdef CALSCAL_HC
            dsalit=(dsalc(ic,jc,kc)+dsalc(ic,jc,km))*0.5d0
#else
            dsalit=0.0d0
#endif

!
!      add buoyangcy force
            dcq3 = -(h31+h32+h33) + byct * densit - bycs * dsalit



! VISCOUS TERMS


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
            dq33=q3(ic,jc,kp3)*app3   &
               +q3(ic,jc,kc)*acc3    &
               +q3(ic,jc,km)*amm3

!

!
!    viscid terms
!
            dcq1=dcq1+nu*(d11q1+d22q1+d33q1)
            dcq2=dcq2+nu*(d22q2+d33q2+d11q2)
            dcq3=dcq3+nu*(dq32+dq33+dq31)
!
!   component of grad(pr) along 2 direction
!
            dpx11=(pr(ic,jc,kc)-pr(im,jc,kc))*udx1
            dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*udx2
            dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*udx3
!
            rhs(ic,jc,kc)=(ga*dcq1+ro*ru1(ic,jc,kc) -dpx11)*dt
            ru1(ic,jc,kc)=dcq1

            rhss(ic,jc,kc)=(ga*dcq2+ro*ru2(ic,jc,kc) -dpx22)*dt
            ru2(ic,jc,kc)=dcq2

            rhsss(ic,jc,kc)=(ga*dcq3+ro*ru3(ic,jc,kc) -dpx33)*dt 
            ru3(ic,jc,kc)=dcq3
!===========================================================
!
          enddo
        enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      enddo

!@cuf istat = cudaDeviceSynchronize() ! JDR TMP


!      usaldto = 1./aldto

!FADLUN TERMS 1
if(infig.eq.1) then
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
                  (1.-ftAS)*dt*ftAS2*(-q1(i,j,k)) + dt*ftAS3
        end do
          
#ifdef USE_CUDA
        !$cuf kernel do(1)
#endif
        do n=1,npunifxAS(1)
          i=indgeoeeAS(1,n,1)
          j=indgeoeeAS(1,n,2)
          k=indgeoeeAS(1,n,3)
          rhs(i,j,k) = (1.-ftAS)*rhs(i,j,k)+ftAS*(-q1(i,j,k)) +  &
                  (1.-ftAS)*dt*ftAS2*(-q1(i,j,k)) + dt*ftAS3
        end do

!@cuf istat = cudaDeviceSynchronize() ! JDR TMP

      ! end if

!FADLUN TERMS 2
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


          rhss(i,j,k) = (1.-ftAS)*rhss(i,j,k) +    &
                  ftAS*(-q2(i,j,k) + q2e*distbAS(2,n))+  &
                  (1.-ftAS)*dt*ftAS2*(-q2(i,j,k))
        end do

#ifdef USE_CUDA
        !$cuf kernel do(1)
#endif
        do n=1,npunifxAS(2)
          i=indgeoeeAS(2,n,1)
          j=indgeoeeAS(2,n,2)
          k=indgeoeeAS(2,n,3)
          rhss(i,j,k) = (1.-ftAS)*rhss(i,j,k)+ftAS*(-q2(i,j,k)) +  &
                  (1.-ftAS)*dt*ftAS2*(-q2(i,j,k))
        end do

!@cuf istat = cudaDeviceSynchronize() ! JDR TMP


!FADLUN TERMS 3
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
          rhsss(i,j,k) = (1.-ftAS)*rhsss(i,j,k) +             &
                  ftAS*(-q3(i,j,k) + q3e*distbAS(3,n))+   &
                  (1.-ftAS)*dt*ftAS2*(-q3(i,j,k)) + dt*ftAS3 
        end do

#ifdef USE_CUDA
        !$cuf kernel do(1)
#endif
        do n=1,npunifxAS(3)
          i=indgeoeeAS(3,n,1)
          j=indgeoeeAS(3,n,2)
          k=indgeoeeAS(3,n,3)
          rhsss(i,j,k) = (1.-ftAS)*rhsss(i,j,k)+ftAS*(-q3(i,j,k)) +  &
                  (1.-ftAS)*dt*ftAS2*(-q3(i,j,k)) + dt*ftAS3 
        end do

!@cuf istat = cudaDeviceSynchronize() ! JDR TMP
endif 


!UPDATE VELOCITY
#ifdef USE_CUDA
      !$cuf kernel do(3)
#endif
      do kc=kstart,kend
#ifndef USE_CUDA
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ic,jc)
#endif
        do jc=1,n2m
          do ic=1,n1m
            q1(ic,jc,kc) = q1(ic,jc,kc) + rhs(ic,jc,kc)
            q2(ic,jc,kc) = q2(ic,jc,kc) + rhss(ic,jc,kc)
            q3(ic,jc,kc) = q3(ic,jc,kc) + rhsss(ic,jc,kc)
          enddo
        enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      enddo

!@cuf istat = cudaDeviceSynchronize() ! JDR TMP


      if (kend .eq. n3m) then !condition n3m, write in n3
#ifdef USE_CUDA
      !$cuf kernel do(2)
#endif
         do jc=1,n2m
            do ic=1,n1m
               q3(ic,jc,n3) = qb3n(ic,jc)
            enddo
         enddo
!@cuf istat = cudaDeviceSynchronize() ! JDR TMP
      endif

      if (kstart .eq. 1) then !condition 1, write in 1
#ifdef USE_CUDA
      !$cuf kernel do(2)
#endif
         do jc=1,n2m
            do ic=1,n1m
               q3(ic,jc,1) = qb3n(ic,jc)
            enddo
         enddo
!@cuf istat = cudaDeviceSynchronize() ! JDR TMP
      endif
     
     
      return
      end


!NOTE that rhs, rhss, rhsss could be replaced by other arrays such as dq dph qcap 
!this replacement would allow gaining some memory but without a speedup, which depends
!on the number of large arrays used by the routine that would still be 3. (i've checked this, 16-12-20)

!NOTE that the Fadlun could be applied directly to q* without writing the rhs as follows
        ! do n=1,npunfxAS(1)
        !   i=indgeoAS(1,n,1)
        !   j=indgeoAS(1,n,2)
        !   k=indgeoAS(1,n,3)
        !   ie=indgeoeAS(1,n,1)
        !   je=indgeoeAS(1,n,2)
        !   ke=indgeoeAS(1,n,3)
        !   q1e = q1(ie,je,ke)
        !   q1(i,j,k) = (1.-ftAS)*q1(i,j,k) +              &
        !                    ftAS*q1e*distbAS(1,n)+    &
        !               (1.-ftAS)*dt*ftAS2*(-q1(i,j,k)) 
        ! end do
        ! do n=1,npunifxAS(1)
        !   i=indgeoeeAS(1,n,1)
        !   j=indgeoeeAS(1,n,2)
        !   k=indgeoeeAS(1,n,3)
        !   q1(i,j,k) = (1.-ftAS)*q1(i,j,k) +              &
        !               (1.-ftAS)*dt*ftAS2*(-q1(i,j,k)) 
        ! end do

