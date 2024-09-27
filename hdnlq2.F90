!***********************************************************************
      subroutine hdnlq2
      use param
      use outflow_vars
      use local_arrays, only: q1,q2,q3,dph,dens
      use mgrd_arrays, only: dsalc
      use mpi_param, only: kstart,kend,kstartp,kendp
      implicit none
      integer :: kc,kp,jp,jm,jc,ic,im,ip
      integer :: kmm,kpp
      real(DP)    :: h22,h23,udx1,udx2,h21
      real(DP)    :: densit,dsalit


      udx1=dx1*0.25d0
      udx2=dx2*0.25d0
!
!     h term for the q2 momentum equation at i+1/2,j,k+1/2
!

#ifdef USE_CUDA
      !$cuf kernel do (3)
      do kc=kstart,kend
        do jc=1,n2m
          do ic=1,n1m
            kmm=kmv(kc)
            kpp=kpv(kc)
            kp=kc+1
            jm=jmv(jc)
            jp=jpv(jc)
            im=imv(ic)
            ip=ipv(ic)
#else
      do kc=kstart,kend
        kmm=kmv(kc)
        kpp=kpv(kc)
        kp=kc+1
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,jm,jp,ic,im,ip) &
!$OMP PRIVATE(h21,h22,h23)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
#endif
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
               )*udx1
      
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
               )*udx2
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
               -(q3(ic,jc,kc)+q3(ic,jm,kc))*(q2(ic,jc,kc)+q2(ic,jc,kmm)) &
               *0.5 )*udx3m(kc)*0.5
          elseif(kc .eq. 1) then
            h23= ( (q3(ic,jc,kp)+q3(ic,jm,kp))*(q2(ic,jc,kc)+q2(ic,jc,kp))*0.5 &
                -qb2s(ic,jc)*(qb3s(ic,jc)+qb3s(ic,jm)))*udx3m(kc)*0.5
          else
            h23=( (q3(ic,jc,kp)+q3(ic,jm,kp)) &
                  *(q2(ic,jc,kpp)+q2(ic,jc,kc)) &
                 -(q3(ic,jc,kc)+q3(ic,jm,kc)) &
                  *(q2(ic,jc,kc)+q2(ic,jc,kmm)) &
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


            dph(ic,jc,kc)=-(h21+h22+h23) &
                             + byct * densit - bycs * dsalit


          enddo
        enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      enddo

      return
      end

