      subroutine hdnlq1
      use param
      use outflow_vars
      use local_arrays, only: q1,q2,q3,dq
      use mpi_param, only: kstart,kend,kstartp,kendp
      implicit none
      integer :: kc,kp,jp,jm,jc,ic,im,ip
      integer :: kmm,kpp
      real(DP)    :: h11,h12,h13,udx1,udx2


      udx1=dx1*0.25d0
      udx2=dx2*0.25d0


#ifdef USE_CUDA
      !$cuf kernel do(3)
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
!$OMP PRIVATE(h11,h12,h13)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
#endif
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
               )*udx1

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
               )*udx2
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
              -(q3(ic,jc,kc)+q3(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jc,kmm)) &
              *0.5)*udx3m(kc)*0.5
            elseif (kc .eq. 1) then
              h13=-(qb1s(ic,jc)*(qb3s(ic,jc)+qb3s(im,jc))  &
                 -(q3(ic,jc,kp)+q3(im,jc,kp))*(q1(ic,jc,kc)+q1(ic,jc,kp)) &
                 *0.5)*udx3m(kc)*0.5
            else
              h13=( (q3(ic,jc,kp)+q3(im,jc,kp)) &
                    *(q1(ic,jc,kpp)+q1(ic,jc,kc)) &
                   -(q3(ic,jc,kc)+q3(im,jc,kc)) &
                    *(q1(ic,jc,kc)+q1(ic,jc,kmm)) &
                  )*udx3m(kc)*0.25d0
            endif
!
            dq(ic,jc,kc)=-(h11+h12+h13)
!
          enddo
        enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      enddo

      return
      end
