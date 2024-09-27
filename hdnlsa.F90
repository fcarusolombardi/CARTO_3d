!***********************************************************************
      subroutine hdnlsa
      use param
      use local_arrays, only: hsal,dsal
      use mgrd_arrays, only: q1lr,q2lr,q3lr
      use mpi_param, only: kstartr,kendr
      implicit none
      integer :: jc,kc,ic
      integer :: kpp,kmm,kp,jp,jmm,jpp,ip,imm,ipp
      real(DP)    :: h31,h32,h33,udx2r,udx1r
!
!     h term for the q3 momentum equation at i+1/2,j+1/2,k
!
      udx1r=dx1r*0.5d0
      udx2r=dx2r*0.5d0
      do kc=kstartr,kendr
        kmm=kmvr(kc)
        kpp=kpvr(kc)
        kp=kc+1
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,jp,jmm,jpp,ic,ip,ipp,imm) &
!$OMP PRIVATE(h31,h32,h33)
        do jc=1,n2mr
          jp=jc+1
          jmm=jmvr(jc)
          jpp=jpvr(jc)
          do ic=1,n1mr
            ip=ic+1
            ipp=ipvr(ic)
            imm=imvr(ic)
!
!
!    rho q1 term
!
!
!                d  rho q_t 
!             -----------
!                d   t      
!
            h31=( q1lr(ip,jc,kc)*(dsal(ipp,jc,kc)+dsal(ic,jc,kc)) &
                -q1lr(ic,jc,kc)*(dsal(ic,jc,kc)+dsal(imm,jc,kc))  &
               )*udx1r
!
!
!    rho q2 term
!
!
!                d  rho q_r 
!             -----------
!                d   r      
!
            h32=( q2lr(ic,jp,kc)*(dsal(ic,jpp,kc)+dsal(ic,jc,kc)) &
                -q2lr(ic,jc,kc)*(dsal(ic,jc,kc)+dsal(ic,jmm,kc))  &
               )*udx2r
!
!    rho q3 term
!
!
!                 d  rho q_x 
!                -----------
!                 d   x      
!
            h33=( (dsal(ic,jc,kpp)*g3rmr(kc)+dsal(ic,jc,kc)*g3rmr(kpp)) &
                   /(g3rmr(kc)+g3rmr(kpp))*q3lr(ic,jc,kp)               &
                -(dsal(ic,jc,kc)*g3rmr(kmm)+dsal(ic,jc,kmm)*g3rmr(kc))  &
                   /(g3rmr(kc)+g3rmr(kmm))*q3lr(ic,jc,kc)               &
               )*udx3mr(kc)

            hsal(ic,jc,kc)= -(h31+h32+h33)
     
          enddo
        enddo
!$OMP  END PARALLEL DO
      enddo

      return
      end

