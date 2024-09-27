
!***********************************************************************
      subroutine hdnlte
      use param
      use local_arrays, only: q1,q2,q3,hro,dens
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,ic
      integer :: kpp,kmm,kp,jp,jmm,jpp,ip,imm,ipp
      real(DP)    :: h32,h33,udx2,udx1,h31
!
!     h term for the q3 momentum equation at i+1/2,j+1/2,k
!
      udx1=dx1*0.5d0
      udx2=dx2*0.5d0
      do kc=kstart,kend
        kmm=kmv(kc)
        kpp=kpv(kc)
        kp=kc+1
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,jp,jmm,jpp,ic,ip,ipp,imm) &
!$OMP PRIVATE(h31,h32,h33)
        do jc=1,n2m
          jp=jc+1
          jmm=jmv(jc)
          jpp=jpv(jc)
          do ic=1,n1m
            ip=ic+1
            ipp=ipv(ic)
            imm=imv(ic)
!
!
!    rho q1 term
!
!
!                d  rho q_t 
!             -----------
!                d   t      
!
            h31=( q1(ip,jc,kc)*(dens(ipp,jc,kc)+dens(ic,jc,kc)) &
                -q1(ic,jc,kc)*(dens(ic,jc,kc)+dens(imm,jc,kc))  &
               )*udx1
!
!
!    rho q2 term
!
!
!                d  rho q_r 
!             -----------
!                d   r      
!
            h32=( q2(ic,jp,kc)*(dens(ic,jpp,kc)+dens(ic,jc,kc)) &
                -q2(ic,jc,kc)*(dens(ic,jc,kc)+dens(ic,jmm,kc))  &
               )*udx2
!
!    rho q3 term
!
!
!                 d  rho q_x 
!                -----------
!                 d   x    
!  
            h33=( (dens(ic,jc,kpp)*g3rm(kc)+dens(ic,jc,kc)*g3rm(kpp)) &
                 *q3(ic,jc,kp)/(g3rm(kc)+g3rm(kpp))                   &
                -(dens(ic,jc,kc)*g3rm(kmm)+dens(ic,jc,kmm)*g3rm(kc))  &
                 *q3(ic,jc,kc)/(g3rm(kc)+g3rm(kmm))                   &
               )*udx3m(kc)

            hro(ic,jc,kc)=-(h31+h32+h33)
     
          enddo
        enddo
!$OMP  END PARALLEL DO
      enddo

      return
      end

