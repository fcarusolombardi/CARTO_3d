!***********************************************************************
      subroutine hdnlq3
      use param
      use local_arrays, only: q1,q2,q3,qcap,dens
      use mgrd_arrays, only: dsalc
      use mpi_param, only: kstart,kend,kstartp,kendp
      implicit none
      integer :: ic,jc,kc
      integer :: km,kp,jm,jp,im,ip
      real(DP)    :: h32,h33,h31
      real(DP)    :: densit,dsalit,udx1,udx2



      udx1=dx1*0.25d0
      udx2=dx2*0.25d0

#ifdef USE_CUDA
      !$cuf kernel do(3)
      do kc=kstartp,kend
        do jc=1,n2m
          do ic=1,n1m
            km=kmv(kc)
            kp=kc+1
            jm=jmv(jc)
            jp=jpv(jc)
            im=imv(ic)
            ip=ipv(ic)
#else
      do kc=kstartp,kend
        km=kmv(kc)
        kp=kc+1
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(jc,jm,jp,ic,im,ip) &
!$OMP PRIVATE(h31,h32,h33) &
!$OMP PRIVATE(densit,dsalit)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
#endif
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
                *(q3(ic,jc,kc)+q3(im,jc,kc)))*udx1
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
                *(q3(ic,jc,kc)+q3(ic,jm,kc)))*udx2
!
!    q3 q3 term
!
!               d  q_x q_x 
!              ------------
!               d   x      
!
            h33=( (q3(ic,jc,kp)+q3(ic,jc,kc)) &
                *(q3(ic,jc,kp)+q3(ic,jc,kc))  &
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
            qcap(ic,jc,kc) = -(h31+h32+h33) &
                            + byct * densit - bycs * dsalit

          enddo
        enddo
#ifndef USE_CUDA
!$OMP  END PARALLEL DO
#endif
      enddo

      return
      end

