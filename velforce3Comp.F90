!------------------------------------------------------
!     force the velocity in single phase domain
!------------------------------------------------------
      subroutine velforce3Comp
      USE param
      USE mpih
      USE mpi_param, only: kstart, kend
      USE local_arrays, only: q1, q2, q3
      USE mls_local, only: for_xc, for_yc, for_zc
      USE mls_param, only: nsstep
!@cuf USE cudafor
      IMPLICIT NONE


      integer ic,jc,kc
      integer im,jm,km
      integer kstartp
      real(DP) nsfac
      real(DP) cksum1,cksum2,cksum3
      real(DP) mck1,mck2,mck3
!@cuf integer istat


!      nsfac = 1.0/dble(nsstep) !to be modified based on nsstep 
      nsfac = 1.0D0
#ifdef USE_CUDA
      call update_add_lower_ghost_3(n1,n2,for_xc,for_yc,for_zc)
#else
      call update_add_lower_ghost(n1,n2,for_xc)
      call update_add_lower_ghost(n1,n2,for_yc)
      call update_add_lower_ghost(n1,n2,for_zc)
#endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

#ifdef USE_CUDA
      call update_add_upper_ghost_3(n1,n2,for_xc,for_yc,for_zc)
#else
      call update_add_upper_ghost(n1,n2,for_xc)
      call update_add_upper_ghost(n1,n2,for_yc)
      call update_add_upper_ghost(n1,n2,for_zc)
#endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!      call update_both_ghosts(n1,n2,for_xc,kstart,kend)
!      call update_both_ghosts(n1,n2,for_yc,kstart,kend)
!      call update_both_ghosts(n1,n2,for_zc,kstart,kend)
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)



      if(kstart.eq.1) then 
       kstartp=2
      else
       kstartp=kstart
      end if

!     Add forcing | Correct velocities
#ifdef USE_CUDA
      !$cuf kernel do(3)
#endif
      do kc=kstartp,kend
       do jc=1,n2m
        do ic=1,n1m
         !km=kmv(kc)
         !jm=jmv(jc)
         !im=imv(ic)
#ifdef SCALINGTEST
         q1(ic,jc,kc)=q1(ic,jc,kc)+0.0D0*for_xc(ic,jc,kc)*nsfac
         q2(ic,jc,kc)=q2(ic,jc,kc)+0.0D0*for_yc(ic,jc,kc)*nsfac
         q3(ic,jc,kc)=q3(ic,jc,kc)+0.0D0*for_zc(ic,jc,kc)*nsfac
#else
         q1(ic,jc,kc)=q1(ic,jc,kc)+for_xc(ic,jc,kc)*nsfac
         q2(ic,jc,kc)=q2(ic,jc,kc)+for_yc(ic,jc,kc)*nsfac
         q3(ic,jc,kc)=q3(ic,jc,kc)+for_zc(ic,jc,kc)*nsfac
#endif
        end do
       end do
      end do
  
!@cuf istat = cudaDeviceSynchronize() !JDR TMP

      return
      end  
