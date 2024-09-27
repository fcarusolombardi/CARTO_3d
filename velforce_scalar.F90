!------------------------------------------------------
!     force the velocity in single phase domain
!------------------------------------------------------
      subroutine velforce_scalar
      USE param
      USE mpih
      USE mpi_param, only: kstartr, kendr
      USE local_arrays, only: dsal
      USE mls_local, only: for_sc
      IMPLICIT NONE


      integer ic,jc,kc
      integer im,jm,km
      integer kstartp
      real(DP) cksum1,cksum2,cksum3
      real(DP) mck1,mck2,mck3


        call update_add_lower_ghost(n1r,n2r,for_sc)
  
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        call update_add_upper_ghost(n1r,n2r,for_sc)
         
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)


      if(kstartr.eq.1) then 
       kstartp=2
      else
       kstartp=kstartr
      end if

      do kc=kstartp,kendr
       km=kmvr(kc)
       do jc=1,n2mr
       jm=jmvr(jc)
        do ic=1,n1mr
         im=imvr(ic)
         dsal(ic,jc,kc)=dsal(ic,jc,kc)+for_sc(ic,jc,kc)
        end do
       end do
      end do
  


      return
      end  
