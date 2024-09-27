#ifdef USE_BUDA 
     module gpu_balance
      contains
! we are  assuming that the block shape is 32 x BlockDim%y with BlockDim%y less than 32
      attributes(global) subroutine balance_kernel(kstart,kend,n1m,n2m,n3m, &
                                          dx1,dx2,nu,kpt, &
                                          dissuc, dissipuc_partial, &
                                          disste, dissipte_partial)


       use local_arrays, only: q1=>q1_d,q2=>q2_d,q3=>q3_d,dens=>dens_d
       use param, only: udx3c=>udx3c_d,g3rm=>g3rm_d,g3rc=>g3rc_d
       use param, only: densbn=>densbn_d
       use param, only : ipv=>ipv_d,jpv=>jpv_d
       use param, only : imv=>imv_d,jmv=>jmv_d
       use param, only : zm=>zm_d,zc=>zc_d

       implicit none

       integer,value:: kstart,kend,n1m,n2m,n3m
       real(DP),value:: dx1,dx2,nu,kpt
       real(DP):: dissuc(kstart:kend), dissipuc_partial(kstart:kend), &
                  disste(kstart:kend), dissipte_partial(kstart:kend)

       integer :: imm,ipp,jmm,jpp,kmm,kpp
       real(DP) :: udx3_m, udx3_p
       real(DP) :: h11,h12,h13,h21,h22,h23,h31,h32,h33
       real(DP) :: h13p,h23p,h31p,h32p
       real(DP) :: nuuc,nute,nusa,nuur,volb,volr
       real(DP) :: udx1,udx2,un12m,unu2
       real(DP) :: dissipur, my_dissipur
       real(DP) :: dissipuc, my_dissipuc
       real(DP) :: dissipte, my_dissipte
       real(DP) :: dissipsa, my_dissipsa
       real(DP) :: loc_var(4),val
       real(DP),shared:: val_s(32)

  
       integer:: tx,ty,i,j,k

       udx1=dx1
       udx2=dx2

       un12m = 1.d0/(dble(n1m)*dble(n2m))
       unu2 = 1.d0/(nu*nu)

       tx = threadIdx%x
       ty = threadIdx%y
       k  = kstart+(blockIdx%x-1)

       kpp=k+1
       kmm=k-1
       if(k.eq.1)then
          udx3_m=1.d0/(zm(k)-zc(k))
       else
          udx3_m=1.d0/(zm(k)-zm(kmm))
       endif
       udx3_p=1.d0/(zm(kpp)-zm(k))

       do i=1,4
        loc_var(i)=0.d0
       end do

       do j=ty,n2m,blockDim%y
        jmm=jmv(j)
        jpp=jpv(j)
        do i=tx,n1m,blockDim%x
            imm= imv(i)
            ipp= ipv(i)

            h11=(q1(ipp,j,k)-q1(i,j,k))*udx1
            h22=(q2(i,jpp,k)-q2(i,j,k))*udx2
            h33=(q3(i,j,kpp)-q3(i,j,k))*udx3c(k)

            h12=(q1(i,j,k)-q1(i,jmm,k))*udx2
            h21=(q2(i,j,k)-q2(imm,j,k))*udx1

            h13=(q1(i,j,k)-q1(i,j,kmm))*udx3_m
            h31=(q3(i,j,k)-q3(imm,j,k))*udx1

            h13p=(q1(i,j,kpp)-q1(i,j,k))*udx3_p
            h31p=(q3(i,j,kpp)-q3(imm,j,kpp))*udx1

            h23=(q2(i,j,k)-q2(i,j,kmm))*udx3_m
            h32=(q3(i,j,k)-q3(i,jmm,k))*udx2

            h23p=(q2(i,j,kpp)-q2(i,j,k))*udx3_p
            h32p=(q3(i,j,kpp)-q3(i,jmm,kpp))*udx2
            
            dissipuc = 2.d0*(h11**2+h22**2+h33**2) &
                + (h21+h12)**2 &
                + ((h31+h13)**2+(h31p+h13p)**2)*0.5d0 &
                + ((h32+h23)**2+(h32p+h23p)**2)*0.5d0  
  

            !my_dissipuc = my_dissipuc+dissipuc*g3rm(k)

            loc_var(1)=loc_var(1)+dissipuc*nu*un12m
            loc_var(2)=loc_var(2)+dissipuc*g3rm(k)

            h31=(dens(ipp,j,k)-dens(i,j,k))*udx1
            h32=(dens(i,jpp,k)-dens(i,j,k))*udx2
            if(k.eq.n3m)then
              h33=(densbn(i,j)-dens(i,j,k))*udx3_p
            else
              h33=(dens(i,j,kpp)-dens(i,j,k))*udx3_p
            endif

            dissipte = h31*h31 + h32*h32 + h33*h33 
            !my_dissipte = my_dissipte+dissipte*g3rc(k)
            loc_var(3)= loc_var(3)+dissipte*kpt*un12m
            loc_var(4)= loc_var(4)+dissipte*g3rc(k)
        end do
       end do

       ! Reduce inside each block
       do i=1,4
        val = __shfl_down(loc_var(i),16)
        loc_var(i) = loc_var(i) + val
        val = __shfl_down(loc_var(i),8)
        loc_var(i) = loc_var(i) + val
        val = __shfl_down(loc_var(i),4)
        loc_var(i) = loc_var(i) + val
        val = __shfl_down(loc_var(i),2)
        loc_var(i) = loc_var(i) + val
        val = __shfl_down(loc_var(i),1)
        loc_var(i) = loc_var(i) + val
         if (tx ==1) then
           val_s(ty)=loc_var(i)
         end if
         call syncthreads()
         ! first warp does final reduction
         if (tx==1 .and. ty ==1) then
            loc_var(i)=0.d0
            do j=1,blockDim%y
            loc_var(i)=loc_var(i)+val_s(j)
            end do
         endif
         call syncthreads()
       end do
       ! First thread in first warp write back 
       if ( tx == 1 .and. ty==1) then
        dissuc(k) = dissuc(k) + loc_var(1)
        dissipuc_partial(k)=  loc_var(2)
        disste(k) = disste(k) + loc_var(3)
        dissipte_partial(k)=  loc_var(4)
       end if

      end subroutine balance_kernel

      attributes(global) subroutine balance_r_kernel(kstartr,kendr,n1mr,n2mr,n3mr, &
                                          dx1r,dx2r,nu,kps, &
                                          dissur, dissipur_partial, &
                                          disssa, dissipsa_partial)


       use local_arrays, only: dsal=>dsal_d
       use mgrd_arrays, only: q1lr=>q1lr_d,q2lr=>q2lr_d,q3lr=>q3lr_d
       use param, only: udx3cr=>udx3cr_d,g3rmr=>g3rmr_d,g3rcr=>g3rcr_d
       use param, only: dsalbn=>dsalbn_d
       use param, only : ipvr=>ipvr_d,jpvr=>jpvr_d
       use param, only : imvr=>imvr_d,jmvr=>jmvr_d
       use param, only : zmr=>zmr_d,zcr=>zcr_d

       implicit none

       integer,value:: kstartr,kendr,n1mr,n2mr,n3mr
       real(DP),value:: dx1r,dx2r,nu,kps
       real(DP):: dissur(kstartr:kendr), dissipur_partial(kstartr:kendr), &
                 disssa(kstartr:kendr), dissipsa_partial(kstartr:kendr)

       integer :: imm,ipp,jmm,jpp,kmm,kpp
       real(DP) :: udx3_m, udx3_p
       real(DP) :: h11,h12,h13,h21,h22,h23,h31,h32,h33
       real(DP) :: h13p,h23p,h31p,h32p
       real(DP) :: nuuc,nute,nusa,nuur,volb,volr
       real(DP) :: udx1,udx2,un12m,unu2
       real(DP) :: dissipur, my_dissipur
       real(DP) :: dissipuc, my_dissipuc
       real(DP) :: dissipte, my_dissipte
       real(DP) :: dissipsa, my_dissipsa
       real(DP) :: loc_var(4),val
       real(DP),shared:: val_s(32)

  
       integer:: tx,ty,i,j,k

       udx1=dx1r
       udx2=dx2r

       un12m = 1.d0/(dble(n1mr)*dble(n2mr))
       unu2 = 1.d0/(nu*nu)

       tx = threadIdx%x
       ty = threadIdx%y
       k  = kstartr+(blockIdx%x-1)

       kpp=k+1
       kmm=k-1
       if(k.eq.1)then
          udx3_m=1.d0/(zmr(k)-zcr(k))
       else
          udx3_m=1.d0/(zmr(k)-zmr(kmm))
       endif
       udx3_p=1.d0/(zmr(kpp)-zmr(k))

       do i=1,4
        loc_var(i)=0.d0
       end do

       do j=ty,n2mr,blockDim%y
        jmm=jmvr(j)
        jpp=jpvr(j)
        do i=tx,n1mr,blockDim%x
            imm= imvr(i)
            ipp= ipvr(i)

            h11=(q1lr(ipp,j,k)-q1lr(i,j,k))*udx1
            h22=(q2lr(i,jpp,k)-q2lr(i,j,k))*udx2
            h33=(q3lr(i,j,kpp)-q3lr(i,j,k))*udx3cr(k)

            h12=(q1lr(i,j,k)-q1lr(i,jmm,k))*udx2
            h21=(q2lr(i,j,k)-q2lr(imm,j,k))*udx1

            h13=(q1lr(i,j,k)-q1lr(i,j,kmm))*udx3_m
            h31=(q3lr(i,j,k)-q3lr(imm,j,k))*udx1

            h13p=(q1lr(i,j,kpp)-q1lr(i,j,k))*udx3_p
            h31p=(q3lr(i,j,kpp)-q3lr(imm,j,kpp))*udx1

            h23=(q2lr(i,j,k)-q2lr(i,j,kmm))*udx3_m
            h32=(q3lr(i,j,k)-q3lr(i,jmm,k))*udx2

            h23p=(q2lr(i,j,kpp)-q2lr(i,j,k))*udx3_p
            h32p=(q3lr(i,j,kpp)-q3lr(i,jmm,kpp))*udx2
            
            dissipur = 2.d0*(h11**2+h22**2+h33**2) &
                + (h21+h12)**2 &
                + ((h31+h13)**2+(h31p+h13p)**2)*0.5d0 &
                + ((h32+h23)**2+(h32p+h23p)**2)*0.5d0  
  

            !my_dissipuc = my_dissipuc+dissipuc*g3rm(k)

            loc_var(1)=loc_var(1)+dissipur*nu*un12m
            loc_var(2)=loc_var(2)+dissipur*g3rmr(k)

            h31=(dsal(ipp,j,k)-dsal(i,j,k))*udx1
            h32=(dsal(i,jpp,k)-dsal(i,j,k))*udx2
            if(k.eq.n3mr)then
              h33=(dsalbn(i,j)-dsal(i,j,k))*udx3_p
            else
              h33=(dsal(i,j,kpp)-dsal(i,j,k))*udx3_p
            endif

            dissipsa = h31*h31 + h32*h32 + h33*h33 
            !my_dissipte = my_dissipte+dissipte*g3rc(k)
            loc_var(3)= loc_var(3)+dissipsa*kps*un12m
            loc_var(4)= loc_var(4)+dissipsa*g3rcr(k)
        end do
       end do

       ! Reduce inside each block
       do i=1,4
        val = __shfl_down(loc_var(i),16)
        loc_var(i) = loc_var(i) + val
        val = __shfl_down(loc_var(i),8)
        loc_var(i) = loc_var(i) + val
        val = __shfl_down(loc_var(i),4)
        loc_var(i) = loc_var(i) + val
        val = __shfl_down(loc_var(i),2)
        loc_var(i) = loc_var(i) + val
        val = __shfl_down(loc_var(i),1)
        loc_var(i) = loc_var(i) + val
         if (tx ==1) then
           val_s(ty)=loc_var(i)
         end if
         call syncthreads()
         ! first warp does final reduction
         if (tx==1 .and. ty ==1) then
            loc_var(i)=0.d0
            do j=1,blockDim%y
            loc_var(i)=loc_var(i)+val_s(j)
            end do
         endif
         call syncthreads()
       end do
       ! First thread in first warp write back 
       if ( tx == 1 .and. ty==1) then
        dissur(k) = dissur(k) + loc_var(1)
        dissipur_partial(k)=  loc_var(2)
        disssa(k) = disssa(k) + loc_var(3)
        dissipsa_partial(k)=  loc_var(4)
       end if

      end subroutine balance_r_kernel

      subroutine balance_gpu
       use cudafor
       use param
       use mpi_param, only: kstart,kend,kstartr,kendr
       use stat_arrays
       implicit none
 
       integer    :: blocks,k
       type(dim3) :: threads
       threads = dim3(32, 8, 1)
       ! Main grid
       blocks =  kend - kstart +1 
       densbn_d=densbn
       call balance_kernel<<<blocks,threads>>>( kstart,kend,n1m,n2m,n3m, &
                                          dx1,dx2,nu,kpt, &
                                          dissuc_d, dissipuc_partial_d, &
                                          disste_d, dissipte_partial_d)
       ! Salinity grid
       blocks =  kendr - kstartr +1 
       dsalbn_d=dsalbn
       call balance_r_kernel<<<blocks,threads>>>( kstartr,kendr,n1mr,n2mr,n3mr, &
                                          dx1r,dx2r,nu,kps, &
                                          dissur_d, dissipur_partial_d, &
                                          disssa_d, dissipsa_partial_d)

       

      end subroutine balance_gpu

      end module gpu_balance
#endif

      subroutine balance
      use mpih
      use param
#ifdef USE_BUDA
      use gpu_balance
#else
      use local_arrays,only: q1,q2,q3,dens,dsal
      use mgrd_arrays,only: q1lr,q2lr,q3lr
#endif
      use mpi_param, only: kstart,kend,kstartr,kendr
      use stat_arrays

      implicit none
      integer :: i,j,k
      integer :: imm,ipp,jmm,jpp,kmm,kpp
      real(DP) :: udx3_m, udx3_p
      real(DP) :: h11,h12,h13,h21,h22,h23,h31,h32,h33
      real(DP) :: h13p,h23p,h31p,h32p
      real(DP) :: nuuc,nute,nusa,nuur,volb,volr
      real(DP) :: udx1,udx2,un12m,unu2
      real(DP) :: dissipur, my_dissipur
      real(DP) :: dissipuc, my_dissipuc
      real(DP) :: dissipte, my_dissipte
      real(DP) :: dissipsa, my_dissipsa
      real(DP) :: dissuc_loc,disste_loc
      real(DP) :: dissur_loc,disssa_loc
      
      my_dissipuc = 0.0d0
      my_dissipte = 0.0d0
      my_dissipsa = 0.0d0
      my_dissipur = 0.0d0

      udx1=dx1
      udx2=dx2

      un12m = 1.d0/(dble(n1m)*dble(n2m))

      unu2 = 1.d0/(nu*nu)
      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Dissipation rates
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       	
!
!                                   1  |         | 2
!                   dissipation:  ---- | nabla  u|
!                                  Re  |         |
!
#ifdef USE_BUDA
         call balance_gpu()
        ! Use cuf kernels to compute final values
        !$cuf kernel do(1) <<<*,*>>>
         do k=kstart,kend
           my_dissipuc= my_dissipuc+ dissipuc_partial_d(k)
           my_dissipte= my_dissipte+ dissipte_partial_d(k)
         end do
        !$cuf kernel do(1) <<<*,*>>>
         do k=kstartr,kendr
           my_dissipur= my_dissipur+ dissipur_partial_d(k)
           my_dissipsa= my_dissipsa+ dissipsa_partial_d(k)
         end do
#else
      do k=kstart,kend
        kpp=k+1
        kmm=k-1
        if(k.eq.1)then
          udx3_m=1.d0/(zm(k)-zc(k))
        else
          udx3_m=1.d0/(zm(k)-zm(kmm))
        endif
        udx3_p=1.d0/(zm(kpp)-zm(k))
!$OMP  PARALLEL DO &
!$OMP  DEFAULT(SHARED) &
!$OMP  PRIVATE(j,i,jmm,jpp,imm,ipp) &
!$OMP  PRIVATE(h11,h12,h13) &
!$OMP  PRIVATE(h21,h22,h23) &
!$OMP  PRIVATE(h31,h32,h33) &
!$OMP  PRIVATE(h13p,h23p,h31p,h32p) &
!$OMP  PRIVATE(dissipuc,dissipte) &
!$OMP  REDUCTION(+: my_dissipuc,my_dissipte) &
!$OMP  REDUCTION(+: dissuc, disste)
        do j=1,n2m       
          jmm=jmv(j)
          jpp=jpv(j)
        
          do i=1,n1m
            imm= imv(i)
            ipp= ipv(i)

            h11=(q1(ipp,j,k)-q1(i,j,k))*udx1
            h22=(q2(i,jpp,k)-q2(i,j,k))*udx2
            h33=(q3(i,j,kpp)-q3(i,j,k))*udx3c(k)

            h12=(q1(i,j,k)-q1(i,jmm,k))*udx2
            h21=(q2(i,j,k)-q2(imm,j,k))*udx1

            h13=(q1(i,j,k)-q1(i,j,kmm))*udx3_m
            h31=(q3(i,j,k)-q3(imm,j,k))*udx1

            h13p=(q1(i,j,kpp)-q1(i,j,k))*udx3_p
            h31p=(q3(i,j,kpp)-q3(imm,j,kpp))*udx1

            h23=(q2(i,j,k)-q2(i,j,kmm))*udx3_m
            h32=(q3(i,j,k)-q3(i,jmm,k))*udx2

            h23p=(q2(i,j,kpp)-q2(i,j,k))*udx3_p
            h32p=(q3(i,j,kpp)-q3(i,jmm,kpp))*udx2

            dissipuc = 2.d0*(h11**2+h22**2+h33**2) &
                + (h21+h12)**2 &
                + ((h31+h13)**2+(h31p+h13p)**2)*0.5d0 &
                + ((h32+h23)**2+(h32p+h23p)**2)*0.5d0  
            my_dissipuc = my_dissipuc+dissipuc*g3rm(k)

            h31=(dens(ipp,j,k)-dens(i,j,k))*udx1
            h32=(dens(i,jpp,k)-dens(i,j,k))*udx2
            if(k.eq.n3m)then
              h33=(densbn(i,j)-dens(i,j,k))*udx3_p
            else
              h33=(dens(i,j,kpp)-dens(i,j,k))*udx3_p
            endif

            dissipte = h31*h31 + h32*h32 + h33*h33 
            my_dissipte = my_dissipte+dissipte*g3rc(k)

            dissuc(k) =  dissuc(k) + dissipuc * nu * un12m
            disste(k) =  disste(k) + dissipte * kpt * un12m

          end do
        end do
!$OMP  END PARALLEL DO
      end do
       
      udx1=dx1r
      udx2=dx2r

      un12m = 1.d0/(dble(n1mr)*dble(n2mr))

      do k=kstartr,kendr
        kpp=k+1
        kmm=k-1
        if(k.eq.1)then
          udx3_m=1.d0/(zmr(k)-zcr(k))
        else
          udx3_m=1.d0/(zmr(k)-zmr(kmm))
        endif
        udx3_p=1.d0/(zmr(kpp)-zmr(k))
!$OMP  PARALLEL DO &
!$OMP  DEFAULT(SHARED) &
!$OMP  PRIVATE(j,i,jmm,jpp,imm,ipp) &
!$OMP  PRIVATE(h11,h12,h13) &
!$OMP  PRIVATE(h21,h22,h23) &
!$OMP  PRIVATE(h31,h32,h33) &
!$OMP  PRIVATE(h13p,h23p,h31p,h32p) &
!$OMP  PRIVATE(dissipur,dissipsa) &
!$OMP  REDUCTION(+: my_dissipsa, my_dissipur) &
!$OMP  REDUCTION(+: dissur, disssa)
        do j=1,n2mr       
          jmm=jmvr(j)
          jpp=jpvr(j)
        
          do i=1,n1mr
            imm= imvr(i)
            ipp= ipvr(i)

            h11=(q1lr(ipp,j,k)-q1lr(i,j,k))*udx1
            h22=(q2lr(i,jpp,k)-q2lr(i,j,k))*udx2
            h33=(q3lr(i,j,kpp)-q3lr(i,j,k))*udx3cr(k)

            h12=(q1lr(i,j,k)-q1lr(i,jmm,k))*udx2
            h21=(q2lr(i,j,k)-q2lr(imm,j,k))*udx1

            h13=(q1lr(i,j,k)-q1lr(i,j,kmm))*udx3_m
            h31=(q3lr(i,j,k)-q3lr(imm,j,k))*udx1

            h13p=(q1lr(i,j,kpp)-q1lr(i,j,k))*udx3_p
            h31p=(q3lr(i,j,kpp)-q3lr(imm,j,kpp))*udx1

            h23=(q2lr(i,j,k)-q2lr(i,j,kmm))*udx3_m
            h32=(q3lr(i,j,k)-q3lr(i,jmm,k))*udx2

            h23p=(q2lr(i,j,kpp)-q2lr(i,j,k))*udx3_p
            h32p=(q3lr(i,j,kpp)-q3lr(i,jmm,kpp))*udx2

            dissipur = 2.d0*(h11**2+h22**2+h33**2) &
                 + (h21+h12)**2 &
                 + ((h31+h13)**2+(h31p+h13p)**2)*0.5d0 &
                 + ((h32+h23)**2+(h32p+h23p)**2)*0.5d0
            my_dissipur = my_dissipur+dissipur*g3rmr(k)

            h31=(dsal(ipp,j,k)-dsal(i,j,k))*udx1
            h32=(dsal(i,jpp,k)-dsal(i,j,k))*udx2
            if(k.eq.n3mr)then
              h33=(dsalbn(i,j)-dsal(i,j,k))*udx3_p
            else
              h33=(dsal(i,j,kpp)-dsal(i,j,k))*udx3_p
            endif 

            dissipsa = h31*h31 + h32*h32 + h33*h33
            my_dissipsa = my_dissipsa + dissipsa*g3rcr(k)

            dissur(k) =  dissur(k) + dissipur * nu * un12m
            disssa(k) =  disssa(k) + dissipsa * kps * un12m

          end do
        end do
!$OMP  END PARALLEL DO
      end do
#endif

      call MPI_REDUCE(my_dissipuc,nuuc,1,MDP,MPI_SUM,0,MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(my_dissipte,nute,1,MDP,MPI_SUM,0,MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(my_dissipur,nuur,1,MDP,MPI_SUM,0,MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(my_dissipsa,nusa,1,MDP,MPI_SUM,0,MPI_COMM_WORLD, ierr)
      
      volb = 1.d0/(alx3*dble(n3m)*dble(n1m)*dble(n2m))
      volr = 1.d0/(alx3*dble(n3mr)*dble(n1mr)*dble(n2mr))
      if(myid.eq.0) then
        nuuc = nuuc*volb
        nute = nute*volb
        nusa = nusa*volr
        nuur = nuur*volr
        write(92,520) time, nute, nusa, nuuc, &
          ( (nusa-1.d0)*Ras/Prs/Prs  &
           -(nute-1.d0)*Rat/Prt/Prt )*nu*nu/(alx3**4)
        write(98,521) time, nuuc, nuur, (nuur-nuuc)/nuuc
      endif
520   format(1x,f10.4,4(1x,ES20.8))
521   format(1x,f10.4,2(1x,ES20.8),f12.5)

      return   
      end
!
