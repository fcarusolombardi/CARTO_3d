module mlsForceTiledOff3Comp_m
  use constants
  contains

      ! Helper functions for computing weights
      real(DP) function get_wt_exp(val, wcon)
        implicit none
        real(DP), intent(in) :: val, wcon
        real(DP) :: rv

        if(val .le. 1.0d0)then
          rv = val / wcon
          get_wt_exp = exp(-rv * rv)
        else
          get_wt_exp = 0.d0
          print*, val
          write(*,*)'Something wrong in support domain-mlsForce'
        end if

        return
      end function get_wt_exp

      real(DP) function get_wt_cub(val)
        implicit none
        real(DP), intent(in) :: val

        if(val .le. 0.5d0) then
          get_wt_cub = (2.d0/3.d0)-4.d0*(val * val)+ &
                          4.d0*(val * val * val)
        elseif(val .le. 1.d0) then
          get_wt_cub = (4.d0/3.d0)*(1.d0-(val * val * val))- &
                    4.d0*(val-(val*val))
        else
          get_wt_cub = 0.d0
        end if

        return
      end function get_wt_cub

#ifdef USE_CUDA
      attributes(global) &
      subroutine get_valid_trianglesTile(trcnt, mytr, pindtOff, nttot, ntOffstart,ntOffend,ntOfftot, kstart, kend)
        implicit none

        integer, value, intent(in) :: nttot, kstart, kend
        integer, value, intent(in) :: ntOffstart,ntOffend,ntOfftot
        integer,  device, intent(inout) :: trcnt
        integer,  device, dimension(nttot), intent(out) :: mytr
        integer,  device, dimension(6, nttot), intent(in) :: pindtOff
        integer ::i, ind_pal, ntr,ioff

        i = (blockIdx%x - 1) * blockDim%x + threadIdx%x

        if (i > ntOfftot) return
        ioff = i + ntOffstart - 1
        ind_pal = pindtOff(6, ioff)

        if(ind_pal .ge. kstart.and. ind_pal .le. kend-1) then
          mytr(atomicAdd(trcnt, 1) + 1) = ioff
        endif

      end subroutine get_valid_trianglesTile


      attributes(device) &
      real(DP) function get_wt_exp_gpu(val, wcon)
        implicit none
        real(DP), intent(in), value :: val, wcon
        real(DP) :: rv

        if(val .le. 1.0d0)then
          rv = val / wcon
          get_wt_exp_gpu = exp(-rv * rv)
        else
          get_wt_exp_gpu = 0.d0
        end if

        return
      end function get_wt_exp_gpu

      attributes(device) &
      real(DP) function get_wt_cub_gpu(val)
        implicit none
        real(DP), intent(in), value :: val

        if(val .le. 0.5d0) then
          get_wt_cub_gpu = (2.d0/3.d0)-4.d0*(val * val)+ &
                          4.d0*(val * val * val)
        elseif(val .le. 1.d0) then
          get_wt_cub_gpu = (4.d0/3.d0)*(1.d0-(val * val * val))- &
                    4.d0*(val-(val*val))
        else
          get_wt_cub_gpu = 0.d0
        end if

        return
      end function get_wt_cub_gpu

      attributes(device) &
      subroutine det4_gpu(a, det_l)
       implicit none
       real(DP), device, intent(in) :: a(4,4)
       real(DP), device, intent(out) :: det_l

       det_l =  a(1,1)*(a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3)) &
                  -a(2,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2)) &
                  +a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2)))&
          -a(1,2)*(a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3)) &
                  -a(2,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1)) &
                  +a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))&
          +a(1,3)*(a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2)) &
                  -a(2,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1)) &
                  +a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))&
          -a(1,4)*(a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2)) &
                  -a(2,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)) &
                  +a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))

      return
      end subroutine det4_gpu

      attributes(device) &
      subroutine invert4_unscaled_gpu(a,ia)
       implicit none
       real(DP), device, intent(in) :: a (4,4)
       real(DP), device, intent(out) :: ia (4,4)

       ia(1,1) =  (a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)* &
      (a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2)))
       ia(2,1) = (-a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(2,3)* &
      (a(3,1)*a(4,4)-a(3,4)*a(4,1))-a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))
       ia(3,1) =  (a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(2,2)* &
      (a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))
       ia(4,1) = (-a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))+a(2,2)* &
      (a(3,1)*a(4,3)-a(3,3)*a(4,1))-a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))

       ia(1,2) = (-a(1,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(1,3)* &
      (a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(1,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2)))
       ia(2,2) =  (a(1,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(1,3)* &
      (a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(1,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))
       ia(3,2) = (-a(1,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(1,2)* &
      (a(3,1)*a(4,4)-a(3,4)*a(4,1))-a(1,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))
       ia(4,2) =  (a(1,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))-a(1,2)* &
      (a(3,1)*a(4,3)-a(3,3)*a(4,1))+a(1,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))

       ia(1,3) =  (a(1,2)*(a(2,3)*a(4,4)-a(2,4)*a(4,3))-a(1,3)* &
      (a(2,2)*a(4,4)-a(2,4)*a(4,2))+a(1,4)*(a(2,2)*a(4,3)-a(2,3)*a(4,2)))
       ia(2,3) = (-a(1,1)*(a(2,3)*a(4,4)-a(2,4)*a(4,3))+a(1,3)* &
      (a(2,1)*a(4,4)-a(2,4)*a(4,1))-a(1,4)*(a(2,1)*a(4,3)-a(2,3)*a(4,1)))
       ia(3,3) =  (a(1,1)*(a(2,2)*a(4,4)-a(2,4)*a(4,2))-a(1,2)* &
      (a(2,1)*a(4,4)-a(2,4)*a(4,1))+a(1,4)*(a(2,1)*a(4,2)-a(2,2)*a(4,1)))
       ia(4,3) = (-a(1,1)*(a(2,2)*a(4,3)-a(2,3)*a(4,2))+a(1,2)* &
      (a(2,1)*a(4,3)-a(2,3)*a(4,1))-a(1,3)*(a(2,1)*a(4,2)-a(2,2)*a(4,1)))

       ia(1,4) = (-a(1,2)*(a(2,3)*a(3,4)-a(2,4)*a(3,3))+a(1,3)* &
      (a(2,2)*a(3,4)-a(2,4)*a(3,2))-a(1,4)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)))
       ia(2,4) =  (a(1,1)*(a(2,3)*a(3,4)-a(2,4)*a(3,3))-a(1,3)* &
      (a(2,1)*a(3,4)-a(2,4)*a(3,1))+a(1,4)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)))
       ia(3,4) = (-a(1,1)*(a(2,2)*a(3,4)-a(2,4)*a(3,2))+a(1,2)* &
      (a(2,1)*a(3,4)-a(2,4)*a(3,1))-a(1,4)*(a(2,1)*a(3,2)-a(2,2)*a(3,1)))
       ia(4,4) =  (a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)* &
      (a(2,1)*a(3,3)-a(2,3)*a(3,1))+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1)))

      return
      end subroutine invert4_unscaled_gpu

      attributes(device) &
      subroutine invert4_gpu(a,ia)
       implicit none
       real(DP), device, intent(in) :: a(4,4)
       real(DP), device, intent(out) :: ia(4,4)
       real(DP) :: det_l, invdet

       det_l =  a(1,1)*(a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3)) &
                  -a(2,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2)) &
                  +a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2)))&
          -a(1,2)*(a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3)) &
                  -a(2,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1)) &
                  +a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))&
          +a(1,3)*(a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2)) &
                  -a(2,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1)) &
                  +a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))&
          -a(1,4)*(a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2)) &
                  -a(2,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)) &
                  +a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))

      invdet = 1.0d0/det_l

       ia(1,1) =  (a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)* &
      (a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))) * invdet
       ia(2,1) = (-a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(2,3)* &
      (a(3,1)*a(4,4)-a(3,4)*a(4,1))-a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))) * invdet
       ia(3,1) =  (a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(2,2)* &
      (a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))) * invdet
       ia(4,1) = (-a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))+a(2,2)* &
      (a(3,1)*a(4,3)-a(3,3)*a(4,1))-a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))) * invdet

       ia(1,2) = (-a(1,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(1,3)* &
      (a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(1,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))) * invdet
       ia(2,2) =  (a(1,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(1,3)* &
      (a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(1,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))) * invdet
       ia(3,2) = (-a(1,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(1,2)* &
      (a(3,1)*a(4,4)-a(3,4)*a(4,1))-a(1,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))) * invdet
       ia(4,2) =  (a(1,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))-a(1,2)* &
      (a(3,1)*a(4,3)-a(3,3)*a(4,1))+a(1,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))) * invdet

       ia(1,3) =  (a(1,2)*(a(2,3)*a(4,4)-a(2,4)*a(4,3))-a(1,3)* &
      (a(2,2)*a(4,4)-a(2,4)*a(4,2))+a(1,4)*(a(2,2)*a(4,3)-a(2,3)*a(4,2))) * invdet
       ia(2,3) = (-a(1,1)*(a(2,3)*a(4,4)-a(2,4)*a(4,3))+a(1,3)* &
      (a(2,1)*a(4,4)-a(2,4)*a(4,1))-a(1,4)*(a(2,1)*a(4,3)-a(2,3)*a(4,1))) * invdet
       ia(3,3) =  (a(1,1)*(a(2,2)*a(4,4)-a(2,4)*a(4,2))-a(1,2)* &
      (a(2,1)*a(4,4)-a(2,4)*a(4,1))+a(1,4)*(a(2,1)*a(4,2)-a(2,2)*a(4,1))) * invdet
       ia(4,3) = (-a(1,1)*(a(2,2)*a(4,3)-a(2,3)*a(4,2))+a(1,2)* &
      (a(2,1)*a(4,3)-a(2,3)*a(4,1))-a(1,3)*(a(2,1)*a(4,2)-a(2,2)*a(4,1))) * invdet

       ia(1,4) = (-a(1,2)*(a(2,3)*a(3,4)-a(2,4)*a(3,3))+a(1,3)* &
      (a(2,2)*a(3,4)-a(2,4)*a(3,2))-a(1,4)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))) * invdet
       ia(2,4) =  (a(1,1)*(a(2,3)*a(3,4)-a(2,4)*a(3,3))-a(1,3)* &
      (a(2,1)*a(3,4)-a(2,4)*a(3,1))+a(1,4)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))) * invdet
       ia(3,4) = (-a(1,1)*(a(2,2)*a(3,4)-a(2,4)*a(3,2))+a(1,2)* &
      (a(2,1)*a(3,4)-a(2,4)*a(3,1))-a(1,4)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))) * invdet
       ia(4,4) =  (a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)* &
      (a(2,1)*a(3,3)-a(2,3)*a(3,1))+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))) * invdet


      return
      end subroutine invert4_gpu


      attributes(device) &
      subroutine inverseLU_gpu(a,c)
       implicit none
       real(DP), device :: a (4,4)
       real(DP), device, intent(out) :: c (4,4)
       real(DP) :: det_l, invdet
       real(DP) :: L(4,4), U(4,4)
       real(DP) :: b(4), d(4), x(4)
       real(DP) :: coeff
       integer i,j,k

      L = 0.0 ; U = 0.0 ; b = 0.0

      !forward elimination
      do k=1,3
       do i=k+1,4
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,4
               a(i,j)=a(i,j)-coeff*a(k,j)
            end do
               
       end do
    end do
          !prepare L U
      do i=1,4
       L(i,i) = 1.0
      end do
      do j=1,4
       do i=1,j
            U(i,j) = a(i,j)
       end do
      end do

      !compute columns of c
      do k=1,4
       b(k)=1.0
       d(1)=b(1)
        do i=2,4
        d(i)=b(i)
         do j=1,i-1
            d(i) = d(i)-L(i,j)*d(j)
!            write(*,*) 'f',d(i)
         end do
        end do
       !solve ux=d with back subs.
       x(4)=d(4)/U(4,4)
       do i=3,1,-1
        x(i) = d(i)
        do j=4,i+1,-1
          x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
       end do
      !fill solns of x(n) to k of C
      do i=1,4
         c(i,k)=x(i)
      end do
            b(k) = 0.0
      end do
      
      return
      end subroutine inverseLU_gpu

      
      attributes(global) &
      subroutine update_forces_kernel_1cmp(cmp, siz, str1, str2, &
        nel, dx1, dx2, dx3, wscl, wcon, wcub, wexp, n1m, &
        n2m, for_xc, for_yc, for_zc, pindtOff, tri_ver, tri_vel,vel_tri, albegaBar, &
        faceid_t, tri_tiling, tstart, mytl, tlcntOff, q1, q2, q3, &
        sur, xm, ym, zm, xc, yc, zc, Hboxx, celvol, nttot, nftot, rigamax, Navamax, n1, n2, n3, ks, ke, lvlhalo, &
        face_to_part,face_to_chamb,body_nsl, body_flux,dOff,tri_nor,dt)

        implicit none

        integer, value, intent(in) :: cmp, siz, str1, str2, n1m, n2m, nftot, nttot, rigamax, Navamax
        integer, value, intent(in) :: nel, wcub, wexp, n1, n2, n3, ks, ke, lvlhalo
        real(DP), value, intent(in) :: dx1, dx2, dx3, wscl, wcon, dOff,dt
        integer, device, dimension(1), intent(in) :: tlcntOff

        integer, device, dimension(*), intent(in) :: body_nsl
        real(DP), device, dimension(*), intent(in) :: body_flux

        integer, device, dimension(*), intent(in) :: face_to_part
        integer, device, dimension(*), intent(in) :: face_to_chamb
        integer, device, dimension(6, nttot), intent(in) :: pindtOff
        real(DP), device, dimension(3, nftot), intent(in) :: vel_tri, tri_nor
        real(DP), device, dimension(9, nftot), intent(in) :: tri_ver,tri_vel
        real(DP), device, dimension(rigamax, Navamax, 3), intent(in) :: albegaBar
        real(DP), device, dimension(nftot), intent(in) :: sur
        real(DP), device, dimension(n1, n2, ks-lvlhalo:ke+lvlhalo), intent(in) :: q1, q2, q3
        real(DP), device, dimension(*), intent(in) :: ym, zm, xc, yc, zc
        real(DP), device, dimension(0:n1), intent(in) :: xm
        real(DP), device, dimension(0:n3), intent(in) :: Hboxx, celvol
        integer, device, dimension(nftot), intent(in) :: tstart
        integer, device, dimension(nftot,2), intent(in) :: tri_tiling
        integer, device, dimension(nttot), intent(in) :: mytl, faceid_t
        real(DP), device, dimension(n1, n2, ks-lvlhalo:ke+lvlhalo-1), intent(inout) :: for_xc, for_yc, for_zc

        integer :: ci, cj, ck
        integer :: i, j, k, inw
        integer :: i1, j1, k1
        real(DP) :: Offx,Offy,Offz
        real(DP), device :: pos_MLS(3), vel_MLS(3), norp(3), el_cx(3), ptx(4), Hbox(3), Wt(3)
        real(DP), device :: ui(3)
        real(DP), device :: B(4)
        real(DP) :: Wtx, rv1, rv2, cfac, epsw
        real(DP) :: Forx, Fory, Forz

        real(DP), shared :: pinvA(16,4), invA(16,4)
        !real(8), device :: pinvA(16), invA(16)
        real(DP), shared :: ptxA(4,4)
        real(DP), shared :: det_l(4)

        integer :: tid, tx, ty, tl, tile, ntr, inp,chamb
        integer :: rowtile, tilex

        !tid = (blockIdx%x - 1) * blockDim%x + threadIdx%x
        !tx = mod(tid-1, 32) + 1
        !ty = (tx - 1)/32 + 1
        tx = threadIdx%x
        ty = threadIdx%y
        tl = (blockIdx%x - 1) * blockDim%y + threadIdx%y

        if (tl > tlcntOff(1)) return

        tile = mytl(tl)
        ntr = faceid_t(tile)

        inp   = face_to_part(ntr)
        chamb = face_to_chamb(ntr)

        Offx = dOff*tri_nor(1,ntr)
        Offy = dOff*tri_nor(2,ntr)
        Offz = dOff*tri_nor(3,ntr)

        rowtile = tri_tiling(ntr,2)
        tilex = tile - tstart(ntr) + 1


        select case (cmp)
        case (1)
          ci = pindtOff(1,tile)
          cj = pindtOff(5,tile)
          ck = pindtOff(6,tile)
        case (2)
          ci = pindtOff(4,tile)
          cj = pindtOff(2,tile)
          ck = pindtOff(6,tile)
        case (3)
          ci = pindtOff(4,tile)
          cj = pindtOff(5,tile)
          ck = pindtOff(3,tile)
        end select

        i1 = mod(tx - 1, str1) + 1
        j1 = mod((tx - 1) / str1, str1) + 1
        k1 = (tx - 1) / str2 + 1

        inw = 1 + (i1 - 1) + (j1 - 1) * str1 + (k1 - 1) * str2
        i = ci - (siz + 1)/2 + i1
        j = cj - (siz + 1)/2 + j1
        k = ck - (siz + 1)/2 + k1


        pos_MLS(1) = albegaBar(rowtile,tilex,1) * tri_ver(1, ntr) + &
                     albegaBar(rowtile,tilex,2) * tri_ver(4, ntr) + &
                     albegaBar(rowtile,tilex,3) * tri_ver(7, ntr) +Offx

        pos_MLS(2) = albegaBar(rowtile,tilex,1) * tri_ver(2, ntr) + &
                     albegaBar(rowtile,tilex,2) * tri_ver(5, ntr) + &
                     albegaBar(rowtile,tilex,3) * tri_ver(8, ntr) +Offy

        pos_MLS(3) = albegaBar(rowtile,tilex,1) * tri_ver(3, ntr) + &
                     albegaBar(rowtile,tilex,2) * tri_ver(6, ntr) + &
                     albegaBar(rowtile,tilex,3) * tri_ver(9, ntr) +Offz

        vel_MLS(1) = albegaBar(rowtile,tilex,1) * tri_vel(1, ntr) + &
                     albegaBar(rowtile,tilex,2) * tri_vel(4, ntr) + &
                     albegaBar(rowtile,tilex,3) * tri_vel(7, ntr)

        vel_MLS(2) = albegaBar(rowtile,tilex,1) * tri_vel(2, ntr) + &
                     albegaBar(rowtile,tilex,2) * tri_vel(5, ntr) + &
                     albegaBar(rowtile,tilex,3) * tri_vel(8, ntr)

        vel_MLS(3) = albegaBar(rowtile,tilex,1) * tri_vel(3, ntr) + &
                     albegaBar(rowtile,tilex,2) * tri_vel(6, ntr) + &
                     albegaBar(rowtile,tilex,3) * tri_vel(9, ntr)


  !     --------------------------------------------------------

        Hbox(1) = dx1
        Hbox(2) = dx2
        Hbox(3) = dx3

  !     --------------------------------------------------------

  !     spline weights for all three components
        Wtx = 0
        epsw = 1e-18

        if (tx <= nel) then
          select case (cmp)
          case(1)
            el_cx(1) = abs(xc(i) - pos_MLS(1))! + epsw)
            el_cx(2) = abs(ym(j) - pos_MLS(2))! + epsw)
            el_cx(3) = abs(zm(k) - pos_MLS(3))! + epsw)
          case(2)
            el_cx(1) = abs(xm(i) - pos_MLS(1))! + epsw)
            el_cx(2) = abs(yc(j) - pos_MLS(2))! + epsw)
            el_cx(3) = abs(zm(k) - pos_MLS(3))! + epsw)
          case(3)
            el_cx(1) = abs(xm(i) - pos_MLS(1))! + epsw)
            el_cx(2) = abs(ym(j) - pos_MLS(2))! + epsw)
            el_cx(3) = abs(zc(k) - pos_MLS(3))! + epsw)
          end select

          norp(1) = el_cx(1) * Hbox(1) / wscl
          norp(2) = el_cx(2) * Hbox(2) / wscl
          norp(3) = el_cx(3) * Hbox(3) / wscl

  !           ----------------CUBIC SPLINES---------------------
          if(wcub .eq. 1) then
            Wt(1) = get_wt_cub_gpu(norp(1))
            Wt(2) = get_wt_cub_gpu(norp(2))
            Wt(3) = get_wt_cub_gpu(norp(3))
  !           ----------------EXPONENTIAL SPLINES------------------
          else if(wexp .eq. 1) then
            Wt(1) = get_wt_exp_gpu(norp(1), wcon)
            Wt(2) = get_wt_exp_gpu(norp(2), wcon)
            Wt(3) = get_wt_exp_gpu(norp(3), wcon)
          end if

  !     ------------------------------------------------

          Wtx = Wt(1)*Wt(2)*Wt(3)
        endif

        ui(1) = 0.d0
        !ui(2) = 0.d0
        !ui(3) = 0.d0

        if (tx <= nel) then
          select case (cmp)
          case (1)
            ui(1) = q1(i,j,k)
          case (2)
            ui(1) = q2(i,j,k)
          case (3)
            ui(1) = q3(i,j,k)
          end select
        endif

        ptx(1) = 0.d0
        ptx(2) = 0.d0
        ptx(3) = 0.d0
        ptx(4) = 0.d0
        B(1) = 0.d0
        B(2) = 0.d0
        B(3) = 0.d0
        B(4) = 0.d0

        ! set pinvA and B
        if (tx <= nel) then
          select case (cmp)
          case (1)
            ptx(1) = 1.0d0
            ptx(2) = xc(i)
            ptx(3) = ym(j)
            ptx(4) = zm(k)
          case (2)
            ptx(1) = 1.0d0
            ptx(2) = xm(i)
            ptx(3) = yc(j)
            ptx(4) = zm(k)
          case (3)
            ptx(1) = 1.0d0
            ptx(2) = xm(i)
            ptx(3) = ym(j)
            ptx(4) = zc(k)
          end select

          B(1) = Wtx * ptx(1)
          B(2) = Wtx * ptx(2)
          B(3) = Wtx * ptx(3)
          B(4) = Wtx * ptx(4)
        endif

        ! each thread does outer product
        do j1 = 1, 4
          do i1 = 1, 4
            rv1 = Wtx*ptx(i1)*ptx(j1)

            ! reduce across warp
            rv2 = __shfl_down(rv1,1)
            rv1 = rv1 +  rv2
            rv2 = __shfl_down(rv1,2)
            rv1 = rv1 +  rv2
            rv2 = __shfl_down(rv1,4)
            rv1 = rv1 +  rv2
            rv2 = __shfl_down(rv1,8)
            rv1 = rv1 +  rv2
            rv2 = __shfl_down(rv1,16)
            rv1 = rv1 +  rv2

            !rv2 = __shfl_xor(rv1,1)
            !rv1 = rv1 +  rv2
            !rv2 = __shfl_xor(rv1,2)
            !rv1 = rv1 +  rv2
            !rv2 = __shfl_xor(rv1,4)
            !rv1 = rv1 +  rv2
            !rv2 = __shfl_xor(rv1,8)
            !rv1 = rv1 +  rv2
            !rv2 = __shfl_xor(rv1,16)
            !rv1 = rv1 +  rv2


            if (tx == 1) then
              !pinvA(i1,j1,ty) = rv1
              pinvA(i1 + (j1-1)*4, ty) = rv1
            endif

            !call syncthreads()
            !pinvA(i1,j1) = rv1
            !pinvA(i1 + (j1-1)*4) = rv1

          enddo
        enddo

        call syncthreads()

        ! Compute ptxA (using 16 threads)
        if (tx == 1) then
          ! call invert4_unscaled_gpu(pinvA(1,ty), invA(1,ty))
          ! call det4_gpu(pinvA(1,ty), det_l(ty))
#ifdef LUINVERSE                      
           call inverseLU_gpu(pinvA(1,ty), invA(1,ty))
#else
           call invert4_gpu(pinvA(1,ty), invA(1,ty))
#endif           
        endif

        call syncthreads()

        if (tx <= 16) then
          !call det4_gpu(pinvA(1,ty), det_l)
          !call invert4_unscaled_gpu(pinvA, invA)

          ptx(1) = 1.d0;
          ptx(2) = pos_MLS(1)
          ptx(3) = pos_MLS(2)
          ptx(4) = pos_MLS(3)

          rv1 = ptx(mod(tx-1, 4) + 1) * invA(tx,ty) !/ det_l(ty)
        endif

        rv2 = __shfl_down(rv1,1)
        rv1 = rv1 +  rv2
        rv2 = __shfl_down(rv1,2)
        rv1 = rv1 +  rv2

        if (tx <= 16) then
          if (mod(tx-1, 4) == 0) then
            ptxA((tx-1)/4 + 1, ty) = rv1
          endif
        endif

        call syncthreads()

        ! Compute ptxAB, store in B(1)
        if (tx <= nel) then
          rv1 = 0.d0
          do i1 = 1, 4
            rv1 = rv1 + ptxA(i1, ty) * B(i1)
          enddo
          B(1) = rv1
        endif

        ! Compute um, store in ui
        do i1 = 1, 1
          ! JR note: if statement unneeded. ui set to zero for tx > nel earlier
          !if (tx <= nel) then
            rv1 = B(1) * ui(i1)
          !else
          !  rv1 = 0.d0
          !endif

          rv2 = __shfl_xor(rv1,1)
          rv1 = rv1 +  rv2
          rv2 = __shfl_xor(rv1,2)
          rv1 = rv1 +  rv2
          rv2 = __shfl_xor(rv1,4)
          rv1 = rv1 +  rv2
          rv2 = __shfl_xor(rv1,8)
          rv1 = rv1 +  rv2
          rv2 = __shfl_xor(rv1,16)
          rv1 = rv1 +  rv2

          ui(i1) = rv1

        enddo

        ! compute mvol_t
        if (tx <= nel) then
          rv1 = B(1) * Hboxx(k)
        else
          rv1 = 0.d0
        endif

        rv2 = __shfl_xor(rv1,1)
        rv1 = rv1 +  rv2
        rv2 = __shfl_xor(rv1,2)
        rv1 = rv1 +  rv2
        rv2 = __shfl_xor(rv1,4)
        rv1 = rv1 +  rv2
        rv2 = __shfl_xor(rv1,8)
        rv1 = rv1 +  rv2
        rv2 = __shfl_xor(rv1,16)
        rv1 = rv1 +  rv2

        cfac = rv1*sur(ntr)

        ! compute tcelvol
        if (tx <= nel) then
          rv1 = B(1) * celvol(k)
        else
          rv1 = 0.d0
        endif

        rv2 = __shfl_xor(rv1,1)
        rv1 = rv1 +  rv2
        rv2 = __shfl_xor(rv1,2)
        rv1 = rv1 +  rv2
        rv2 = __shfl_xor(rv1,4)
        rv1 = rv1 +  rv2
        rv2 = __shfl_xor(rv1,8)
        rv1 = rv1 +  rv2
        rv2 = __shfl_xor(rv1,16)
        rv1 = rv1 +  rv2

        cfac = cfac / (rv1 * dble(tri_tiling(ntr,1)))

    !   For moving particles
!        Forx = B(1) * cfac * blend * (vel_tri(cmp,ntr) - Ui(1))
        Forx = B(1) * cfac * (vel_MLS(cmp) - Ui(1)) !FV
        !Fory = B(1) * cfac * blend * (vel_tri(2,ntr) - Ui(2))
        !Forz = B(1) * cfac * blend * (vel_tri(3,ntr) - Ui(3))

        if ( (chamb.EQ.5).OR.(chamb.EQ.2) ) then
        ! Forx = cfac * tri_nor(cmp,ntr)
           Forx = 0.0d0
        endif!condizione inp e chamb 

        select case (cmp)
        case (1)
          rv1 = atomicadd(for_xc(i,j,k), Forx)
        case (2)
          rv1 = atomicadd(for_yc(i,j,k), Forx)
        case (3)
          rv1 = atomicadd(for_zc(i,j,k), Forx)
        end select

      end subroutine
#endif

end module


!-----------------------------------------------------------------------
!     read in position of a individual trial marker for MLS and compute
!     compute support domain, shape function and interpolate
!------------------------------------------------------------------------

      subroutine mlsForceTiledOff3Comp
      USE mpih
      USE param
      USE mls_param
      USE local_arrays, only: q1,q2,q3
      USE mpi_param, only: kstart, kend
      USE mls_local, only: for_xc, for_yc, for_zc, coll
      USE tile_arrays
      USE ibm_param
      USE ieee_arithmetic
      USE mlsForceTiledOff3Comp_m
!@cuf USE cudafor
      IMPLICIT NONE

      integer :: inp,seed,merr,ntr,chamb
      real(DP) :: pos_MLS(3),vel_MLS(3)
      integer :: siz, stride(2),tlcntOff(1)

      integer :: inw,i,j,k,ii,jj,kk,cmp,ind_pal
      integer :: i1, j1, k1,ist,jst,kst
      integer :: isto,isti,jsto,jsti,ksto,ksti
      real(DP) :: norp(3),el_cx(3),Offxyz(3)
      real(DP) :: elmag,norpd,normd,epsw
      real(DP) :: pinvA(4,4),invA(4,4),B(4,nel),pxk(4)
      real(DP) :: ui(nel,3),Wt(3),Wtx(nel),Hbox(3)
      real(DP) :: ptx(4)
      real(DP) :: ptxA(4),ptxAB(nel)
!     --------------------------------------------------
      real(DP) :: Um(3)
      real(DP) :: tcelvol,cfac,cel_fx,cel_fy,cel_fz
      real(DP) :: Forx,Fory,Forz
      real(DP) :: post1,post2,post3
      real(DP) :: Offx,Offy,Offz
      real(DP) :: P1(3), P2(3), P3(3)
      real(DP) :: vP1(3), vP2(3), vP3(3)

      integer :: ci,cj,ck
      integer :: tile, rowtile, tilex
      real(DP) :: mvol_t

      integer :: tr, tl
#ifdef USE_CUDA
      type(dim3) :: blocks, threads
      attributes(managed) :: tlcntOff
#endif
!@cuf integer :: istat


!FV This is has been moved to the routine findindicesTiledOff
! !FINDINDICETILEDOFF
! #ifdef USE_CUDA
!       !$cuf kernel do (1)
! #endif
!       do tile = ntOffstart,ntOffend

!         ntr = faceid_t(tile)

!         Offx = dOff*tri_nor(1,ntr)
!         Offy = dOff*tri_nor(2,ntr)
!         Offz = dOff*tri_nor(3,ntr)

!         rowtile = tri_tiling(ntr,2)

!         tilex = tile - tstart(ntr) + 1

!         post1 = albegaBar(rowtile,tilex,1)*tri_ver(1,ntr) + &
!                 albegaBar(rowtile,tilex,2)*tri_ver(4,ntr) + &
!                 albegaBar(rowtile,tilex,3)*tri_ver(7,ntr) + Offx
!         post2 = albegaBar(rowtile,tilex,1)*tri_ver(2,ntr) + &
!                 albegaBar(rowtile,tilex,2)*tri_ver(5,ntr) + &
!                 albegaBar(rowtile,tilex,3)*tri_ver(8,ntr) + Offy
!         post3 = albegaBar(rowtile,tilex,1)*tri_ver(3,ntr) + &
!                 albegaBar(rowtile,tilex,2)*tri_ver(6,ntr) + &
!                 albegaBar(rowtile,tilex,3)*tri_ver(9,ntr) + Offz



! !       ++++++++Indices of the marker++++++++++++++++++++++++
! !       X - indices
!         i1=NINT((post1+xyz_tr(1))*dx1) +1
!         ist=FLOOR((post1+xyz_tr(1))*dx1) +1      
!         if(ist.eq.0)ist=n1m
! !       Y - indices
!         j1=NINT((post2+xyz_tr(2))*dx2) +1
!         jst=FLOOR((post2+xyz_tr(2))*dx2) + 1
!         if(jst.eq.0)jst=n2m
! !       Z - indices
! !       uniform grid : wall-normal direction
!         if(istr3.eq.0)then
!         k1=NINT((post3+xyz_tr(3))*dx3) +1
!         kst=FLOOR((post3+xyz_tr(3))*dx3) + 1
!         endif
! !       OUTER and INNER INDICES
!         isto = ist+1 ; isti = ist-1
!         jsto = jst+1 ; jsti = jst-1
!         ksto = kst+1 ; ksti = kst-1
! !       -------------------------------------------------------------
!         pindtOff(1,tile)=i1 ; pindtOff(2,tile)=j1 ; pindtOff(3,tile)=k1
!         pindtOff(4,tile)=ist; pindtOff(5,tile)=jst; pindtOff(6,tile)=kst
! !       -------------------------------------------------------------
!       end do

! !@cuf istat = cudaDeviceSynchronize !JDR TMP
! !ENDFINDINDICESTILED
!     --------------------------------------------------------

      ! Setup shift and strides based on nel
      if (nel .eq. 27) then
        siz = 3
        stride(1) = 3
        stride(2) = 9
      elseif (nel .eq. 125) then
        siz = 5
        stride(1) = 5
        stride(2) = 25
      endif

!      Main loop
      ! precompute valid tiles
#ifdef USE_CUDA
      tlcntOff = 0
      threads = dim3(128, 1, 1)
      blocks = dim3(ceiling(real(ntOfftot)/128), 1, 1)
!      ! Can reuse get_valid_triangles kernel
      call get_valid_trianglesTile<<<blocks, threads>>>(tlcntOff(1), mytl, pindtOff, &
                                                    nttot,ntOffstart,ntOffend,ntOfftot, &
                                                    kstart, kend)
#else
      tl = 1
      tlcntOff = 0
      do tile = ntOffstart, ntOffend
        ntr = faceid_t(tile)
        ind_pal = pindtOff(6,tile)
        if(ind_pal .ge. kstart .and. ind_pal .le. kend-1) then
          mytl(tl) = tile
          tl = tl + 1
          tlcntOff(1) = tlcntOff(1) + 1
        endif
      enddo
#endif

!@cuf istat = cudaDeviceSynchronize !JDR TMP
#ifdef USE_CUDA

      if (nel .ne. 27) then
        print*, "GPU implementation of mlsForces only supports nel == 27 currently!"
        flush(6); stop
      endif

      threads = dim3(32, 4, 1)
      blocks = dim3(ceiling(real(ntOfftot)/4), 1, 1) !controlla FV

      do cmp = 1, 3
        call update_forces_kernel_1cmp<<<blocks, threads>>>(cmp, siz, stride(1), stride(2), &
        nel, dx1, dx2, dx3, wscl, wcon, wcub, wexp, &
        n1m, n2m, for_xc, for_yc, for_zc, pindtOff, tri_ver, tri_vel,vel_tri, albegaBar, &
        faceid_t, tri_tiling, tstart, mytl, tlcntOff, &
        q1, q2, q3, sur, xm, ym, zm, xc, yc, zc, Hboxx, celvol, &
        nttot, nftot, size(albegaBar, 1), size(albegaBar, 2), n1, n2, n3, kstart, kend, lvlhalo, & 
        face_to_part,face_to_chamb,body_nsl, body_flux,dOff,tri_nor,dt)
      end do

#else

      do tl = 1, tlcntOff(1)
        tile = mytl(tl)
        ntr = faceid_t(tile)

        inp   = face_to_part(ntr)
        chamb = face_to_chamb(ntr)

        Offxyz(1:3) = dOff*tri_nor(1:3,ntr)

        rowtile = tri_tiling(ntr,2)
        tilex = tile - tstart(ntr) + 1


        P1(1:3) = tri_ver(1:3,ntr)
        P2(1:3) = tri_ver(4:6,ntr)
        P3(1:3) = tri_ver(7:9,ntr)

        pos_MLS(1:3) = albegaBar(rowtile,tilex,1) * P1(1:3) + &
                       albegaBar(rowtile,tilex,2) * P2(1:3) + &
                       albegaBar(rowtile,tilex,3) * P3(1:3) + &
                       Offxyz(1:3)

        ! vP1(1:3) = tri_vel(l,1:3)
        ! vP2(1:3) = tri_vel(l,4:6)
        ! vP3(1:3) = tri_vel(l,7:9)
        vP1(1:3) = tri_vel(1:3,ntr)
        vP2(1:3) = tri_vel(4:6,ntr)
        vP3(1:3) = tri_vel(7:9,ntr)

        vel_MLS(1:3) = albegaBar(rowtile,tilex,1) * vP1(1:3) + &
                       albegaBar(rowtile,tilex,2) * vP2(1:3) + &
                       albegaBar(rowtile,tilex,3) * vP3(1:3)

    !   initialise pre-factor matrix

        ptx(1) = 1.d0;
        ptx(2) = pos_MLS(1)
        ptx(3) = pos_MLS(2)
        ptx(4) = pos_MLS(3)

    !   --------------------------------------------------------

        Hbox(1) = dx1
        Hbox(2) = dx2
        Hbox(3) = dx3

        ! Loop over 3 components separately
        do cmp = 1,3
          select case (cmp)
          case (1)
            ci = pindtOff(1,tile)
            cj = pindtOff(5,tile)
            ck = pindtOff(6,tile)
          case (2)
            ci = pindtOff(4,tile)
            cj = pindtOff(2,tile)
            ck = pindtOff(6,tile)
          case (3)
            ci = pindtOff(4,tile)
            cj = pindtOff(5,tile)
            ck = pindtOff(3,tile)
          end select
    !     --------------------------------------------------------

    !     spline weights for all three components
          epsw = 1e-18
          do k1 = 1, siz
            do j1 = 1, siz
              do i1 = 1, siz
                inw = 1 + (i1 - 1) + (j1 - 1) * stride(1) + (k1 - 1) * stride(2)
                i = ci - (siz + 1)/2 + i1
                j = cj - (siz + 1)/2 + j1
                k = ck - (siz + 1)/2 + k1

                select case (cmp)
                case (1)
                  el_cx(1) = abs(xc(i) - pos_MLS(1))! + epsw)
                  el_cx(2) = abs(ym(j) - pos_MLS(2))! + epsw)
                  el_cx(3) = abs(zm(k) - pos_MLS(3))! + epsw)
                case (2)
                  el_cx(1) = abs(xm(i) - pos_MLS(1))! + epsw)
                  el_cx(2) = abs(yc(j) - pos_MLS(2))! + epsw)
                  el_cx(3) = abs(zm(k) - pos_MLS(3))! + epsw)
                case (3)
                  el_cx(1) = abs(xm(i) - pos_MLS(1))! + epsw)
                  el_cx(2) = abs(ym(j) - pos_MLS(2))! + epsw)
                  el_cx(3) = abs(zc(k) - pos_MLS(3))! + epsw)
                end select

                norp(1) = el_cx(1) * Hbox(1) / wscl
                norp(2) = el_cx(2) * Hbox(2) / wscl
                norp(3) = el_cx(3) * Hbox(3) / wscl

    !           ----------------CUBIC SPLINES---------------------
                if(wcub .eq. 1) then
                  Wt(1) = get_wt_cub(norp(1))
                  Wt(2) = get_wt_cub(norp(2))
                  Wt(3) = get_wt_cub(norp(3))
    !           ----------------EXPONENTIAL SPLINES------------------
                else if(wexp .eq. 1) then
                  Wt(1) = get_wt_exp(norp(1), wcon)
                  Wt(2) = get_wt_exp(norp(2), wcon)
                  Wt(3) = get_wt_exp(norp(3), wcon)
                end if

    !     ------------------------------------------------

                Wtx(inw) = Wt(1)*Wt(2)*Wt(3)

              end do
            end do
          end do

    !     -------------------------------------------------
          do k1 = 1, siz
            do j1 = 1, siz
              do i1 = 1, siz
                inw = 1 + (i1 - 1) + (j1 - 1) * stride(1) + (k1 - 1) * stride(2)
                i = ci - (siz + 1)/2 + i1
                j = cj - (siz + 1)/2 + j1
                k = ck - (siz + 1)/2 + k1

                select case (cmp)
                case (1)
                  !ui(inw,1) = 0.5*(q1(i,j,k)+q1(i+1,j,k))
                  ui(inw,1) = q1(i,j,k)
                case (2)
                  !ui(inw,1) = 0.5*(q2(i,j,k)+q2(i,j+1,k))
                  ui(inw,1) = q2(i,j,k)
                case (3)
                  !ui(inw,1) = 0.5*(q3(i,j,k)+q3(i,j,k+1))
                  ui(inw,1) = q3(i,j,k)
                end select

              end do
            end do
          end do

    !     --------------------------------------------------------

          pinvA(:,:) = 0.0d0

    !     pre-inverse matrix A(x) and B(x)
          do k1 = 1, siz
            do j1 = 1, siz
              do i1 = 1, siz
                inw = 1 + (i1 - 1) + (j1 - 1) * stride(1) + (k1 - 1) * stride(2)
                i = ci - (siz + 1)/2 + i1
                j = cj - (siz + 1)/2 + j1
                k = ck - (siz + 1)/2 + k1

                select case (cmp)
                case (1)
                  pxk(1) = 1.0d0
                  pxk(2) = xc(i)
                  pxk(3) = ym(j)
                  pxk(4) = zm(k)
                case (2)
                  pxk(1) = 1.0d0
                  pxk(2) = xm(i)
                  pxk(3) = yc(j)
                  pxk(4) = zm(k)
                case (3)
                  pxk(1) = 1.0d0
                  pxk(2) = xm(i)
                  pxk(3) = ym(j)
                  pxk(4) = zc(k)
                end select

                call DGEMM('N','T',4,4,1,Wtx(inw),pxk,4,pxk,4, &
                                   1.0d0,pinvA,4)

    !     -------------------------------------------------------------------

                B(:,inw) = Wtx(inw)*pxk(:)

    !     -----------------------------------------------------------------

              end do
            end do
          end do

    !     calling routine to compute inverse
          !call invert4(pinvA,invA)
#ifdef LUINVERSE                     
          call inverseLU(pinvA,invA)
#else          
          call invert4(pinvA, invA)
#endif
    !     ------------------------------------------------------
    !     matrix multiplications for final interpolation
    !     DGEMM(transA,transB,m,n,k,alpha,A,LDA,b,LDB,beta,c,LDC)
    !     C = alpha * A * B + beta * C

    !     ---------------Shape function calculation---------------
          call DGEMM('N','N',1,4,4,1.0d0,ptx,1,invA,4,0.0d0,ptxA,1)
          call DGEMM('N','N',1,nel,4,1.0d0,ptxA,1,B,4,0.0d0,ptxAB,1)

    !     --------------Velocity Interpolation from \phi---------
          call DGEMM('N','N',1,1,nel,1.0d0,ptxAB,1,ui,nel,0.0d0,Um,1)


    !     ------------------------------------------------------

    !     For moving particles
          Forx = (vel_tri(cmp,ntr) - Um(1))
          !Fory = blend * (vel_tri(2,ntr) - Um(2))
          !Forz = blend * (vel_tri(3,ntr) - Um(3))

    !!     For stationary particles
    !      Forx = -Um(1)
    !      Fory = -Um(2)
    !      Forz = -Um(3)

    !!     For zero forcing
    !      Forx = 0.0
    !      Fory = 0.0
    !      Forz = 0.0

    !     -------------FORCING FUNCTION------------------------
    !     volume of a face with a specific marker - thickness taken as average
    !     of grid spacing

          mvol_t = 0.d0
          tcelvol = 0.d0
          do k1 = 1, siz
            do j1 = 1, siz
              do i1 = 1, siz
                inw = 1 + (i1 - 1) + (j1 - 1) * stride(1) + (k1 - 1) * stride(2)
                i = ci - (siz + 1)/2 + i1
                j = cj - (siz + 1)/2 + j1
                k = ck - (siz + 1)/2 + k1

                mvol_t = mvol_t + ptxAB(inw)*Hboxx(k)
                tcelvol = tcelvol + ptxAB(inw)*celvol(k)

              end do
            end do
          end do

          cfac = (mvol_t*sur(ntr)) / (tcelvol * dble(tri_tiling(ntr,1)))
    !     ------------------------------------------------------------------
    !     Resolutions need to be such that cfac is close to 1
    !     cfac is ratio of volume of face to the volume of cell it resides in
    !     For test case - can be artificially set to 1

          cel_fx = cfac*Forx
          !cel_fy = cfac*Fory
          !cel_fz = cfac*Forz


          ! if ( (chamb.EQ.5).OR.(chamb.EQ.2).OR.(LabelH(ntr).NE.0)  ) then
          !    cel_fx = 0.0D0
          ! endif

    !     ------------------------------------------------------

          do k1 = 1, siz
            do j1 = 1, siz
              do i1 = 1, siz
                inw = 1 + (i1 - 1) + (j1 - 1) * stride(1) + (k1 - 1) * stride(2)
                i = ci - (siz + 1)/2 + i1
                j = cj - (siz + 1)/2 + j1
                k = ck - (siz + 1)/2 + k1

                select case (cmp)
                case (1)
                  for_xc(i,j,k) = for_xc(i,j,k)+ptxAB(inw)*cel_fx
                case (2)
                  for_yc(i,j,k) = for_yc(i,j,k)+ptxAB(inw)*cel_fx
                case (3)
                  for_zc(i,j,k) = for_zc(i,j,k)+ptxAB(inw)*cel_fx
                end select

              end do
            end do
          end do

        end do

    !   -------------------------------------------------------
      end do
#endif

!@cuf istat = cudaDeviceSynchronize !JDR TMP
      !print*, "for_xc", sum(for_xc), minval(for_xc), maxval(for_xc)
      !print*, "for_yc", sum(for_yc), minval(for_yc), maxval(for_yc)
      !print*, "for_zc", sum(for_zc), minval(for_zc), maxval(for_zc)
!     --------------------------------------------------------

      return
      end
