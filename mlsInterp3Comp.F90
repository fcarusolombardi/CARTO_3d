module mlsInterp3Comp_m
  use constants
  contains

      ! JR TODO: consolidate weight functions into single module
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
        end if
        return
      end function get_wt_exp

      real(DP) function get_dwt_exp(val, el, H, wcon, wscl)
        implicit none
        real(DP), intent(in) :: val, el, H, wcon, wscl
        real(DP) :: rv

        if(val .le. 1.0d0)then
          rv = val / wcon
          get_dwt_exp = (-2.d0/(wcon * wcon))* val * &
                            exp(-rv * rv) * (el/abs(el)) * &
                            H/wscl
        else
          get_dwt_exp = 0.d0
        end if

        return
      end function get_dwt_exp

      
      subroutine compute_velp(siz, stride, pind_probe)
        use mpih
        use param
        use mls_param
        use local_arrays, only: q1,q2,q3,pr
        implicit none

        integer, intent(in) :: siz, stride(*)
        integer, dimension(6, *), intent(in) :: pind_probe
        integer :: inw, i, j, k, i1, j1, k1, ci, cj, ck
        integer :: tr, ntr, inp, cmp, cmp2
        integer :: v1, v2, v3
        real(DP) :: press, pr_p, epsw, facsign
        real(DP) :: du(3,3)
        real(DP) :: pos_pro(3), fsur_xyz(3), tau_f(3)
        real(DP) :: Wt(3), Wtx(nel), dWt(3), dWti(3, nel)
        real(DP) :: Hbox(3), el_cx(3), norp(3)
        real(DP) :: prp(nel), ui(nel, 3), B(4,nel), Bi(4,nel)
        real(DP) :: pinvA(4,4), invA(4,4)
        real(DP) :: Gmat(4), Gmati(4), Gmat1(4), Gmat2(4), Gmat3(4)
        real(DP) :: Gtib(nel), Gtbi(nel), PhiTi(nel)
        real(DP) :: ptx(4), pxk(4), ptxA(4), ptxAB(nel), ptxABu(1,1)
        real(DP) :: fnca(3)
        real(DP) :: sr_xx, sr_xy, sr_xz
        real(DP) :: sr_yy, sr_yz, sr_yx
        real(DP) :: sr_zz, sr_zx, sr_zy
        real(DP) :: Um1,Um2,Um3,Um4

        do tr = 1, tpcnt(1)
          ntr = mypr(tr)

          !  probe location
          pos_pro(1) = xyzMLSprobe(1,ntr)  
          pos_pro(2) = xyzMLSprobe(2,ntr)  
          pos_pro(3) = xyzMLSprobe(3,ntr)  


          do cmp = 1, 4 ! cmp = 4 is for pressure
!           Support domain for probe
            !call partindicesMLS(pos_pro,pind_i,pind_o,dismaxpro)
            ! JR Note: Inlining what is needed from partindicesMLS
            select case (cmp)
            case (1) ! 1,5,6
              ci = pind_probe(1,ntr)
              cj = pind_probe(5,ntr)
              ck = pind_probe(6,ntr)
            case (2) ! 4,2,6
              ci = pind_probe(4,ntr)
              cj = pind_probe(2,ntr)
              ck = pind_probe(6,ntr)
            case(3) ! 4,5,3
              ci = pind_probe(4,ntr)
              cj = pind_probe(5,ntr)
              ck = pind_probe(3,ntr)
            case (4) ! 4,5,6
              ci = pind_probe(4,ntr)
              cj = pind_probe(5,ntr)
              ck = pind_probe(6,ntr)
            end select

!           initialise pre-factor matrix

            ptx(1) = 1.d0
            ptx(2) = pos_pro(1)
            ptx(3) = pos_pro(2)
            ptx(4) = pos_pro(3)

!           --------------------------------------------------------

            Hbox(1) = dx1
            Hbox(2) = dx2
            Hbox(3) = dx3

!           --------------------------------------------------------

            epsw = 1e-18

!           spline weights for all three components and initialise ui matrix
            do k1 = 1, siz
              do j1 = 1, siz
                do i1 = 1, siz
                  inw = 1 + (i1 - 1) + (j1 - 1) * stride(1) + (k1 - 1) * stride(2)
                  i = ci - (siz + 1)/2 + i1
                  j = cj - (siz + 1)/2 + j1
                  k = ck - (siz + 1)/2 + k1

                  select case (cmp)
                  case (1)
                    el_cx(1) = xc(i) - pos_pro(1) !+ epsw
                    el_cx(2) = ym(j) - pos_pro(2) !+ epsw
                    el_cx(3) = zm(k) - pos_pro(3) !+ epsw
                  case (2)
                    el_cx(1) = xm(i) - pos_pro(1) !+ epsw
                    el_cx(2) = yc(j) - pos_pro(2) !+ epsw
                    el_cx(3) = zm(k) - pos_pro(3) !+ epsw
                  case (3)
                    el_cx(1) = xm(i) - pos_pro(1) !+ epsw
                    el_cx(2) = ym(j) - pos_pro(2) !+ epsw
                    el_cx(3) = zc(k) - pos_pro(3) !+ epsw
                  case (4)
                    el_cx(1) = xm(i) - pos_pro(1) !+ epsw
                    el_cx(2) = ym(j) - pos_pro(2) !+ epsw
                    el_cx(3) = zm(k) - pos_pro(3) !+ epsw
                  end select

                  norp(1) = abs(el_cx(1)) * Hbox(1) / wscl
                  norp(2) = abs(el_cx(2)) * Hbox(2) / wscl
                  norp(3) = abs(el_cx(3)) * Hbox(3) / wscl

    !             ----------------EXPONENTIAL SPLINES------------------
                  if(wexp .eq. 1) then
                    Wt(1) = get_wt_exp(norp(1), wcon)
                    Wt(2) = get_wt_exp(norp(2), wcon)
                    Wt(3) = get_wt_exp(norp(3), wcon)

                    dWt(1) = get_dwt_exp(norp(1), el_cx(1), Hbox(1), wcon, wscl)
                    dWt(2) = get_dwt_exp(norp(2), el_cx(2), Hbox(2), wcon, wscl)
                    dWt(3) = get_dwt_exp(norp(3), el_cx(3), Hbox(3), wcon, wscl)
                  end if

                  Wtx(inw) = Wt(1)*Wt(2)*Wt(3)

                  dWti(1, inw) = dWt(1)*Wt(2)*Wt(3)
                  dWti(2, inw) = Wt(1)*dWt(2)*Wt(3)
                  dWti(3, inw) = Wt(1)*Wt(2)*dWt(3)
!                 ------------------------------------------------

!                 Allocating variables for support domain
                  select case (cmp)
                  case (1)
                    !ui(inw,1) = 0.5d0*(q1(i,j,k)+q1(i+1,j,k))
                    ui(inw,1) = q1(i,j,k)
                  case (2)
                    !ui(inw,1) = 0.5d0*(q2(i,j,k)+q2(i,j+1,k))
                    ui(inw,1) = q2(i,j,k)
                  case (3)
                    !ui(inw,1) = 0.5d0*(q3(i,j,k)+q3(i,j,k+1))
                    ui(inw,1) = q3(i,j,k)
                  case (4)
!                     prp(inw) = pr(i,j,k)
                     ui(inw,1) = pr(i,j,k) !FV
                  end select

                end do
              end do
            end do

!           --------------------------------------------------------

            pinvA(:,:) = 0.0d0

!           pre-inverse matrix A(x) and B(x)
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
                  case (4)
                    pxk(1) = 1.0d0
                    pxk(2) = xm(i)
                    pxk(3) = ym(j)
                    pxk(4) = zm(k)
                  end select

                  call DGEMM('N','T',4,4,1,Wtx(inw),pxk,4,pxk,4,1.0d0,pinvA,4)

                  B(:,inw) = Wtx(inw)*pxk(:)

                end do
              end do
            end do

!           calling routine to compute inverse
! #ifdef LUINVERSE                       
            call inverseLU(pinvA,invA) !FV NOTA!this is more accurate than the analytical one
! #else
!             call invert4(pinvA,invA)
! #endif

            
! ! !           Compute Gamma=inva*p for derivatives of shape functions
!             call DGEMM('N','T',4,1,4,1.0d0,invA,4,ptx,1,0.0d0,Gmat,4)

!             ------------------------------------------------------
!             matrix multiplications for final interpolation
!             DGEMM(transA,transB,m,n,k,alpha,A,LDA,b,LDB,beta,c,LDC)
!             C = alpha * A * B + beta * C

!             ---------------Shape function calculation---------------
              call DGEMM('N','N',1,4,4,1.0d0,ptx,1,invA,4,0.0d0,ptxA,1)
              call DGEMM('N','N',1,nel,4,1.0d0,ptxA,1,B,4,0.0d0,ptxAB,1)
!             --------------Pressure Interpolation at probe from \phi---------
!              call DGEMM('N','N',1,1,nel,1.0d0,ptxAB,1,prp(:),nel,0.0d0,ptxABu,1)
              call DGEMM('N','N',1,1,nel,1.0d0,ptxAB,1,ui,nel,0.0d0,ptxABu,1) !FV


              select CASE(cmp)
              CASE(1)
                 Um1 = ptxABu(1,1)
              CASE(2)
                 Um2 = ptxABu(1,1)
              CASE(3)
                 Um3 = ptxABu(1,1)
              CASE(4)
                 Um4 = ptxABu(1,1)
              end select              
           end do !end loop cmp

           !Calculate forces on triangle faces - pressure+viscous
          outMLSprobe(1,ntr) = Um1
          outMLSprobe(2,ntr) = Um2
          outMLSprobe(3,ntr) = Um3
          outMLSprobe(4,ntr) = Um4

        end do !end loop probes
      end subroutine compute_velp

#ifdef USE_CUDA
      attributes(global) &
      subroutine get_valid_triangles(trcnt, mytr, pind, nftot, kstart, kend)
        implicit none

        integer, value, intent(in) :: nftot, kstart, kend
        integer, device, intent(inout) :: trcnt
        integer, device, dimension(nftot), intent(out) :: mytr
        integer, device, dimension(6, nftot), intent(in) :: pind

        integer ::i, ind_pal

        i = (blockIdx%x - 1) * blockDim%x + threadIdx%x

        if (i > nftot) return

        ind_pal = pind(6, i)

        if(ind_pal .ge. kstart .and. ind_pal .le. kend) then
          mytr(atomicAdd(trcnt, 1) + 1) = i
        endif

      end subroutine get_valid_triangles
  

       attributes(device) &     
       real(DP) function get_wt_exp_gpu(val, wcon)
        implicit none
        real(DP), intent(in), device :: val, wcon
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
      real(DP) function get_dwt_exp_gpu(val, el, H, wcon, wscl)
        implicit none
        real(DP), intent(in), device :: val, el, H, wcon, wscl
        real(DP) :: rv

        if(val .le. 1.0d0)then
          rv = val / wcon
          get_dwt_exp_gpu = (-2.d0/(wcon * wcon))* val * &
                            exp(-rv * rv) * (el/abs(el)) * &
                            H/wscl
        else
          get_dwt_exp_gpu = 0.d0
        end if

        return
      end function get_dwt_exp_gpu

      attributes(device) subroutine det4_gpu(a, det_l)
       implicit none
       real(DP), device, intent(in) :: a (4,4)
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

      attributes(device) subroutine invert4_unscaled_gpu(a,ia)
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
       real(DP), device, intent(in) :: a (4,4)
       real(DP), device, intent(out) :: ia (4,4)
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

      subroutine compute_velp_kernel(siz, str1, str2, &
        nel, dx1, dx2, dx3, wscl, wcon, wcub, wexp, n1m, &
        n2m, nu, outprobe, pind_probe, xyzprobe, q1, q2, q3, pr, &
        mypr, tpcnt, xm, ym, zm, xc, yc, zc, &
        n1, n2, n3, ks, ke, ksr, ker, lvlhalo, xyz_tr)
        implicit none

        integer, value, intent(in) :: siz, str1, str2, n1m, n2m
        integer, value, intent(in) ::n1, n2, n3, ks, ke, ksr, ker, lvlhalo
        integer, value, intent(in) :: nel, wcub, wexp
        real(DP), value, intent(in) :: dx1, dx2, dx3, wscl, wcon, nu
        integer, device, intent(in), dimension(1) :: tpcnt

        integer, device, dimension(6, *), intent(in) :: pind_probe
        real(DP), device, dimension(n1, n2, ks-lvlhalo:ke+lvlhalo), intent(in) :: q1, q2, q3, pr
        real(DP), device, dimension(3, *), intent(in) :: xyzprobe        
        real(DP), device, dimension(*), intent(in) :: ym, zm, xc, yc, zc
        real(DP), device, dimension(0:n3), intent(in) :: xm
        integer, device, dimension(*), intent(in) :: mypr
        real(DP), device, dimension(3), intent(in) :: xyz_tr

        real(DP), device, dimension(4, *), intent(inout) :: outprobe

        integer :: ci, cj, ck
        integer :: i, j, k, ii, jj, inw
        integer :: i1, j1, k1
        integer :: i2, j2
        real(DP), device :: pos_pro(3), norp(3), el_cx(3), ptx(4), Hbox(3), Wt(3), dWt(3)
        real(DP), device :: ui(4)
        real(DP), device :: fnca(3)
        real(DP), device :: B(4), Bi(4)
        real(DP) :: Wtx, dWti(3), rv1, rv2, cfac, epsw, Um1,Um2,Um3,Um4
        real(DP) :: For, tau, press

        !real(DP), shared :: pinvA(4,4,4)
        !real(DP), device :: pinvA(4,4), invA(4,4)
        !real(DP), device :: pinvA(16), invA(16)
        real(DP), shared :: pinvA(16,4), invA(16,4), det_l(4)
        real(DP), shared :: ptxA(4, 4)
        real(DP), shared :: gmat(4, 4)
        real(DP), shared :: gmat2(4, 4)
        real(DP), shared :: gmat3(4, 4)
        real(DP), shared :: gmati(4, 4)
        real(DP), shared :: du(9, 4)
        real(DP) :: phiTi
        real(DP) :: rcxx, rcyy, rczz, appr_ch, rcsqi

        integer :: tid, tx, ty, tr, ntr, inp, cmp, cmp2, v1, v2, v3

        !tid = (blockIdx%x - 1) * blockDim%x + threadIdx%x
        !tx = mod(tid-1, 32) + 1
        !ty = (tx - 1)/32 + 1
        tx = threadIdx%x
        ty = threadIdx%y
        tr = (blockIdx%x - 1) * blockDim%y + threadIdx%y

        if (tr > tpcnt(1)) return
        ntr = mypr(tr)

        i1 = mod(tx - 1, str1) + 1
        j1 = mod((tx - 1) / str1, str1) + 1
        k1 = (tx - 1) / str2 + 1
        inw = 1 + (i1 - 1) + (j1 - 1) * str1 + (k1 - 1) * str2

!       Compute positions of probes from centroids
        pos_pro(1) = xyzprobe(1,ntr)
        pos_pro(2) = xyzprobe(2,ntr)
        pos_pro(3) = xyzprobe(3,ntr)


        do cmp = 1,4 ! cmp = 4 is for pressure
!         Support domain for probe
          !call partindicesMLS(pos_pro,pind_i,pind_o,dismaxpro)
           ! JR Note: Inlining what is needed from partindicesMLS
            select case (cmp)
            case (1) ! 1,5,6
              ci = pind_probe(1,ntr)
              cj = pind_probe(5,ntr)
              ck = pind_probe(6,ntr)
            case (2) ! 4,2,6
              ci = pind_probe(4,ntr)
              cj = pind_probe(2,ntr)
              ck = pind_probe(6,ntr)
            case(3) ! 4,5,3
              ci = pind_probe(4,ntr)
              cj = pind_probe(5,ntr)
              ck = pind_probe(3,ntr)
            case (4) ! 4,5,6
              ci = pind_probe(4,ntr)
              cj = pind_probe(5,ntr)
              ck = pind_probe(6,ntr)
            end select
          i = ci - (siz + 1)/2 + i1
          j = cj - (siz + 1)/2 + j1
          k = ck - (siz + 1)/2 + k1

  !       --------------------------------------------------------

          Hbox(1) = dx1
          Hbox(2) = dx2
          Hbox(3) = dx3

  !       --------------------------------------------------------

  !       spline weights for all three components
          Wtx = 0
          epsw = 1e-18
          dWti(1) = 0
          dWti(2) = 0
          dWti(3) = 0

          if (tx <= nel) then

            select case (cmp)
            case (1)
              el_cx(1) = xc(i) - pos_pro(1) !+ epsw
              el_cx(2) = ym(j) - pos_pro(2) !+ epsw
              el_cx(3) = zm(k) - pos_pro(3) !+ epsw
            case (2)
              el_cx(1) = xm(i) - pos_pro(1) !+ epsw
              el_cx(2) = yc(j) - pos_pro(2) !+ epsw
              el_cx(3) = zm(k) - pos_pro(3) !+ epsw
            case (3)
              el_cx(1) = xm(i) - pos_pro(1) !+ epsw
              el_cx(2) = ym(j) - pos_pro(2) !+ epsw
              el_cx(3) = zc(k) - pos_pro(3) !+ epsw
            case (4)
              el_cx(1) = xm(i) - pos_pro(1) !+ epsw
              el_cx(2) = ym(j) - pos_pro(2) !+ epsw
              el_cx(3) = zm(k) - pos_pro(3) !+ epsw
            end select

            norp(1) = abs(el_cx(1)) * Hbox(1) / wscl
            norp(2) = abs(el_cx(2)) * Hbox(2) / wscl
            norp(3) = abs(el_cx(3)) * Hbox(3) / wscl

  !             ----------------CUBIC SPLINES---------------------
            if(wexp .eq. 1) then
              Wt(1) = get_wt_exp_gpu(norp(1), wcon)
              Wt(2) = get_wt_exp_gpu(norp(2), wcon)
              Wt(3) = get_wt_exp_gpu(norp(3), wcon)

              dWt(1) = get_dwt_exp_gpu(norp(1), el_cx(1), Hbox(1), wcon, wscl)
              dWt(2) = get_dwt_exp_gpu(norp(2), el_cx(2), Hbox(2), wcon, wscl)
              dWt(3) = get_dwt_exp_gpu(norp(3), el_cx(3), Hbox(3), wcon, wscl)
            end if

  !       ------------------------------------------------

            Wtx = Wt(1)*Wt(2)*Wt(3)

            dWti(1) = dWt(1)*Wt(2)*Wt(3)
            dWti(2) = Wt(1)*dWt(2)*Wt(3)
            dWti(3) = Wt(1)*Wt(2)*dWt(3)
          endif

          ui(1) = 0.d0
          ! ui(2) = 0.d0
          ! ui(3) = 0.d0
          ! ui(4) = 0.d0
          
          if (tx <= nel) then
            select case (cmp)
            case (1)
              !ui(1) = 0.5*(q1(i,j,k)+q1(i+1,j,k))
              ui(1) = q1(i,j,k)
            case (2)
              !ui(1) = 0.5*(q2(i,j,k)+q2(i,j+1,k))
              ui(1) = q2(i,j,k)
            case (3)
              !ui(1) = 0.5*(q3(i,j,k)+q3(i,j,k+1))
              ui(1) = q3(i,j,k)
            case (4)
              ui(1) = pr(i,j,k)
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
            case (4)
              ptx(1) = 1.0d0
              ptx(2) = xm(i)
              ptx(3) = ym(j)
              ptx(4) = zm(k)
            end select

            B(1) = Wtx * ptx(1)
            B(2) = Wtx * ptx(2)
            B(3) = Wtx * ptx(3)
            B(4) = Wtx * ptx(4)
          endif

          ! each thread does outer product
          do j2 = 1, 4
            do i2 = 1, 4
              rv1 = Wtx*ptx(i2)*ptx(j2)

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
                pinvA(i2 + (j2-1)*4, ty) = rv1
              endif

              !call syncthreads()
              !pinvA(i1 + (j1-1)*4) = rv1

            enddo
          enddo

          call syncthreads()

          !call invert4_gpu(pinvA, invA)
          if (tx == 1) then
             ! call det4_gpu(pinvA(1,ty), det_l(ty))
             ! call invert4_unscaled_gpu(pinvA(1,ty), invA(1,ty))             
! #ifdef LUINVERSE                        
             call inverseLU_gpu(pinvA(1,ty), invA(1,ty))  !FV NOTA!this is more accurate than the analytical one
! #else
!              call invert4_gpu(pinvA(1,ty), invA(1,ty))
! #endif        


          endif
          call syncthreads()

          ! ! Compute gmat (using 16 threads)
          ! if (tx <= 16) then
          !   ptx(1) = 1.d0;
          !   ptx(2) = pos_pro(1)
          !   ptx(3) = pos_pro(2)
          !   ptx(4) = pos_pro(3)

          !   rv1 = ptx(mod(tx-1, 4) + 1) * invA(4*(mod(tx-1, 4)) + (tx-1)/4 + 1, ty) / det_l(ty)
          ! endif

          ! rv2 = __shfl_down(rv1,1)
          ! rv1 = rv1 +  rv2
          ! rv2 = __shfl_down(rv1,2)
          ! rv1 = rv1 +  rv2

          ! if (tx <= 16) then
          !   if (mod(tx-1, 4) == 0) then
          !     gmat((tx-1)/4 + 1, ty) = rv1
          !   endif
          ! endif

          ! call syncthreads()

          ! ! Compute gradient terms componentwise
          ! if (tx <= nel) then
          !   select case (cmp)
          !   case (1)
          !     ptx(1) = 1.0d0
          !     ptx(2) = xc(i)
          !     ptx(3) = ym(j)
          !     ptx(4) = zm(k)
          !   case (2)
          !     ptx(1) = 1.0d0
          !     ptx(2) = xm(i)
          !     ptx(3) = yc(j)
          !     ptx(4) = zm(k)
          !   case (3)
          !     ptx(1) = 1.0d0
          !     ptx(2) = xm(i)
          !     ptx(3) = ym(j)
          !     ptx(4) = zc(k)
          !   case (4)
          !     ptx(1) = 1.0d0
          !     ptx(2) = xm(i)
          !     ptx(3) = ym(j)
          !     ptx(4) = zm(k)
          !   end select
          ! endif

!vel p components          
            ! Compute ptxA (using 16 threads)
            if (tx <= 16) then
              ptx(1) = 1.d0;
              ptx(2) = pos_pro(1)
              ptx(3) = pos_pro(2)
              ptx(4) = pos_pro(3)

              rv1 = ptx(mod(tx-1, 4) + 1) * invA(tx, ty) !/ det_l(ty)
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
              do i2 = 1, 4
                rv1 = rv1 + ptxA(i2, ty) * B(i2)
              enddo
              B(1) = rv1
            endif

            ! Compute prp
! JR note: if statement unneeded. ui set to zero for tx > nel earlier
            ! if (tx <= nel) then
              rv1 = B(1) * ui(1)
            ! else
            !   rv1 = 0.d0
            ! endif

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

            select CASE(cmp)
            CASE(1)
               Um1 = rv1
            CASE(2)
               Um2 = rv1
            CASE(3)
               Um3 = rv1
            CASE(4)
               Um4 = rv1
            end select

        end do
        !Calculate forces on triangle faces - pressure+viscous

        outprobe(1,ntr) = Um1
        outprobe(2,ntr) = Um2
        outprobe(3,ntr) = Um3
        outprobe(4,ntr) = Um4

        

      end subroutine compute_velp_kernel

#endif

end module mlsInterp3Comp_m
!-----------------------------------------------------------------------
!     read in position of a individual trial marker for MLS and compute
!     compute support domain, shape function and interpolate
!
!     compute internal forces on the nodes
!------------------------------------------------------------------------

      subroutine mlsInterp3Comp
      USE mpih
      USE param
      USE mls_param
      USE local_arrays, only: q1,q2,q3
      USE mpi_param, only: kstart, kend, kstartr, kendr
      USE mls_local, only: coll
      USE local_arrays, only: pr
      USE ieee_arithmetic
      use mlsInterp3Comp_m
      use nvtx
!@cuf use cudafor
      IMPLICIT NONE

      integer :: inp, chamb, ntr
      integer :: siz, stride(2)
      integer :: i, j, k, tr

      real(DP) :: pos_MLS(3)
      integer i1,j1,k1,ist,jst,kst
      integer :: f1, f2, v1, v2, v3, v4, iv, ie
      real(DP), dimension(3) :: tvec1, tvec2, tvec3, tvec4, csi, zet, Vij
      real(DP), dimension(3) :: a32, a13, a34, a21, a42, a23, a31, a24
      real(DP), dimension(3) :: csi1, zet1,tcv
      real(DP) :: modcsi, modzet, betab, betav
      real(DP) :: alphat, alphal, alphae, b11, b12, b22
      real(DP) :: tdum,VV, Vij1, Vij2, Vij3

      integer :: ci, cj, ck
      real(DP), dimension(3) :: fcoll
      real(DP) :: rcxx, rcyy, rczz, rcsqi
      real(DP) :: appr_ch, dtr
      real(DP) :: d, d0, kei

      integer :: imasZ
      integer :: imasXYZ(3)
      real(DP) :: angleHelix, angleSpiral, masZ
      real(DP) :: xCC, yCC, zCC, surT
      real(DP) :: alphaES, betaES, gammaES, smoothf, distMEAN, CosePhi, SenoPhi
      real(DP) :: visc_coeff, visc_coeff2
      real(DP) :: fsmth, fsmth2, fsmth3, xV, yV, zV, xV1, yV1, zV1, rV1
      real(DP) :: dinvm,pos1,pos2,pos3,ztmp
      integer kok,kini,kfin,kmid,indprobe
!     --------------------------------------------------------

#ifdef USE_CUDA
      type(dim3) :: blocks, threads
#endif
!@cuf integer :: istat

      if (istr3.EQ.0) then      
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do indprobe = 1,Nmlsprobe
         pos1=xyzMLSprobe(1,indprobe)
         pos2=xyzMLSprobe(2,indprobe)
         pos3=xyzMLSprobe(3,indprobe)
         
!       ++++++++Indices of the marker++++++++++++++++++++++++	
!       X - indices
        i1=NINT((pos1+xyz_tr(1))*dx1) +1
        ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m

!       Y - indices
        j1=NINT((pos2+xyz_tr(2))*dx2) +1
        jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m

!       Z - indices
!
        k1=NINT((pos3+xyz_tr(3))*dx3) +1
        kst=FLOOR((pos3+xyz_tr(3))*dx3) + 1

!       -------------------------------------------------------------
        pindMLSprobe(1,indprobe)=i1 ; pindMLSprobe(2,indprobe)=j1 ; pindMLSprobe(3,indprobe)=k1
        pindMLSprobe(4,indprobe)=ist; pindMLSprobe(5,indprobe)=jst; pindMLSprobe(6,indprobe)=kst
!       -------------------------------------------------------------
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP


      else !GENERAL CASE
         
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do indprobe = 1,Nmlsprobe
         pos1=xyzMLSprobe(1,indprobe)
         pos2=xyzMLSprobe(2,indprobe)
         pos3=xyzMLSprobe(3,indprobe)
         
!       ++++++++Indices of the marker++++++++++++++++++++++++	
!       X - indices
        i1=NINT((pos1+xyz_tr(1))*dx1) +1
        ist=FLOOR((pos1+xyz_tr(1))*dx1) +1      
        if(ist.eq.0)ist=n1m

!       Y - indices
        j1=NINT((pos2+xyz_tr(2))*dx2) +1
        jst=FLOOR((pos2+xyz_tr(2))*dx2) + 1
        if(jst.eq.0)jst=n2m

!       Z - indices
        k1=0
        kini=1
        kfin=n3        
        do while (k1.eq.0)
           if ((kfin-kini).eq.1) then
              ztmp=0.5D0*(zc(kfin)+zc(kini))
              if (pos3.le.ztmp) then
                 k1=kini
              else
                 k1=kfin
              endif
           endif
           kmid=(kini+kfin)/2
           if (pos3.gt.zc(kini).and.pos3.le.zc(kmid)) then
              kini=kini
              kfin=kmid
           else
              kini=kmid
              kfin=kfin
           endif
        enddo
        kst=k1
!        if((pos3).gt.zc(k1))kst=k1   !fv
        if((pos3).le.zc(k1))kst=k1-1 !fv

!       -------------------------------------------------------------
        pindMLSprobe(1,indprobe)=i1 ; pindMLSprobe(2,indprobe)=j1 ; pindMLSprobe(3,indprobe)=k1
        pindMLSprobe(4,indprobe)=ist; pindMLSprobe(5,indprobe)=jst; pindMLSprobe(6,indprobe)=kst
!       -------------------------------------------------------------
      end do

!@cuf istat = cudaDeviceSynchronize !JDR TMP

   endif
   
!      Get valid triangles precompute valid triangles
#ifdef USE_CUDA
      tpcnt = 0
      threads = dim3(128, 1, 1)
      blocks = dim3(ceiling(real(Nmlsprobe)/128), 1, 1)
      call get_valid_triangles<<<blocks, threads>>>(tpcnt(1), mypr, pindMLSprobe, &
           Nmlsprobe, kstart, kend) 
!@cuf istat = cudaDeviceSynchronize !JDR TMP     
#else
!      precompute valid triangles
      tr = 1
      tpcnt = 0
      do ntr = 1,Nmlsprobe 
        indprobe = pindMLSprobe(6,ntr)
        if(indprobe .ge. kstart .and. indprobe .le. kend) then
          mypr(tr) = ntr
          tr = tr + 1
          tpcnt(1) = tpcnt(1) + 1
        endif
      enddo
#endif

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

      outMLSprobe(:,:)=0.D0
      
      call nvtxStartRange("compute_velp", 9)
#ifdef USE_CUDA
      threads = dim3(32, 4, 1)
      blocks = dim3(ceiling(real(Nmlsprobe)/4), 1, 1)

      call compute_velp_kernel<<<blocks, threads>>>(siz, stride(1), stride(2), &
      nel, dx1, dx2, dx3, wscl, wcon, wcub, wexp, n1m, n2m, &
      nu, outMLSprobe, pindMLSprobe, xyzMLSprobe, &
      q1, q2, q3, pr, mypr, tpcnt, xm, ym, zm, xc, yc, zc, &
      n1, n2, n3, kstart, kend, kstartr, kendr, lvlhalo, &
      xyz_tr)
!@cuf istat = cudaDeviceSynchronize !JDR TMP
#else
      call compute_velp(siz, stride, pindMLSprobe)
#endif
      call nvtxEndRange


!     --------------------------------------------------------
!     Reduce the forces from all processors over each particle
      if (numtasks > 1) then
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,outMLSprobe,4*Nmlsprobe,MPI_DOUBLE_PRECISION, &
              MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      endif


      return
      end
