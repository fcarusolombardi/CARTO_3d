module mlsStruc3CompStrong_m
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

      subroutine compute_fpxyz(fac, iopen, siz, stride, pind_probe)
        use mpih
        use param
        use mls_param
        use local_arrays, only: q1,q2,q3,pr
        implicit none

        real(DP), intent(in) :: fac
        integer, intent(in) :: siz, stride(*), iopen(*)
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

        facsign =  sign(1.D0,fac) !modificaFV verificare lo vedano tutti

        do tr = 1, trcnt(1)

          ntr = mytr(tr)
          inp = face_to_part(ntr)

          ! Skip if contribution from this probe not needed
          if (iopen(inp) .ne. 1) cycle

          fsur_xyz(:) = 0.d0

          !  probe location
          pos_pro(1) = tri_bar(1,ntr)  + fac * tri_nor(1,ntr) * Hboxx(pind(6,ntr))
          pos_pro(2) = tri_bar(2,ntr)  + fac * tri_nor(2,ntr) * Hboxx(pind(6,ntr))
          pos_pro(3) = tri_bar(3,ntr)  + fac * tri_nor(3,ntr) * Hboxx(pind(6,ntr))


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

            epsw = 1e-8

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
                    el_cx(1) = xc(i) - pos_pro(1) + epsw
                    el_cx(2) = ym(j) - pos_pro(2) + epsw
                    el_cx(3) = zm(k) - pos_pro(3) + epsw
                  case (2)
                    el_cx(1) = xm(i) - pos_pro(1) + epsw
                    el_cx(2) = yc(j) - pos_pro(2) + epsw
                    el_cx(3) = zm(k) - pos_pro(3) + epsw
                  case (3)
                    el_cx(1) = xm(i) - pos_pro(1) + epsw
                    el_cx(2) = ym(j) - pos_pro(2) + epsw
                    el_cx(3) = zc(k) - pos_pro(3) + epsw
                  case (4)
                    el_cx(1) = xm(i) - pos_pro(1) + epsw
                    el_cx(2) = ym(j) - pos_pro(2) + epsw
                    el_cx(3) = zm(k) - pos_pro(3) + epsw
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
                    prp(inw) = pr(i,j,k)
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
#ifdef LUINVERSE                       
            call inverseLU(pinvA,invA)
#else
            call invert4(pinvA,invA)            
#endif
!           Compute Gamma=inva*p for derivatives of shape functions
            call DGEMM('N','T',4,1,4,1.0d0,invA,4,ptx,1,0.0d0,Gmat,4)

            if (cmp <= 3) then
              ! Compute gradient terms componentwise
              do cmp2 = 1, 3
                pinvA = 0.d0

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

  !                   ------ derivative of A matrix - req for Gam derivative ------------
                      call DGEMM('N','T',4,4,1,dWti(cmp2, inw),pxk,4,pxk,4,1.0d0,pinvA,4)

  !                   -------derivative of B matrix - req for shape derivative---------
                      Bi(:,inw) = dWti(cmp2,inw)*pxk(:)

                    end do
                  end do
                end do

!               Compute Gamma=inva*p for derivatives of shape functions
                !call DGEMM('N','N',4,1,4,1.0d0,invA,4,pxx,4,0.0d0,Gmatx1,4)
                call DGEMM('N','N',4,1,4,1.0d0,pinvA,4,Gmat,4,0.0d0,Gmat2,4)
                call DGEMM('N','N',4,1,4,1.0d0,invA,4,Gmat2,4,0.0d0,Gmat3,4)
                Gmati = invA(:,cmp2+1) - Gmat3

!               Compute shape function derivatives
                call DGEMM('T','N',1,nel,4,1.0d0,Gmati,4,B,4,0.0d0,Gtib,1)
                call DGEMM('T','N',1,nel,4,1.0d0,Gmat,4,Bi,4,0.0d0,Gtbi,1)
                PhiTi = Gtib + Gtbi

!               ------Velocity gradient interpolation-------------------
                call DGEMM('N','N',1,1,nel,1.0d0,PhiTi,1,ui,nel,0.0d0,du(cmp,cmp2),1)

              enddo
            else

!             ------------------------------------------------------
!             matrix multiplications for final interpolation
!             DGEMM(transA,transB,m,n,k,alpha,A,LDA,b,LDB,beta,c,LDC)
!             C = alpha * A * B + beta * C

!             ---------------Shape function calculation---------------
              call DGEMM('N','N',1,4,4,1.0d0,ptx,1,invA,4,0.0d0,ptxA,1)
              call DGEMM('N','N',1,nel,4,1.0d0,ptxA,1,B,4,0.0d0,ptxAB,1)

!             --------------Pressure Interpolation at probe from \phi---------
              call DGEMM('N','N',1,1,nel,1.0d0,ptxAB,1,prp(:),nel,0.0d0,ptxABu,1)

              pr_p = ptxABu(1,1)

            endif
          end do

!         -----------------------------------------------------------------
!         Face normals
          fnca(1:3) = tri_nor(1:3,ntr)

!         Compute pressure at marker from probe
          press = pr_p + fac * Hboxx(pind(6,ntr))*(acc_tri(1,ntr)*fnca(1)+ &
                 acc_tri(2,ntr)*fnca(2)+acc_tri(3,ntr)*fnca(3))

!         Compute shear stress from vel.gradients
          sr_xx = du(1,1)
          sr_yy = du(2,2)
          sr_zz = du(3,3)

          sr_xy = 0.5d0*(du(1,2)+du(2,1))
          sr_yz = 0.5d0*(du(2,3)+du(3,2))
          sr_zx = 0.5d0*(du(3,1)+du(1,3))

          sr_yx = sr_xy
          sr_zy = sr_yz
          sr_xz = sr_zx

          tau_f(1) = 2.d0*nu*(sr_xx*fnca(1)+sr_xy*fnca(2)+sr_xz*fnca(3))
          tau_f(2) = 2.d0*nu*(sr_yx*fnca(1)+sr_yy*fnca(2)+sr_yz*fnca(3))
          tau_f(3) = 2.d0*nu*(sr_zx*fnca(1)+sr_zy*fnca(2)+sr_zz*fnca(3))

!         -----------------------------------------------------------------
          fsur_xyz(1:3) = ( tau_f(1:3) - press  * fnca(1:3)) * (sur(ntr)/3.0)

          v1 = vert_of_face(1,ntr)
          v2 = vert_of_face(2,ntr)
          v3 = vert_of_face(3,ntr)

!         Calculate forces on triangle faces - pressure+viscous
          fpxyz(1:3,v1) = fpxyz(1:3,v1) + facsign * fsur_xyz(1:3)
          fpxyz(1:3,v2) = fpxyz(1:3,v2) + facsign * fsur_xyz(1:3)
          fpxyz(1:3,v3) = fpxyz(1:3,v3) + facsign * fsur_xyz(1:3)

        end do

      end subroutine compute_fpxyz

#ifdef USE_CUDA
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
      subroutine compute_fpxyz_kernel(siz, str1, str2, &
        nel, dx1, dx2, dx3, wscl, wcon, wcub, wexp, n1m, &
        n2m, nu, fpxyz, pind, pind_probe, tri_bar, tri_nor, vel_tri, acc_tri, sur, coll, q1, q2, q3, pr, &
        vert_of_face, face_to_part, mytr, trcnt, xm, ym, zm, xc, yc, zc, Hboxx, nftot, nvtot, &
        n1, n2, n3, ks, ke, ksr, ker, lvlhalo, fac, iopen, xyz_tr)
        implicit none

        integer, value, intent(in) :: siz, str1, str2, n1m, n2m
        integer, value, intent(in) :: nvtot, nftot, n1, n2, n3, ks, ke, ksr, ker, lvlhalo
        integer, value, intent(in) :: nel, wcub, wexp
        real(DP), value, intent(in) :: dx1, dx2, dx3, wscl, wcon, nu
        real(DP), value, intent(in) :: fac
        integer, device, intent(in), dimension(1) :: trcnt

        integer, device, dimension(6, *), intent(in) :: pind
        integer, device, dimension(6, *), intent(in) :: pind_probe        
        integer, device, dimension(3, *), intent(in) :: vert_of_face
        integer, device, dimension(*), intent(in) :: face_to_part, iopen
        real(DP), device, dimension(3, *), intent(in) :: tri_bar, tri_nor, vel_tri, acc_tri
        real(DP), device, dimension(*), intent(in) :: sur
        real(DP), device, dimension(n1, n2, ks-lvlhalo:ke+lvlhalo), intent(in) :: q1, q2, q3, pr
        integer, device, dimension(n1, n2, ksr-1:ker), intent(in) :: coll
        real(DP), device, dimension(*), intent(in) :: ym, zm, Hboxx, xc, yc, zc
        real(DP), device, dimension(0:n3), intent(in) :: xm
        integer, device, dimension(*), intent(in) :: mytr
        real(DP), device, dimension(3), intent(in) :: xyz_tr

        real(DP), device, dimension(3, *), intent(inout) :: fpxyz

        integer :: ci, cj, ck
        integer :: i, j, k, ii, jj, inw
        integer :: i1, j1, k1
        integer :: i2, j2
        real(DP), device :: pos_pro(3), norp(3), el_cx(3), ptx(4), Hbox(3), Wt(3), dWt(3)
        real(DP), device :: ui(3)
        real(DP), device :: fnca(3)
        real(DP), device :: B(4), Bi(4)
        real(DP) :: Wtx, dWti(3), rv1, rv2, cfac, prp, epsw
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
        real(DP) :: phiTi, facsign
        real(DP) :: rcxx, rcyy, rczz, appr_ch, rcsqi

        integer :: tid, tx, ty, tr, ntr, inp, cmp, cmp2, v1, v2, v3

        !tid = (blockIdx%x - 1) * blockDim%x + threadIdx%x
        !tx = mod(tid-1, 32) + 1
        !ty = (tx - 1)/32 + 1
        tx = threadIdx%x
        ty = threadIdx%y
        tr = (blockIdx%x - 1) * blockDim%y + threadIdx%y

        facsign =  sign(1.D0,fac) !modificaFV verificare lo vedano tutti
        
        if (tr > trcnt(1)) return
        ntr = mytr(tr)
        inp = face_to_part(ntr)
        if (iopen(inp) .ne. 1) return

        i1 = mod(tx - 1, str1) + 1
        j1 = mod((tx - 1) / str1, str1) + 1
        k1 = (tx - 1) / str2 + 1
        inw = 1 + (i1 - 1) + (j1 - 1) * str1 + (k1 - 1) * str2

!       Compute positions of probes from centroids
        pos_pro(1) = tri_bar(1,ntr) + fac * tri_nor(1,ntr) * Hboxx(pind(6,ntr))
        pos_pro(2) = tri_bar(2,ntr) + fac * tri_nor(2,ntr) * Hboxx(pind(6,ntr))
        pos_pro(3) = tri_bar(3,ntr) + fac * tri_nor(3,ntr) * Hboxx(pind(6,ntr))


        do cmp = 1, 4 ! cmp = 4 is for pressure
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
           case(3) ! 4, 5, 3
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
          epsw = 1e-8
          dWti(1) = 0
          dWti(2) = 0
          dWti(3) = 0

          if (tx <= nel) then

            select case (cmp)
            case (1)
              el_cx(1) = xc(i) - pos_pro(1) + epsw
              el_cx(2) = ym(j) - pos_pro(2) + epsw
              el_cx(3) = zm(k) - pos_pro(3) + epsw
            case (2)
              el_cx(1) = xm(i) - pos_pro(1) + epsw
              el_cx(2) = yc(j) - pos_pro(2) + epsw
              el_cx(3) = zm(k) - pos_pro(3) + epsw
            case (3)
              el_cx(1) = xm(i) - pos_pro(1) + epsw
              el_cx(2) = ym(j) - pos_pro(2) + epsw
              el_cx(3) = zc(k) - pos_pro(3) + epsw
            case (4)
              el_cx(1) = xm(i) - pos_pro(1) + epsw
              el_cx(2) = ym(j) - pos_pro(2) + epsw
              el_cx(3) = zm(k) - pos_pro(3) + epsw
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
          ui(2) = 0.d0
          ui(3) = 0.d0
          prp = 0.d0

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
              prp = pr(i,j,k)
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
#ifdef LUINVERSE             
             call inverseLU_gpu(pinvA(1,ty), invA(1,ty))
#else
             call invert4_gpu(pinvA(1,ty), invA(1,ty))
#endif

             
          endif

          call syncthreads()

          ! Compute gmat (using 16 threads)
          if (tx <= 16) then
            ptx(1) = 1.d0;
            ptx(2) = pos_pro(1)
            ptx(3) = pos_pro(2)
            ptx(4) = pos_pro(3)

            rv1 = ptx(mod(tx-1, 4) + 1) * invA(4*(mod(tx-1, 4)) + (tx-1)/4 + 1, ty) !/ det_l(ty)
          endif

          rv2 = __shfl_down(rv1,1)
          rv1 = rv1 +  rv2
          rv2 = __shfl_down(rv1,2)
          rv1 = rv1 +  rv2

          if (tx <= 16) then
            if (mod(tx-1, 4) == 0) then
              gmat((tx-1)/4 + 1, ty) = rv1
            endif
          endif

          call syncthreads()

          ! Compute gradient terms componentwise
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
          endif

          if (cmp <= 3) then
            do cmp2 = 1, 3
              if (tx <= nel) then
                Bi(1) = dWti(cmp2) * ptx(1)
                Bi(2) = dWti(cmp2) * ptx(2)
                Bi(3) = dWti(cmp2) * ptx(3)
                Bi(4) = dWti(cmp2) * ptx(4)
              endif

              ! each thread does outer product
              do j2 = 1, 4
                do i2 = 1, 4
                  rv1 = dWti(cmp2)*ptx(i2)*ptx(j2)

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

              ! Compute gmat2 (using 16 threads)
              if (tx <= 16) then
                rv1 = gmat(mod(tx-1, 4) + 1, ty) * pinvA(4*(mod(tx-1, 4)) + (tx-1)/4 + 1, ty)
              endif

              rv2 = __shfl_down(rv1,1)
              rv1 = rv1 +  rv2
              rv2 = __shfl_down(rv1,2)
              rv1 = rv1 +  rv2

              if (tx <= 16) then
                if (mod(tx-1, 4) == 0) then
                  gmat2((tx-1)/4 + 1, ty) = rv1
                endif
              endif

              call syncthreads()

              ! Compute gmat3 (using 16 threads)
              if (tx <= 16) then
                rv1 = gmat2(mod(tx-1, 4) + 1, ty) * invA(4*(mod(tx-1, 4)) + (tx-1)/4 + 1, ty) !/ det_l(ty)
              endif

              rv2 = __shfl_down(rv1,1)
              rv1 = rv1 +  rv2
              rv2 = __shfl_down(rv1,2)
              rv1 = rv1 +  rv2

              if (tx <= 16) then
                if (mod(tx-1, 4) == 0) then
                  gmat3((tx-1)/4 + 1, ty) = rv1
                endif
              endif

              call syncthreads()

              ! Compute gmati (using 4 threads)
              if (tx <= 4) then
                 ! gmati(tx,ty) = invA(tx + 4 * cmp2, ty) / det_l(ty) - gmat3(tx,ty)
                 gmati(tx,ty) = invA(tx + 4 * cmp2, ty) - gmat3(tx,ty)                
              endif

              call syncthreads()

              ! Compute phiTi
              phiTi = 0.d0
              if (tx <= nel) then
                do i2 = 1, 4
                  phiTi = phiTi + gmati(i2, ty) * B(i2) + gmat(i2, ty) * Bi(i2)
                enddo
              endif

              ! Compute gradients
              do i2 = 1, 1
                if (tx <= nel) then
                  rv1 = phiTi * ui(i2)
                else
                  rv1 = 0.d0
                endif

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

                if (tx == 1) then
                  !du(i1 + (cmp2-1)*3, ty) = rv1
                  du(cmp + (cmp2-1)*3, ty) = rv1
                endif

                call syncthreads()

              enddo

            enddo

          else 

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
            if (tx <= nel) then
              rv1 = B(1) * prp
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

            prp = rv1
          endif
        end do

        ! complete computation with 3 threads per warp
        if (tx <= 3) then
          for = 0

          ! Face normals
          fnca(1) = tri_nor(1, ntr)
          fnca(2) = tri_nor(2, ntr)
          fnca(3) = tri_nor(3, ntr)

          press = prp + fac * Hboxx(pind(6,ntr))*(acc_tri(1,ntr)*fnca(1) + &
              acc_tri(2,ntr)*fnca(2)+acc_tri(3,ntr)*fnca(3))

          i2 = tx
          j2 = mod(tx, 3) + 1
          du(i2 + (j2-1)*3, ty) = 0.5d0*(du(i2 + (j2-1)*3, ty) + du(j2 + (i2-1)*3, ty))
          call syncthreads()

          du(j2 + (i2-1)*3, ty) = du(i2 + (j2-1)*3, ty)
          call syncthreads()

          tau = 2.d0*nu*(du(tx, ty)*fnca(1) + du(tx+3, ty)*fnca(2) + du(tx+6, ty)*fnca(3))

          for = for + facsign * (tau - press*fnca(tx)) * (sur(ntr)/3.0)

          v1 = vert_of_face(1,ntr)
          v2 = vert_of_face(2,ntr)
          v3 = vert_of_face(3,ntr)

!         Calculate forces on triangle faces - pressure+viscous
          rv1 = atomicadd(fpxyz(tx,v1),for)
          rv1 = atomicadd(fpxyz(tx,v2),for)
          rv1 = atomicadd(fpxyz(tx,v3),for)
        endif

      end subroutine compute_fpxyz_kernel

#endif

end module mlsStruc3CompStrong_m
!-----------------------------------------------------------------------
!     read in position of a individual trial marker for MLS and compute
!     compute support domain, shape function and interpolate
!
!     compute internal forces on the nodes
!------------------------------------------------------------------------

      subroutine mlsStruc3CompStrong
      USE mpih
      USE param
      USE mls_param
      USE local_arrays, only: q1,q2,q3
      USE mpi_param, only: kstart, kend, kstartr, kendr
      USE mls_local, only: coll
      USE local_arrays, only: pr
      USE ieee_arithmetic
      use mlsStruc3CompStrong_m
      use nvtx
#ifdef NCCLAVAIL
      use nccl
#endif            
!@cuf use cudafor
      IMPLICIT NONE

      integer :: inp, ntr
      integer :: siz, stride(2)
      integer :: i, j, k, tr

      real(DP) :: pos_MLS(3)

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
      real(DP) :: fsmth, fsmth2, xV, yV, zV, xV1, yV1, zV1, rV1
      real(DP) :: dinvm,ciao1,ciao2,ciao3

!     --------------------------------------------------------

#ifdef USE_CUDA
      type(dim3) :: blocks, threads
#endif
!@cuf integer :: istat


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

      fsmth  = 1.0
      fsmth2 = 0.0

      !fpxyz(:,:) = 0.d0
      fpxyz = 0.d0
      fpxyz_3d = 0.d0

      dtr = dt/dt_o

      call nvtxStartRange("compute_fpxyz", 9)
#ifdef USE_CUDA
      threads = dim3(32, 4, 1)
      blocks = dim3(ceiling(real(nftot)/4), 1, 1)

      ! Surface forces
      ! outer (positive) probe component
      call compute_fpxyz_kernel<<<blocks, threads>>>(siz, stride(1), stride(2), &
      nel, dx1, dx2, dx3, wscl, wcon, wcub, wexp, n1m, n2m, &
      nu, fpxyz, pind, pind_probeP,tri_bar, tri_nor, vel_tri, acc_tri, sur, coll, &
      q1, q2, q3, pr, vert_of_face, face_to_part, mytr, trcnt, xm, ym, zm, xc, yc, zc, &
      Hboxx, nftot, nvtot, n1, n2, n3, kstart, kend, kstartr, kendr, lvlhalo, &
      h_probe, iopen_pos, xyz_tr)

      ! inner (negative) probe component
      call compute_fpxyz_kernel<<<blocks, threads>>>(siz, stride(1), stride(2), &
      nel, dx1, dx2, dx3, wscl, wcon, wcub, wexp, n1m, n2m, &
      nu, fpxyz, pind, pind_probeN, tri_bar, tri_nor, vel_tri, acc_tri, sur, coll, &
      q1, q2, q3, pr, vert_of_face, face_to_part, mytr, trcnt, xm, ym, zm, xc, yc, zc, &
      Hboxx, nftot, nvtot, n1, n2, n3, kstart, kend, kstartr, kendr, lvlhalo, &
      -h_probe, iopen_neg, xyz_tr)

#else

      ! Surface forces
      ! outer (positive) probe component
      call compute_fpxyz(h_probe, iopen_pos, siz, stride, pind_probeP)
      ! inner (negative) probe component
      call compute_fpxyz(-h_probe, iopen_neg, siz, stride, pind_probeN)
      !print*, "fpxyz", sum(fpxyz), minval(fpxyz), maxval(fpxyz)
#endif
      call nvtxEndRange

!@cuf istat = cudaDeviceSynchronize !JDR TMP

!     --------------------------------------------------------
!     Reduce the forces from all processors over each particle
      if (numtasks > 1) then
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#ifndef NCCLAVAIL                  
        call MPI_ALLREDUCE(MPI_IN_PLACE,fpxyz,3*nvtot,MPI_DOUBLE_PRECISION, &
              MPI_SUM,MPI_COMM_WORLD,ierr)
#else
        nccl_result = ncclAllReduce(fpxyz, fpxyz, 3*nvtot, ncclDouble, ncclSum, nccl_comm, 0)
#endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      endif

!===================================
!=============3D====================
!===================================

!@cuf istat = cudaDeviceSynchronize !JDR TMP

!Overwrite tagged wet boundaries 
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do ntr = nvstart_2dwet,nvend_2dwet
         v1 = tag_2dwet(ntr)
         fpxyz_3d(1,v1) = fpxyz(1,ntr)
         fpxyz_3d(2,v1) = fpxyz(2,ntr)
         fpxyz_3d(3,v1) = fpxyz(3,ntr)
      enddo
!@cuf istat = cudaDeviceSynchronize !JDR TMP

      return
      end
