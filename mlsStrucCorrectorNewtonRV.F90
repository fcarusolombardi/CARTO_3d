      subroutine mlsStrucCorrectorNewtonRV
      USE mpih
      USE param
      USE mls_param
      USE local_arrays, only: q1,q2,q3
      USE mpi_param, only: kstart, kend, kstartr, kendr
      USE mls_local, only: coll
      USE local_arrays, only: pr
      USE local_arrays, only: q1g,q2g,q3g,prg
      USE ieee_arithmetic
      use nvtx
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

      integer :: imasZ,ic,jc,kc
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


      dtr = dt/dt_o



!store stuff
#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = 1, nvtot_3d
         dinvm = 1.d0/mass_of_vert_3d(ntr)
         xyzakp1_3d(1,ntr)=(fpxyz_3d(1,ntr)+fxyz_3d(1,ntr))*dinvm
         xyzakp1_3d(2,ntr)=(fpxyz_3d(2,ntr)+fxyz_3d(2,ntr))*dinvm
         xyzakp1_3d(3,ntr)=(fpxyz_3d(3,ntr)+fxyz_3d(3,ntr))*dinvm
      enddo

#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = nvstart_2dwet,nvend_2dwet
         v1 = tag_2dwet(ntr)
         xyzakp1(1,ntr)=xyzakp1_3d(1,v1)
         xyzakp1(2,ntr)=xyzakp1_3d(2,v1)
         xyzakp1(3,ntr)=xyzakp1_3d(3,v1)
      enddo

#ifdef USE_CUDA
      !$cuf kernel do (1)                                                                   
#endif
      do ntr = nvstart_2dstr,nvend_2dstr
         dinvm = 1.d0/mass_of_vert(ntr)
         xyzakp1(1,ntr)=(fpxyz(1,ntr)+fxyz(1,ntr))*dinvm
         xyzakp1(2,ntr)=(fpxyz(2,ntr)+fxyz(2,ntr))*dinvm
         xyzakp1(3,ntr)=(fpxyz(3,ntr)+fxyz(3,ntr))*dinvm
      enddo
!end store stuff





      !====== Corrector ======================================================== 
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

!     --------------------------------------------------------
!     Modify pres+vis forces with relevant prefactors
      ! store old positions before update
      !xyzold = xyz
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do ntr = 1, nvtot_3d
        xyzold_3d(1, ntr) = xyz_3d(1, ntr)
        xyzold_3d(2, ntr) = xyz_3d(2, ntr)
        xyzold_3d(3, ntr) = xyz_3d(3, ntr)
      enddo

!@cuf istat = cudaDeviceSynchronize !JDR TMP

!     hard-coded now; should be cleaned later;
!     apex position + offset
      xV = 1.129
      yV = 0.0
      zV = 3.321

      call nvtxStartRange("update_xyz", 10)
      ! JR TODO: Transpose xyz arrays so that first dim is last for better memory accesses

#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do ntr = 1, nvtot_3d
        inp = vert_to_part_3d(ntr)

        if (inp .eq. 1) then
          xV1 = xyz_3d(1,ntr)-xV
          yV1 = xyz_3d(2,ntr)-yV
          zV1 = xyz_3d(3,ntr)-zV

          rV1 = sqrt(xV1**2 + zV1**2 + yV1**2)
          fsmth = -0.5d0 * (tanh(5.d0 * (rV1 - 3.20d0)) + 1.0d0) + 1.0d0
!pianovalvmobile fsmooth = -0.5d0*( tanh(10.d0*(xyz(3,i)-2.50d0))+1.d0)+1.d0 !FV cambia quando sblocco piano valvolare 
        else
          fsmth = 1.d0
        end if


        xyzauk_3d(1,ntr)=xyzakp1_3d(1,ntr)
        xyzauk_3d(2,ntr)=xyzakp1_3d(2,ntr)
        xyzauk_3d(3,ntr)=xyzakp1_3d(3,ntr)


        xyzvkp1_3d(1,ntr) = fsmth * (xyzv_3d(1,ntr) + ( (xyzauk_3d(1,ntr) + xyza_3d(1,ntr) )* (0.5d0*dtr)) * dt)
        xyzvkp1_3d(2,ntr) = fsmth * (xyzv_3d(2,ntr) + ( (xyzauk_3d(2,ntr) + xyza_3d(2,ntr) )* (0.5d0*dtr)) * dt)
        xyzvkp1_3d(3,ntr) = SMZ*fsmth * (xyzv_3d(3,ntr) + ( (xyzauk_3d(3,ntr) + xyza_3d(3,ntr) )* (0.5d0*dtr)) * dt)


        xyzvuk_3d(1,ntr)=xyzvkp1_3d(1,ntr)
        xyzvuk_3d(2,ntr)=xyzvkp1_3d(2,ntr)
        xyzvuk_3d(3,ntr)=xyzvkp1_3d(3,ntr)


        xyzkp1_3d(1,ntr) = xyz_3d(1,ntr) + (1.0-undr)*( (xyzvuk_3d(1,ntr) + xyzv_3d(1,ntr) )* (0.5d0*dtr)) * dt + undr*xyzk_3d(1,ntr)
        xyzkp1_3d(2,ntr) = xyz_3d(2,ntr) + (1.0-undr)*( (xyzvuk_3d(2,ntr) + xyzv_3d(2,ntr) )* (0.5d0*dtr)) * dt + undr*xyzk_3d(2,ntr)
        xyzkp1_3d(3,ntr) = xyz_3d(3,ntr) + (1.0-undr)*( (xyzvuk_3d(3,ntr) + xyzv_3d(3,ntr) )* (0.5d0*dtr)) * dt + undr*xyzk_3d(3,ntr)



      enddo

!===================================
!===========BC 3D===================
!===================================
!boundary condition 3d
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do i = 1,count2_3d
        v1 = boundary2_3d(1,i)
        v2 = boundary2_3d(2,i) 
        if (v2.LT.0) then
           xyzauk_3d(1,v1) = 0.d0
           xyzauk_3d(2,v1) = 0.d0
           xyzauk_3d(3,v1) = 0.d0
           xyzvuk_3d(1,v1) = 0.d0
           xyzvuk_3d(2,v1) = 0.d0
           xyzvuk_3d(3,v1) = 0.d0
           xyzvkp1_3d(1,v1) = 0.d0
           xyzvkp1_3d(2,v1) = 0.d0
           xyzvkp1_3d(3,v1) = 0.d0
           xyzkp1_3d(1,v1) = xyz0_3d(1,v1)
           xyzkp1_3d(2,v1) = xyz0_3d(2,v1)
           xyzkp1_3d(3,v1) = xyz0_3d(3,v1)
        else !it shouldn't enter here
           xyzauk_3d(1,v1) = xyzauk_3d(1,v2)
           xyzauk_3d(2,v1) = xyzauk_3d(2,v2)
           xyzauk_3d(3,v1) = xyzauk_3d(3,v2)
           xyzvuk_3d(1,v1) = xyzvuk_3d(1,v2)
           xyzvuk_3d(2,v1) = xyzvuk_3d(2,v2)
           xyzvuk_3d(3,v1) = xyzvuk_3d(3,v2)
           xyzvkp1_3d(1,v1) = xyzvkp1_3d(1,v2)
           xyzvkp1_3d(2,v1) = xyzvkp1_3d(2,v2)
           xyzvkp1_3d(3,v1) = xyzvkp1_3d(3,v2)
           xyzkp1_3d(1,v1) = xyzkp1_3d(1,v2)+BCoffset_3d(1,i)   
           xyzkp1_3d(2,v1) = xyzkp1_3d(2,v2)+BCoffset_3d(2,i)   
           xyzkp1_3d(3,v1) = xyzkp1_3d(3,v2)+BCoffset_3d(3,i)   
        endif
      end do
!--------------------------------------------------------------------

!@cuf istat = cudaDeviceSynchronize !JDR TMP
 
!===================================
!=============2D====================
!===================================
!    --------------------------------------------------------
!    Modify pres+vis forces with relevant prefactors
!    store old positions before update
      !xyzold = xyz
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do ntr = nvstart_2dstr,nvend_2dstr
        xyzold(1, ntr) = xyz(1, ntr)
        xyzold(2, ntr) = xyz(2, ntr)
        xyzold(3, ntr) = xyz(3, ntr)
      enddo

!@cuf istat = cudaDeviceSynchronize !JDR TMP

!     hard-coded now; should be cleaned later;
!     apex position + offset
      xV = 1.129
      yV = 0.0
      zV = 3.321

      call nvtxStartRange("update_xyz", 10)
      ! JR TODO: Transpose xyz arrays so that first dim is last for better memory accesses

#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do ntr = nvstart_2dstr,nvend_2dstr
        inp = vert_to_part(ntr)

         if (inp .eq. 1) then
           xV1 = xyz(1,ntr)-xV
           yV1 = xyz(2,ntr)-yV
           zV1 = xyz(3,ntr)-zV

!           rV1 = sqrt(xV1**2 + zV1**2)
           rV1 = sqrt(xV1**2 + zV1**2 +yV1**2)
           fsmth = -0.5d0 * (tanh(5.d0 * (rV1 - 3.20d0)) + 1.0d0) + 1.0d0
! ! !pianovalvmobile fsmooth = -0.5d0*( tanh(10.d0*(xyz(3,i)-2.50d0))+1.d0)+1.d0 !FV cambia quando sblocco piano valvolare 
         elseif (inp .eq. 8) then
           fsmth = 0.5d0 * (tanh(5.d0 * (-xyz(3,ntr) - 0.1d0)) + 1.d0)
         else
          fsmth = 1.d0
         end if

        xyzauk(1,ntr)=xyzakp1(1,ntr)
        xyzauk(2,ntr)=xyzakp1(2,ntr)
        xyzauk(3,ntr)=xyzakp1(3,ntr)


        xyzvkp1(1,ntr) = fsmth * (xyzv(1,ntr) + ( (xyzauk(1,ntr) + xyza(1,ntr) )* (0.5d0*dtr)) * dt)
        xyzvkp1(2,ntr) = fsmth * (xyzv(2,ntr) + ( (xyzauk(2,ntr) + xyza(2,ntr) )* (0.5d0*dtr)) * dt)
        xyzvkp1(3,ntr) = fsmth * (xyzv(3,ntr) + ( (xyzauk(3,ntr) + xyza(3,ntr) )* (0.5d0*dtr)) * dt)


        xyzvuk(1,ntr)=xyzvkp1(1,ntr)
        xyzvuk(2,ntr)=xyzvkp1(2,ntr)
        xyzvuk(3,ntr)=xyzvkp1(3,ntr)


        xyzkp1(1,ntr) = xyz(1,ntr) + (1.0-undr)*( (xyzvuk(1,ntr) + xyzv(1,ntr) )* (0.5d0*dtr)) * dt + undr*xyzk(1,ntr)
        xyzkp1(2,ntr) = xyz(2,ntr) + (1.0-undr)*( (xyzvuk(2,ntr) + xyzv(2,ntr) )* (0.5d0*dtr)) * dt + undr*xyzk(2,ntr)
        xyzkp1(3,ntr) = xyz(3,ntr) + (1.0-undr)*( (xyzvuk(3,ntr) + xyzv(3,ntr) )* (0.5d0*dtr)) * dt + undr*xyzk(3,ntr)

      enddo

!@cuf istat = cudaDeviceSynchronize !JDR TMP

!===================================
!===========BC 2D===================
!===================================
!boundary condition 2d
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do i = 1,count2
        v1 = boundary2(1,i)
        v2 = boundary2(2,i) 
        if (v2.LT.0) then
           if (v2.EQ.-1) then
           xyzauk(1,v1) = 0.d0
           xyzauk(2,v1) = 0.d0
           xyzauk(3,v1) = 0.d0
           xyzvuk(1,v1) = 0.d0
           xyzvuk(2,v1) = 0.d0
           xyzvuk(3,v1) = 0.d0
           xyzvkp1(1,v1) = 0.d0
           xyzvkp1(2,v1) = 0.d0
           xyzvkp1(3,v1) = 0.d0
           xyzkp1(1,v1) = xyz0(1,v1)
           xyzkp1(2,v1) = xyz0(2,v1)
           xyzkp1(3,v1) = xyz0(3,v1)
           elseif (v2.EQ.-2) then
           xyzauk(3,v1) = 0.d0
           xyzvkp1(3,v1) = 0.d0
           xyzkp1(3,v1) = xyz0(3,v1)
           endif
        else
           xyzauk(1,v1) = xyzauk(1,v2) !now slaved again to 2d
           xyzauk(2,v1) = xyzauk(2,v2)
           xyzauk(3,v1) = xyzauk(3,v2)
           xyzvuk(1,v1) = xyzvuk(1,v2) !now slaved again to 2d
           xyzvuk(2,v1) = xyzvuk(2,v2)
           xyzvuk(3,v1) = xyzvuk(3,v2)
           xyzvkp1(1,v1) = xyzvkp1(1,v2) !now slaved again to 2d
           xyzvkp1(2,v1) = xyzvkp1(2,v2)
           xyzvkp1(3,v1) = xyzvkp1(3,v2)
           xyzkp1(1,v1) = xyzkp1(1,v2)+BCoffset(1,i)   
           xyzkp1(2,v1) = xyzkp1(2,v2)+BCoffset(2,i)   
           xyzkp1(3,v1) = xyzkp1(3,v2)+BCoffset(3,i)   
        endif
      end do
!--------------------------------------------------------------------

!@cuf istat = cudaDeviceSynchronize !JDR TMP

!Overwrite tagged 2dwet
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do ntr = nvstart_2dwet,nvend_2dwet
         v1 = tag_2dwet(ntr)
         xyzkp1(1,ntr) = xyzkp1_3d(1,v1)
         xyzkp1(2,ntr) = xyzkp1_3d(2,v1)
         xyzkp1(3,ntr) = xyzkp1_3d(3,v1)
         xyzvkp1(1,ntr) = xyzvkp1_3d(1,v1)
         xyzvkp1(2,ntr) = xyzvkp1_3d(2,v1)
         xyzvkp1(3,ntr) = xyzvkp1_3d(3,v1)
         xyzvuk(1,ntr) = xyzvuk_3d(1,v1) !ci va questo?
         xyzvuk(2,ntr) = xyzvuk_3d(2,v1) !ci va questo!!!!
         xyzvuk(3,ntr) = xyzvuk_3d(3,v1) !ci va questo!!!!
         xyzauk(1,ntr) = xyzauk_3d(1,v1) !ci va questo!!!!
         xyzauk(2,ntr) = xyzauk_3d(2,v1) !ci va questo!!!!
         xyzauk(3,ntr) = xyzauk_3d(3,v1) !ci va questo!!!!
      enddo

!@cuf istat = cudaDeviceSynchronize !JDR TMP

!===================================
!========CONTACT MODEL==============
!===================================      
      call nvtxEndRange

      ! Collision operations
      call nvtxStartRange("contact", 16)
#ifdef CONTACT
      call ContactModel(1,nvtot,xyzkp1,xyzvkp1)
#endif
#ifdef CONTACTV
      call ContactModelV(1,nvtot,xyzkp1,xyzvkp1)
#endif
      call nvtxEndRange


!@cuf istat = cudaDeviceSynchronize !JDR TMP

 


      return
      end
