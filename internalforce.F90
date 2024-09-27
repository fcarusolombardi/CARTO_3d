module internalforce_m
  use constants
  contains

#ifdef USE_CUDA

      attributes(device) &
      subroutine  cross_gpu(xcy,x,y)
        implicit none
        real(DP), intent(in), device ::  x(3),y(3)
        real(DP), intent(out), device ::  xcy(3)

        xcy(1) = x(2)*y(3)-y(2)*x(3)
        xcy(2) = x(3)*y(1)-y(3)*x(1)
        xcy(3) = x(1)*y(2)-y(1)*x(2)
        return
      end subroutine cross_gpu

      attributes(global) &
      subroutine compute_favxyz_kernel(fstart, fend, kat, kal, kv, kae, Surface, Surface0, Volume, Volume0, &
        xyzL, vert_of_face, face_to_part, tri_nor, sur, sur0, aalpha, aalpha0, nfi, fxyz)
        use mls_param, only: StiffV
        use param, only: nbeat
        implicit none

        integer, value, intent(in) :: fstart, fend
        integer, device, dimension(3,*), intent(in) :: vert_of_face
        integer, device, dimension(*), intent(in) :: face_to_part
        real(DP), device, dimension(*), intent(in) :: Surface, Surface0, Volume, Volume0
        real(DP), device, dimension(*), intent(in) :: kat, kal, kv, kae
        integer, device, dimension(*), intent(in) :: nfi
        real(DP), device, dimension(*), intent(in) :: sur, sur0
        real(DP), device, dimension(3,*), intent(in) :: xyzL
        real(DP), device, dimension(3,*), intent(in) :: tri_nor, aalpha, aalpha0
        real(DP), device, dimension(3,*), intent(inout) :: fxyz
        integer :: v1, v2, v3
        real(DP), device :: a32(3), a13(3), a21(3), csi(3), tcv(3)
        real(DP), device :: tvec1(3), tvec2(3), tvec3(3), val(3)
        real(DP) :: alphat, alphal, betav, rv1, rv2, rv3, surT
        real(DP) :: alphae1, alphae2, alphae3
        real(DP) :: smoothv1,smoothv2,smoothv3

        integer :: i, inp

        i = (blockIdx%x - 1) * blockDim%x + threadIdx%x + fstart - 1

        if (i > fend) return

        inp = face_to_part(i)

        v1 = vert_of_face(1,i)
        v2 = vert_of_face(2,i)
        v3 = vert_of_face(3,i)

        a32(1:3) = xyzL(1:3,v3)-xyzL(1:3,v2)
        a13(1:3) = xyzL(1:3,v1)-xyzL(1:3,v3)
        a21(1:3) = xyzL(1:3,v2)-xyzL(1:3,v1)

        ! local, area, and angle potential
        csi(1:3) = tri_nor(1:3,i)

        alphat = -kat(inp) * (Surface(inp)-Surface0(inp)) / Surface0(inp)
        alphal = -kal(inp) * (sur(i)-sur0(i)) / sur0(i)

        alphae1 = kae(inp) * (aalpha(1,i) - aalpha0(1,i))
        alphae2 = kae(inp) * (aalpha(2,i) - aalpha0(2,i))
        alphae3 = kae(inp) * (aalpha(3,i) - aalpha0(3,i))

        call cross_gpu(tvec1,csi,a32)
        rv1 = alphat*tvec1(1) +  alphal*tvec1(1) + alphae2*a13(1) - alphae3*a21(1)
        rv1 = atomicadd(fxyz(1,v1),rv1)

        rv2 = alphat*tvec1(2) +  alphal*tvec1(2) + alphae2*a13(2) - alphae3*a21(2)
        rv2 = atomicadd(fxyz(2,v1),rv2)

        rv3 = alphat*tvec1(3) +  alphal*tvec1(3) + alphae2*a13(3) - alphae3*a21(3)
        rv3 = atomicadd(fxyz(3,v1),rv3)

        call cross_gpu(tvec2,csi,a13)
        rv1 = alphat*tvec2(1) + alphal*tvec2(1) - alphae1*a32(1) + alphae3*a21(1)
        rv1 = atomicadd(fxyz(1,v2), rv1)

        rv2 = alphat*tvec2(2) + alphal*tvec2(2) - alphae1*a32(2) + alphae3*a21(2)
        rv2 = atomicadd(fxyz(2,v2), rv2)

        rv3 = alphat*tvec2(3) + alphal*tvec2(3) - alphae1*a32(3) + alphae3*a21(3)
        rv3 = atomicadd(fxyz(3,v2), rv3)

        call cross_gpu(tvec3,csi,a21)
        rv1 = alphat*tvec3(1) + alphal*tvec3(1) + alphae1*a32(1) - alphae2*a13(1)
        rv1 = atomicadd(fxyz(1,v3),rv1)

        rv2 = alphat*tvec3(2) + alphal*tvec3(2) + alphae1*a32(2) - alphae2*a13(2)
        rv2 = atomicadd(fxyz(2,v3),rv2)

        rv3 = alphat*tvec3(3) + alphal*tvec3(3) + alphae1*a32(3) - alphae2*a13(3)
        rv3 = atomicadd(fxyz(3,v3),rv3)


!         ! volume potential
!         if (inp .eq. 1 .or. inp .eq. 8) then ! test
!           csi(1:3) = 0.d0 !csi(1:3) / 3.d0 !FV se serve parte normale va cambiato con csi/3

!           tcv(1:3) = (xyzL(1:3,v1)+xyzL(1:3,v2)+ &
!                       xyzL(1:3,v3))/3.d0

!           !betav = -kv*(Volume(inp)-Volume0(inp)) / Volume0(inp)
! !          betav = -kv(inp)*(pot(v1)+pot(v2)+pot(v3))/3.
! !          betav = betav / 6.d0

!           smoothv1 = (1.d0-StiffV(v1))
!           smoothv2 = (1.d0-StiffV(v2))
!           smoothv3 = (1.d0-StiffV(v3))

!           call cross_gpu(tvec1,tcv,a32)
!           rv1 = atomicadd(fxyz(1,v1),smoothv1*betav*(csi(1) + tvec1(1)))
!           rv1 = atomicadd(fxyz(2,v1),smoothv1*betav*(csi(2) + tvec1(2)))
!           rv1 = atomicadd(fxyz(3,v1),smoothv1*betav*(csi(3) + tvec1(3)))

!           call cross_gpu(tvec2,tcv,a13)
!           rv1 = atomicadd(fxyz(1,v2),smoothv2*betav*(csi(1) + tvec2(1)))
!           rv1 = atomicadd(fxyz(2,v2),smoothv2*betav*(csi(2) + tvec2(2)))
!           rv1 = atomicadd(fxyz(3,v2),smoothv2*betav*(csi(3) + tvec2(3)))

!           call cross_gpu(tvec3,tcv,a21)
!           rv1 = atomicadd(fxyz(1,v3),smoothv3*betav*(csi(1) + tvec3(1)))
!           rv1 = atomicadd(fxyz(2,v3),smoothv3*betav*(csi(2) + tvec3(2)))
!           rv1 = atomicadd(fxyz(3,v3),smoothv3*betav*(csi(3) + tvec3(3)))
!         endif


      end subroutine compute_favxyz_kernel

      attributes(global) &
      subroutine compute_fbxyz_kernel(estart, eend, kb, face_of_edge, &
        v1234, xyzL, theta, theta0, fxyz)
        implicit none

        integer, value, intent(in) :: estart, eend
        integer, device, dimension(2,*), intent(in) :: face_of_edge
        integer, device, dimension(4,*), intent(in) :: v1234
        real(DP), device, dimension(3,*), intent(in) :: xyzL
        real(DP), device, dimension(*), intent(in) :: theta, theta0, kb
        real(DP), device, dimension(3,*), intent(inout) :: fxyz
        integer :: v1, v2, v3, v4, f1, f2
        real(DP), device :: a21(3), a31(3), a34(3), a24(3)
        real(DP), device :: a32(3), a13(3), a42(3), a23(3)
        real(DP), device :: csi(3), zet(3)
        real(DP), device :: tvec1(3), tvec2(3), tvec3(3), tvec4(3)
        real(DP) :: rv1, modcsi, modzet, betab, b11, b12, b22

        integer :: i

        i = (blockIdx%x - 1) * blockDim%x + threadIdx%x + estart - 1


        if (i > eend) return

        f1 = face_of_edge(1,i)
        f2 = face_of_edge(2,i)

        if (f1 .ne. 0.and. f2 .ne. 0) then

          v1 = v1234(1,i)
          v2 = v1234(2,i)
          v3 = v1234(3,i)
          v4 = v1234(4,i)

          a21(1:3) = xyzL(1:3,v2) - xyzL(1:3,v1)
          a31(1:3) = xyzL(1:3,v3) - xyzL(1:3,v1)
          a34(1:3) = xyzL(1:3,v3) - xyzL(1:3,v4)
          a24(1:3) = xyzL(1:3,v2) - xyzL(1:3,v4)

          a32(1:3) = xyzL(1:3,v3) - xyzL(1:3,v2)
          a13(1:3) = xyzL(1:3,v1) - xyzL(1:3,v3)
          a42(1:3) = xyzL(1:3,v4) - xyzL(1:3,v2)
          a23(1:3) = xyzL(1:3,v2) - xyzL(1:3,v3)

          call cross_gpu(csi,a21,a31)
          call cross_gpu(zet,a34,a24)

          modcsi = sqrt(csi(1)*csi(1) + csi(2)*csi(2) + csi(3)*csi(3))
          modzet = sqrt(zet(1)*zet(1) + zet(2)*zet(2) + zet(3)*zet(3))


          if ((1. - cos(theta(i))**2) .le. 1.e-8) then
            betab = 0.0
          else
            betab = kb(i)*(sin(theta(i))*cos(theta0(i))- &
                    cos(theta(i))*sin(theta0(i)))/       &
                    (sqrt(1.-(cos(theta(i))**2)))
          endif

          b11 = -betab*cos(theta(i))/(modcsi**2)
          b12 = betab/(modcsi*modzet)
          b22 = -betab*cos(theta(i))/(modzet**2)

          call cross_gpu(tvec1,csi,a32)
          call cross_gpu(tvec2,zet,a32)

          !fxyz(1:3,v1,inp) = fxyz(1:3,v1,inp) + b11*tvec1(1:3) + b12*tvec2(1:3)
          tvec1(1:3) = b11*tvec1(1:3) + b12*tvec2(1:3)
          rv1 = atomicadd(fxyz(1,v1), tvec1(1))
          rv1 = atomicadd(fxyz(2,v1), tvec1(2))
          rv1 = atomicadd(fxyz(3,v1), tvec1(3))

          call cross_gpu(tvec1,csi,a13)
          call cross_gpu(tvec2,csi,a34)
          call cross_gpu(tvec3,zet,a13)
          call cross_gpu(tvec4,zet,a34)

          !fxyz(1:3,v2,inp) = fxyz(1:3,v2,inp) + b11*tvec1(1:3)+ b12*(tvec2(1:3)+tvec3(1:3))+ b22*tvec4(1:3)
          tvec1(1:3) = b11*tvec1(1:3)+ b12*(tvec2(1:3)+tvec3(1:3))+ b22*tvec4(1:3)
          rv1 = atomicadd(fxyz(1,v2), tvec1(1))
          rv1 = atomicadd(fxyz(2,v2), tvec1(2))
          rv1 = atomicadd(fxyz(3,v2), tvec1(3))

          call cross_gpu(tvec1,csi,a21)
          call cross_gpu(tvec2,csi,a42)
          call cross_gpu(tvec3,zet,a21)
          call cross_gpu(tvec4,zet,a42)

          !fxyz(1:3,v3,inp) = fxyz(1:3,v3,inp) + b11*tvec1(1:3) + b12*(tvec2(1:3)+tvec3(1:3)) + b22*tvec4(1:3)
          tvec1(1:3) =  b11*tvec1(1:3) + b12*(tvec2(1:3)+tvec3(1:3)) + b22*tvec4(1:3)
          rv1 = atomicadd(fxyz(1,v3), tvec1(1))
          rv1 = atomicadd(fxyz(2,v3), tvec1(2))
          rv1 = atomicadd(fxyz(3,v3), tvec1(3))


          call cross_gpu(tvec1,csi,a23)
          call cross_gpu(tvec2,zet,a23)

          !fxyz(1:3,v4,inp) = fxyz(1:3,v4,inp) + b12*tvec1(1:3) + b22*tvec2(1:3)
          tvec1(1:3) = b12*tvec1(1:3) + b22*tvec2(1:3)
          rv1 = atomicadd(fxyz(1,v4), tvec1(1))
          rv1 = atomicadd(fxyz(2,v4), tvec1(2))
          rv1 = atomicadd(fxyz(3,v4), tvec1(3))
        endif


      end subroutine compute_fbxyz_kernel
#endif

end module internalforce_m

      subroutine internalforce(snv,env,xyzL,xyzvL)
      use mpih
      use param
      use mls_param
      use internalforce_m
      use nvtx
!@cuf use cudafor
      implicit none

      integer :: inp,i,j,k
      integer :: imasZ, snv,env
      integer :: v1,v2,v3,v4,f1,f2,f3,iv,ie
      integer :: nverti_2dstr,nedges_2dstr,nfaces_2dstr      
      real(DP),dimension(3,snv:env),intent(in) :: xyzL,xyzvL
      real(DP),dimension(3)::tvec1,tvec2,tvec3,tvec4,csi,zet,Vij
      real(DP),dimension(3)::a32,a13,a34,a21,a42,a23,a31,a24,csi1,zet1,tcv
      real(DP)::modcsi,modzet,betab,betav,alphat,alphal,b11,b12,b22,tdum,VV,alphae
 
      real(DP) :: angleHelix, angleSpiral, masZ
      real(DP) :: xCC, yCC, zCC, surT
      real(DP) :: zMax
      real(DP) :: alphaES,betaES,gammaES,smoothf,distMEAN,CosePhi,SenoPhi
      real(DP) :: visc_coeff, visc_coeff2, epsG, FungNL, cstiff
      real(DP) :: d, d0, kei, NUPOIS
      real(DP) :: Vij1, Vij2, Vij3

      integer :: ne_per_rank, my_es, my_ee
      integer :: nv_per_rank, my_vs, my_ve
      integer :: nf_per_rank, my_fs, my_fe

      real(DP), dimension(3) :: xyzCC
#ifdef USE_CUDA
      attributes(managed) :: xyzCC
      attributes(managed) :: xyzL
      attributes(managed) :: xyzvL
      type(dim3) :: blocks, threads
#endif
!@cuf integer :: istat



      call nvtxStartRange("geo_update", 11)

      ! Just calculate all geometry properties first
      call calculate_distance(dist(1:nvend_2dstr), 1, nvend_2dstr, 1, neend_2dstr, &
               xyzL(:,1:nvend_2dstr), vert_of_edge(:,1:neend_2dstr))

      !potrei farle calcolare per tutti se serve per pproc
      call calculate_normal(nvstart_2dstr,nvend_2dstr, nfstart_2dstr,nfend_2dstr, &
               xyzL(:,nvstart_2dstr:nvend_2dstr),  &
               vert_of_face(:,nfstart_2dstr:nfend_2dstr), &
               tri_nor(:,nfstart_2dstr:nfend_2dstr))

      call calculate_angle(theta(nestart_2dstr:neend_2dstr), & 
               nestart_2dstr,neend_2dstr, nfstart_2dstr,nfend_2dstr, &
               face_of_edge(:,nestart_2dstr:neend_2dstr),  &
               tri_nor(:,nfstart_2dstr:nfend_2dstr))

      call calculate_anglealp(aalpha(:,nfstart_2dstr:nfend_2dstr),  &
                  nvstart_2dstr,nvend_2dstr, nfstart_2dstr,nfend_2dstr, &
                  xyzL(:,nvstart_2dstr:nvend_2dstr),  &
                  vert_of_face(:,nfstart_2dstr:nfend_2dstr))

      ! Only surface area needs to be computed per particle
      do inp = inpstart_2dstr, inpend_2dstr
        call calculate_area(Surface(inp),vstart(inp),vend(inp), &
                 fstart(inp),fend(inp),xyzL(:,vstart(inp):vend(inp)), &
                 vert_of_face(:,fstart(inp):fend(inp)), &
                 sur(fstart(inp):fend(inp)))
      end do

      call nvtxEndRange

!     --------------------------------------------------------
      call nvtxStartRange("compute_s_forces", 12)
      fxyz = 0.d0

!     Computing internal force for particles in domain
      ! For parallelization, split edge, vertex, and face ranges across processes

      nverti_2dstr = nvend_2dstr-nvstart_2dstr+1
      nedges_2dstr = neend_2dstr-nestart_2dstr+1
      nfaces_2dstr = nfend_2dstr-nfstart_2dstr+1

      nv_per_rank = ceiling(real(nverti_2dstr)/numtasks)
      my_vs = myid * nv_per_rank + nvstart_2dstr
      my_ve = min((myid + 1) * nv_per_rank + nvstart_2dstr -1, nvend_2dstr)
      ne_per_rank = ceiling(real(nedges_2dstr)/numtasks)
      my_es = myid * ne_per_rank + nestart_2dstr
      my_ee = min((myid + 1) * ne_per_rank + nestart_2dstr -1, neend_2dstr)
      nf_per_rank = ceiling(real(nfaces_2dstr)/numtasks)
      my_fs = myid * nf_per_rank + nfstart_2dstr
      my_fe = min((myid + 1) * nf_per_rank + nfstart_2dstr -1, nfend_2dstr)

      ! Van Gelder model for elastic constant. rke is the 2d Young modulus
      ! Note: All ranks compute ke and kb for all edges
      NUPOIS = 0.333
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do i = nestart_2dstr,neend_2dstr
        inp = edge_to_part(i)
        f1 = face_of_edge(1,i)
        f2 = face_of_edge(2,i)
#ifdef HYPERELASTIC
!hyperelastic        
        if(f1 .ne. 0 .and. f2 .ne. 0) then
          ke(i) = rkc(inp)*AFung(i)*(sur(f1)+sur(f2))/dist(i)**2
        elseif (f1 .ne. 0 .and. f2 .eq. 0) then
          ke(i) = rkc(inp)*AFung(i)*sur(f1)/dist(i)**2
        elseif (f2 .ne. 0 .and. f1 .eq. 0) then
          ke(i)=rkc(inp)*AFung(i)*sur(f2)/dist(i)**2
        end if
#else
!linear                                                                                                                                                                   
        if(f1 .ne. 0 .and. f2 .ne. 0) then
          ke(i) = rke(inp) * (sur(f1)+sur(f2)) / dist(i)**2
        elseif (f1 .ne. 0 .and. f2 .eq. 0) then
	  ke(i) = rke(inp) * sur(f1) / dist(i)**2
	elseif (f2 .ne. 0 .and. f1 .eq. 0) then
	  ke(i) = rke(inp) * sur(f2) / dist(i)**2
        end if
#endif

        ke(i) = ke(i) * Surface0(inp) / Surface(inp)
        kb(i) = rke(inp) * thck(inp)**2 / &
         (12.*(1-NUPOIS**2)) * rkb(inp)*2 * sqrt(3.)  !check if o / sqrt(3)
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!     elastic potential AND viscous damping
#ifdef USE_CUDA
!$cuf kernel do (1)
#endif
      do i = my_vs, my_ve
        inp = vert_to_part(i)
        visc_coeff = cv1(inp)
        visc_coeff2 = cv2(inp)
        cstiff = (3.-1.)*StiffV(i)+1.
        do j = 1, n_edge_of_vert(i)
          iv = vert_of_vert(j,i)
          ie = edge_of_vert(j,i)
          d = dist(ie)
          d0 = dist0(ie)
          kei = ke(ie)

          epsG = (d-d0) / d
#ifdef HYPERELASTIC
          FungNL = exp(AFung(ie)*epsG**2) !hyperelastic model
#else
          FungNL = 1.0D0
#endif
          fxyz(1,i) = fxyz(1,i)-cstiff*kei*FungNL*epsG * &
                           (xyzL(1,i)-xyzL(1,iv))
          fxyz(2,i) = fxyz(2,i)-cstiff*kei*FungNL*epsG * &
                           (xyzL(2,i)-xyzL(2,iv))
          fxyz(3,i) = fxyz(3,i)-cstiff*kei*FungNL*epsG * &
                           (xyzL(3,i)-xyzL(3,iv))

          ! fxyz(1,i) = fxyz(1,i)-kei*(d-d0) / d * &
          !                  (xyzL(1,i)-xyzL(1,iv))
          ! fxyz(2,i) = fxyz(2,i)-kei*(d-d0) / d * &
          !                  (xyzL(2,i)-xyzL(2,iv))
          ! fxyz(3,i) = fxyz(3,i)-kei*(d-d0) / d * &
          !                  (xyzL(3,i)-xyzL(3,iv))

          Vij1 = (xyzvL(1,i)-xyzvL(1,iv))
          Vij2 = (xyzvL(2,i)-xyzvL(2,iv))
          Vij3 = (xyzvL(3,i)-xyzvL(3,iv))

          VV = sqrt(Vij1*Vij1 + Vij2*Vij2 + Vij3*Vij3)

!         F = F - \mu V   (Stokes damping - linear)
          if (VV .gt. 0.0) then
            fxyz(1,i) = fxyz(1,i) - visc_coeff*Vij1 / VV
            fxyz(2,i) = fxyz(2,i) - visc_coeff*Vij2 / VV
            fxyz(3,i) = fxyz(3,i) - visc_coeff*Vij3 / VV
          endif
        enddo

        fxyz(1,i) = fxyz(1,i) - visc_coeff2*xyzvL(1,i)
        fxyz(2,i) = fxyz(2,i) - visc_coeff2*xyzvL(2,i)
        fxyz(3,i) = fxyz(3,i) - visc_coeff2*xyzvL(3,i)

      enddo

!@cuf istat = cudaDeviceSynchronize !JDR TMP

!       bending forces
      ! NOTE: This routine is causes variation in GPU results at startup.
#ifdef USE_CUDA
      threads = dim3(128, 1, 1)
      blocks = dim3(ceiling(real(ne_per_rank)/128), 1, 1)
      call compute_fbxyz_kernel<<<blocks, threads>>>(my_es, my_ee, kb, &
                                                     face_of_edge, v1234, xyzL, theta, theta0, fxyz)
#else
      do i = my_es, my_ee

        f1 = face_of_edge(1,i)
        f2 = face_of_edge(2,i)

        if (f1 .ne. 0.and. f2 .ne. 0) then

          v1 = v1234(1,i)
          v2 = v1234(2,i)
          v3 = v1234(3,i)
          v4 = v1234(4,i)

          a21(1:3) = xyzL(1:3,v2) - xyzL(1:3,v1)
          a31(1:3) = xyzL(1:3,v3) - xyzL(1:3,v1)
          a32(1:3) = xyzL(1:3,v3) - xyzL(1:3,v2)
          a13(1:3) = xyzL(1:3,v1) - xyzL(1:3,v3)
          a34(1:3) = xyzL(1:3,v3) - xyzL(1:3,v4)
          a42(1:3) = xyzL(1:3,v4) - xyzL(1:3,v2)
          a23(1:3) = xyzL(1:3,v2) - xyzL(1:3,v3)
          a24(1:3) = xyzL(1:3,v2) - xyzL(1:3,v4)

          call cross(csi,a21,a31)
          call cross(zet,a34,a24)

          modcsi = sqrt(csi(1)**2+csi(2)**2+csi(3)**2)
          modzet = sqrt(zet(1)**2+zet(2)**2+zet(3)**2)

          if ((1. - cos(theta(i))**2) .le. 1.e-8) then
            betab = 0.0
          else
            betab = kb(i)*(sin(theta(i))*cos(theta0(i))- &
                    cos(theta(i))*sin(theta0(i)))/       &
                    (sqrt(1.-(cos(theta(i))**2)))
          endif

          b11 = -betab*cos(theta(i))/(modcsi**2)
          b12 = betab/(modcsi*modzet)
          b22 = -betab*cos(theta(i))/(modzet**2)

          call cross(tvec1,csi,a32)
          call cross(tvec2,zet,a32)
          fxyz(1:3,v1) = fxyz(1:3,v1)+b11*tvec1(1:3)+ &
                                        b12*tvec2(1:3)

          call cross(tvec1,csi,a13)
          call cross(tvec2,csi,a34)
          call cross(tvec3,zet,a13)
          call cross(tvec4,zet,a34)
          fxyz(1:3,v2) = fxyz(1:3,v2)+b11*tvec1(1:3)+ &
                            b12*(tvec2(1:3)+tvec3(1:3))+ &
                                        b22*tvec4(1:3)

          call cross(tvec1,csi,a21)
          call cross(tvec2,csi,a42)
          call cross(tvec3,zet,a21)
          call cross(tvec4,zet,a42)
          fxyz(1:3,v3) = fxyz(1:3,v3)+b11*tvec1(1:3)+ &
                            b12*(tvec2(1:3)+tvec3(1:3))+ &
                                        b22*tvec4(1:3)

          call cross(tvec1,csi,a23)
          call cross(tvec2,zet,a23)
          fxyz(1:3,v4) = fxyz(1:3,v4)+b12*tvec1(1:3)+ &
                                        b22*tvec2(1:3)

        endif

      enddo
#endif

!@cuf istat = cudaDeviceSynchronize !JDR TMP

!     Volume constraint, local, and total area potential
#ifdef USE_CUDA
      threads = dim3(128, 1, 1)
      blocks = dim3(ceiling(real(nf_per_rank)/128), 1, 1)
      call compute_favxyz_kernel<<<blocks, threads>>>(my_fs, my_fe, kat, kal, kv, kae,&
                                                      Surface, Surface0, Volume, Volume0, &
                                                      xyzL, vert_of_face, face_to_part, &
                                                      tri_nor, sur, sur0, aalpha, aalpha0, &
                                                      nfi, fxyz)
#else
      do i = my_fs, my_fe
        inp = face_to_part(i)
        v1 = vert_of_face(1,i)
        v2 = vert_of_face(2,i)
        v3 = vert_of_face(3,i)

        a32(1:3) = xyzL(1:3,v3)-xyzL(1:3,v2)
        a13(1:3) = xyzL(1:3,v1)-xyzL(1:3,v3)
        a21(1:3) = xyzL(1:3,v2)-xyzL(1:3,v1)


        csi(1:3) = tri_nor(1:3,i)

        tcv(1:3) = (xyzL(1:3,v1)+xyzL(1:3,v2)+ &
                    xyzL(1:3,v3))/3.


        alphat = -kat(inp)*(Surface(inp)-Surface0(inp))/Surface0(inp)
        alphal = -kal(inp)*(sur(i)-sur0(i))/sur0(i)

        call cross(tvec1,csi,a32)
        fxyz(1:3,v1) = fxyz(1:3,v1)+(alphat+alphal)*tvec1(1:3)
        call cross(tvec1,csi,a13)
        fxyz(1:3,v2) = fxyz(1:3,v2)+(alphat+alphal)*tvec1(1:3)
        call cross(tvec1,csi,a21)
        fxyz(1:3,v3) = fxyz(1:3,v3)+(alphat+alphal)*tvec1(1:3)

        ! angle potential
        alphae = kae(inp) * (aalpha(1,i) - aalpha0(1,i))
        fxyz(1:3,v2) = fxyz(1:3,v2) - alphae*a32(1:3)
        fxyz(1:3,v3) = fxyz(1:3,v3) + alphae*a32(1:3)

        alphae = kae(inp) * (aalpha(2,i) - aalpha0(2,i))
        fxyz(1:3,v3) = fxyz(1:3,v3) - alphae*a13(1:3)
        fxyz(1:3,v1) = fxyz(1:3,v1) + alphae*a13(1:3)

        alphae = kae(inp)*(aalpha(3,i) - aalpha0(3,i))
        fxyz(1:3,v1)=fxyz(1:3,v1) - alphae*a21(1:3)
        fxyz(1:3,v2)=fxyz(1:3,v2) + alphae*a21(1:3)
      enddo
#endif

!@cuf istat = cudaDeviceSynchronize !JDR TMP

      call nvtxEndRange

!     --------------------------------------------------------


      return
      end 
