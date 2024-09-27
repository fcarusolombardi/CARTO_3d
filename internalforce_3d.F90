module internalforce_3d_m
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
      subroutine compute_favxyz_3d_kernel(fstart, fend, kae_3d, &
        xyzL_3d, vert_of_face_3d, face_to_part_3d, aalpha_3d, aalpha0_3d, fxyz_3d)
!        use mls_param, only: StiffV
        implicit none

        integer, value, intent(in) :: fstart, fend
        integer, device, dimension(3,*), intent(in) :: vert_of_face_3d
        integer, device, dimension(*), intent(in) :: face_to_part_3d
        real(DP), device, dimension(*), intent(in) :: kae_3d
        real(DP), device, dimension(3,*), intent(in) :: xyzL_3d
        real(DP), device, dimension(3,*), intent(in) :: aalpha_3d, aalpha0_3d
        real(DP), device, dimension(3,*), intent(inout) :: fxyz_3d
        integer :: v1, v2, v3
        real(DP), device :: a32(3), a13(3), a21(3), csi(3), tcv(3)
        real(DP), device :: tvec1(3), tvec2(3), tvec3(3), val(3)
        real(DP) :: rv1, rv2, rv3
        real(DP) :: alphae1, alphae2, alphae3

        integer :: i, inp

        i = (blockIdx%x - 1) * blockDim%x + threadIdx%x + fstart - 1

        if (i > fend) return

        inp = face_to_part_3d(i)

        v1 = vert_of_face_3d(1,i)
        v2 = vert_of_face_3d(2,i)
        v3 = vert_of_face_3d(3,i)

        a32(1:3) = xyzL_3d(1:3,v3)-xyzL_3d(1:3,v2)
        a13(1:3) = xyzL_3d(1:3,v1)-xyzL_3d(1:3,v3)
        a21(1:3) = xyzL_3d(1:3,v2)-xyzL_3d(1:3,v1)

        alphae1 = kae_3d(inp) * (aalpha_3d(1,i) - aalpha0_3d(1,i))
        alphae2 = kae_3d(inp) * (aalpha_3d(2,i) - aalpha0_3d(2,i))
        alphae3 = kae_3d(inp) * (aalpha_3d(3,i) - aalpha0_3d(3,i))

        rv1 = alphae2*a13(1) - alphae3*a21(1)
        rv1 = atomicadd(fxyz_3d(1,v1),rv1)

        rv2 = alphae2*a13(2) - alphae3*a21(2)
        rv2 = atomicadd(fxyz_3d(2,v1),rv2)

        rv3 = alphae2*a13(3) - alphae3*a21(3)
        rv3 = atomicadd(fxyz_3d(3,v1),rv3)


        rv1 =- alphae1*a32(1) + alphae3*a21(1)
        rv1 = atomicadd(fxyz_3d(1,v2), rv1)

        rv2 =- alphae1*a32(2) + alphae3*a21(2)
        rv2 = atomicadd(fxyz_3d(2,v2), rv2)

        rv3 =- alphae1*a32(3) + alphae3*a21(3)
        rv3 = atomicadd(fxyz_3d(3,v2), rv3)


        rv1 = alphae1*a32(1) - alphae2*a13(1)
        rv1 = atomicadd(fxyz_3d(1,v3),rv1)

        rv2 = alphae1*a32(2) - alphae2*a13(2)
        rv2 = atomicadd(fxyz_3d(2,v3),rv2)

        rv3 = alphae1*a32(3) - alphae2*a13(3)
        rv3 = atomicadd(fxyz_3d(3,v3),rv3)
      end subroutine compute_favxyz_3d_kernel


      attributes(global) &
      subroutine compute_fcvxyz_3d_kernel(cstart, cend, kv_3d, &
        xyzL_3d, vert_of_cell_3d, cell_to_part_3d, cell_to_chamb_3d, fxyz_3d)
        use mls_param, only: StiffV_3d
        use param, only: nbeat
        implicit none

        integer, value, intent(in) :: cstart, cend
        integer, device, dimension(4,*), intent(in) :: vert_of_cell_3d
        integer, device, dimension(*), intent(in) :: cell_to_part_3d
        integer, device, dimension(*), intent(in) :: cell_to_chamb_3d
        real(DP), device, dimension(*), intent(in) :: kv_3d
        real(DP), device, dimension(3,*), intent(in) :: xyzL_3d
        real(DP), device, dimension(3,*), intent(inout) :: fxyz_3d
        integer :: v1, v2, v3, v4
        real(DP), device :: tcv(3),tvec1(3), tvec2(3), tvec3(3), tvec4(3)
        real(DP) :: rv1, rv2, rv3, rv4
        real(DP) :: smoothv1,smoothv2,smoothv3,smoothv4
        real(DP) :: betav

        integer :: i, inp,chamb

        i = (blockIdx%x - 1) * blockDim%x + threadIdx%x + cstart - 1

        if (i > cend) return

        inp   = cell_to_part_3d(i)
        
        if (inp .eq. 1) then 
        v1 = vert_of_cell_3d(1,i)
        v2 = vert_of_cell_3d(2,i)
        v3 = vert_of_cell_3d(3,i)
        v4 = vert_of_cell_3d(4,i)

        tvec1(1:3) = xyzL_3d(1:3,v1)
        tvec2(1:3) = xyzL_3d(1:3,v2)
        tvec3(1:3) = xyzL_3d(1:3,v3)
        tvec4(1:3) = xyzL_3d(1:3,v4)

        tcv(1:3) = (tvec1(1:3)+tvec2(1:3)+tvec3(1:3)+tvec4(1:3))/4.d0

        betav = -kv_3d(inp)
        betav = 0.D0
        
        smoothv1 = (1.d0)!-StiffV_3d(v1))
        smoothv2 = (1.d0)!-StiffV_3d(v2))
        smoothv3 = (1.d0)!-StiffV_3d(v3))
        smoothv4 = (1.d0)!-StiffV_3d(v4))

        rv1 = atomicadd(fxyz_3d(1,v1),smoothv1*betav*(tvec1(1)-tcv(1)))
        rv1 = atomicadd(fxyz_3d(2,v1),smoothv1*betav*(tvec1(2)-tcv(2)))
        rv1 = atomicadd(fxyz_3d(3,v1),smoothv1*betav*(tvec1(3)-tcv(3)))

        rv2 = atomicadd(fxyz_3d(1,v2),smoothv2*betav*(tvec2(1)-tcv(1)))
        rv2 = atomicadd(fxyz_3d(2,v2),smoothv2*betav*(tvec2(2)-tcv(2)))
        rv2 = atomicadd(fxyz_3d(3,v2),smoothv2*betav*(tvec2(3)-tcv(3)))

        rv3 = atomicadd(fxyz_3d(1,v3),smoothv3*betav*(tvec3(1)-tcv(1)))
        rv3 = atomicadd(fxyz_3d(2,v3),smoothv3*betav*(tvec3(2)-tcv(2)))
        rv3 = atomicadd(fxyz_3d(3,v3),smoothv3*betav*(tvec3(3)-tcv(3)))

        rv4 = atomicadd(fxyz_3d(1,v4),smoothv4*betav*(tvec4(1)-tcv(1)))
        rv4 = atomicadd(fxyz_3d(2,v4),smoothv4*betav*(tvec4(2)-tcv(2)))
        rv4 = atomicadd(fxyz_3d(3,v4),smoothv4*betav*(tvec4(3)-tcv(3)))

        ! rv1 = atomicadd(fxyz_3d(1,v1),-betav*xyzL_3d(1,v1))
        ! rv1 = atomicadd(fxyz_3d(2,v1),-betav*xyzL_3d(2,v1))
        ! rv1 = atomicadd(fxyz_3d(3,v1),-betav*xyzL_3d(3,v1))

        ! rv2 = atomicadd(fxyz_3d(1,v2),-betav*xyzL_3d(1,v2))
        ! rv2 = atomicadd(fxyz_3d(2,v2),-betav*xyzL_3d(2,v2))
        ! rv2 = atomicadd(fxyz_3d(3,v2),-betav*xyzL_3d(3,v2))

        ! rv3 = atomicadd(fxyz_3d(1,v3),-betav*xyzL_3d(1,v3))
        ! rv3 = atomicadd(fxyz_3d(2,v3),-betav*xyzL_3d(2,v3))
        ! rv3 = atomicadd(fxyz_3d(3,v3),-betav*xyzL_3d(3,v3))

        ! rv4 = atomicadd(fxyz_3d(1,v4),-betav*xyzL_3d(1,v4))
        ! rv4 = atomicadd(fxyz_3d(2,v4),-betav*xyzL_3d(2,v4))
        ! rv4 = atomicadd(fxyz_3d(3,v4),-betav*xyzL_3d(3,v4))

        endif

      end subroutine compute_fcvxyz_3d_kernel

#endif

end module internalforce_3d_m

      subroutine internalforce_3d(snv,env,xyzL_3d,xyzvL_3d)
      use mpih
      use param
      use mls_param
      use internalforce_3d_m
      use nvtx
!@cuf use cudafor
      implicit none

      integer :: inp,i,j,k,chamb
      integer :: imasZ, snv,env
      integer :: v1,v2,v3,v4,f1,f2,f3,iv,ie,c1
      
      real(DP),dimension(3,snv:env),intent(in) :: xyzL_3d,xyzvL_3d
      real(DP),dimension(3)::tvec1,tvec2,tvec3,tvec4,csi,zet,Vij
      real(DP),dimension(3)::a32,a13,a34,a21,a42,a23,a31,a24,csi1,zet1,tcv
      real(DP)::modcsi,modzet,betab,betav,alphat,alphal,b11,b12,b22,tdum,VV,alphae
 
      real(DP) :: angleHelix, angleSpiral, masZ,cosAngFE,astressnondim
      real(DP) :: xCC, yCC, zCC, surT
      real(DP) :: zMax
      real(DP) :: alphaES,betaES,gammaES,smoothf,distMEAN,CosePhi,SenoPhi
      real(DP) :: visc_coeff, visc_coeff2, epsG, FungNL, cstiff
      real(DP) :: d, d0, kei,volcells,facti
      real(DP) :: Vij1, Vij2, Vij3
      real(DP) :: x_inf_loc,y_inf_loc,z_inf_loc,size_inf
      real(DP) :: xedge,yedge,zedge,campo_inf,num,den

      integer :: ne_per_rank, my_es, my_ee
      integer :: nv_per_rank, my_vs, my_ve
      integer :: nf_per_rank, my_fs, my_fe
      integer :: nc_per_rank, my_cs, my_ce

      real(DP), dimension(3) :: xyzCC
#ifdef USE_CUDA
      attributes(managed) :: xyzCC
      attributes(managed) :: xyzL_3d
      attributes(managed) :: xyzvL_3d
      type(dim3) :: blocks, threads
#endif
!@cuf integer :: istat


      call nvtxStartRange("geo_update_3d", 110)

      ! Just calculate all geometry properties first
      call calculate_distance(dist_3d, 1, nvtot_3d, 1, netot_3d, &
               xyzL_3d, vert_of_edge_3d)

      call calculate_anglealp(aalpha_3d, 1, nvtot_3d, 1, nftot_3d, &
               xyzL_3d, vert_of_face_3d)

      call calculate_volume_cells(1,nvtot_3d,1,nctot_3d,xyzL_3d(:,1:nvtot_3d),  &
               vert_of_cell_3d(:,1:nctot_3d),vol_3d(1:nctot_3d))

      call nvtxEndRange

!     --------------------------------------------------------
      call nvtxStartRange("compute_s_forces", 120)
      fxyz_3d = 0.d0

!     Computing internal force for particles in domain
      ! For parallelization, split edge, vertex, and face ranges across processes
      ne_per_rank = ceiling(real(netot_3d)/numtasks)
      my_es = myid * ne_per_rank + 1
      my_ee = min((myid + 1) * ne_per_rank, netot_3d)
      nv_per_rank = ceiling(real(nvtot_3d)/numtasks)
      my_vs = myid * nv_per_rank + 1
      my_ve = min((myid + 1) * nv_per_rank, nvtot_3d)
      nf_per_rank = ceiling(real(nftot_3d)/numtasks)
      my_fs = myid * nf_per_rank + 1
      my_fe = min((myid + 1) * nf_per_rank, nftot_3d)
      nc_per_rank = ceiling(real(nctot_3d)/numtasks)
      my_cs = myid * nc_per_rank + 1
      my_ce = min((myid + 1) * nc_per_rank, nctot_3d)


      size_inf=xyz_inf(4)
      x_inf_loc=xyz_inf(1)
      y_inf_loc=xyz_inf(2)
      z_inf_loc=xyz_inf(3)           
      ! Van Gelder model for elastic constant. rke is the 2d Young modulus
      ! Note: All ranks compute ke and kb for all edges
#ifdef USE_CUDA
      !$cuf kernel do (1)
#endif
      do i = 1, netot_3d 
         inp = edge_to_part_3d(i)
         chamb = edge_to_chamb_3d(i)
         volcells = 0.D0
         do j=1, n_cell_of_edge_3d(i)
            c1 = cell_of_edge_3d(j,i)          
            volcells =  volcells + vol_3d(c1)
         enddo
#ifdef HYPERELASTIC
        !hyperelastic for each chamber
         ke_3d(i) = rkc_3d(chamb)*AFung_3d(i)*volcells/dist_3d(i)**2
#else
        !linear elastic for the whole body
         ke_3d(i) = rke_3d(inp)*volcells/dist_3d(i)**2
#endif    
         if (minf.GE.1) then   !for infarction
            v1=vert_of_edge_3d(1,i)
            v2=vert_of_edge_3d(2,i)
            xedge=0.5D0*(xyz0_3d(1,v1)+xyz0_3d(1,v2))
            yedge=0.5D0*(xyz0_3d(2,v1)+xyz0_3d(2,v2))
            zedge=0.5D0*(xyz0_3d(3,v1)+xyz0_3d(3,v2))
            campo_inf=1.0D0 + (2.0D0-1.0D0)*exp( -( ( (xedge-x_inf_loc)**2+(yedge-y_inf_loc)**2+(zedge-z_inf_loc)**2)**4/size_inf**8 ) )  !per raddoppiare stiffness          
            ke_3d(i) = ke_3d(i)*campo_inf
         endif
                  
         !for active tension
         cosAngFE   = edgefiber_cosangle_3d(i)
         astressnondim=astressEFedge_3d(i)
         f_act_3d(i) = cosAngFE*astressnondim*volcells/dist_3d(i)**2
      end do
!@cuf istat = cudaDeviceSynchronize !JDR TMP

!     elastic potential AND viscous damping
#ifdef USE_CUDA
!$cuf kernel do (1)
#endif
      do i = my_vs, my_ve  
         inp = vert_to_part_3d(i)
         chamb = vert_to_chamb_3d(i)
         visc_coeff  = cv1_3d(inp)         
         visc_coeff2 = cv2_3d(inp)
         cstiff = 1.d0!(3.-1.)*StiffV_3d(i)+1.

        do j = 1, n_edge_of_vert_3d(i)
          iv = vert_of_vert_3d(j,i)
          ie = edge_of_vert_3d(j,i)
          d = dist_3d(ie)
          d0 = dist0_3d(ie)
          kei = ke_3d(ie)
          facti= f_act_3d(ie)

          epsG = (d-d0) / d
#ifdef HYPERELASTIC          
          FungNL = exp(AFung_3d(ie)*epsG**2)
#else
          FungNL = 1.
#endif
          fxyz_3d(1,i) = fxyz_3d(1,i)-cstiff*kei*FungNL*epsG * &
                           (xyzL_3d(1,i)-xyzL_3d(1,iv))
          fxyz_3d(2,i) = fxyz_3d(2,i)-cstiff*kei*FungNL*epsG * &
                           (xyzL_3d(2,i)-xyzL_3d(2,iv))
          fxyz_3d(3,i) = fxyz_3d(3,i)-cstiff*kei*FungNL*epsG * &
                           (xyzL_3d(3,i)-xyzL_3d(3,iv))

          ! fxyz_3d(1,i) = fxyz_3d(1,i)-kei*(d-d0) / d * &
          !                  (xyzL_3d(1,i)-xyzL_3d(1,iv))
          ! fxyz_3d(2,i) = fxyz_3d(2,i)-kei*(d-d0) / d * &
          !                  (xyzL_3d(2,i)-xyzL_3d(2,iv))
          ! fxyz_3d(3,i) = fxyz_3d(3,i)-kei*(d-d0) / d * &
          !                  (xyzL_3d(3,i)-xyzL_3d(3,iv))

!active tension                                                                                                                                                                  
          fxyz_3d(1,i) = fxyz_3d(1,i)-facti*(xyzL_3d(1,i)-xyzL_3d(1,iv))
          fxyz_3d(2,i) = fxyz_3d(2,i)-facti*(xyzL_3d(2,i)-xyzL_3d(2,iv))
          fxyz_3d(3,i) = fxyz_3d(3,i)-facti*(xyzL_3d(3,i)-xyzL_3d(3,iv))

          
          Vij1 = (xyzvL_3d(1,i)-xyzvL_3d(1,iv))
          Vij2 = (xyzvL_3d(2,i)-xyzvL_3d(2,iv))
          Vij3 = (xyzvL_3d(3,i)-xyzvL_3d(3,iv))

          VV = sqrt(Vij1*Vij1 + Vij2*Vij2 + Vij3*Vij3)

!         F = F - \mu V   (Stokes damping - linear)
          if (VV .gt. 0.0) then
            fxyz_3d(1,i) = fxyz_3d(1,i) - visc_coeff*Vij1 / VV
            fxyz_3d(2,i) = fxyz_3d(2,i) - visc_coeff*Vij2 / VV
            fxyz_3d(3,i) = fxyz_3d(3,i) - visc_coeff*Vij3 / VV
         endif

        enddo

        fxyz_3d(1,i) = fxyz_3d(1,i) - visc_coeff2*xyzvL_3d(1,i)
        fxyz_3d(2,i) = fxyz_3d(2,i) - visc_coeff2*xyzvL_3d(2,i)
        fxyz_3d(3,i) = fxyz_3d(3,i) - visc_coeff2*xyzvL_3d(3,i)
      enddo

!@cuf istat = cudaDeviceSynchronize !JDR TMP

!     in plane angle 
#ifdef USE_CUDA
      threads = dim3(128, 1, 1)
      blocks = dim3(ceiling(real(nf_per_rank)/128), 1, 1)
      call compute_favxyz_3d_kernel<<<blocks, threads>>>(my_fs, my_fe, kae_3d,&
                                                      xyzL_3d, vert_of_face_3d, face_to_part_3d, &
                                                      aalpha_3d, aalpha0_3d, fxyz_3d)
#else
      do i = my_fs, my_fe
        inp = face_to_part_3d(i)
        v1 = vert_of_face_3d(1,i)
        v2 = vert_of_face_3d(2,i)
        v3 = vert_of_face_3d(3,i)

        a32(1:3) = xyzL_3d(1:3,v3)-xyzL_3d(1:3,v2)
        a13(1:3) = xyzL_3d(1:3,v1)-xyzL_3d(1:3,v3)
        a21(1:3) = xyzL_3d(1:3,v2)-xyzL_3d(1:3,v1)

        ! angle potential
        alphae = kae_3d(inp) * (aalpha_3d(1,i) - aalpha0_3d(1,i))
        fxyz_3d(1:3,v2) = fxyz_3d(1:3,v2) - alphae*a32(1:3)
        fxyz_3d(1:3,v3) = fxyz_3d(1:3,v3) + alphae*a32(1:3)

        alphae = kae_3d(inp) * (aalpha_3d(2,i) - aalpha0_3d(2,i))
        fxyz_3d(1:3,v3) = fxyz_3d(1:3,v3) - alphae*a13(1:3)
        fxyz_3d(1:3,v1) = fxyz_3d(1:3,v1) + alphae*a13(1:3)

        alphae = kae_3d(inp)*(aalpha_3d(3,i) - aalpha0_3d(3,i))
        fxyz_3d(1:3,v1)=fxyz_3d(1:3,v1) - alphae*a21(1:3)
        fxyz_3d(1:3,v2)=fxyz_3d(1:3,v2) + alphae*a21(1:3)
      enddo
#endif



!@cuf istat = cudaDeviceSynchronize !JDR TMP

!     cell shrinking  
#ifdef USE_CUDA
      threads = dim3(128, 1, 1)
      blocks = dim3(ceiling(real(nc_per_rank)/128), 1, 1)
      call compute_fcvxyz_3d_kernel<<<blocks, threads>>>(my_cs, my_ce, kv_3d,&
                                                      xyzL_3d, vert_of_cell_3d, cell_to_part_3d, &
                                                       cell_to_chamb_3d, fxyz_3d)
#else
      do i = my_cs, my_ce
        inp = cell_to_part_3d(i)
        if (inp .eq. 1) then
        v1 = vert_of_cell_3d(1,i)
        v2 = vert_of_cell_3d(2,i)
        v3 = vert_of_cell_3d(3,i)
        v4 = vert_of_cell_3d(4,i)

        tcv(1:3) = (xyzL_3d(1:3,v1)+xyzL_3d(1:3,v2)+xyzL_3d(1:3,v3)+xyzL_3d(1:3,v4))/4.d0
        
        betav = -kv_3d(inp)!*(pot_3d(v1)+pot_3d(v2)+pot_3d(v3)+pot_3d(v4))/4.d0
        betav = 0.D0

        ! fxyz_3d(1:3,v1) = fxyz_3d(1:3,v1)-(1.-StiffV_3d(v1))*betav*(xyzL_3d(1:3,v1)-tcv(1:3))
        ! fxyz_3d(1:3,v2) = fxyz_3d(1:3,v2)-(1.-StiffV_3d(v2))*betav*(xyzL_3d(1:3,v2)-tcv(1:3))
        ! fxyz_3d(1:3,v3) = fxyz_3d(1:3,v3)-(1.-StiffV_3d(v3))*betav*(xyzL_3d(1:3,v3)-tcv(1:3))
        ! fxyz_3d(1:3,v4) = fxyz_3d(1:3,v4)-(1.-StiffV_3d(v4))*betav*(xyzL_3d(1:3,v4)-tcv(1:3))
        ! fxyz_3d(1:3,v1) = fxyz_3d(1:3,v1)-betav*(xyzL_3d(1:3,v1)-tcv(1:3))
        ! fxyz_3d(1:3,v2) = fxyz_3d(1:3,v2)-betav*(xyzL_3d(1:3,v2)-tcv(1:3))
        ! fxyz_3d(1:3,v3) = fxyz_3d(1:3,v3)-betav*(xyzL_3d(1:3,v3)-tcv(1:3))
        ! fxyz_3d(1:3,v4) = fxyz_3d(1:3,v4)-betav*(xyzL_3d(1:3,v4)-tcv(1:3))
        fxyz_3d(1:3,v1) = fxyz_3d(1:3,v1)+betav*(xyzL_3d(1:3,v1)-tcv(1:3))
        fxyz_3d(1:3,v2) = fxyz_3d(1:3,v2)+betav*(xyzL_3d(1:3,v2)-tcv(1:3))
        fxyz_3d(1:3,v3) = fxyz_3d(1:3,v3)+betav*(xyzL_3d(1:3,v3)-tcv(1:3))
        fxyz_3d(1:3,v4) = fxyz_3d(1:3,v4)+betav*(xyzL_3d(1:3,v4)-tcv(1:3))
        endif

      enddo
#endif


!@cuf istat = cudaDeviceSynchronize !JDR TMP

      call nvtxEndRange

!     --------------------------------------------------------


      return
      end 
