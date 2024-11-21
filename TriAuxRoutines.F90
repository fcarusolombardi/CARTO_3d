!----------------------------------------------------------------------
!     Auxillary routines for triangulated geometry
!----------------------------------------------------------------------
        subroutine calculate_normal(snv,env,snf,enf,xyz,vert_of_face, tri_nor)
        use constants
!@cuf   use cudafor

        implicit none
        integer::snv,env,snf,enf,v1,v2,v3,i
        integer,dimension(3,snf:enf),intent(in) :: vert_of_face
        real(DP),dimension(3,snv:env),intent(in) :: xyz
        real(DP),dimension(3,snf:enf),intent(out) :: tri_nor
        real(DP) :: ve11, ve12, ve13
        real(DP) :: ve21, ve22, ve23
        real(DP) :: fn1, fn2, fn3, fnmag
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: vert_of_face, xyz, tri_nor
#endif

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i = snf,enf
          v1 = vert_of_face(1,i)
          v2 = vert_of_face(2,i)
          v3 = vert_of_face(3,i)

          ve11 = xyz(1,v2)-xyz(1,v1)
          ve12 = xyz(2,v2)-xyz(2,v1)
          ve13 = xyz(3,v2)-xyz(3,v1)

          ve21 = xyz(1,v3)-xyz(1,v1)
          ve22 = xyz(2,v3)-xyz(2,v1)
          ve23 = xyz(3,v3)-xyz(3,v1)

          fn1 = ve12*ve23 - ve13*ve22
          fn2 = ve13*ve21 - ve11*ve23
          fn3 = ve11*ve22 - ve12*ve21

          fnmag = sqrt(fn1*fn1 + fn2*fn2 + fn3*fn3)

          tri_nor(1,i) = fn1 / fnmag
          tri_nor(2,i) = fn2 / fnmag
          tri_nor(3,i) = fn3 / fnmag

        enddo

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_normal
!------------------------------------------------------
        subroutine calculate_volume(Volume,snv,env,snf,enf,xyz,vert_of_face)
        use constants
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snf,enf,v1,v2,v3,i
        integer, dimension (3,snf:enf), intent(in) :: vert_of_face
        real(DP), dimension (3,snv:env), intent(in) ::xyz
        real(DP), intent(out) :: Volume
        real(DP) :: d12,d23,d31,sp
        real(DP) :: x1, x2, x3
        real(DP) :: y1, y2, y3
        real(DP) :: z1, z2, z3
        real(DP) :: vl
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: vert_of_face, xyz!, vl
        attributes(device) :: Volume
#endif

        Volume=0.0
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snf,enf
           v1=vert_of_face(1,i)
           v2=vert_of_face(2,i)
           v3=vert_of_face(3,i)

           x1 = xyz(1,v1)
           x2 = xyz(1,v2)
           x3 = xyz(1,v3)

           y1 = xyz(2,v1)
           y2 = xyz(2,v2)
           y3 = xyz(2,v3)

           z1 = xyz(3,v1)
           z2 = xyz(3,v2)
           z3 = xyz(3,v3)

           d12 = x1*(y2*z3-z2*y3)
           d23 = x2*(y3*z1-z3*y1)
           d31 = x3*(y1*z2-z1*y2)

           vl = (d12+d23+d31)/6.d0

           Volume = Volume + vl
        enddo

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_volume
!------------------------------------------------------
!------------------------------------------------------
        subroutine calculate_volumeCPU(Volume,snv,env,snf,enf,xyz,vert_of_face)
        use constants

        implicit none
        integer :: snv,env,snf,enf,v1,v2,v3,i
        integer, dimension (3,snf:enf) :: vert_of_face
        real(DP), dimension (3,snv:env) ::xyz
        real(DP), dimension (3,3) ::dxyz
        real(DP), dimension (3) ::xyz_tr
        real(DP) :: Volume

        Volume=0.0
        do i=snf,enf
           v1=vert_of_face(1,i)
           v2=vert_of_face(2,i)
           v3=vert_of_face(3,i)

           dxyz(1,1)=xyz(1,v1);dxyz(2,1)=xyz(2,v1);dxyz(3,1)=xyz(3,v1)
           dxyz(1,2)=xyz(1,v2);dxyz(2,2)=xyz(2,v2);dxyz(3,2)=xyz(3,v2)
           dxyz(1,3)=xyz(1,v3);dxyz(2,3)=xyz(2,v3);dxyz(3,3)=xyz(3,v3)

           Volume = Volume + (dxyz(1,1)*(dxyz(2,2)*dxyz(3,3)-&  
                          dxyz(3,2)*dxyz(2,3)) +         & 
                         dxyz(1,2)*(dxyz(2,3)*dxyz(3,1)- &
                          dxyz(3,3)*dxyz(2,1)) +         &
                         dxyz(1,3)*(dxyz(2,1)*dxyz(3,2)- &
                          dxyz(3,1)*dxyz(2,2)))
        enddo
        Volume=Volume/6.

        return
        end subroutine calculate_volumeCPU
!------------------------------------------------------
        subroutine calculate_volume_chamb(Volume_chamb,snv,env,snf,enf,xyz,vert_of_face,face_to_chamb)
          use constants
          use mls_param, only:Volume_chambShift
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snf,enf,v1,v2,v3,i,chamb
        integer, dimension (3,snf:enf), intent(in) :: vert_of_face
        integer, dimension (snf:enf), intent(in) :: face_to_chamb
        real(DP), dimension (3,snv:env), intent(in) ::xyz
        real(DP), dimension (1:8), intent(out) :: Volume_chamb
!        real(DP), dimension (1) :: sp_chamb
        real(DP) :: sp_chamb
        real(DP), dimension(1) :: sp_chambvet
        real(DP) :: d12,d23,d31
        real(DP) :: x1, x2, x3
        real(DP) :: y1, y2, y3
        real(DP) :: z1, z2, z3
        real(DP) :: volLV,volLA,volRV,volRA,volAO,volVP,volAP,volVC
        real(DP) :: vl
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: face_to_chamb,vert_of_face, xyz,Volume_chamb
!        attributes(device) :: sp_chamb, sp_chambvet !Volume_chamb
#endif

! %1 LV, 2 LA, 3 RV, 4 RA, 5 AO, 6 VP, 7 AP, 8 VC 

        volLV = 0.0
        volLA = 0.0
        volRV = 0.0
        volRA = 0.0
        volAO = 0.0
        volVP = 0.0
        volAP = 0.0
        volVC = 0.0
        
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snf,enf
           chamb = face_to_chamb(i)

           v1=vert_of_face(1,i)
           v2=vert_of_face(2,i)
           v3=vert_of_face(3,i)

           x1 = xyz(1,v1)
           x2 = xyz(1,v2)
           x3 = xyz(1,v3)

           y1 = xyz(2,v1)
           y2 = xyz(2,v2)
           y3 = xyz(2,v3)

           z1 = xyz(3,v1)
           z2 = xyz(3,v2)
           z3 = xyz(3,v3)

           d12 = x1*(y2*z3-z2*y3)
           d23 = x2*(y3*z1-z3*y1)
           d31 = x3*(y1*z2-z1*y2)

           vl = (d12+d23+d31)/6.d0

           select case (chamb)
           case (1)
              volLV = volLV + vl
           case (2)
              volLA = volLA + vl
           case (3)
              volRV = volRV + vl
           case (4)
              volRA = volRA + vl
           case (5)
              volAO = volAO + vl
           case (6)
              volVP = volVP + vl
           case (7)
              volAP = volAP + vl
           case (8)
              volVC = volVC + vl
           end select
! %1 LV, 2 LA, 3 RV, 4 RA, 5 AO, 6 VP, 7 AP, 8 VC 
        enddo
        Volume_chamb(1) = volLV+Volume_chambShift(1)
        Volume_chamb(2) = volLA+Volume_chambShift(2)
        Volume_chamb(3) = volRV+Volume_chambShift(3)
        Volume_chamb(4) = volRA+Volume_chambShift(4)
        Volume_chamb(5) = volAO+Volume_chambShift(5)
        Volume_chamb(6) = volVP+Volume_chambShift(6)
        Volume_chamb(7) = volAP+Volume_chambShift(7)
        Volume_chamb(8) = volVC+Volume_chambShift(8)
        

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_volume_chamb
!------------------------------------------------------
        subroutine calculate_area_chamb(Surface_chamb,snv,env,snf,enf,xyz,vert_of_face,face_to_chamb)
        use constants
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snf,enf,v1,v2,v3,i,chamb
        integer, dimension (3,snf:enf), intent(in) :: vert_of_face
        integer, dimension (snf:enf), intent(in) :: face_to_chamb
        real(DP), dimension (3,snv:env), intent(in) ::xyz
        real(DP), dimension (1:8), intent(out) :: Surface_chamb
        real(DP) :: surLV,surLA,surRV,surRA,surAO,surVP,surAP,surVC
        real(DP) :: d12,d23,d31,sp
        real(DP) :: x1, x2, x3
        real(DP) :: y1, y2, y3
        real(DP) :: z1, z2, z3
        real(DP) :: s
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: vert_of_face, xyz, face_to_chamb,Surface_chamb
!        attributes(device) :: Surface
#endif

 ! %1 LV, 2 LA, 3 RV, 4 RA, 5 AO, 6 VP, 7 AP, 8 VC        
        surLV = 0.0
        surLA = 0.0
        surRV = 0.0
        surRA = 0.0
        surAO = 0.0
        surVP = 0.0
        surAP = 0.0
        surVC = 0.0
! #ifdef USE_CUDA
!         !$cuf kernel do (1)
! #endif
        do i=snf,enf
          chamb = face_to_chamb(i)

          v1=vert_of_face(1,i)
          v2=vert_of_face(2,i)
          v3=vert_of_face(3,i)

          x1 = xyz(1,v1)
          x2 = xyz(1,v2)
          x3 = xyz(1,v3)

          y1 = xyz(2,v1)
          y2 = xyz(2,v2)
          y3 = xyz(2,v3)

          z1 = xyz(3,v1)
          z2 = xyz(3,v2)
          z3 = xyz(3,v3)

          d12 = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))
          d23 = sqrt( (x2-x3)*(x2-x3) + (y2-y3)*(y2-y3) + (z2-z3)*(z2-z3))
          d31 = sqrt( (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) + (z3-z1)*(z3-z1))

          sp = 0.5d0 * (d12+d23+d31)
          s = sqrt(sp*(sp-d12)*(sp-d23)*(sp-d31))

           select case (chamb)
           case (1)
              surLV = surLV + s
           case (2)
              surLA = surLA + s
           case (3)
              surRV = surRV + s
           case (4)
              surRA = surRA + s
           case (5)
              surAO = surAO + s
           case (6)
              surVP = surVP + s
           case (7)
              surAP = surAP + s
           case (8)
              surVC = surVC + s
           end select
        enddo
        Surface_chamb(1) = surLV
        Surface_chamb(2) = surLA
        Surface_chamb(3) = surRV
        Surface_chamb(4) = surRA
        Surface_chamb(5) = surAO
        Surface_chamb(6) = surVP
        Surface_chamb(7) = surAP
        Surface_chamb(8) = surVC

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_area_chamb
!------------------------------------------------------
        subroutine calculate_area (Surface,snv,env,snf,enf,xyz,vert_of_face,sur)
        use constants
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snf,enf,v1,v2,v3,i
        integer, dimension (3,snf:enf), intent(in) :: vert_of_face
        real(DP), dimension (3,snv:env), intent(in) ::xyz
        real(DP), dimension (snf:enf), intent(out) :: sur
        real(DP), intent(out) :: Surface
        real(DP) :: d12,d23,d31,sp
        real(DP) :: x1, x2, x3
        real(DP) :: y1, y2, y3
        real(DP) :: z1, z2, z3
        real(DP) :: s
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: vert_of_face, xyz, sur
        attributes(device) :: Surface
#endif

        Surface=0.0
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snf,enf
          v1=vert_of_face(1,i)
          v2=vert_of_face(2,i)
          v3=vert_of_face(3,i)

          x1 = xyz(1,v1)
          x2 = xyz(1,v2)
          x3 = xyz(1,v3)

          y1 = xyz(2,v1)
          y2 = xyz(2,v2)
          y3 = xyz(2,v3)

          z1 = xyz(3,v1)
          z2 = xyz(3,v2)
          z3 = xyz(3,v3)

          d12 = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))
          d23 = sqrt( (x2-x3)*(x2-x3) + (y2-y3)*(y2-y3) + (z2-z3)*(z2-z3))
          d31 = sqrt( (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) + (z3-z1)*(z3-z1))

          sp = 0.5d0 * (d12+d23+d31)
          s = sqrt(sp*(sp-d12)*(sp-d23)*(sp-d31))
          sur(i) = s
          Surface = Surface + s
        enddo

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_area
!------------------------------------------------------
        subroutine calculate_anglealp(aalpha,snv,env,snf,enf,xyz,vert_of_face)
        use constants
!@cuf   use cudafor
        
        implicit none
        integer :: snv,env,snf,enf,i,v1,v2,v3
        integer,dimension(3,snf:enf), intent(in) :: vert_of_face
        real(DP),dimension(3,snv:env), intent(in) :: xyz
        real(DP),dimension(3,snf:enf), intent(out) :: aalpha
        real(DP) :: a32_1, a32_2, a32_3
        real(DP) :: a13_1, a13_2, a13_3
        real(DP) :: a21_1, a21_2, a21_3
        real(DP) :: a32n, a13n, a21n
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: vert_of_face, xyz, aalpha
#endif

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snf,enf
          v1 = vert_of_face(1,i)
          v2 = vert_of_face(2,i)
          v3 = vert_of_face(3,i)

          a32_1 = xyz(1,v3) - xyz(1,v2)
          a32_2 = xyz(2,v3) - xyz(2,v2)
          a32_3 = xyz(3,v3) - xyz(3,v2)
          a32n = sqrt(a32_1 * a32_1 + a32_2 * a32_2 + a32_3 * a32_3)

          a13_1 = xyz(1,v1) - xyz(1,v3)
          a13_2 = xyz(2,v1) - xyz(2,v3)
          a13_3 = xyz(3,v1) - xyz(3,v3)
          a13n = sqrt(a13_1 * a13_1 + a13_2 * a13_2 + a13_3 * a13_3)

          a21_1 = xyz(1,v2) - xyz(1,v1)
          a21_2 = xyz(2,v2) - xyz(2,v1)
          a21_3 = xyz(3,v2) - xyz(3,v1)
          a21n = sqrt(a21_1 * a21_1 + a21_2 * a21_2 + a21_3 * a21_3)

          aalpha(1,i) = acos((a13_1*a21_1 + a13_2*a21_2 + a13_3*a21_3) / (a13n * a21n))
          aalpha(2,i) = acos((a21_1*a32_1 + a21_2*a32_2 + a21_3*a32_3) / (a21n * a32n))
          aalpha(3,i) = acos((a32_1*a13_1 + a32_2*a13_2 + a32_3*a13_3) / (a32n * a13n))
        end do

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        end subroutine calculate_anglealp
!------------------------------------------------------
        ! NOTE: This routine is causes variation in GPU results at startup.
        subroutine calculate_angle(theta,sne,ene,snf,enf,face_of_edge,face_normal)
        use constants
!@cuf   use cudafor

        implicit none
        integer :: sne,ene,snf,enf,f1,f2,i
        integer, dimension (2,sne:ene), intent(in) :: face_of_edge
        real(DP), dimension (3,snf:enf), intent(in) :: face_normal
        real(DP), dimension (sne:ene), intent(out) :: theta
        real(DP) :: n11, n12, n13, n21, n22, n23
        real(DP) :: pvx,pvy,pvz,tmp,th,tmp1,tmp2
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: face_of_edge, face_normal, theta
#endif

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=sne,ene
          f1 = face_of_edge(1,i)
          f2 = face_of_edge(2,i)

          if (f1 .ne. 0 .and. f2 .ne. 0) then
            n11 = face_normal(1,f1)
            n12 = face_normal(2,f1)
            n13 = face_normal(3,f1)

            n21 = face_normal(1,f2)
            n22 = face_normal(2,f2)
            n23 = face_normal(3,f2)

            pvx = n12*n23 - n13*n22
            pvy = n13*n21 - n11*n23
            pvz = n11*n22 - n12*n21

            tmp1 = sqrt(pvx*pvx + pvy*pvy + pvz*pvz)
            tmp2 = (n11*n21 + n12*n22 + n13*n23)
            th = atan2(tmp1,tmp2)
          else
            th = 0.d0
          endif

          theta(i) = th
        enddo

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_angle
!------------------------------------------------------
        subroutine calculate_ginterp_face_cells(snv,env,snf,enf,snc,enc,xyz,vert_of_face,cell_of_face,cell_bar,versCFface,distCFface,g1interpface)
        use constants
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snf,enf,snc,enc,v1,v2,v3,i
        integer :: ce1,ce2
        integer, dimension (3,snf:enf), intent(in) :: vert_of_face
        integer, dimension (2,snf:enf), intent(in) :: cell_of_face
        real(DP), dimension (3,snv:env), intent(in) ::xyz
        real(DP), dimension (3,snc:enc), intent(in) ::cell_bar
        real(DP), dimension (3,snf:enf), intent(out) :: versCFface
        real(DP), dimension (snf:enf), intent(out) :: distCFface
        real(DP), dimension (snf:enf), intent(out) :: g1interpface
        real(DP) :: xBface,yBface,zBface
        real(DP) :: xM_ce1,yM_ce1,zM_ce1,xM_ce2,yM_ce2,zM_ce2
        real(DP) :: distce1,distce2,normve
        real(DP) :: xvettore,yvettore,zvettore
        real(DP) :: x1, x2, x3
        real(DP) :: y1, y2, y3
        real(DP) :: z1, z2, z3
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: cell_bar,vert_of_face,cell_of_face,xyz, versCFface, distCFface, g1interpface
#endif

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snf,enf
          v1=vert_of_face(1,i)
          v2=vert_of_face(2,i)
          v3=vert_of_face(3,i)

          x1 = xyz(1,v1)
          x2 = xyz(1,v2)
          x3 = xyz(1,v3)

          y1 = xyz(2,v1)
          y2 = xyz(2,v2)
          y3 = xyz(2,v3)

          z1 = xyz(3,v1)
          z2 = xyz(3,v2)
          z3 = xyz(3,v3)

          xBface = (x1+x2+x3)/3.D0
          yBface = (y1+y2+y3)/3.D0
          zBface = (z1+z2+z3)/3.D0
          
          ce1=cell_of_face(1,i)
          ce2=cell_of_face(2,i)
             
          if ((ce1.NE.0).AND.(ce2.NE.0)) then
             xM_ce1 = cell_bar(1,ce1)
             yM_ce1 = cell_bar(2,ce1)
             zM_ce1 = cell_bar(3,ce1)
             xM_ce2 = cell_bar(1,ce2)
             yM_ce2 = cell_bar(2,ce2)
             zM_ce2 = cell_bar(3,ce2)

             distce2 = sqrt( (xBface-xM_ce2)**2+(yBface-yM_ce2)**2+(zBface-zM_ce2)**2 )
             distce1 = sqrt( (xM_ce2-xM_ce1)**2+(yM_ce2-yM_ce1)**2+(zM_ce2-zM_ce1)**2 )

             g1interpface(i)=distce2/distce1
             xvettore=xM_ce2-xM_ce1    
             yvettore=yM_ce2-yM_ce1    
             zvettore=zM_ce2-zM_ce1    
          elseif ((ce1.NE.0).AND.(ce2.EQ.0)) then
             xM_ce1 = cell_bar(1,ce1)
             yM_ce1 = cell_bar(2,ce1)
             zM_ce1 = cell_bar(3,ce1)
             g1interpface(i)=1
             xvettore=xBface-xM_ce1   
             yvettore=yBface-yM_ce1    
             zvettore=zBface-zM_ce1    
          elseif ((ce1.EQ.0).AND.(ce2.NE.0)) then
             xM_ce2 = cell_bar(1,ce2)
             yM_ce2 = cell_bar(2,ce2)
             zM_ce2 = cell_bar(3,ce2)
             g1interpface(i)=0;
             xvettore=xM_ce2-xBface;      
             yvettore=yM_ce2-yBface;      
             zvettore=zM_ce2-zBface;      
          endif
          normve = sqrt(xvettore**2+yvettore**2+zvettore**2)
          versCFface(1,i)=xvettore/normve;
          versCFface(2,i)=yvettore/normve;
          versCFface(3,i)=zvettore/normve;
          distCFface(i)=normve;

        enddo

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_ginterp_face_cells
!------------------------------------------------------
        subroutine calculate_distance (dist,snv,env,sne,ene,xyz,vert_of_edge)
        use constants
!@cuf   use cudafor

        implicit none
        integer :: snv,env,sne,ene,v1,v2,i
        integer, dimension (2,sne:ene), intent(in) :: vert_of_edge
        real(DP), dimension (3,snv:env), intent(in) ::xyz
        real(DP), dimension (sne:ene), intent(out) :: dist
        real(DP) :: x1, x2, y1, y2, z1, z2
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: vert_of_edge, xyz, dist
#endif

#ifdef USE_CUDA
        !$cuf kernel do(1)
#endif
        do i=sne,ene
          v1 = vert_of_edge(1,i)
          v2 = vert_of_edge(2,i)

          x1 = xyz(1,v1)
          y1 = xyz(2,v1)
          z1 = xyz(3,v1)

          x2 = xyz(1,v2)
          y2 = xyz(2,v2)
          z2 = xyz(3,v2)

          dist(i) = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))
        enddo

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_distance
!------------------------------------------------------
        subroutine find_quartet (v1234,snv,env,sne,ene,snf,enf,xyz,face_of_edge, &
                             edge_of_face,vert_of_edge,face_normal)
        use constants
        !@cuf use cudafor
     
        implicit none

        !@cuf integer :: istat
        integer :: snv,env,sne,ene,snf,enf
        integer :: i,f1,f2,e1,e2,e3,v,v1,v2,v3,v4,vtmp
        integer, dimension (2,sne:ene) :: face_of_edge,vert_of_edge
        integer, dimension (3,snf:enf) :: edge_of_face
        real(DP), dimension (3,snf:enf) :: face_normal
        real(DP), dimension (3,snv:env) ::xyz
        integer, dimension (4,sne:ene) :: v1234
        real(DP), dimension (3) :: csi,zet,csi1,zet1,a21,a31,a24,a34
        real(DP) :: res
#ifdef USE_CUDA
        attributes(managed) :: face_of_edge,vert_of_edge
        attributes(managed) :: edge_of_face
        attributes(managed) :: face_normal
        attributes(managed) :: xyz
        attributes(managed) :: v1234
        attributes(managed) :: csi,zet,csi1,zet1,a21,a31,a24,a34
#endif
        do i=sne,ene
           f1=face_of_edge(1,i)
           f2=face_of_edge(2,i)
           if (f1.ne.0.and.f2.ne.0) then
              ! get opposite vertex of edge
              e1=edge_of_face(1,f1)
              e2=edge_of_face(2,f1)
              e3=edge_of_face(3,f1)
              if (e1.eq.i) then
                 v=vert_of_edge(1,e2)
              if(v.ne.vert_of_edge(1,i).and.v.ne.vert_of_edge(2,i))then
                    v1=v
                 else
                    v1=vert_of_edge(2,e2)
              endif
              endif
              if (e2.eq.i) then
                 v=vert_of_edge(1,e1)
              if(v.ne.vert_of_edge(1,i).and.v.ne.vert_of_edge(2,i))then
                    v1=v
                 else
                    v1=vert_of_edge(2,e1)
              endif
              endif
              if (e3.eq.i) then
                 v=vert_of_edge(1,e2)
              if(v.ne.vert_of_edge(1,i).and.v.ne.vert_of_edge(2,i))then
                    v1=v
                 else
                    v1=vert_of_edge(2,e2)
              endif
              endif
              ! get opposite vertex of edge
              e1=edge_of_face(1,f2)
              e2=edge_of_face(2,f2)
              e3=edge_of_face(3,f2)
              if (e1.eq.i) then
                 v=vert_of_edge(1,e2)
              if(v.ne.vert_of_edge(1,i).and.v.ne.vert_of_edge(2,i))then
                    v4=v
                 else
                    v4=vert_of_edge(2,e2)
              endif
              endif
              if (e2.eq.i) then
                 v=vert_of_edge(1,e1)
              if(v.ne.vert_of_edge(1,i).and.v.ne.vert_of_edge(2,i))then
                    v4=v
                 else
                    v4=vert_of_edge(2,e1)
              endif
              endif
              if (e3.eq.i) then
                 v=vert_of_edge(1,e2)
              if(v.ne.vert_of_edge(1,i).and.v.ne.vert_of_edge(2,i))then
                    v4=v
                 else
                 v4=vert_of_edge(2,e2)
              endif
              endif
              ! Other two vertices, on edge
              v2=vert_of_edge(1,i)
              v3=vert_of_edge(2,i)

              csi(1:3)=face_normal(1:3,f1)
              zet(1:3)=face_normal(1:3,f2)

              a21(1:3)=xyz(1:3,v2)-xyz(1:3,v1)
              a31(1:3)=xyz(1:3,v3)-xyz(1:3,v1)
              a34(1:3)=xyz(1:3,v3)-xyz(1:3,v4)
              a24(1:3)=xyz(1:3,v2)-xyz(1:3,v4)

              call cross(csi1,a21,a31)
              call cross(zet1,a34,a24)

              call dot(res,csi,csi1)
              if (res.lt.0.0) then
                 vtmp=v2
                 v2=v3
                 v3=vtmp
              endif

              v1234(1,i)=v1
              v1234(2,i)=v2
              v1234(3,i)=v3
              v1234(4,i)=v4
           endif
        enddo
        !@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine find_quartet     
!     ---------------------------------------------------------------------
        subroutine convert_geo(snv,env,sne,ene,snf,enf,xyz,xyzv,xyza, &
                  vert_of_face,tri_ver,tri_vel,tri_bar,vel_tri,acc_tri)
        use constants
!@cuf   use cudafor

        implicit none
        integer :: snv,env,sne,ene,snf,enf,v1,v2,v3,i
        real(DP), dimension (3,snv:env), intent(in) :: xyz,xyzv,xyza
        integer, dimension (3,snf:enf), intent(in) :: vert_of_face
        real(DP), dimension(9,snf:enf), intent(out) :: tri_ver,tri_vel
        real(DP), dimension(3,snf:enf), intent(out) :: tri_bar,vel_tri,acc_tri
        real(DP) :: v11, v12, v13, v21, v22, v23, v31, v32, v33
        real(DP) :: vv11, vv12, vv13, vv21, vv22, vv23, vv31, vv32, vv33
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: xyz,xyzv,xyza
        attributes(managed) :: vert_of_face
        attributes(managed) :: tri_ver,tri_vel,tri_bar,vel_tri,acc_tri
#endif
!        write(*,*)  "DENTRO CONVERT"
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snf,enf
          v1=vert_of_face(1,i)
          v2=vert_of_face(2,i)
          v3=vert_of_face(3,i)

          v11 = xyz(1,v1)
          v12 = xyz(2,v1)
          v13 = xyz(3,v1)

          v21 = xyz(1,v2)
          v22 = xyz(2,v2)
          v23 = xyz(3,v2)

          v31 = xyz(1,v3)
          v32 = xyz(2,v3)
          v33 = xyz(3,v3)

          tri_ver(1,i) = v11
          tri_ver(2,i) = v12
          tri_ver(3,i) = v13
          tri_ver(4,i) = v21
          tri_ver(5,i) = v22
          tri_ver(6,i) = v23
          tri_ver(7,i) = v31
          tri_ver(8,i) = v32
          tri_ver(9,i) = v33

!        write(*,*)  "DENTRO CONVERT1"
!FV
          vv11 = xyzv(1,v1)
          vv12 = xyzv(2,v1)
          vv13 = xyzv(3,v1)

          vv21 = xyzv(1,v2)
          vv22 = xyzv(2,v2)
          vv23 = xyzv(3,v2)

          vv31 = xyzv(1,v3)
          vv32 = xyzv(2,v3)
          vv33 = xyzv(3,v3)

          tri_vel(1,i) = vv11
          tri_vel(2,i) = vv12
          tri_vel(3,i) = vv13
          tri_vel(4,i) = vv21
          tri_vel(5,i) = vv22
          tri_vel(6,i) = vv23
          tri_vel(7,i) = vv31
          tri_vel(8,i) = vv32
          tri_vel(9,i) = vv33


          ! Find triangles' baricentre
          tri_bar(1,i) = (v11+v21+v31) / 3.d0
          tri_bar(2,i) = (v12+v22+v32) / 3.d0
          tri_bar(3,i) = (v13+v23+v33) / 3.d0

          ! Find triangle's velocity and acceleration 
          vel_tri(1,i)=(xyzv(1,v1)+xyzv(1,v2)+xyzv(1,v3)) / 3.d0
          vel_tri(2,i)=(xyzv(2,v1)+xyzv(2,v2)+xyzv(2,v3)) / 3.d0
          vel_tri(3,i)=(xyzv(3,v1)+xyzv(3,v2)+xyzv(3,v3)) / 3.d0

          acc_tri(1,i)=(xyza(1,v1)+xyza(1,v2)+xyza(1,v3)) / 3.d0
          acc_tri(2,i)=(xyza(2,v1)+xyza(2,v2)+xyza(2,v3)) / 3.d0
          acc_tri(3,i)=(xyza(3,v1)+xyza(3,v2)+xyza(3,v3)) / 3.d0
        enddo

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        end subroutine convert_geo
!     ----------------------------------------------------------------
        subroutine dot(xy,x,y)
        use constants

          implicit none
          real(DP) x(3),y(3),xy

          xy = x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
          return
        end subroutine dot
!     ----------------------------------------------------------------
        subroutine  sub(xmy,x,y)
        use constants

          implicit none
          real(DP) x(3),y(3),xmy(3)

          xmy = x-y
          return
        end subroutine sub
!     ----------------------------------------------------------------
        subroutine  cross(xcy,x,y)
        use constants
          implicit none
          real(DP) x(3),y(3),xcy(3)

          xcy(1) = x(2)*y(3)-y(2)*x(3)
          xcy(2) = x(3)*y(1)-y(3)*x(1)
          xcy(3) = x(1)*y(2)-y(1)*x(2)
          return
        end subroutine cross
!     ----------------------------------------------------------------
        subroutine dotangle(alpha,x,y)
        use constants

          implicit none
          real(DP) x(3),y(3),xy,alpha,xn,yn

          call dot(xy,x,y)
          xn = sqrt(x(1)**2+x(2)**2+x(3)**2)
          yn = sqrt(y(1)**2+y(2)**2+y(3)**2)
!          call norm2(xn,x)
!          call norm2(yn,y)
          
          alpha = acos( xy / (xn*yn) ) 

          return
        end subroutine dotangle
!------------------------------------------------------
        subroutine calculate_volume_cells(snv,env,snc,enc,xyz,vert_of_cell,vol)
        use constants
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snc,enc,v1,v2,v3,v4,i
        integer, dimension (4,snc:enc), intent(in) :: vert_of_cell
        real(DP), dimension (3,snv:env), intent(in) ::xyz
        real(DP), dimension (snc:enc), intent(out) :: vol
        real(DP) :: d12,d23
        real(DP) :: x1, x2, x3, x4
        real(DP) :: y1, y2, y3, y4
        real(DP) :: z1, z2, z3, z4
        real(DP) :: v
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: vert_of_cell, xyz, vol
#endif

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snc,enc
          v1=vert_of_cell(1,i)
          v2=vert_of_cell(2,i)
          v3=vert_of_cell(3,i)
          v4=vert_of_cell(4,i)


          x1 = xyz(1,v1)
          x2 = xyz(1,v2)
          x3 = xyz(1,v3)
          x4 = xyz(1,v4)

          y1 = xyz(2,v1)
          y2 = xyz(2,v2)
          y3 = xyz(2,v3)
          y4 = xyz(2,v4)

          z1 = xyz(3,v1)
          z2 = xyz(3,v2)
          z3 = xyz(3,v3)
          z4 = xyz(3,v4)

          d12 = (x1-x4)*(y2-y4)*(z3-z4) + (x2-x4)*(y3-y4)*(z1-z4) + (x3-x4)*(y1-y4)*(z2-z4)  
          d23 = (x2-x4)*(y1-y4)*(z3-z4) + (x1-x4)*(y3-y4)*(z2-z4) + (x3-x4)*(y2-y4)*(z1-z4)  
          
          v = abs(d12-d23)/6.d0
          vol(i) = v
        enddo

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_volume_cells
!------------------------------------------------------
        subroutine calculate_fiberdirection(snv,env,snc,enc,xyz,cell_bar,AmatrFibers)
        use constants
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snc,enc,v1,v2,v3,v4,i
        real(DP), dimension (3,snc:enc), intent(in) :: cell_bar
        real(DP), dimension (3,snv:env), intent(in) ::xyz
        real(DP), dimension (3,3,snc:enc), intent(out) :: AmatrFibers
        real(DP) :: x1, x2, x3, x4
        real(DP) :: y1, y2, y3, y4
        real(DP) :: z1, z2, z3, z4
        real(DP) :: xversZ,yversZ,zversZ,normxvett1
        real(DP) :: xM_ce1,yM_ce1,zM_ce1,xvett1,yvett1,zvett1
        real(DP) :: xfv,yfv,zfv,xsv,ysv,zsv,xnv,ynv,znv
        real(DP) :: phif
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: cell_bar, xyz, AmatrFibers
#endif

        phif = PI/4.D0
        xversZ=0.0
        yversZ=0.0
        zversZ=1.0

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snc,enc
           xM_ce1 = cell_bar(1,i)
           yM_ce1 = cell_bar(2,i)
           zM_ce1 = cell_bar(3,i)
           
           xvett1 = yM_ce1
           yvett1 =-xM_ce1
           zvett1 = 0.0D0
           normxvett1 = sqrt(xvett1**2+yvett1**2+zvett1**2)
           xvett1 = yM_ce1/normxvett1
           yvett1 =-xM_ce1/normxvett1
           zvett1 = 0.0D0/normxvett1
           
           xfv = cos(phif)*xvett1 + sin(phif)*xversZ;
           yfv = cos(phif)*yvett1 + sin(phif)*yversZ;
           zfv = cos(phif)*zvett1 + sin(phif)*zversZ;

           xsv =-sin(phif)*xvett1 + cos(phif)*xversZ;
           ysv =-sin(phif)*yvett1 + cos(phif)*yversZ;
           zsv =-sin(phif)*zvett1 + cos(phif)*zversZ;

           xnv = yfv*zsv - zfv*ysv
           ynv = zfv*xsv - xfv*zsv
           znv = xfv*ysv - yfv*xsv

           !fill the matrix by rows
          ! if (scar_cell(i).EQ.0) then
              AmatrFibers(1,1,i) = xfv
              AmatrFibers(2,1,i) = yfv
              AmatrFibers(3,1,i) = zfv
              
              AmatrFibers(1,2,i) = xsv
              AmatrFibers(2,2,i) = ysv
              AmatrFibers(3,2,i) = zsv
           
              AmatrFibers(1,3,i) = xnv
              AmatrFibers(2,3,i) = ynv
              AmatrFibers(3,3,i) = znv
           ! else
           !    AmatrFibers(1,1,i) = 1.0D0
           !    AmatrFibers(2,1,i) = 0.0D0
           !    AmatrFibers(3,1,i) = 0.0D0
              
           !    AmatrFibers(1,2,i) = 0.0D0
           !    AmatrFibers(2,2,i) = 1.0D0
           !    AmatrFibers(3,2,i) = 0.0D0
           
           !    AmatrFibers(1,3,i) = 0.0D0
           !    AmatrFibers(2,3,i) = 0.0D0
           !    AmatrFibers(3,3,i) = 1.0D0
           ! endif
              
              
        enddo

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_fiberdirection
!------------------------------------------------------
        subroutine calculate_normals_face_cells(snv,env,snf,enf,snc,enc,xyz,vert_of_face,face_of_cell,cell_bar,normalfaceofcells_3d)
        use constants
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snf,enf,snc,enc,v1,v2,v3,v4,i,f1,jj
        integer, dimension (3,snf:enf), intent(in) :: vert_of_face
        integer, dimension (4,snc:enc), intent(in) :: face_of_cell
        real, dimension (3,snc:enc), intent(in) :: cell_bar
        real(DP), dimension (3,snv:env), intent(in) ::xyz
        real(DP), dimension (3,4,snc:enc), intent(out) :: normalfaceofcells_3d
        real(DP) :: normve,psnormalecong
        real(DP) :: xM_ce1,yM_ce1,zM_ce1
        real(DP) :: xf1M,yf1M,zf1M
        real(DP) :: xnormalefaccia,ynormalefaccia,znormalefaccia
        real(DP) :: x1, x2, x3, x4
        real(DP) :: y1, y2, y3, y4
        real(DP) :: z1, z2, z3, z4
        real(DP) :: v
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: xyz,vert_of_face, face_of_cell,cell_bar, normalfaceofcells_3d
#endif

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snc,enc
           xM_ce1 = cell_bar(1,i)
           yM_ce1 = cell_bar(2,i)
           zM_ce1 = cell_bar(3,i)
           do jj=1,4 !%loop sulle facce
              f1=face_of_cell(jj,i);
              v1=vert_of_face(1,f1);
              v2=vert_of_face(2,f1);
              v3=vert_of_face(3,f1);                
              x1=xyz(1,v1);
              y1=xyz(2,v1);
              z1=xyz(3,v1);
              x2=xyz(1,v2);
              y2=xyz(2,v2);
              z2=xyz(3,v2);
              x3=xyz(1,v3);
              y3=xyz(2,v3);
              z3=xyz(3,v3);
              xf1M=(x1+x2+x3)/3.0D0
              yf1M=(y1+y2+y3)/3.0D0
              zf1M=(z1+z2+z3)/3.0D0
! !cross([x2-x1;y2-y1;z2-z1],[x3-x1;y3-y1;z3-z1]);
              xnormalefaccia = (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1)
              ynormalefaccia = (x3-x1)*(z2-z1) - (x2-x1)*(z3-z1)
              znormalefaccia = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
              normve = sqrt(xnormalefaccia**2+ynormalefaccia**2+znormalefaccia**2)

              xnormalefaccia = xnormalefaccia/normve
              ynormalefaccia = ynormalefaccia/normve
              znormalefaccia = znormalefaccia/normve

              psnormalecong=xnormalefaccia*(xf1M-xM_ce1)+ynormalefaccia*(yf1M-yM_ce1)+znormalefaccia*(zf1M-zM_ce1)
              if (psnormalecong.LT.0) then  
                 xnormalefaccia = -xnormalefaccia
                 ynormalefaccia = -ynormalefaccia
                 znormalefaccia = -znormalefaccia
              endif
              normalfaceofcells_3d(1,jj,i) = xnormalefaccia
              normalfaceofcells_3d(2,jj,i) = ynormalefaccia
              normalfaceofcells_3d(3,jj,i) = znormalefaccia
           end do !jj
        enddo !i cells

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_normals_face_cells
!------------------------------------------------------
        subroutine calculate_bar_cells(snv,env,snc,enc,xyz,vert_of_cell,cell_bar)
        use constants
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snc,enc,v1,v2,v3,v4,i
        integer, dimension (4,snc:enc), intent(in) :: vert_of_cell
        real(DP), dimension (3,snv:env), intent(in) ::xyz
        real(DP), dimension (3,snc:enc), intent(out) :: cell_bar
        real(DP) :: x1, x2, x3, x4
        real(DP) :: y1, y2, y3, y4
        real(DP) :: z1, z2, z3, z4
        real(DP) :: xB, yB, zB
        real(DP) :: v
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: vert_of_cell, xyz, cell_bar
#endif

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snc,enc
          v1=vert_of_cell(1,i)
          v2=vert_of_cell(2,i)
          v3=vert_of_cell(3,i)
          v4=vert_of_cell(4,i)


          x1 = xyz(1,v1)
          x2 = xyz(1,v2)
          x3 = xyz(1,v3)
          x4 = xyz(1,v4)

          y1 = xyz(2,v1)
          y2 = xyz(2,v2)
          y3 = xyz(2,v3)
          y4 = xyz(2,v4)

          z1 = xyz(3,v1)
          z2 = xyz(3,v2)
          z3 = xyz(3,v3)
          z4 = xyz(3,v4)

          xB = (x1+x2+x3+x4)/4.D0
          yB = (y1+y2+y3+y4)/4.D0
          zB = (z1+z2+z3+z4)/4.D0

          cell_bar(1,i) = xB
          cell_bar(2,i) = yB
          cell_bar(3,i) = zB
        enddo

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_bar_cells
!------------------------------------------------------
        subroutine convert_geoCPU(snv,env,sne,ene,snf,enf,xyz,xyzv,xyza, &
                  vert_of_face,tri_ver,tri_vel,tri_bar,vel_tri,acc_tri)
        use constants
!!@cuf   use cudafor

        implicit none
        integer :: snv,env,sne,ene,snf,enf,v1,v2,v3,i
        real(DP), dimension (3,snv:env), intent(in) :: xyz,xyzv,xyza
        integer, dimension (3,snf:enf), intent(in) :: vert_of_face
        real(DP), dimension(9,snf:enf), intent(out) :: tri_ver,tri_vel
        real(DP), dimension(3,snf:enf), intent(out) :: tri_bar,vel_tri,acc_tri
        real(DP) :: v11, v12, v13, v21, v22, v23, v31, v32, v33
        real(DP) :: vv11, vv12, vv13, vv21, vv22, vv23, vv31, vv32, vv33
!!@cuf   integer :: istat
! #ifdef USE_CUDA
!         attributes(managed) :: xyz,xyzv,xyza
!         attributes(managed) :: vert_of_face
!         attributes(managed) :: tri_ver,tri_vel,tri_bar,vel_tri,acc_tri
! #endif

! #ifdef USE_CUDA
!         !$cuf kernel do (1)
!         write(*,*) "A",i,snf,enf
! #endif
        do i=snf,enf
          v1=vert_of_face(1,i)
          v2=vert_of_face(2,i)
          v3=vert_of_face(3,i)

          v11 = xyz(1,v1)
          v12 = xyz(2,v1)
          v13 = xyz(3,v1)

          v21 = xyz(1,v2)
          v22 = xyz(2,v2)
          v23 = xyz(3,v2)

          v31 = xyz(1,v3)
          v32 = xyz(2,v3)
          v33 = xyz(3,v3)

          tri_ver(1,i) = v11
          tri_ver(2,i) = v12
          tri_ver(3,i) = v13
          tri_ver(4,i) = v21
          tri_ver(5,i) = v22
          tri_ver(6,i) = v23
          tri_ver(7,i) = v31
          tri_ver(8,i) = v32
          tri_ver(9,i) = v33

          vv11 = xyzv(1,v1)
          vv12 = xyzv(2,v1)
          vv13 = xyzv(3,v1)

          vv21 = xyzv(1,v2)
          vv22 = xyzv(2,v2)
          vv23 = xyzv(3,v2)

          vv31 = xyzv(1,v3)
          vv32 = xyzv(2,v3)
          vv33 = xyzv(3,v3)

          tri_vel(1,i) = vv11
          tri_vel(2,i) = vv12
          tri_vel(3,i) = vv13
          tri_vel(4,i) = vv21
          tri_vel(5,i) = vv22
          tri_vel(6,i) = vv23
          tri_vel(7,i) = vv31
          tri_vel(8,i) = vv32
          tri_vel(9,i) = vv33


          ! Find triangles' baricentre
          tri_bar(1,i) = (v11+v21+v31) / 3.d0
          tri_bar(2,i) = (v12+v22+v32) / 3.d0
          tri_bar(3,i) = (v13+v23+v33) / 3.d0

          ! Find triangle's velocity and acceleration 
          vel_tri(1,i)=(xyzv(1,v1)+xyzv(1,v2)+xyzv(1,v3)) / 3.d0
          vel_tri(2,i)=(xyzv(2,v1)+xyzv(2,v2)+xyzv(2,v3)) / 3.d0
          vel_tri(3,i)=(xyzv(3,v1)+xyzv(3,v2)+xyzv(3,v3)) / 3.d0

          acc_tri(1,i)=(xyza(1,v1)+xyza(1,v2)+xyza(1,v3)) / 3.d0
          acc_tri(2,i)=(xyza(2,v1)+xyza(2,v2)+xyza(2,v3)) / 3.d0
          acc_tri(3,i)=(xyza(3,v1)+xyza(3,v2)+xyza(3,v3)) / 3.d0
        enddo

!!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        end subroutine convert_geoCPU
!     ----------------------------------------------------------------
        subroutine calculate_normalCPU(snv,env,snf,enf,xyz,vert_of_face, tri_nor)
        use constants
!!@cuf   use cudafor

        implicit none
        integer::snv,env,snf,enf,v1,v2,v3,i
        integer,dimension(3,snf:enf),intent(in) :: vert_of_face
        real(DP),dimension(3,snv:env),intent(in) :: xyz
        real(DP),dimension(3,snf:enf),intent(out) :: tri_nor
        real(DP) :: ve11, ve12, ve13
        real(DP) :: ve21, ve22, ve23
        real(DP) :: fn1, fn2, fn3, fnmag
! !@cuf   integer :: istat
! #ifdef USE_CUDA
!         attributes(managed) :: vert_of_face, xyz, tri_nor
! #endif

! #ifdef USE_CUDA
!         !$cuf kernel do (1)
! #endif
        do i = snf,enf
          v1 = vert_of_face(1,i)
          v2 = vert_of_face(2,i)
          v3 = vert_of_face(3,i)

          ve11 = xyz(1,v2)-xyz(1,v1)
          ve12 = xyz(2,v2)-xyz(2,v1)
          ve13 = xyz(3,v2)-xyz(3,v1)

          ve21 = xyz(1,v3)-xyz(1,v1)
          ve22 = xyz(2,v3)-xyz(2,v1)
          ve23 = xyz(3,v3)-xyz(3,v1)

          fn1 = ve12*ve23 - ve13*ve22
          fn2 = ve13*ve21 - ve11*ve23
          fn3 = ve11*ve22 - ve12*ve21

          fnmag = sqrt(fn1*fn1 + fn2*fn2 + fn3*fn3)

          tri_nor(1,i) = fn1 / fnmag
          tri_nor(2,i) = fn2 / fnmag
          tri_nor(3,i) = fn3 / fnmag

        enddo

!!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_normalCPU
!------------------------------------------------------
        subroutine calculate_areaCPU (Surface,snv,env,snf,enf,xyz,vert_of_face,sur)
        use constants
!!@cuf   use cudafor

        implicit none
        integer :: snv,env,snf,enf,v1,v2,v3,i
        integer, dimension (3,snf:enf), intent(in) :: vert_of_face
        real(DP), dimension (3,snv:env), intent(in) ::xyz
        real(DP), dimension (snf:enf), intent(out) :: sur
        real(DP), intent(out) :: Surface
        real(DP) :: d12,d23,d31,sp
        real(DP) :: x1, x2, x3
        real(DP) :: y1, y2, y3
        real(DP) :: z1, z2, z3
        real(DP) :: s
! !@cuf   integer :: istat
! #ifdef USE_CUDA
!         attributes(managed) :: vert_of_face, xyz, sur
!         attributes(device) :: Surface
! #endif

        Surface=0.0
! #ifdef USE_CUDA
!         !$cuf kernel do (1)
! #endif
        do i=snf,enf
          v1=vert_of_face(1,i)
          v2=vert_of_face(2,i)
          v3=vert_of_face(3,i)

          x1 = xyz(1,v1)
          x2 = xyz(1,v2)
          x3 = xyz(1,v3)

          y1 = xyz(2,v1)
          y2 = xyz(2,v2)
          y3 = xyz(2,v3)

          z1 = xyz(3,v1)
          z2 = xyz(3,v2)
          z3 = xyz(3,v3)

          d12 = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))
          d23 = sqrt( (x2-x3)*(x2-x3) + (y2-y3)*(y2-y3) + (z2-z3)*(z2-z3))
          d31 = sqrt( (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) + (z3-z1)*(z3-z1))

          sp = 0.5d0 * (d12+d23+d31)
          s = sqrt(sp*(sp-d12)*(sp-d23)*(sp-d31))
          sur(i) = s
          Surface = Surface + s
        enddo

!!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_areaCPU
!------------------------------------------------------




!     ---------------------------------------------------------------------
        subroutine convert_geoEF_2d(snv,env,sne,ene,snf,enf,xyz, &
                  vert_of_face,tri_bar)
        use constants
!@cuf   use cudafor

        implicit none
        integer :: snv,env,sne,ene,snf,enf,v1,v2,v3,i
        real(DP), dimension (3,snv:env), intent(in) :: xyz
        integer, dimension (3,snf:enf), intent(in) :: vert_of_face
        real(DP), dimension(3,snf:enf), intent(out) :: tri_bar
        real(DP) :: v11, v12, v13, v21, v22, v23, v31, v32, v33
        real(DP) :: vv11, vv12, vv13, vv21, vv22, vv23, vv31, vv32, vv33
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: xyz
        attributes(managed) :: vert_of_face
        attributes(managed) :: tri_bar
#endif
!        write(*,*)  "DENTRO CONVERT"
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snf,enf
          v1=vert_of_face(1,i)
          v2=vert_of_face(2,i)
          v3=vert_of_face(3,i)

          v11 = xyz(1,v1)
          v12 = xyz(2,v1)
          v13 = xyz(3,v1)

          v21 = xyz(1,v2)
          v22 = xyz(2,v2)
          v23 = xyz(3,v2)

          v31 = xyz(1,v3)
          v32 = xyz(2,v3)
          v33 = xyz(3,v3)

          ! Find triangles' baricentre
          tri_bar(1,i) = (v11+v21+v31) / 3.d0
          tri_bar(2,i) = (v12+v22+v32) / 3.d0
          tri_bar(3,i) = (v13+v23+v33) / 3.d0
        enddo

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        end subroutine convert_geoEF_2d
!     ---------------------------------------------------------------------
        subroutine convert_geoEF_1d(snv,env,sne,ene,xyz, &
                  vert_of_edge,edg_bar)
        use constants
!@cuf   use cudafor

        implicit none
        integer :: snv,env,sne,ene,v1,v2,v3,i
        real(DP), dimension (3,snv:env), intent(in) :: xyz
        integer, dimension (2,sne:ene), intent(in) :: vert_of_edge
        real(DP), dimension(3,sne:ene), intent(out) :: edg_bar
        real(DP) :: v11, v12, v13, v21, v22, v23, v31, v32, v33
        real(DP) :: vv11, vv12, vv13, vv21, vv22, vv23, vv31, vv32, vv33
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: xyz
        attributes(managed) :: vert_of_edge
        attributes(managed) :: edg_bar
#endif
!        write(*,*)  "DENTRO CONVERT"
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=sne,ene
          v1=vert_of_edge(1,i)
          v2=vert_of_edge(2,i)

          v11 = xyz(1,v1)
          v12 = xyz(2,v1)
          v13 = xyz(3,v1)

          v21 = xyz(1,v2)
          v22 = xyz(2,v2)
          v23 = xyz(3,v2)

          ! Find edges' baricentre
          edg_bar(1,i) = (v11+v21) / 2.d0
          edg_bar(2,i) = (v12+v22) / 2.d0
          edg_bar(3,i) = (v13+v23) / 2.d0
        enddo

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        end subroutine convert_geoEF_1d
!     ----------------------------------------------------------------

! !     ----------------------------------------------------------------
!         subroutine dot(xy,x,y)
!         use constants

!           implicit none
!           real(DP) x(3),y(3),xy

!           xy = x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
!           return
!         end subroutine dot
! !     ----------------------------------------------------------------
!         subroutine  sub(xmy,x,y)
!         use constants

!           implicit none
!           real(DP) x(3),y(3),xmy(3)

!           xmy = x-y
!           return
!         end subroutine sub
! !     ----------------------------------------------------------------
!         subroutine  cross(xcy,x,y)
!         use constants
!           implicit none
!           real(DP) x(3),y(3),xcy(3)

!           xcy(1) = x(2)*y(3)-y(2)*x(3)
!           xcy(2) = x(3)*y(1)-y(3)*x(1)
!           xcy(3) = x(1)*y(2)-y(1)*x(2)
!           return
!         end subroutine cross
! !     ----------------------------------------------------------------
!         subroutine dotangle(alpha,x,y)
!         use constants

!           implicit none
!           real(DP) x(3),y(3),xy,alpha,xn,yn

!           call dot(xy,x,y)
!           xn = sqrt(x(1)**2+x(2)**2+x(3)**2)
!           yn = sqrt(y(1)**2+y(2)**2+y(3)**2)
! !          call norm2(xn,x)
! !          call norm2(yn,y)
          
!           alpha = acos( xy / (xn*yn) ) 

!           return
!         end subroutine dotangle
! !------------------------------------------------------

!------------------------------------------------------
        subroutine calculate_ginterp_edge_facesEF_2d(snv,env,sne,ene,snf,enf,xyz,vert_of_edge,face_of_edge,tri_bar,versCFedge,distCFedge,g1interpedge)
        use constants
!@cuf   use cudafor

        implicit none
        integer :: snv,env,sne,ene,snf,enf,v1,v2,v3,i
        integer :: f1,f2
        integer, dimension (2,sne:ene), intent(in) :: vert_of_edge
        integer, dimension (2,sne:ene), intent(in) :: face_of_edge
        real(DP), dimension (3,snv:env), intent(in) ::xyz
        real(DP), dimension (3,snf:enf), intent(in) ::tri_bar
        real(DP), dimension (3,sne:ene), intent(out) :: versCFedge
        real(DP), dimension (sne:ene), intent(out) :: distCFedge
        real(DP), dimension (sne:ene), intent(out) :: g1interpedge
        real(DP) :: xBedge,yBedge,zBedge
        real(DP) :: xM_f1,yM_f1,zM_f1,xM_f2,yM_f2,zM_f2
        real(DP) :: distf1,distf2,normve
        real(DP) :: xvettore,yvettore,zvettore
        real(DP) :: x1, x2, x3
        real(DP) :: y1, y2, y3
        real(DP) :: z1, z2, z3
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: tri_bar,vert_of_edge,face_of_edge,xyz, versCFedge, distCFedge, g1interpedge
#endif

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=sne,ene
          v1=vert_of_edge(1,i)
          v2=vert_of_edge(2,i)
          x1 = xyz(1,v1)
          y1 = xyz(2,v1)
          z1 = xyz(3,v1)
          x2 = xyz(1,v2)
          y2 = xyz(2,v2)
          z2 = xyz(3,v2)

          xBedge = (x1+x2)/2.D0
          yBedge = (y1+y2)/2.D0
          zBedge = (z1+z2)/2.D0
          
          f1=face_of_edge(1,i)
          f2=face_of_edge(2,i)
             
          if ((f1.NE.0).AND.(f2.NE.0)) then
             xM_f1 = tri_bar(1,f1)
             yM_f1 = tri_bar(2,f1)
             zM_f1 = tri_bar(3,f1)
             xM_f2 = tri_bar(1,f2)
             yM_f2 = tri_bar(2,f2)
             zM_f2 = tri_bar(3,f2)

             distf2 = sqrt( (xBedge-xM_f2)**2+(yBedge-yM_f2)**2+(zBedge-zM_f2)**2 )
             distf1 = sqrt( (xM_f2-xM_f1)**2+(yM_f2-yM_f1)**2+(zM_f2-zM_f1)**2 )

             g1interpedge(i)=distf2/distf1
             xvettore=xM_f2-xM_f1    
             yvettore=yM_f2-yM_f1    
             zvettore=zM_f2-zM_f1    
          elseif ((f1.NE.0).AND.(f2.EQ.0)) then
             xM_f1 = tri_bar(1,f1)
             yM_f1 = tri_bar(2,f1)
             zM_f1 = tri_bar(3,f1)
             g1interpedge(i)=1
             xvettore=xBedge-xM_f1    
             yvettore=yBedge-yM_f1    
             zvettore=zBedge-zM_f1    
          elseif ((f1.EQ.0).AND.(f2.NE.0)) then
             xM_f2 = tri_bar(1,f2)
             yM_f2 = tri_bar(2,f2)
             zM_f2 = tri_bar(3,f2)
             g1interpedge(i)=0;
             xvettore=xM_f1-xBedge
             yvettore=yM_f2-yBedge      
             zvettore=zM_f2-zBedge     
          endif
          normve = sqrt(xvettore**2+yvettore**2+zvettore**2)
          versCFedge(1,i)=xvettore/normve;
          versCFedge(2,i)=yvettore/normve;
          versCFedge(3,i)=zvettore/normve;
          distCFedge(i)=normve;
        enddo

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_ginterp_edge_facesEF_2d
!------------------------------------------------------
!------------------------------------------------------
        subroutine calculate_normals_edge_facesEF_2d(snv,env,sne,ene,snf,enf,xyz,vert_of_edge,edge_of_face,tri_bar,tri_nor,normaledgeoffacesEF_2d)
        use constants
!@cuf   use cudafor

        implicit none
        integer :: snv,env,sne,ene,snf,enf,v1,v2,v3,v4,i,e1,jj
        integer, dimension (2,sne:ene), intent(in) :: vert_of_edge
        integer, dimension (3,snf:enf), intent(in) :: edge_of_face
        real, dimension (3,snf:enf), intent(in) :: tri_bar,tri_nor
        real(DP), dimension (3,snv:env), intent(in) ::xyz
        real(DP), dimension (3,3,snf:enf), intent(out) :: normaledgeoffacesEF_2d
        real(DP) :: normve,psnormalecong
        real(DP) :: xM_f1,yM_f1,zM_f1
        real(DP) :: xe1M,ye1M,ze1M
        real(DP) :: xnormalefaccia,ynormalefaccia,znormalefaccia
        real(DP) :: xnormaleedge,ynormaleedge,znormaleedge
        real(DP) :: x1, x2, x3, x4
        real(DP) :: y1, y2, y3, y4
        real(DP) :: z1, z2, z3, z4
        real(DP) :: v
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: xyz, vert_of_edge,edge_of_face,tri_bar, tri_nor,normaledgeoffacesEF_2d
#endif

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snf,enf
           xM_f1 = tri_bar(1,i)
           yM_f1 = tri_bar(2,i)
           zM_f1 = tri_bar(3,i)
           xnormalefaccia = tri_nor(1,i)
           ynormalefaccia = tri_nor(2,i)
           znormalefaccia = tri_nor(3,i)
           do jj=1,3 !%loop sugli edges
              e1=edge_of_face(jj,i);
              v1=vert_of_edge(1,e1);
              v2=vert_of_edge(2,e1);
              x1=xyz(1,v1);
              y1=xyz(2,v1);
              z1=xyz(3,v1);
              x2=xyz(1,v2);
              y2=xyz(2,v2);
              z2=xyz(3,v2);
              xe1M=(x1+x2)/2.0D0
              ye1M=(y1+y2)/2.0D0
              ze1M=(z1+z2)/2.0D0
! !cross([x2-x1;y2-y1;z2-z1],[x3-x1;y3-y1;z3-z1]);
              xnormaleedge = (y2-y1)*znormalefaccia - ynormalefaccia*(z2-z1)
              ynormaleedge = xnormalefaccia*(z2-z1) - (x2-x1)*znormalefaccia
              znormaleedge = (x2-x1)*ynormalefaccia - xnormalefaccia*(y2-y1)
              normve = sqrt(xnormaleedge**2+ynormaleedge**2+znormaleedge**2)

              xnormaleedge = xnormaleedge/normve
              ynormaleedge = ynormaleedge/normve
              znormaleedge = znormaleedge/normve

              psnormalecong=xnormaleedge*(xe1M-xM_f1)+ynormaleedge*(ye1M-yM_f1)+znormaleedge*(ze1M-zM_f1)
              if (psnormalecong.LT.0) then  
                 xnormaleedge = -xnormaleedge
                 ynormaleedge = -ynormaleedge
                 znormaleedge = -znormaleedge
              endif
              normaledgeoffacesEF_2d(1,jj,i) = xnormaleedge
              normaledgeoffacesEF_2d(2,jj,i) = ynormaleedge
              normaledgeoffacesEF_2d(3,jj,i) = znormaleedge
           end do !jj
        enddo !i cells

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_normals_edge_facesEF_2d
!------------------------------------------------------
        subroutine calculate_fiberdirectionEF_2d(snv,env,snf,enf,xyz,tri_bar,tri_nor,AmatrFibers)
        use constants
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snf,enf,v1,v2,v3,v4,i
        real(DP), dimension (3,snf:enf), intent(in) :: tri_bar
        real(DP), dimension (3,snf:enf), intent(in) :: tri_nor
        real(DP), dimension (3,snv:env), intent(in) :: xyz
        real(DP), dimension (3,3,snf:enf), intent(out) :: AmatrFibers
        real(DP) :: x1, x2, x3, x4
        real(DP) :: y1, y2, y3, y4
        real(DP) :: z1, z2, z3, z4
        real(DP) :: xversZ,yversZ,zversZ,normxvett1,normxvett2
        real(DP) :: xM_f1,yM_f1,zM_f1,xvett1,yvett1,zvett1
        real(DP) :: xBface,yBface,zBface,xvett2,yvett2,zvett2
        real(DP) :: xnormalefaccia,ynormalefaccia,znormalefaccia
        real(DP) :: xfv,yfv,zfv,xsv,ysv,zsv,xnv,ynv,znv
        real(DP) :: phif
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: tri_bar, tri_nor, xyz, AmatrFibers
#endif

        phif =-PI/4.D0
        xversZ=0.0
        yversZ=0.0
        zversZ=1.0

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snf,enf
           xBface = tri_bar(1,i)
           yBface = tri_bar(2,i)
           zBface = tri_bar(3,i)
           
           xnormalefaccia=tri_nor(1,i)
           ynormalefaccia=tri_nor(2,i)
           znormalefaccia=tri_nor(3,i)
            
           ! vx1 = cross(normalefaccia,vettz0);
           xvett1= ynormalefaccia
           yvett1=-xnormalefaccia
           zvett1= 0.0D0
           normxvett1 = sqrt(xvett1**2+yvett1**2+zvett1**2)
           xvett1= xvett1/normxvett1
           yvett1= yvett1/normxvett1
           zvett1= zvett1/normxvett1

           !vett2=cross(normalefaccia,vett1);
           xvett2=ynormalefaccia*zvett1-znormalefaccia*yvett1
           yvett2=znormalefaccia*xvett1-xnormalefaccia*zvett1
           zvett2=xnormalefaccia*yvett1-ynormalefaccia*xvett1
           normxvett2 = sqrt(xvett2**2+yvett2**2+zvett2**2)
           xvett2= xvett2/normxvett2
           yvett2= yvett2/normxvett2
           zvett2= zvett2/normxvett2

           xfv = cos(phif)*xvett1 + sin(phif)*xvett2;
           yfv = cos(phif)*yvett1 + sin(phif)*yvett2;
           zfv = cos(phif)*zvett1 + sin(phif)*zvett2;

           xsv =-sin(phif)*xvett1 + cos(phif)*xvett2;
           ysv =-sin(phif)*yvett1 + cos(phif)*yvett2;
           zsv =-sin(phif)*zvett1 + cos(phif)*zvett2;

           xnv =-xnormalefaccia
           ynv =-ynormalefaccia
           znv =-znormalefaccia
!fill the matrix by rows FV
           AmatrFibers(1,1,i) = xfv
           AmatrFibers(2,1,i) = yfv
           AmatrFibers(3,1,i) = zfv

           AmatrFibers(1,2,i) = xsv
           AmatrFibers(2,2,i) = ysv
           AmatrFibers(3,2,i) = zsv

           AmatrFibers(1,3,i) = xnv
           AmatrFibers(2,3,i) = ynv
           AmatrFibers(3,3,i) = znv
        enddo

!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        ! open(109,file='../d17eCUBOGIULIOATRIOh/fiberdirection_h13.txt')
        ! do i=snc,enc
        !    read(109,*)xfv
        !    read(109,*)yfv
        !    read(109,*)zfv
        !    read(109,*)xsv
        !    read(109,*)ysv
        !    read(109,*)zsv
        !    read(109,*)xnv
        !    read(109,*)ynv
        !    read(109,*)znv

        !    AmatrFibers(1,1,i) = xfv
        !    AmatrFibers(2,1,i) = yfv
        !    AmatrFibers(3,1,i) = zfv

        !    AmatrFibers(1,2,i) = xsv
        !    AmatrFibers(2,2,i) = ysv
        !    AmatrFibers(3,2,i) = zsv

        !    AmatrFibers(1,3,i) = xnv
        !    AmatrFibers(2,3,i) = ynv
        !    AmatrFibers(3,3,i) = znv



        !    ! AmatrFibers(1,1,i) = 1.d0
        !    ! AmatrFibers(2,1,i) = 0.d0
        !    ! AmatrFibers(3,1,i) = 0.d0

        !    ! AmatrFibers(1,2,i) = 0.d0
        !    ! AmatrFibers(2,2,i) = 1.d0
        !    ! AmatrFibers(3,2,i) = 0.d0

        !    ! AmatrFibers(1,3,i) = 0.d0
        !    ! AmatrFibers(2,3,i) = 0.d0
        !    ! AmatrFibers(3,3,i) = 1.d0

        ! enddo


        return
        end subroutine calculate_fiberdirectionEF_2d
!------------------------------------------------------
!-----------------------------------------------------
        subroutine calculate_norm (normve,snv,env,xyz)
        use constants
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snf,enf,v1,v2,v3,i
        real(DP), dimension (snv:env), intent(in) ::xyz
        real(DP), intent(out) :: normve
        real(DP) :: d12,d23,d31,sp
        real(DP) :: x1, x2, x3
        real(DP) :: y1, y2, y3
        real(DP) :: z1, z2, z3
        real(DP) :: s
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: xyz
        attributes(device) :: normve
#endif

        normve=0.0
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snv,env
          x1 = xyz(i)
          x2 = x1**2

          normve = normve + x2
        enddo

!@cuf   istat = cudaDeviceSynchronize !JDR TMP
        
        normve = sqrt(normve)

        return
        end subroutine calculate_norm
!------------------------------------------------------
!-----------------------------------------------------
        subroutine calculate_normdiff (normve,snv,env,xyz1,xyz2)
        use constants
!@cuf   use cudafor

        !questa non l'ho validata
        implicit none
        integer :: snv,env,snf,enf,v1,v2,v3,i
        real(DP), dimension (snv:env), intent(in) ::xyz1,xyz2
        real(DP), intent(out) :: normve
        real(DP) :: d12,d23,d31,sp
        real(DP) :: x1, x2, x3
        real(DP) :: y1, y2, y3
        real(DP) :: z1, z2, z3
        real(DP) :: s
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: xyz1,xyz2
        attributes(device) :: normve
#endif

        normve=0.0
#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snv,env
          x1 = xyz1(i)-xyz2(i)
          x2 = x1**2

          normve = normve + x2
        enddo

!@cuf   istat = cudaDeviceSynchronize !JDR TMP
        
        normve = sqrt(normve)

        return
        end subroutine calculate_normdiff
!------------------------------------------------------
!------------------------------------------------------
        subroutine calculate_ECG(ECG_Pot,snc,enc,ECG_E,cell_bar,vol_3d,grad,LabelSmooth_cell,cell_to_chamb_3d,ECG_D)
        use constants
        use mls_param, only:LabelStenosi,LabelSkipConn
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snf,enf,v1,v2,v3,i,snc,enc
        real(DP), intent(in) :: ECG_D
        real(DP), dimension (3), intent(in) :: ECG_E
        real(DP), dimension (3,snc:enc), intent(in) ::cell_bar
        real(DP), dimension (snc:enc), intent(in) ::vol_3d
        real(DP), dimension (snc:enc), intent(in) ::LabelSmooth_cell
        real(DP), dimension (3,snc:enc), intent(in) ::grad
        integer, dimension (snc:enc), intent(in) ::cell_to_chamb_3d
        real(DP), intent(out) :: ECG_Pot
        real(DP) :: d12,d23,d31,sp
        real(DP) :: xEc,yEc,zEc
        real(DP) :: xBc,yBc,zBc
        real(DP) :: z1, z2, z3
        real(DP) :: vl,Rm1,ECG_Pot_tmp
        real(DP) :: ECG_gradRm1_x,ECG_gradRm1_y,ECG_gradRm1_z,chamb
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: cell_bar,vol_3d,ECG_E,grad,LabelSmooth_cell,cell_to_chamb_3d
        attributes(device) :: ECG_Pot
#endif

        ECG_Pot=0.0
!        ECG_Pot_tmp=0.0

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snc,enc
!           if (LabelStenosi(i).EQ.0) then
        if ((LabelStenosi(i).EQ.0).AND.(LabelSkipConn(i).EQ.0)) then
        ! if ((LabelSkipConn(i).EQ.0)) then
           chamb=cell_to_chamb_3d(i)
           if (chamb.LE.4) then
           
           xEc=ECG_E(1)
           yEc=ECG_E(2)
           zEc=ECG_E(3)
           
           xBc=cell_bar(1,i)
           yBc=cell_bar(2,i)
           zBc=cell_bar(3,i)

          vl = vol_3d(i)

           Rm1=1./ ( sqrt( (xBc-xEc)**2 +(yBc-yEc)**2 +(zBc-zEc)**2  ) )
           ECG_gradRm1_x=-(xBc-xEc)*(Rm1**3);
           ECG_gradRm1_y=-(yBc-yEc)*(Rm1**3);
           ECG_gradRm1_z=-(zBc-zEc)*(Rm1**3);

          ECG_Pot=ECG_Pot + (grad(1,i)*ECG_gradRm1_x + grad(2,i)*ECG_gradRm1_y + grad(3,i)*ECG_gradRm1_z)*vl*ECG_D*LabelSmooth_cell(i)**6
          endif

        endif
        enddo
!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_ECG
!------------------------------------------------------
!------------------------------------------------------
        subroutine calculate_ECG_Atr(ECG_Pot,snc,enc,ECG_E,cell_bar,vol_3d,grad,LabelSmooth_cell,cell_to_chamb_3d,ECG_D)
        use constants
        use mls_param, only:LabelStenosi,LabelSkipConn
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snf,enf,v1,v2,v3,i,snc,enc
        real(DP), intent(in) :: ECG_D
        real(DP), dimension (3), intent(in) :: ECG_E
        real(DP), dimension (3,snc:enc), intent(in) ::cell_bar
        real(DP), dimension (snc:enc), intent(in) ::vol_3d
        real(DP), dimension (snc:enc), intent(in) ::LabelSmooth_cell
        real(DP), dimension (3,snc:enc), intent(in) ::grad
        integer, dimension (snc:enc), intent(in) ::cell_to_chamb_3d
        real(DP), intent(out) :: ECG_Pot
        real(DP) :: d12,d23,d31,sp
        real(DP) :: xEc,yEc,zEc
        real(DP) :: xBc,yBc,zBc
        real(DP) :: z1, z2, z3
        real(DP) :: vl,Rm1,ECG_Pot_tmp
        real(DP) :: ECG_gradRm1_x,ECG_gradRm1_y,ECG_gradRm1_z,chamb
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: cell_bar,vol_3d,ECG_E,grad,LabelSmooth_cell,cell_to_chamb_3d
        attributes(device) :: ECG_Pot
#endif

        ECG_Pot=0.0
!        ECG_Pot_tmp=0.0

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snc,enc
        if ((LabelStenosi(i).EQ.0).AND.(LabelSkipConn(i).EQ.0)) then
        ! if ((LabelSkipConn(i).EQ.0)) then
           chamb=cell_to_chamb_3d(i)
           if ((chamb.EQ.2).OR.(chamb.EQ.4)) then
           
           xEc=ECG_E(1)
           yEc=ECG_E(2)
           zEc=ECG_E(3)
           
           xBc=cell_bar(1,i)
           yBc=cell_bar(2,i)
           zBc=cell_bar(3,i)

           vl = vol_3d(i)

           Rm1=1./ ( sqrt( (xBc-xEc)**2 +(yBc-yEc)**2 +(zBc-zEc)**2  ) )
           ECG_gradRm1_x=-(xBc-xEc)*(Rm1**3);
           ECG_gradRm1_y=-(yBc-yEc)*(Rm1**3);
           ECG_gradRm1_z=-(zBc-zEc)*(Rm1**3);

          ECG_Pot=ECG_Pot + (grad(1,i)*ECG_gradRm1_x + grad(2,i)*ECG_gradRm1_y + grad(3,i)*ECG_gradRm1_z)*vl*ECG_D*LabelSmooth_cell(i)**6
          endif
        endif
        enddo
!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_ECG_Atr
!------------------------------------------------------
!------------------------------------------------------
        subroutine calculate_ECG_Ventr(ECG_Pot,snc,enc,ECG_E,cell_bar,vol_3d,grad,LabelSmooth_cell,cell_to_chamb_3d,ECG_D)
        use constants
        use mls_param, only:LabelStenosi,LabelSkipConn
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snf,enf,v1,v2,v3,i,snc,enc
        real(DP), intent(in) :: ECG_D
        real(DP), dimension (3), intent(in) :: ECG_E
        real(DP), dimension (3,snc:enc), intent(in) ::cell_bar
        real(DP), dimension (snc:enc), intent(in) ::vol_3d
        real(DP), dimension (snc:enc), intent(in) ::LabelSmooth_cell
        real(DP), dimension (3,snc:enc), intent(in) ::grad
        integer, dimension (snc:enc), intent(in) ::cell_to_chamb_3d
        real(DP), intent(out) :: ECG_Pot
        real(DP) :: d12,d23,d31,sp
        real(DP) :: xEc,yEc,zEc
        real(DP) :: xBc,yBc,zBc
        real(DP) :: z1, z2, z3
        real(DP) :: vl,Rm1,ECG_Pot_tmp
        real(DP) :: ECG_gradRm1_x,ECG_gradRm1_y,ECG_gradRm1_z,chamb
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: cell_bar,vol_3d,ECG_E,grad,LabelSmooth_cell,cell_to_chamb_3d
        attributes(device) :: ECG_Pot
#endif

        ECG_Pot=0.0
!        ECG_Pot_tmp=0.0

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snc,enc
        if ((LabelStenosi(i).EQ.0).AND.(LabelSkipConn(i).EQ.0)) then
        ! if ((LabelSkipConn(i).EQ.0)) then
           chamb=cell_to_chamb_3d(i)
!           if (chamb.LE.4) then
           if ((chamb.EQ.1).OR.(chamb.EQ.3)) then           
           xEc=ECG_E(1)
           yEc=ECG_E(2)
           zEc=ECG_E(3)
           
           xBc=cell_bar(1,i)
           yBc=cell_bar(2,i)
           zBc=cell_bar(3,i)

           vl = vol_3d(i)

           Rm1=1./ ( sqrt( (xBc-xEc)**2 +(yBc-yEc)**2 +(zBc-zEc)**2  ) )
           ECG_gradRm1_x=-(xBc-xEc)*(Rm1**3);
           ECG_gradRm1_y=-(yBc-yEc)*(Rm1**3);
           ECG_gradRm1_z=-(zBc-zEc)*(Rm1**3);

          ECG_Pot=ECG_Pot + (grad(1,i)*ECG_gradRm1_x + grad(2,i)*ECG_gradRm1_y + grad(3,i)*ECG_gradRm1_z)*vl*ECG_D*LabelSmooth_cell(i)**6
          endif
        endif  
        enddo
!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_ECG_Ventr
!------------------------------------------------------
!------------------------------------------------------
        subroutine calculate_ECG_2D(ECG_Pot_2D,snf,enf,ECG_E,tri_barEF_2d,surEF_2d,grad,ECG_D)
        use constants
        use param, only:LSTAR
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snf,enf,v1,v2,v3,i,snc,enc
        real(DP), intent(in) :: ECG_D
        real(DP), dimension (3), intent(in) :: ECG_E
        real(DP), dimension (3,snf:enf), intent(in) ::tri_barEF_2d
        real(DP), dimension (snf:enf), intent(in) ::surEF_2d
        real(DP), dimension (3,snf:enf), intent(in) ::grad
        real(DP), intent(out) :: ECG_Pot_2D
        real(DP) :: d12,d23,d31,sp
        real(DP) :: xEc,yEc,zEc
        real(DP) :: xBc,yBc,zBc
        real(DP) :: z1, z2, z3
        real(DP) :: vl,Rm1,ECG_Pot_tmp,Purk_Thick
        real(DP) :: ECG_gradRm1_x,ECG_gradRm1_y,ECG_gradRm1_z,chamb
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: tri_barEF_2d,surEF_2d,ECG_E,grad
        attributes(device) :: ECG_Pot_2D
#endif

        ECG_Pot_2D=0.0
!        ECG_Pot_tmp=0.0
        ! Purk_Thick=1!mm
        Purk_Thick=1/(LSTAR*1000.D0)!mm

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snf,enf

           xEc=ECG_E(1)
           yEc=ECG_E(2)
           zEc=ECG_E(3)
           
           xBc=tri_barEF_2d(1,i)
           yBc=tri_barEF_2d(2,i)
           zBc=tri_barEF_2d(3,i)

           vl = surEF_2d(i)*Purk_Thick

           Rm1=1./ ( sqrt( (xBc-xEc)**2 +(yBc-yEc)**2 +(zBc-zEc)**2  ) )
           ECG_gradRm1_x=-(xBc-xEc)*(Rm1**3);
           ECG_gradRm1_y=-(yBc-yEc)*(Rm1**3);
           ECG_gradRm1_z=-(zBc-zEc)*(Rm1**3);

          ECG_Pot_2D=ECG_Pot_2D + (grad(1,i)*ECG_gradRm1_x + grad(2,i)*ECG_gradRm1_y + grad(3,i)*ECG_gradRm1_z)*vl*ECG_D

        enddo
!@cuf   istat = cudaDeviceSynchronize !JDR TMP

        return
        end subroutine calculate_ECG_2D
!------------------------------------------------------
!------------------------------------------------------
        subroutine calculate_EGM(ECG_Pot,snc,enc,ECG_E,cell_bar,vol_3d,grad,cell_to_chamb_3d)
        use constants
        use mls_param, only:LabelStenosi,LabelSkipConn,Mintcells_3d,CARTO_Dcell3d
        use param, only:LSTAR
!@cuf   use cudafor

        implicit none
        integer :: snv,env,snf,enf,v1,v2,v3,i,snc,enc
        real(DP), dimension (3), intent(in) :: ECG_E
        real(DP), dimension (3,snc:enc), intent(in) ::cell_bar
        real(DP), dimension (snc:enc), intent(in) ::vol_3d
        real(DP), dimension (3,snc:enc), intent(in) ::grad
        integer, dimension (snc:enc), intent(in) ::cell_to_chamb_3d
        real(DP)::mint11,mint12,mint13
        real(DP)::mint21,mint22,mint23
        real(DP)::mint31,mint32,mint33 
        real(DP), intent(out) :: ECG_Pot
        real(DP) :: d12,d23,d31,sp
        real(DP) :: xEc,yEc,zEc
        real(DP) :: xBc,yBc,zBc
        real(DP) :: z1, z2, z3
        real(DP) :: vl,Rm1,ECG_Pot_tmp,EGMc1,EGMc2,EGMc3
        real(DP) :: ECG_gradRm1_x,ECG_gradRm1_y,ECG_gradRm1_z,chamb
!@cuf   integer :: istat
#ifdef USE_CUDA
        attributes(managed) :: cell_bar,vol_3d,ECG_E,grad,cell_to_chamb_3d
        attributes(device) :: ECG_Pot
#endif

        ECG_Pot=0.0d0
!        ECG_Pot_tmp=0.0

#ifdef USE_CUDA
        !$cuf kernel do (1)
#endif
        do i=snc,enc
!           if (LabelStenosi(i).EQ.0) then
        !if ((LabelStenosi(i).EQ.0).AND.(LabelSkipConn(i).EQ.0)) then
        ! if ((LabelSkipConn(i).EQ.0)) then
           chamb=cell_to_chamb_3d(i)
           if (chamb.gt.4) cycle
           
           xEc=ECG_E(1)*1000.0D0*LSTAR
           yEc=ECG_E(2)*1000.0D0*LSTAR
           zEc=ECG_E(3)*1000.0D0*LSTAR
           
           xBc=cell_bar(1,i)*1000.0D0*LSTAR
           yBc=cell_bar(2,i)*1000.0D0*LSTAR
           zBc=cell_bar(3,i)*1000.0D0*LSTAR

           vl = vol_3d(i)*(1000.0D0*LSTAR)**3

           mint11 = Mintcells_3d(1,1,i)
           mint12 = Mintcells_3d(1,2,i)
           mint13 = Mintcells_3d(1,3,i)
           mint21 = Mintcells_3d(2,1,i)
           mint22 = Mintcells_3d(2,2,i)
           mint23 = Mintcells_3d(2,3,i)
           mint31 = Mintcells_3d(3,1,i)
           mint32 = Mintcells_3d(3,2,i)
           mint33 = Mintcells_3d(3,3,i)

           Rm1=1./ ( sqrt( (xBc-xEc)**2 +(yBc-yEc)**2 +(zBc-zEc)**2  ) )
           ECG_gradRm1_x=-(xBc-xEc)*(Rm1**3);
           ECG_gradRm1_y=-(yBc-yEc)*(Rm1**3);
           ECG_gradRm1_z=-(zBc-zEc)*(Rm1**3);

           EGMc1 = CARTO_Dcell3d(i)*grad(1,i)!mint11*grad(1,i)+mint12*grad(2,i)+mint13*grad(3,i)
           EGMc2 = CARTO_Dcell3d(i)*grad(2,i)!mint21*grad(1,i)+mint22*grad(2,i)+mint23*grad(3,i)
           EGMc3 = CARTO_Dcell3d(i)*grad(3,i)!mint31*grad(1,i)+mint32*grad(2,i)+mint33*grad(3,i)
           
           ECG_Pot=ECG_Pot + (EGMc1*ECG_gradRm1_x + EGMc2*ECG_gradRm1_y + EGMc3*ECG_gradRm1_z)*vl*0.25D0/(pi*0.667D0)

        enddo
!@cuf   istat = cudaDeviceSynchronize !JDR TMP
 
        return
        end subroutine calculate_EGM
!------------------------------------------------------


