module scar_tagging
contains

  subroutine scaray(dirx, diry, dirz)
    ! This routine performs a tagging of the Eulerian grid points, generating a scalar field define dover the Euleria points
    ! that is 1 if the point lies in the first layer of external points, -1 if it lies on the first layer of internal points, 0 otherwise
    ! Ray tracing based on the Möller–Trumbore intersection algorithm

    !#--------------------------------------------------------------#
    !# version 2.0: Parallel, seems to work only with pressure grid #
    !#--------------------------------------------------------------#

    !-------!
    ! NOTES !
    !-------!
    ! in gcurv.F90 it must be activated tagging_2dwet to perform the full heart computation 

    !-------!
    ! TO DO !
    !-------!
    ! - Implement list of geometries as input 
    ! 

    !---- INPUTS --------------------------------------------------------
    ! gridflag: if 0 performs ray-tracing only for the pressure grid;
    !           if 1 performs ray-tracing for all the grids
    ! dirx, diry, dirz: flags for performing ray tracing in x, y and z directions

    !---- OUTPUTS --------------------------------------------------------
    ! color: scalar map defined over the pressure grid. It has the value:
    ! +1 - if the point lies on the first external layer
    ! -1 - if the point lies on the first internal layer
    !  0 - elsewhere
    !
    !---------------------------------------------------------------------
    use constants
    use param
    use mpih
    use rt_arrays
    use mls_param, only: scar_cell,cell_bar,xyz_gz,xyz_sc, vert_of_face_gz,vert_of_face_sc,nctot_3d, &
         vert_of_edge_gz,edge_of_face_gz,vert_of_edge_sc,edge_of_face_sc
    use mpi_param, only: kstart,kend,kstartp,kendp
    !@cuf use cudafor
    !@cuf use cudadevice

    implicit none
    !-------------!
    ! DECLARATION !
    !-------------!
    !    integer,dimension(nctot_3d), managed ::intersect_IZ
    integer, intent(in) :: dirx, diry, dirz
    integer ::negz,nvgz,nfgz,nesc,nvsc,nfsc
    character(len=4)   :: dummy
    !integer,dimension(nctot_3d),intent(out) :: scarpoint,idscar ! IDs for identifying scar regions in the myiocardium
    integer :: e1,e2,ic,i, j ,k, n, nn, i1, i2, im, ip, jm, jp, km, kp, found, isave!, ncross
    integer :: tmp
    integer :: loop, ngrid, iz
    !    integer :: nout, nint, nintx, ninty, nintz, noutx, nouty, noutz
    integer :: ibv0, ibv1, ibv2, ibv3, ibv4, ibv5, ibv6, ibsum
    integer :: v1, v2, v3,nf_gz,nf_sc
    real(DP) :: xcc, ycc, zcc, square_max_x, square_min_x, square_max_y, square_min_y, square_max_z, square_min_z
    real(DP) :: b_min_x,b_min_y,b_min_z
    real(DP) :: b_max_x,b_max_y,b_max_z
    ! rays and triangles
    real(DP), parameter :: epsilon =  1.0e-12
    real(DP) :: u, v, t, det, inv_det, xh, yh , xtmin, ytmin, ztmin, zh
    !integer, dimension(mpugeo,mpugeo), device :: intersect
    ! real(DP),  dimension(3), managed :: q, r, pvec, tvec, qvec, edge1, edge2
    ! real(DP), dimension(9), managed :: va
    real(DP) :: va1, va2, va3, va4, va5, va6, va7, va8, va9
    real(DP) :: q1, q2, q3, r1, r2, r3, pvec1, pvec2, pvec3, tvec1, tvec2, tvec3, qvec1, qvec2, qvec3
    real(DP) :: edge11, edge12, edge13, edge21, edge22, edge23
    !    real(DP),  dimension(100), managed :: zcross, zcrosso, ycross, ycrosso, xcross, xcrosso
    integer :: flag, res, Nparticle_IZ
    !integer :: intersect_IZ

    !Scar voxelization
    real(DP) :: dvox
    integer  :: nvx,nvy,nvz
    integer  :: nvmax,nemax,nfmax

    integer,dimension(nctot_3d) :: intersect_IZ
    real(DP),dimension(:,:,:),allocatable :: zcross,zcrosso
    integer,allocatable, dimension(:,:) :: icross,ncross
    real(DP):: xb,yb,zb
    real(DP), allocatable, dimension(:) :: vsx,vsy,vsz !Scar voxels baricentric coordinates
    integer, allocatable, dimension(:,:,:) :: scar_id
#ifdef USE_CUDA
!    !@cuf use cudafor
!    !@cuf use cudadevice
    attributes(device) :: tmp
    attributes(managed) :: xb,yb,zb
    attributes(managed) :: icross,ncross
    attributes(managed) :: scar_id,vsx,vsy,vsz
    attributes(managed) :: ncross,Intersect_IZ,zcross,zcrosso,icross
    !@cuf integer :: istat

#endif

    Nparticle_IZ = 1
    dvox = 0.1

    !Loop over the Gray Zone geometries
    do loop = 1,Nparticle_IZ

       !Load and allocate IZ defect geometry
       open(109,file='meshes/sur_grayzone.gts')
       read(109,*)nvgz,negz,nfgz
       close(109)
       open(110,file='meshes/sur_scarcore.gts')
       read(110,*)nvsc,nesc,nfsc
       close(110)

       !Allocate temporary array for the defect
       nvmax = max(nvgz,nvsc)
       nemax = max(negz,nesc)
       nfmax = max(nfgz,nfsc)
   
       allocate(xyz_gz(3,nvgz))
       allocate(vert_of_face_gz(3,nfgz))
       allocate(edge_of_face_gz(3,nfgz))
       allocate(vert_of_edge_gz(2,negz))
       allocate(xyz_sc(3,nvsc))
       allocate(vert_of_face_sc(3,nfsc))
       allocate(edge_of_face_sc(3,nfsc))
       allocate(vert_of_edge_sc(2,nesc))

       ! GRAY ZONE
       open(109,file='meshes/sur_grayzone.gts')
       read(109,*) dummy
       do i=1,nvgz
          read(109,*)xyz_gz(1,i),xyz_gz(2,i),xyz_gz(3,i)
       end do
       do i=1,negz
          read(109,*) vert_of_edge_gz(1,i),vert_of_edge_gz(2,i)
       end do
       do i=1,nfgz
          read(109,*) edge_of_face_gz(1,i),edge_of_face_gz(2,i),edge_of_face_gz(3,i)
       end do
       do i=1,nfgz                                                                                                                  
         e1=edge_of_face_gz(1,i)
         e2=edge_of_face_gz(2,i)
         if (vert_of_edge_gz(2,e1).eq.vert_of_edge_gz(1,e2)) then
            v1=vert_of_edge_gz(1,e1)                          
            v2=vert_of_edge_gz(2,e1)                          
            v3=vert_of_edge_gz(2,e2)
         elseif(vert_of_edge_gz(2,e1).eq.vert_of_edge_gz(2,e2)) then
            v1=vert_of_edge_gz(1,e1)                             
            v2=vert_of_edge_gz(2,e1)                             
            v3=vert_of_edge_gz(1,e2)                            
         elseif(vert_of_edge_gz(1,e1).eq.vert_of_edge_gz(1,e2)) then
            v1=vert_of_edge_gz(2,e1)                             
            v2=vert_of_edge_gz(1,e1)                             
            v3=vert_of_edge_gz(2,e2)                             
         else                    
            v1=vert_of_edge_gz(2,e1)
            v2=vert_of_edge_gz(1,e1)
            v3=vert_of_edge_gz(1,e2)
         endif               

         vert_of_face_gz(1,i)=v1
         vert_of_face_gz(2,i)=v2
         vert_of_face_gz(3,i)=v3
       end do
       close(109)
       ! SCAR CORE
       open(110,file='meshes/sur_scarcore.gts')
       read(110,*) dummy
       do i=1,nvsc
          read(110,*)xyz_sc(1,i),xyz_sc(2,i),xyz_sc(3,i)
       end do
       do i=1,nesc
          read(110,*) vert_of_edge_sc(1,i),vert_of_edge_sc(2,i)
       end do
       do i=1,nfsc
          read(110,*)edge_of_face_sc(1,i), edge_of_face_sc(2,i),edge_of_face_sc(3,i)
       end do
       do i=1,nfgz
         e1=edge_of_face_sc(1,i)
         e2=edge_of_face_sc(2,i)
         if (vert_of_edge_sc(2,e1).eq.vert_of_edge_sc(1,e2)) then
            v1=vert_of_edge_sc(1,e1)                          
            v2=vert_of_edge_sc(2,e1)                          
            v3=vert_of_edge_sc(2,e2)
         elseif(vert_of_edge_sc(2,e1).eq.vert_of_edge_sc(2,e2)) then
            v1=vert_of_edge_sc(1,e1)                             
            v2=vert_of_edge_sc(2,e1)                             
            v3=vert_of_edge_sc(1,e2)                            
         elseif(vert_of_edge_sc(1,e1).eq.vert_of_edge_sc(1,e2)) then
            v1=vert_of_edge_sc(2,e1)                             
            v2=vert_of_edge_sc(1,e1)                             
            v3=vert_of_edge_sc(2,e2)                             
         else                    
            v1=vert_of_edge_sc(2,e1)
            v2=vert_of_edge_sc(1,e1)
            v3=vert_of_edge_sc(1,e2)
         endif               

         vert_of_face_sc(1,i)=v1
         vert_of_face_sc(2,i)=v2
         vert_of_face_sc(3,i)=v3
      end do
      
       close(110)
       write(*,*) "Defect readen"
       
       b_min_x = minval(xyz_gz(1,:))
       b_max_x = maxval(xyz_gz(1,:))
       b_min_y = minval(xyz_gz(2,:))
       b_max_y = maxval(xyz_gz(2,:))
       b_min_z = minval(xyz_gz(3,:))
       b_max_z = maxval(xyz_gz(3,:))
       
       write(*,*) 'x min',b_min_x
       write(*,*) 'x max',b_max_x
       write(*,*) 'y min',b_min_y
       write(*,*) 'y max',b_max_y
       write(*,*) 'z min',b_min_z
       write(*,*) 'z max',b_max_z
       !Voxelize the scar region
       nvx = ceiling((abs(b_max_x-b_min_x)+0.1)/dvox)+2.
       nvy = ceiling((abs(b_max_y-b_min_y)+0.1)/dvox)+2.
       nvz = ceiling((abs(b_max_z-b_min_z)+0.1)/dvox)+2.
       write(*,*) "Cubic voxel edge length:",dvox
       write(*,*) "Number of voxels in x,y,z direction :",nvx,nvy,nvz
       write(*,*) "Bounding box Lx:",((b_max_x-b_min_x)+0.1)
       write(*,*) "Bounding box Ly:",((b_max_y-b_min_y)+0.1)
       write(*,*) "Bounding box Lz:",((b_max_z-b_min_z)+0.1)

       
       !Baricenters of the voxels in the bounding box
       allocate(vsx(nvx))
       allocate(vsy(nvy))
       allocate(vsz(nvz))
       allocate(scar_id(nvx,nvy,nvz))
       allocate(zcross(nvx,nvy,nfmax))
       allocate(zcrosso(nvx,nvy,nfmax))
       allocate(icross(nvx,nvy))
       allocate(ncross(nvx,nvy))
       
!       icross(:,:)=0.
       scar_id(:,:,:)=1.
       zcross(:,:,:)=0.
       icross(:,:)=0.
       ncross(:,:)=0.
       vsx(:)=0.
       vsy(:)=0.
       vsz(:)=0.
       vsx(1) = b_min_x-0.5*dvox
       vsy(1) = b_min_y-0.5*dvox
       vsz(1) = b_min_z-0.5*dvox
       do i = 2,nvx
          vsx(i) = vsx(i-1) + dvox
       end do !i
       do i = 2,nvy
          vsy(i) = vsy(i-1) + dvox
       end do !i
       do i = 2,nvz
          vsz(i) = vsz(i-1) + dvox
       end do !i

       write(*,*)'Start ray-tracing'
       !-----------------------------------!
       !-----------------------------------!
       ! RAY TRACING IN Z, X, Y DIRECTIONS !
       !-----------------------------------!
       !-----------------------------------!
       ! The rays are parallel to coordinate axis and parametrized as:
       ! ray(t) = q + t * r

       !------------------------------------------------!
       ! Start the ray-tracing algorithm in z direction !
       !------------------------------------------------!
       ! performs the ray tracing algorithm only if dirz is 1
       if (dirz.eq.1) then
          write(*,*) 'Z'

          ! direction of the ray (it is fixed, does not depend on i, j)
          r1 = 0
          r2 = 0
          r3 = vsz(1) - vsz(nvz)

          ! the z coordinate of the ray's origin is fixed
          q3 = vsz(1)        

          ! loop all over the x-y plane's points and all the triangles
#ifdef USE_CUDA
          !$cuf kernel do(2)
#endif
          write(*,*)'Loops i-j GZ'
          do j = 1, nvy
             do i = 1, nvx

                ! Fix an Eulerian point on the x-y plane as origin of the ray
                q1 = vsx(i)
                q2 = vsy(j)

                do n = 1, nfgz

                   ! Vertices of the n-th triangle
                   v1 = vert_of_face_gz(1,n)
                   v2 = vert_of_face_gz(2,n)
                   v3 = vert_of_face_gz(3,n)

                   ! fix the vertices' coordinates of the triangle 
                   va1 = xyz_gz(1,v1)
                   va2 = xyz_gz(2,v1)
                   va3 = xyz_gz(3,v1)
                   va4 = xyz_gz(1,v2)
                   va5 = xyz_gz(2,v2)
                   va6 = xyz_gz(3,v2)
                   va7 = xyz_gz(1,v3)
                   va8 = xyz_gz(2,v3)
                   va9 = xyz_gz(3,v3)


                   square_max_x = max(va1,va4,va7)
                   square_min_x = min(va1,va4,va7)
                   square_max_y = max(va2,va5,va8)
                   square_min_y = min(va2,va5,va8)              

                   if ((q2.ge.square_min_y.and.q2.le.square_max_y).and.       &
                        (q1.ge.square_min_x.and.q1.le.square_max_x)) then

                      ! INTERSECTION ROUTINE - Z

                      ! initialize the parameters
                      t = 0. 
                      u = 0.
                      v = 0.

                      ! vectors of the edges
                      edge11 = va4 - va1
                      edge12 = va5 - va2
                      edge13 = va6 - va3

                      edge21 = va7 - va1
                      edge22 = va8 - va2
                      edge23 = va9 - va3

                      !----------------------------!
                      ! INTERSECTION ALGORITHM - Z !
                      !----------------------------!

                      ! calculate the volumes
                      pvec1 = r2*edge23-edge22*r3
                      pvec2 = r3*edge21-edge23*r1
                      pvec3 = r1*edge22-edge21*r2

                      det = edge11*pvec1+edge12*pvec2+edge13*pvec3
                      inv_det = 1.0 / det

                      ! calculate auxiliary vector
                      tvec1 = q1 - va1
                      tvec2 = q2 - va2
                      tvec3 = q3 - va3

                      ! calculate u parameter
                      u = (tvec1*pvec1 + tvec2*pvec2 + tvec3*pvec3)*inv_det

                      ! if u is not in [0, 1] then no intersection occurs
                      if (u.lt.0.0.or.u.gt.1.0) cycle !go to 112 ! skip to the next triangle

                      ! another auxiliary vector
                      qvec1 = tvec2*edge13 - edge12*tvec3
                      qvec2 = tvec3*edge11 - edge13*tvec1
                      qvec3 = tvec1*edge12 - edge11*tvec2

                      ! calculate v parameter
                      v = (r1*qvec1 + r2*qvec2 + r3*qvec3)*inv_det

                      ! if v is below 0 or u+v is above 1 then no intersection occurs
                      if (v.lt.0.0.or.(u+v).gt.1.0) cycle !go to 112 ! skip to the next triangle

                      ! if we arrive here the intersection occurs

                      ! parameter of intersection
                      t = (edge21*qvec1 + edge22*qvec2 + edge23*qvec3)*inv_det

                      ! z coordinate of intersection
                      zh = q3+t*r3

                      ! end do
                      icross(i,j) = icross(i,j)+1
                      zcross(i,j,icross(i,j)) = zh

                      ! ------------------------------------!
                      ! END OF THE Z-INTERSECTION ALGORITHM !
                      ! ------------------------------------!
                   end if !q1 q2 check square_max/min

                   !112                continue ! jump here if no intersection occurs

                end do !n
                
             end do !i
          end do!j

          !@cuf istat = cudaDeviceSynchronize

                
       end if ! dirz
    end do !loop

    
#ifdef USE_CUDA
    !$cuf kernel do(2)
#endif    
    do j=1,nvy
       do i=1,nvx
          ! total number of triangles that intersect the ray 
          ncross(i,j) = icross(i,j)

          ! Eliminate duplicate intersection
          if (ncross(i,j).gt.2) then
             ! fix the firts point of intersection and look for other points that have the same z-coordinate
             icross(i,j)  = 1
             zcrosso(i,j,1)=zcross(i,j,1)

             ! loop on the other intersection points
             do i1=2,ncross(i,j)
                found = 0
                ! compare the icross point with all the others
                do i2=1,icross(i,j)
                   if (abs(zcrosso(i,j,i2)-zcross(i,j,i1)).lt.1.0E-12) found = 1 ! an intersection point is duplicated
                enddo!i2

                if (found.eq.0) then ! if the point is not duplicated, check the next points
                   icross(i,j) = icross(i,j)+1
                   zcrosso(i,j,icross(i,j)) = zcross(i,j,i1)
                endif
             enddo!i1
             ncross(i,j) = icross(i,j) ! real number of intersections
             zcross(i,j,1:ncross(i,j)) = zcrosso(i,j,1:ncross(i,j)) ! z-coordinates of the real intersection
          endif
          ! sort the intersection points with increasing z: zcrosso = sort(zcross)
          do i1=1,ncross(i,j)
             ztmin = 1.0e10
             do i2=1,ncross(i,j)
                if (zcross(i,j,i2).le.ztmin) then
                   ztmin = zcross(i,j,i2)
                   isave = i2;
                endif
             enddo!i2
             zcrosso(i,j,i1) = ztmin
             zcross (i,j,isave) = 2.0e10
          enddo!i1

          do k=1,nvz
             zcc = vsz(k) ! we fix a point in the space (i, j are already fixed)

             ! if zcc lies within an odd number of intersection is ibs = -1, else ibs = 1
             do ic=1,ncross(i,j)
                if (zcrosso(i,j,ic) .le. zcc) scar_id(i,j,k) = -scar_id(i,j,k)
             enddo
          enddo

          
       enddo
    enddo
    
    !@cuf istat = cudaDeviceSynchronize

    write(*,*)'Transferring voxels informations to unstructured grid (GRAY ZONE)'

    !XEF_3d(22,:)=0.
    scar_cell(:)=0
    do i2 = 1,nctot_3d

       if ((cell_bar(1,i2).le.b_min_x).or.(cell_bar(1,i2).ge.b_max_x)) cycle
       if ((cell_bar(2,i2).le.b_min_y).or.(cell_bar(2,i2).ge.b_max_y)) cycle
       if ((cell_bar(3,i2).le.b_min_z).or.(cell_bar(3,i2).ge.b_max_z)) cycle

       do k = 1,nvz
          if (cell_bar(3,i2).ge.(vsz(k)-0.5*dvox).and.cell_bar(3,i2).le.(vsz(k)+0.5*dvox)) then

             do j = 1,nvy
                if (cell_bar(2,i2).ge.(vsy(j)-0.5*dvox).and.cell_bar(2,i2).le.(vsy(j)+0.5*dvox)) then

                   do i =1,nvx
                      if (cell_bar(1,i2).ge.(vsx(i)-0.5*dvox).and.cell_bar(1,i2).le.(vsx(i)+0.5*dvox)) then
                         if (scar_id(i,j,k).eq.-1.) then
                            !XEF_3d(22,i2)=1.
                            scar_cell(i2)=1
                         end if
                         
                      end if
                      !write(*,*)'cell',i2,'flag',XEF_3d(22,i2)
                   enddo!i
                end if
             enddo!j
          end if
       enddo!k
    enddo!i2

    write(*,*) 'Gray Zone tagged'
    
    !      !@cuf istat = cudaDeviceSynchronize
    
    ! #########################################
    !       Start Ray tracing scar core
    ! #########################################

    ! Reinitialize colormap for second intersection loop
    scar_id(:,:,:)=1.
    zcross(:,:,:)=0.
    zcrosso(:,:,:)=0.
    icross(:,:)=0.
    ncross(:,:)=0.
    
#ifdef USE_CUDA
          !$cuf kernel do(2)
#endif
          write(*,*)'Loops i-j SC'
          do j = 1, nvy
             do i = 1, nvx

                ! Fix an Eulerian point on the x-y plane as origin of the ray
                q1 = vsx(i)
                q2 = vsy(j)

                do n = 1, nfsc

                   ! Vertices of the n-th triangle
                   v1 = vert_of_face_sc(1,n)
                   v2 = vert_of_face_sc(2,n)
                   v3 = vert_of_face_sc(3,n)

                   ! fix the vertices' coordinates of the triangle 
                   va1 = xyz_sc(1,v1)
                   va2 = xyz_sc(2,v1)
                   va3 = xyz_sc(3,v1)
                   va4 = xyz_sc(1,v2)
                   va5 = xyz_sc(2,v2)
                   va6 = xyz_sc(3,v2)
                   va7 = xyz_sc(1,v3)
                   va8 = xyz_sc(2,v3)
                   va9 = xyz_sc(3,v3)


                   square_max_x = max(va1,va4,va7)
                   square_min_x = min(va1,va4,va7)
                   square_max_y = max(va2,va5,va8)
                   square_min_y = min(va2,va5,va8)              

                   if ((q2.ge.square_min_y.and.q2.le.square_max_y).and.       &
                        (q1.ge.square_min_x.and.q1.le.square_max_x)) then

                      ! INTERSECTION ROUTINE - Z

                      ! initialize the parameters
                      t = 0. 
                      u = 0.
                      v = 0.

                      ! vectors of the edges
                      edge11 = va4 - va1
                      edge12 = va5 - va2
                      edge13 = va6 - va3

                      edge21 = va7 - va1
                      edge22 = va8 - va2
                      edge23 = va9 - va3

                      !----------------------------!
                      ! INTERSECTION ALGORITHM - Z !
                      !----------------------------!

                      ! calculate the volumes
                      pvec1 = r2*edge23-edge22*r3
                      pvec2 = r3*edge21-edge23*r1
                      pvec3 = r1*edge22-edge21*r2

                      det = edge11*pvec1+edge12*pvec2+edge13*pvec3
                      inv_det = 1.0 / det

                      ! calculate auxiliary vector
                      tvec1 = q1 - va1
                      tvec2 = q2 - va2
                      tvec3 = q3 - va3

                      ! calculate u parameter
                      u = (tvec1*pvec1 + tvec2*pvec2 + tvec3*pvec3)*inv_det

                      ! if u is not in [0, 1] then no intersection occurs
                      if (u.lt.0.0.or.u.gt.1.0) cycle !go to 112 ! skip to the next triangle

                      ! another auxiliary vector
                      qvec1 = tvec2*edge13 - edge12*tvec3
                      qvec2 = tvec3*edge11 - edge13*tvec1
                      qvec3 = tvec1*edge12 - edge11*tvec2

                      ! calculate v parameter
                      v = (r1*qvec1 + r2*qvec2 + r3*qvec3)*inv_det

                      ! if v is below 0 or u+v is above 1 then no intersection occurs
                      if (v.lt.0.0.or.(u+v).gt.1.0) cycle !go to 112 ! skip to the next triangle

                      ! if we arrive here the intersection occurs

                      ! parameter of intersection
                      t = (edge21*qvec1 + edge22*qvec2 + edge23*qvec3)*inv_det

                      ! z coordinate of intersection
                      zh = q3+t*r3

                      ! end do
                      icross(i,j) = icross(i,j)+1
                      zcross(i,j,icross(i,j)) = zh

                      ! ------------------------------------!
                      ! END OF THE Z-INTERSECTION ALGORITHM !
                      ! ------------------------------------!
                   end if !q1 q2 check square_max/min

                   !112                continue ! jump here if no intersection occurs

                end do !n
                
             end do !i
          end do!j

          !@cuf istat = cudaDeviceSynchronize

    
#ifdef USE_CUDA
    !$cuf kernel do(2)
#endif    
    do j=1,nvy
       do i=1,nvx
          ! total number of triangles that intersect the ray 
          ncross(i,j) = icross(i,j)

          ! Eliminate duplicate intersection
          if (ncross(i,j).gt.2) then
             ! fix the firts point of intersection and look for other points that have the same z-coordinate
             icross(i,j)  = 1
             zcrosso(i,j,1)=zcross(i,j,1)

             ! loop on the other intersection points
             do i1=2,ncross(i,j)
                found = 0
                ! compare the icross point with all the others
                do i2=1,icross(i,j)
                   if (abs(zcrosso(i,j,i2)-zcross(i,j,i1)).lt.1.0E-12) found = 1 ! an intersection point is duplicated
                enddo!i2

                if (found.eq.0) then ! if the point is not duplicated, check the next points
                   icross(i,j) = icross(i,j)+1
                   zcrosso(i,j,icross(i,j)) = zcross(i,j,i1)
                endif
             enddo!i1
             ncross(i,j) = icross(i,j) ! real number of intersections
             zcross(i,j,1:ncross(i,j)) = zcrosso(i,j,1:ncross(i,j)) ! z-coordinates of the real intersection
          endif
          ! sort the intersection points with increasing z: zcrosso = sort(zcross)
          do i1=1,ncross(i,j)
             ztmin = 1.0e10
             do i2=1,ncross(i,j)
                if (zcross(i,j,i2).le.ztmin) then
                   ztmin = zcross(i,j,i2)
                   isave = i2;
                endif
             enddo!i2
             zcrosso(i,j,i1) = ztmin
             zcross (i,j,isave) = 2.0e10
          enddo!i1

          do k=1,nvz
             zcc = vsz(k) ! we fix a point in the space (i, j are already fixed)

             ! if zcc lies within an odd number of intersection is ibs = -1, else ibs = 1
             do ic=1,ncross(i,j)
                if (zcrosso(i,j,ic) .le. zcc) scar_id(i,j,k) = -scar_id(i,j,k)
             enddo
          enddo

          
       enddo
    enddo
    
    !@cuf istat = cudaDeviceSynchronize

    write(*,*)'Transferring voxels informations to unstructured grid (SCAR CORE)'
   
    do i2 = 1,nctot_3d

       if ((cell_bar(1,i2).le.b_min_x).or.(cell_bar(1,i2).ge.b_max_x)) cycle
       if ((cell_bar(2,i2).le.b_min_y).or.(cell_bar(2,i2).ge.b_max_y)) cycle
       if ((cell_bar(3,i2).le.b_min_z).or.(cell_bar(3,i2).ge.b_max_z)) cycle

       do k = 1,nvz
          if (cell_bar(3,i2).ge.(vsz(k)-0.5*dvox).and.cell_bar(3,i2).le.(vsz(k)+0.5*dvox)) then

             do j = 1,nvy
                if (cell_bar(2,i2).ge.(vsy(j)-0.5*dvox).and.cell_bar(2,i2).le.(vsy(j)+0.5*dvox)) then

                   do i =1,nvx
                      if (cell_bar(1,i2).ge.(vsx(i)-0.5*dvox).and.cell_bar(1,i2).le.(vsx(i)+0.5*dvox)) then
                         if (scar_id(i,j,k).eq.-1.) then
                            !XEF_3d(22,i2)=2.
                            scar_cell(i2)=2
                         end if
                         
                      end if
                      !write(*,*)'cell',i2,'flag',XEF_3d(22,i2)
                   enddo!i
                end if
             enddo!j
          end if
       enddo!k
    enddo!i2

    write(*,*) 'Scar Core tagged'
    
    deallocate(xyz_gz)
    deallocate(vert_of_face_gz)
    deallocate(edge_of_face_gz)
    deallocate(vert_of_edge_gz)
    deallocate(xyz_sc)
    deallocate(vert_of_face_sc)
    deallocate(edge_of_face_sc)
    deallocate(vert_of_edge_sc)

    write(*,*) 'Exiting Scaray'

    return

  end subroutine scaray

  ! -------------------------------------------------------------------------------------
  ! =====================================================================================

  subroutine write_scar_vtk

    use mpi_param
    use mls_param
    use mpih
    use param

    implicit none

    integer :: i,inp,itime,chamb,j,ie
    integer :: tnv, tnc,optionvtk,cstart,cend
    integer, dimension(nctot_3d) :: cell_type      
    integer, dimension(4,nctot_3d) :: vert_cell_dum
    real(DP) :: x_inf_loc,y_inf_loc,z_inf_loc,size_inf
    real(DP) :: campo_inf,num,den,d,d0,kei,epsG
    real(DP) :: tprfi,stimleads
    character*70 filname
    character*100 ipfi,ipfip




104 format((2x,f12.8))

    cell_type(:) = 10 !10 is for tretrahedral elements


    do inp=1,Nparticle_3d
       !       scalar_face(1:maxnf) = 0.0
       vert_cell_dum(1:4,cstart_3d(inp):cend_3d(inp))=vert_of_cell_3d(1:4,cstart_3d(inp):cend_3d(inp))-vstart_3d(inp)
    end do

    do inp=1,Nparticle_3d
       tnv = vend_3d(inp) - vstart_3d(inp) + 1
       tnc = cend_3d(inp) - cstart_3d(inp) + 1

       write(ipfip,94)inp
94     format(i2.2)

       filname = 'vtkfiles/scar.vtk'

       open(121,file = filname)
       !Header     
       write(121,'(a)')adjustl('# vtk DataFile Version 3.1')
       write(121,'(a)')adjustl('Stores tetrahdral mesh')
       write(121,'(a)')adjustl('ASCII')
       write(121,'(a)')adjustl('DATASET UNSTRUCTURED_GRID')
       !        write(121,*)''
       write(121,*)'POINTS ',tnv,' FLOAT'
       do i=vstart_3d(inp),vend_3d(inp)
          write(121,*)xyz_3d(1:3,i)
       end do
       write(121,*)''
       write(121,*)'CELLS ',tnc, 5*tnc
       do i=cstart_3d(inp),cend_3d(inp)
          write(121,*)'4 ',vert_cell_dum(1:4,i)
       end do
       write(121,*)''
       write(121,*)'CELL_TYPES ',tnc
       !write(121,*)cell_type(cstart_3d(inp):cend_3d(inp))
       do i=cstart_3d(inp),cend_3d(inp)
          write(121,*)cell_type(i)
       end do
       write(121,*)''
       write(121,*)'CELL_DATA ',nci_3d(inp)
       write(121,*)'Scalars SCAR  FLOAT'
       write(121,*)'LOOKUP_TABLE default'
       do i=cstart_3d(inp),cend_3d(inp)
          !write(121,*) XEF_3d(22,i)
          write(121,*) scar_cell(i)
       enddo

      ! write(121,*)''
      ! write(121,*)'CELL_DATA ',nci_3d(inp)
      ! write(121,*)'Scalars StimS1  FLOAT'
      ! write(121,*)'LOOKUP_TABLE default'
      !do i=cstart_3d(inp),cend_3d(inp)
      !    write(121,*) IstimEF_3dS1(i)
      ! enddo

     !  write(121,*)''
     !  write(121,*)'CELL_DATA ',nci_3d(inp)
     !  write(121,*)'Scalars StimS2  FLOAT'
     !  write(121,*)'LOOKUP_TABLE default'
     !  do i=cstart_3d(inp),cend_3d(inp)
     !     write(121,*) IstimEF_3dS2(i)
     !  enddo
    end do
  endsubroutine write_scar_vtk

end module scar_tagging
