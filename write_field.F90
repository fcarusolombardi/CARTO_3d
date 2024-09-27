!============================================================
      subroutine write_grid_info0
      use mpih
      use param
      use hdf5

      IMPLICIT none

      integer hdf_error

      integer(HID_T) :: file_grid
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid

      integer(HSIZE_T) :: dims_grid(1)

      character namfile*70

      call h5open_f(hdf_error)

      IF(myid.eq.0)THEN      

      namfile='field_gridc.h5'

      call h5fcreate_f(namfile,H5F_ACC_TRUNC_F,file_grid,hdf_error)

      dims_grid(1)=n1

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'xf', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xc(1:n1), &
             dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'xc', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xm(1:n1), &
             dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n2
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'yf', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, yc(1:n2), &
            dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'yc', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ym(1:n2), &
            dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n3m
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'zf', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zc(1:n3m), &
             dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'zc', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zm(1:n3m), &
             dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5fclose_f(file_grid, hdf_error)

      
      namfile='field_gridd.h5'

      call h5fcreate_f(namfile,H5F_ACC_TRUNC_F,file_grid,hdf_error)

      dims_grid(1)=n1r
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'xf', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xcr(1:n1r), &
             dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'xc', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xmr(1:n1r), &
             dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n2r
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'yf', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ycr(1:n2r), &
            dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'yc', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ymr(1:n2r), &
            dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n3mr
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'zf', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zcr(1:n3mr), &
             dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'zc', H5T_NATIVE_DOUBLE, &
                     dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zmr(1:n3mr), &
             dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5fclose_f(file_grid, hdf_error)


!==================================================================
!  xmf file for continua data files.

      write(namfile,'(a)')'field_ms.xmf'
      open(45,file=namfile,status='unknown')
      rewind(45)
      write(45,'("<?xml version=""1.0"" ?>")')
      write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
      write(45,'("<Xdmf Version=""2.0"">")')
      write(45,'("<Domain>")')

      write(45,'("<Grid Name=""DDC_T"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""3DRectMesh""  &
      NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2,n1
      write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
      write(45,'("<DataItem Dimensions=""",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1
      write(45,'("field_gridc.h5:/xc")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2
      write(45,'("field_gridc.h5:/yc")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m
      write(45,'("field_gridc.h5:/zc")')
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""T"" &
      AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")') n3m,n2,n1
      write(45,'("continua_dens.h5:/dens")')
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Time Value=""",e12.5,""" />")') time
      write(45,'("</Grid>")')

      write(45,'("<Grid Name=""DDC_u"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""3DRectMesh""  &
      NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2,n1
      write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
      write(45,'("<DataItem Dimensions=""",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1
      write(45,'("field_gridc.h5:/xf")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2
      write(45,'("field_gridc.h5:/yc")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m
      write(45,'("field_gridc.h5:/zc")')
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""u"" & 
      AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")') n3m,n2,n1
      write(45,'("continua_q1.h5:/Vth")')
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Time Value=""",e12.5,""" />")') time
      write(45,'("</Grid>")')

      write(45,'("<Grid Name=""DDC_v"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""3DRectMesh""  &
      NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2,n1
      write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
      write(45,'("<DataItem Dimensions=""",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1
      write(45,'("field_gridc.h5:/xc")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2
      write(45,'("field_gridc.h5:/yf")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m
      write(45,'("field_gridc.h5:/zc")')
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""v"" &
      AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")') n3m,n2,n1
      write(45,'("continua_q2.h5:/Vr")')
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Time Value=""",e12.5,""" />")') time
      write(45,'("</Grid>")')

      write(45,'("<Grid Name=""DDC_w"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""3DRectMesh""  &
      NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2,n1
      write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
      write(45,'("<DataItem Dimensions=""",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1
      write(45,'("field_gridc.h5:/xc")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2
      write(45,'("field_gridc.h5:/yc")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m
      write(45,'("field_gridc.h5:/zf")')
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""w"" &
      AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")') n3m,n2,n1
      write(45,'("continua_q3.h5:/Vz")')
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Time Value=""",e12.5,""" />")') time
      write(45,'("</Grid>")')

      write(45,'("<Grid Name=""DDC_S"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""3DRectMesh"" &
      NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3mr,n2r,n1r
      write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
      write(45,'("<DataItem Dimensions=""",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1r
      write(45,'("field_gridd.h5:/xc")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2r
      write(45,'("field_gridd.h5:/yc")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3mr
      write(45,'("field_gridd.h5:/zc")')
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""S"" &
      AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" &
      NumberType=""Float"" Precision=""4"" Format=""HDF"">")')  &
      n3mr,n2r,n1r
      write(45,'("continua_dsal.h5:/dsal")')
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Time Value=""",e12.5,""" />")') time
      write(45,'("</Grid>")')
      write(45,'("</Domain>")')
      write(45,'("</Xdmf>")')
      close(45) 

      ENDIF

      call h5close_f(hdf_error)

      return
      end subroutine write_grid_info0
!============================================================
! !  xmf file for continua data files.
!         open(45,file='restart/field_'//trim(ipfi)//'.xmf',status='unknown')
!         rewind(45)
!         write(45,'("<?xml version=""1.0"" ?>")')
!         write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
!         write(45,'("<Xdmf Version=""2.0"">")')
!         write(45,'("<Domain>")')

! #ifdef CALTEMP        
!         write(45,'("<Grid Name=""DDC_T"" GridType=""Uniform"">")')
!         write(45,'("<Topology TopologyType=""3DRectMesh""  &
!              NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2,n1
!         write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1
!         write(45,'("field_gridc.h5:/xc")')
!         write(45,'("</DataItem>")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2
!         write(45,'("field_gridc.h5:/yc")')
!         write(45,'("</DataItem>")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m
!         write(45,'("field_gridc.h5:/zc")')
!         write(45,'("</DataItem>")')
!         write(45,'("</Geometry>")')
!         write(45,'("<Attribute Name=""T"" &
!              AttributeType=""Scalar"" Center=""Node"">")')
!         write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")') n3m,n2,n1
!         write(45,'("continua_dens'//trim(ipfi)//'.h5:/dens")')
!         write(45,'("</DataItem>")')
!         write(45,'("</Attribute>")')
!         write(45,'("<Time Value=""",e12.5,""" />")') time
!         write(45,'("</Grid>")')
! #endif        

!         write(45,'("<Grid Name=""DDC_u"" GridType=""Uniform"">")')
!         write(45,'("<Topology TopologyType=""3DRectMesh""  &
!              NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2,n1
!         write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1
!         write(45,'("field_gridc.h5:/xf")')
!         write(45,'("</DataItem>")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2
!         write(45,'("field_gridc.h5:/yc")')
!         write(45,'("</DataItem>")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m
!         write(45,'("field_gridc.h5:/zc")')
!         write(45,'("</DataItem>")')
!         write(45,'("</Geometry>")')
!         write(45,'("<Attribute Name=""u"" & 
!              AttributeType=""Scalar"" Center=""Node"">")')
!         write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")') n3m,n2,n1
!         write(45,'("continua_q1'//trim(ipfi)//'.h5:/Vth")')
!         write(45,'("</DataItem>")')
!         write(45,'("</Attribute>")')
!         write(45,'("<Time Value=""",e12.5,""" />")') time
!         write(45,'("</Grid>")')
        
!         write(45,'("<Grid Name=""DDC_v"" GridType=""Uniform"">")')
!         write(45,'("<Topology TopologyType=""3DRectMesh""  &
!              NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2,n1
!         write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1
!         write(45,'("field_gridc.h5:/xc")')
!         write(45,'("</DataItem>")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2
!         write(45,'("field_gridc.h5:/yf")')
!         write(45,'("</DataItem>")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m
!         write(45,'("field_gridc.h5:/zc")')
!         write(45,'("</DataItem>")')
!         write(45,'("</Geometry>")')
!         write(45,'("<Attribute Name=""v"" &
!              AttributeType=""Scalar"" Center=""Node"">")')
!         write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")') n3m,n2,n1
!         write(45,'("continua_q2'//trim(ipfi)//'.h5:/Vr")')
!         write(45,'("</DataItem>")')
!         write(45,'("</Attribute>")')
!         write(45,'("<Time Value=""",e12.5,""" />")') time
!         write(45,'("</Grid>")')
        
!         write(45,'("<Grid Name=""DDC_w"" GridType=""Uniform"">")')
!         write(45,'("<Topology TopologyType=""3DRectMesh""  &
!              NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2,n1
!         write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1
!         write(45,'("field_gridc.h5:/xc")')
!         write(45,'("</DataItem>")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2
!         write(45,'("field_gridc.h5:/yc")')
!         write(45,'("</DataItem>")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m
!         write(45,'("field_gridc.h5:/zf")')
!         write(45,'("</DataItem>")')
!         write(45,'("</Geometry>")')
!         write(45,'("<Attribute Name=""w"" &
!              AttributeType=""Scalar"" Center=""Node"">")')
!         write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")') n3m,n2,n1
!         write(45,'("continua_q3'//trim(ipfi)//'.h5:/Vz")')
!         write(45,'("</DataItem>")')
!         write(45,'("</Attribute>")')
!         write(45,'("<Time Value=""",e12.5,""" />")') time
!         write(45,'("</Grid>")')
                
!         write(45,'("<Grid Name=""DDC_pr"" GridType=""Uniform"">")')
!         write(45,'("<Topology TopologyType=""3DRectMesh""  &
!              NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2,n1
!         write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1
!         write(45,'("field_gridc.h5:/xf")')
!         write(45,'("</DataItem>")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2
!         write(45,'("field_gridc.h5:/yf")')
!         write(45,'("</DataItem>")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m
!         write(45,'("field_gridc.h5:/zf")')
!         write(45,'("</DataItem>")')
!         write(45,'("</Geometry>")')
!         write(45,'("<Attribute Name=""pr"" &
!              AttributeType=""Scalar"" Center=""Node"">")')
!         write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")') n3m,n2,n1
!         write(45,'("continua_pr'//trim(ipfi)//'.h5:/Pr")')
!         write(45,'("</DataItem>")')
!         write(45,'("</Attribute>")')
!         write(45,'("<Time Value=""",e12.5,""" />")') time
!         write(45,'("</Grid>")')
        
! #ifdef CALSCAL
!         write(45,'("<Grid Name=""DDC_S"" GridType=""Uniform"">")')
!         write(45,'("<Topology TopologyType=""3DRectMesh"" &
!              NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3mr,n2r,n1r
!         write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1r
!         write(45,'("field_gridd.h5:/xc")')
!         write(45,'("</DataItem>")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2r
!         write(45,'("field_gridd.h5:/yc")')
!         write(45,'("</DataItem>")')
!         write(45,'("<DataItem Dimensions=""",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3mr
!         write(45,'("field_gridd.h5:/zc")')
!         write(45,'("</DataItem>")')
!         write(45,'("</Geometry>")')
!         write(45,'("<Attribute Name=""S"" &
!              AttributeType=""Scalar"" Center=""Node"">")')
!         write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" &
!              NumberType=""Float"" Precision=""4"" Format=""HDF"">")')  &
!              n3mr,n2r,n1r
!         write(45,'("continua_dsal'//trim(ipfi)//'.h5:/dsal")')
!         write(45,'("</DataItem>")')
!         write(45,'("</Attribute>")')
!         write(45,'("<Time Value=""",e12.5,""" />")') time
!         write(45,'("</Grid>")')
! #endif        
!         write(45,'("</Domain>")')
!         write(45,'("</Xdmf>")')
!         close(45) 

!================================================
      subroutine write_3dfield 
      use param
      use mpih
      use mpi_param, only: kstart,kend,kstartr,kendr
      use local_arrays, only: dens,q2,q3,q1,dsal
      use hdf5
      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_q1
      integer(HID_T) :: dset_q2
      integer(HID_T) :: dset_q3
      integer(HID_T) :: dset_dens
      integer(HID_T) :: dset_dsal

      integer(HID_T) :: plist_id

      integer(HSIZE_T), dimension(3) :: dims
      integer(HSIZE_T), dimension(3) :: data_count
      integer(HSSIZE_T), dimension(3) :: data_offset

      integer :: comm, info
      integer :: ndims

      character filnam1*30
      character filnam2*30
      character filnam3*30
      character filnam4*30
      character filnam5*30

      call h5open_f(hdf_error)

!   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!   Form the name of the file

      write(filnam1,'(a,i6.6,a)')'field/field3d_dens_',int(ntime),'.h5'
      write(filnam2,'(a,i6.6,a)')'field/field3d_velx_',int(ntime),'.h5'
      write(filnam3,'(a,i6.6,a)')'field/field3d_vely_',int(ntime),'.h5'
      write(filnam4,'(a,i6.6,a)')'field/field3d_velz_',int(ntime),'.h5'
      write(filnam5,'(a,i6.6,a)')'field/field3d_dsal_',int(ntime),'.h5'

!   Set offsets and element counts

      ndims = 3

      dims(1)=n1
      dims(2)=n2
      dims(3)=n3m

      data_count(1) = n1
      data_count(2) = n2
      data_count(3) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstart-1

!   dens

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id, &
                      hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Te', H5T_NATIVE_DOUBLE, &
                     filespace, dset_dens, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_dens, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)
      call h5dwrite_f(dset_dens, H5T_NATIVE_DOUBLE,         &
                     dens(1:n1,1:n2,kstart:kend), dims,     &
                     hdf_error, file_space_id = slabspace,  &
                     mem_space_id = memspace,               &
                     xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_dens, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!   q1

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id, &
                      hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vx', H5T_NATIVE_DOUBLE, &
                     filespace, dset_q1, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q1, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)
      call h5dwrite_f(dset_q1, H5T_NATIVE_DOUBLE,          &
                     q1(1:n1,1:n2,kstart:kend), dims,      &
                     hdf_error, file_space_id = slabspace, &
                     mem_space_id = memspace,              &
                     xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q1, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!   q2

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam3, H5F_ACC_TRUNC_F, file_id, &
                      hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vy', H5T_NATIVE_DOUBLE, &
                     filespace, dset_q2, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q2, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)
      call h5dwrite_f(dset_q2, H5T_NATIVE_DOUBLE,            &
                      q2(1:n1,1:n2,kstart:kend), dims,       &
                      hdf_error, file_space_id = slabspace,  &
                      mem_space_id = memspace,               &
                      xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q2, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!   q3

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam4, H5F_ACC_TRUNC_F, file_id, &
                      hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE, &
                     filespace, dset_q3, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_q3, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)
      call h5dwrite_f(dset_q3, H5T_NATIVE_DOUBLE,             &
                     q3(1:n1,1:n2,kstart:kend), dims,         &
                     hdf_error, file_space_id = slabspace,    &
                     mem_space_id = memspace,                 &
                     xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q3, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!   Set offsets and element counts

      ndims = 3
      dims(1)=n1r
      dims(2)=n2r
      dims(3)=n3mr

      data_count(1) = n1r
      data_count(2) = n2r
      data_count(3) = kendr-kstartr+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstartr-1

!   dsal

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam5, H5F_ACC_TRUNC_F, file_id, &
                      hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'sa', H5T_NATIVE_DOUBLE, &
                      filespace, dset_dsal, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_dsal, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)
      call h5dwrite_f(dset_dsal, H5T_NATIVE_DOUBLE,         &
                     dsal(1:n1r,1:n2r,kstartr:kendr), dims, &
                     hdf_error, file_space_id = slabspace,  &
                     mem_space_id = memspace,               &
                     xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_dsal, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      call h5close_f(hdf_error)

      return
      end subroutine write_3dfield

!================================================
      subroutine write_base
      use param
      use mpih
      use mpi_param, only: kstart,kend
      use local_arrays, only:q1, q2, q3
      use hdf5
      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_dens

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      character filnam1*30
      character filnam2*30
      character filnam3*30

      call h5open_f(hdf_error)

!   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!   Form the name of the file

      write(filnam1,'(a,i5.5,a)')'field_q1c_',int(ntime),'.h5'
      write(filnam2,'(a,i5.5,a)')'field_q2c_',int(ntime),'.h5'
      write(filnam3,'(a,i5.5,a)')'field_q3c_',int(ntime),'.h5'

!   Set offsets and element counts
   
      ndims = 3

      dims(1)=n1
      dims(2)=n2
      dims(3)=n3m

      call h5screate_simple_f(ndims, dims,  &
                             filespace, hdf_error)

      data_count(1) = n1
      data_count(2) = n2
      data_count(3) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstart-1

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id, &
                      hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'uc', H5T_NATIVE_DOUBLE, &
                     filespace, dset_dens, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_dens, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)
      call h5dwrite_f(dset_dens, H5T_NATIVE_DOUBLE,          &
                     q1(:,:,kstart:kend), dims,              & 
                     hdf_error, file_space_id = slabspace,   &
                     mem_space_id = memspace,                &
                     xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_dens, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id, &
                      hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'vc', H5T_NATIVE_DOUBLE, &
                     filespace, dset_dens, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

      call h5dget_space_f(dset_dens, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)
      call h5dwrite_f(dset_dens, H5T_NATIVE_DOUBLE,           &
                     q2(:,:,kstart:kend), dims,               &
                     hdf_error, file_space_id = slabspace,    &
                     mem_space_id = memspace,                 &
                     xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_dens, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fcreate_f(filnam3, H5F_ACC_TRUNC_F, file_id,    &
                      hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'wc', H5T_NATIVE_DOUBLE,     &
                     filespace, dset_dens, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

      call h5dget_space_f(dset_dens, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,  &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)
      call h5dwrite_f(dset_dens, H5T_NATIVE_DOUBLE,           &
                     q3(:,:,kstart:kend), dims,               &
                     hdf_error, file_space_id = slabspace,    &
                     mem_space_id = memspace,                 &
                     xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_dens, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      call h5sclose_f(filespace, hdf_error)
      call h5close_f(hdf_error)

      return
      end subroutine write_base


!================================================
      subroutine write_fine
      use param
      use mpih
      use mpi_param, only: kstartr,kendr
      use mgrd_arrays, only:q1lr, q2lr, q3lr
      use hdf5
      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_dens

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      character filnam1*30
      character filnam2*30
      character filnam3*30

      call h5open_f(hdf_error)

!   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!   Form the name of the file

      write(filnam1,'(a,i5.5,a)')'field_q1r_',int(ntime),'.h5'
      write(filnam2,'(a,i5.5,a)')'field_q2r_',int(ntime),'.h5'
      write(filnam3,'(a,i5.5,a)')'field_q3r_',int(ntime),'.h5'


!   Set offsets and element counts
   
      ndims = 3

      dims(1)=n1r
      dims(2)=n2r
      dims(3)=n3mr

      call h5screate_simple_f(ndims, dims, &
                             filespace, hdf_error)

      data_count(1) = n1r
      data_count(2) = n2r
      data_count(3) = kendr-kstartr+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstartr-1

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id, &
                      hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'ur', H5T_NATIVE_DOUBLE, &
                     filespace, dset_dens, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

      call h5dget_space_f(dset_dens, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)
      call h5dwrite_f(dset_dens, H5T_NATIVE_DOUBLE,           &
                     q1lr(:,:,kstartr:kendr), dims,           &
                     hdf_error, file_space_id = slabspace,    &
                     mem_space_id = memspace,                 &
                     xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_dens, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)


      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id, &
                      hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'vr', H5T_NATIVE_DOUBLE, &
                     filespace, dset_dens, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

      call h5dget_space_f(dset_dens, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)
      call h5dwrite_f(dset_dens, H5T_NATIVE_DOUBLE,          &
                     q2lr(:,:,kstartr:kendr), dims,          &
                     hdf_error, file_space_id = slabspace,   &
                     mem_space_id = memspace,                &
                     xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_dens, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)


      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam3, H5F_ACC_TRUNC_F, file_id,    &
                      hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'wr', H5T_NATIVE_DOUBLE,     &
                     filespace, dset_dens, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

      call h5dget_space_f(dset_dens, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)
      call h5dwrite_f(dset_dens, H5T_NATIVE_DOUBLE,          &
                     q3lr(:,:,kstartr:kendr), dims,          &
                     hdf_error, file_space_id = slabspace,   &
                     mem_space_id = memspace,                &
                     xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_dens, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)


      call h5sclose_f(filespace, hdf_error)
      call h5close_f(hdf_error)

      return
      end subroutine write_fine

