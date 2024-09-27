module trans
#ifdef USE_CUDA
  use constants
  use cudafor
contains
  attributes(global) subroutine trans_xy_kernel(idata,odata,n1,n2,n3)
    implicit none
    integer, value, intent(IN) :: n1,n2,n3
    real(DP), dimension(n1,n2,n3), device, intent(IN)  :: idata
    real(DP), dimension(n2,n1,n3), device, intent(OUT) :: odata
    integer :: i,j,k
    i = threadIdx%x + (blockIdx%x-1)*blockDim%x
    j = threadIdx%y + (blockIdx%y-1)*blockDim%y
    k = blockIdx%z
    if(i<=n2 .and. j<=n1) then
      odata(i,j,k) = idata(j,i,k)
    end if
  end subroutine trans_xy_kernel
  subroutine trans_xy(idata,odata,n1,n2,n3)
    implicit none
    integer, intent(IN) :: n1,n2,n3
    real(DP), dimension(n1,n2,n3), device, intent(IN)  :: idata
    real(DP), dimension(n2,n1,n3), device, intent(OUT) :: odata
    type(dim3) :: tBlock, grid
    tBlock = dim3(16,16,1)
    grid = dim3(ceiling(real(n2)/tBlock%x), ceiling(real(n1)/tBlock%y), n3)
    call trans_xy_kernel<<<grid,tBlock>>>(idata,odata,n1,n2,n3)
  end subroutine trans_xy
  attributes(global) subroutine trans_yx_kernel(idata,odata,n1,n2,n3)
    implicit none
    integer, value, intent(IN) :: n1,n2,n3
    complex(DP), dimension(n1,n2,n3), device, intent(IN)  :: idata
    complex(DP), dimension(n2,n1,n3), device, intent(OUT) :: odata
    integer :: i,j,k
    i = threadIdx%x + (blockIdx%x-1)*blockDim%x
    j = threadIdx%y + (blockIdx%y-1)*blockDim%y
    k = blockIdx%z
    if(i<=n2 .and. j<=n1) then
      odata(i,j,k) = idata(j,i,k)
    end if
  end subroutine trans_yx_kernel
  subroutine trans_yx(idata,odata,n1,n2,n3)
    implicit none
    integer, intent(IN) :: n1,n2,n3
    complex(DP), dimension(n1,n2,n3), device, intent(IN)  :: idata
    complex(DP), dimension(n2,n1,n3), device, intent(OUT) :: odata
    type(dim3) :: tBlock, grid
    tBlock = dim3(16,16,1)
    grid = dim3(ceiling(real(n2)/tBlock%x), ceiling(real(n1)/tBlock%y), n3)
    call trans_yx_kernel<<<grid,tBlock>>>(idata,odata,n1,n2,n3)
  end subroutine trans_yx
  attributes(global) subroutine trans_yz_kernel(idata,odata,n1,n2,n3)
    implicit none
    integer, value, intent(IN) :: n1,n2,n3
    complex(DP), dimension(n1,n2,n3), device, intent(IN)  :: idata
    complex(DP), dimension(n3,n1,n2), device, intent(OUT) :: odata
    integer :: i,j,k
    i = threadIdx%x + (blockIdx%x-1)*blockDim%x
    j = threadIdx%y + (blockIdx%y-1)*blockDim%y
    k = blockIdx%z
    if(i<=n3 .and. j<=n1) then
      odata(i,j,k) = idata(j,k,i)
    end if
  end subroutine trans_yz_kernel
  subroutine trans_yz(idata,odata,n1,n2,n3)
    implicit none
    integer, intent(IN) :: n1,n2,n3
    complex(DP), dimension(n1,n2,n3), device, intent(IN)  :: idata
    complex(DP), dimension(n3,n1,n2), device, intent(OUT) :: odata
    type(dim3) :: tBlock, grid
    tBlock = dim3(16,16,1)
    grid = dim3(ceiling(real(n3)/tBlock%x), ceiling(real(n1)/tBlock%y), n2)
    call trans_yz_kernel<<<grid,tBlock>>>(idata,odata,n1,n2,n3)
  end subroutine trans_yz
  attributes(global) subroutine trans_zy_kernel(idata,odata,n1,n2,n3)
    implicit none
    integer, value, intent(IN) :: n1,n2,n3
    complex(DP), dimension(n1,n2,n3), device, intent(IN)  :: idata
    complex(DP), dimension(n2,n3,n1), device, intent(OUT) :: odata
    integer :: i,j,k
    i = threadIdx%x + (blockIdx%x-1)*blockDim%x
    j = threadIdx%y + (blockIdx%y-1)*blockDim%y
    k = blockIdx%z
    if(i<=n2 .and. j<=n3) then
      odata(i,j,k) = idata(k,i,j)
    end if
  end subroutine trans_zy_kernel
  subroutine trans_zy(idata,odata,n1,n2,n3)
    implicit none
    integer, intent(IN) :: n1,n2,n3
    complex(DP), dimension(n1,n2,n3), device, intent(IN)  :: idata
    complex(DP), dimension(n2,n3,n1), device, intent(OUT) :: odata
    type(dim3) :: tBlock, grid
    tBlock = dim3(16,16,1)
    grid = dim3(ceiling(real(n2)/tBlock%x), ceiling(real(n3)/tBlock%y), n1)
    call trans_zy_kernel<<<grid,tBlock>>>(idata,odata,n1,n2,n3)
  end subroutine trans_zy
#endif
end module trans


