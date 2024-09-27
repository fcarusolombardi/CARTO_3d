
module tridiag
  use constants
#ifdef USE_CUDA
  use cudafor
#endif

contains

  subroutine trisolve_cpu(n,nrhs,a,b,c,x,ldx,mask,tmp)
    implicit none
    integer, value, intent(in) :: n,nrhs,ldx
    real(DP), dimension(*), intent(in) :: a,b,c
    real(DP), dimension(ldx,nrhs) :: x,mask,tmp

    real(DP) :: m
    integer irhs,i

    do irhs = 1,nrhs
      tmp(1,irhs) = c(1) * mask(1,irhs)     /b(1)
        x(1,irhs) = x(1,irhs)/b(1)

      do i=2,N
        m = b(i) - tmp(i-1,irhs)*a(i)*mask(i,irhs)
        x(i,irhs) = (x(i,irhs) - x(i-1,irhs)*a(i)*mask(i,irhs))/m

        if (i < N) tmp(i,irhs) = c(i)*mask(i,irhs)/m
      end do

      do i=N-1,1,-1
        x(i,irhs) = x(i,irhs) - tmp(i,irhs)*x(i+1,irhs)
      end do
    end do

  end subroutine trisolve_cpu

  subroutine trisolve_periodic_cpu(n,nrhs,a,b,c,x,ldx,mask,q,s,qe)
    implicit none
    integer, value, intent(in) :: n,nrhs,ldx
    real(DP), dimension(*), intent(in) :: a,b,c
    real(DP), dimension(ldx,nrhs) :: x,mask,q,s,qe

    real(DP) :: fn, p
    integer ::  irhs,i

    do irhs = 1,nrhs
      q(1,irhs) = -c(1)*mask(1,irhs) / b(1)
      s(1,irhs) = -a(1)*mask(1,irhs) / b(1)
      fn = x(n,irhs)
      x(1,irhs) = x(1,irhs) / b(1)

      ! forward elimination sweep
      do i = 2,n
        p = 1.0_DP / (b(i) + a(i)*mask(i,irhs)*q(i-1,irhs))
        q(i,irhs) = -c(i)*mask(i,irhs) * p
        s(i,irhs) = -a(i)*mask(i,irhs) * s(i-1,irhs)*p
        x(i,irhs) = (x(i,irhs) - a(i)*mask(i,irhs)*x(i-1,irhs)) * p
      end do

      s(n,irhs) = 1.0_DP
      qe(n,irhs) = 0.0_DP

      ! backward pass
      do i = n-1,1,-1
        s(i,irhs) = s(i,irhs) + q(i,irhs)*s(i+1,irhs)
        qe(i,irhs) = x(i,irhs) + q(i,irhs)*qe(i+1,irhs)
      end do

      x(n,irhs) = (fn - c(1)*qe(1,irhs) -           &
             a(1)*mask(1,irhs)*qe(n-1,irhs)) / (c(1)*mask(1,irhs)*s(1,irhs) + &
             a(1)*mask(1,irhs)*s(n-1,irhs) + b(1))

      ! backward elimination pass
      do i = n-1,1,-1
        x(i,irhs) = x(n,irhs)*s(i,irhs) + qe(i,irhs)
      end do
    end do

    return

  end subroutine trisolve_periodic_cpu

#ifdef USE_CUDA
  attributes(global) subroutine trisolve_periodic_3D_kernel(n,nj,nk,a,b,c,x,mask,q,s,qe)
    implicit none
    integer, value, intent(in) :: n,nj,nk
    real(DP), dimension(*), device, intent(in) :: a,b,c
    real(DP), dimension(nj,n,nk), device :: x,mask,q,s,qe

    real(DP) :: fn, p, m, sp, xp, qp, qep
    integer ::  tid, i, j, k
    tid = (blockIdx%x-1) * blockDim%x + threadIdx%x
    j = mod(tid-1, nj) + 1
    k = (tid-1) / nj + 1

    if (k > nk) return

    fn = x(j,n,k)

    ! forward elimination sweep
    m = mask(j,1,k) / b(1)
    qp= -c(1) * m
    sp = -a(1) * m
    xp = x(j,1,k) / b(1)

    s(j,1,k) = sp
    x(j,1,k) = xp
    q(j,1,k) = qp

    do i = 2,n
      m = mask(j,i,k)
      p = 1.0_DP / (b(i) + a(i)*m*qp)
      qp= -c(i)*m * p
      sp = -a(i)*m * sp*p
      xp = (x(j,i,k) - a(i)*m*xp) * p
      q(j,i,k) = qp
      s(j,i,k) = sp
      x(j,i,k) = xp
    end do

    ! backward pass
    sp = 1.0_DP
    qep = 0.0_DP
    s(j,n,k) = sp
    qe(j,n,k) = qep
    do i = n-1,1,-1
      sp = s(j,i,k) + q(j,i,k) * sp
      qep = x(j,i,k) + q(j,i,k) * qep
      s(j,i,k) = sp
      qe(j,i,k) = qep
    end do

    m = mask(j,1,k)
    x(j,n,k) = (fn - c(1)*qep -           &
           a(1)*m*qe(j,n-1,k)) / (c(1)*m*sp + &
           a(1)*m*s(j,n-1,k) + b(1))

    ! backward elimination pass
    do i = n-1,1,-1
      x(j,i,k) = x(j,n,k)*s(j,i,k) + qe(j,i,k)
    end do

    return

  end subroutine trisolve_periodic_3D_kernel

  subroutine trisolve_periodic_3D_gpu(n,nj,nk,a,b,c,x,mask,xt,maskt,rwork,rwork1)
    use trans
    implicit none
    integer, intent(in) :: n,nj,nk
    real(DP), dimension(*), device, intent(IN) :: a,b,c
    real(DP), dimension(n,nj,nk), device :: x, mask, rwork, rwork1
    real(DP), dimension(nj,n,nk), device :: xt, maskt
    integer :: tBlock, grid

    call trans_xy(x,xt,n,nj,nk)
    call trans_xy(mask,maskt,n,nj,nk)

    tBlock = 128
    grid = ceiling(real(nj*nk)/tBlock)

    call trisolve_periodic_3D_kernel<<<grid,tBlock>>>(n,nj,nk,a,b,c,xt,maskt,x,rwork,rwork1)

    call trans_xy(xt,x,nj,n,nk)

  end subroutine trisolve_periodic_3D_gpu

  subroutine trisolve_periodic_3D_notrans_gpu(n,nj,nk,a,b,c,x,mask,rwork,rwork1,rwork2)
    use trans
    implicit none
    integer, intent(in) :: n,nj,nk
    real(DP), dimension(*), device, intent(IN) :: a,b,c
    real(DP), dimension(nj,n,nk), device :: x, mask, rwork, rwork1, rwork2
    integer :: tBlock, grid

    tBlock = 128
    grid = ceiling(real(nj*nk)/tBlock)

    call trisolve_periodic_3D_kernel<<<grid,tBlock>>>(n,nj,nk,a,b,c,x,mask,rwork,rwork1,rwork2)

  end subroutine trisolve_periodic_3D_notrans_gpu

  attributes(global) subroutine trisolve_kernel(n,nrhs,a,b,c,x,ldx,mask,tmp)
    implicit none
    integer, value, intent(in) :: n,nrhs,ldx
    real(DP), dimension(*), device, intent(in) :: a,b,c
    real(DP), dimension(nrhs,ldx), device :: x,mask,tmp

    real(DP) :: m, xp, tmpp, div
    integer irhs,i

    irhs = threadIdx%x + (blockIdx%x-1)*blockDim%x

    if (irhs <= nrhs) then
      tmpp = c(1) * mask(irhs,1)/b(1)
      xp = x(irhs,1)/b(1)

      tmp(irhs,1) = tmpp
      x(irhs,1) = xp

      do i = 2,n
        m = mask(irhs,i)
        div = b(i) - tmpp*a(i)*m
        xp = (x(irhs,i) - xp*a(i)*m)/div
        x(irhs,i) = xp

        if (i < n) then
          tmpp = c(i)*m/div
          tmp(irhs,i) = tmpp
        endif
      end do

      xp = x(irhs, n)
      do i = n-1,1,-1
        xp = x(irhs,i) - tmp(irhs,i)*xp
        x(irhs,i) = xp
      end do
    end if

  end subroutine trisolve_kernel

  attributes(global) subroutine tran_rhs_kernel(idata,odata,n1,n2)
    implicit none
    integer, value, intent(IN) :: n1,n2
    real(DP), dimension(n2,n1), device, intent(IN)  :: idata
    real(DP), dimension(n1,n2), device, intent(OUT) :: odata
    integer :: i,j
    i = threadIdx%x + (blockIdx%x-1)*blockDim%x
    j = threadIdx%y + (blockIdx%y-1)*blockDim%y
    if(i<=n1 .and. j<=n2) then
      odata(i,j) = idata(j,i)
    end if
  end subroutine tran_rhs_kernel

  subroutine trisolve_gpu(n,nrhs,a,b,c,x,ldx,mask,xt,maskt)
    implicit none
    integer, intent(in) :: n,nrhs,ldx
    real(DP), dimension(*), device, intent(IN) :: a,b,c
    real(DP), dimension(LDX,NRHS), device :: x, mask
    real(DP), dimension(NRHS,LDX), device :: xt, maskt
    type(dim3) :: tBlock, grid, tBlock_tran, grid_tran

    tBlock_tran = dim3(16,16,1)
    grid_tran = dim3(ceiling(real(nrhs)/tBlock_tran%x), ceiling(real(n)/tBlock_tran%y), 1)

    call tran_rhs_kernel<<<grid_tran,tBlock_tran>>>(x,xt,nrhs,ldx)
    call tran_rhs_kernel<<<grid_tran,tBlock_tran>>>(mask,maskt,nrhs,ldx)

    tBlock = dim3(128,1,1)
    grid = dim3(ceiling(real(nrhs)/tBlock%x), 1, 1)

    call trisolve_kernel<<<grid,tBlock>>>(n,nrhs,a,b,c,xt,ldx,maskt,x)

    tBlock_tran = dim3(16,16,1)
    grid_tran = dim3(ceiling(real(n)/tBlock_tran%x),ceiling(real(nrhs)/tBlock_tran%y), 1)

    call tran_rhs_kernel<<<grid_tran,tBlock_tran>>>(xt,x,ldx,nrhs)

  end subroutine trisolve_gpu

  attributes(global) subroutine tepDgtsv_pres_nopivot_kernel(n,nrhs,n1,n2,amphk,acphk,apphk,ak1,ak2,jmhv,x,ldx,tmp,joff)
    implicit none
    integer, value, intent(in) :: n,nrhs,n1,n2,ldx,joff
    integer, dimension(*), device, intent(in) :: jmhv
    real(DP), dimension(*), device, intent(in) :: amphk,acphk,apphk,ak1,ak2
    real(DP), dimension(nrhs,ldx), device :: x,tmp

    real(DP) :: m,acphT_b, a, eps
    integer irhs,i,j,k

    irhs = threadIdx%x + (blockIdx%x-1)*blockDim%x

    k = mod(irhs-1,n1)+1
    j = (irhs+n1-1)/n1 + joff - 1
    j = jmhv(j)


    if ( irhs <= nrhs) then
      eps = epsilon(0.d0)

      acphT_b = 1.d0/(acphk(1) - ak2(j) - ak1(k))
      tmp(irhs,1) = apphk(1)*acphT_b
      x(irhs,1) = x(irhs,1) * acphT_b

      do i=2,N 
        acphT_b = 1.d0/(acphk(i) - ak2(j) - ak1(k))
        a = amphk(i) * acphT_b

        ! For certain wavenumbers, system can be singular (result is ignored later). Need
        ! to clamp m to avoid true NaNs from appearing.
        m = 1.d0 - tmp(irhs,i-1)*a
        if (m == 0.d0) m = eps

        x(irhs,i) = (x(irhs,i)*acphT_b - x(irhs,i-1)*a)/m

        if (i < N) tmp(irhs,i) = apphk(i) * acphT_b/m
      end do

      do i=N-1,1,-1
        x(irhs,i) = x(irhs,i) - tmp(irhs,i)*x(irhs,i+1)
      end do
    end if

  end subroutine tepDgtsv_pres_nopivot_kernel

  subroutine tepDgtsv_pres_nopivot(n,nrhs,n1,n2,amphk,acphk,apphk,ak1,ak2,jmhv,x,ldx,xt,joff)
    !use local_arrays, ONLY: rhs_t_d
    implicit none
    integer, intent(in) :: n,nrhs,n1,n2,ldx,joff
    integer, dimension(*), device, intent(IN) :: jmhv
    real(DP), dimension(*), device, intent(IN) :: amphk, acphk, apphk, ak1, ak2
    real(DP), dimension(LDX,NRHS), device :: x
    real(DP), dimension(NRHS,LDX), device :: xt
    type(dim3) :: tBlock, grid, tBlock_tran, grid_tran
 
    tBlock_tran = dim3(16,16,1)
    grid_tran = dim3(ceiling(real(nrhs)/tBlock_tran%x), ceiling(real(n)/tBlock_tran%y), 1)

    call tran_rhs_kernel<<<grid_tran,tBlock_tran>>>(x,xt,nrhs,ldx)

    tBlock = dim3(128,1,1)
    grid = dim3(ceiling(real(nrhs)/tBlock%x), 1, 1)

    call tepDgtsv_pres_nopivot_kernel<<<grid,tBlock>>>(n,nrhs,n1,n2,amphk,acphk,apphk,ak1,ak2,jmhv,xt,ldx,x,joff)

    tBlock_tran = dim3(16,16,1)
    grid_tran = dim3(ceiling(real(n)/tBlock_tran%x),ceiling(real(nrhs)/tBlock_tran%y), 1)
 
    call tran_rhs_kernel<<<grid_tran,tBlock_tran>>>(xt,x,ldx,nrhs)

  end subroutine tepDgtsv_pres_nopivot


#endif
end module tridiag
