!================================================
!   old grids
      subroutine iniitp_grid
      use param
      use mpih
      use init_itp
      implicit none
      integer :: i,j,k,l
      integer :: nclipr, n3mallr
      real(DP) :: tstr3, z2dp
      real(DP) :: x1,x2,x3,delet
      real(DP) stfn, wdth, dxi

      real(DP),allocatable,dimension(:) :: etaz,etazm,dzr

      allocate(etaz(1:n3or),etazm(0:n3or*2))
      allocate(dzr(1:n3or*2))

!   x (theta) direction

      do i=1,n1or
        x1=dble(i-1)/dble(n1omr)
        xcrold(i) = rext1*x1
      enddo
      do i=1,n1o
        xcold(i) = xcrold((i-1)*mref1o+1)
      enddo      

      do i=1,n1om
        xmold(i)=(xcold(i)+xcold(i+1))*0.5d0
      end do
      xmold(n1o) = 2.d0*xmold(n1om)-xmold(n1om-1)

      do i=1,n1omr
       xmrold(i) = (xcrold(i+1)+xcrold(i))*0.5d0
      end do
      xmrold(n1or) = 2.d0*xmrold(n1omr)-xmrold(n1omr-1)

!   y (radial) direction

      do j=1,n2or
        x2=dble(j-1)/dble(n2omr)
        ycrold(j) = rext2*x2
      enddo
      do j=1,n2o
        ycold(j) = ycrold((j-1)*mref2o+1)
      enddo

      do j=1,n2om
        ymold(j)=(ycold(j)+ycold(j+1))*0.5d0
      end do
      ymold(n2o) = 2.d0*ymold(n2om)-ymold(n2om-1)

      do j=1,n2omr
       ymrold(j) = (ycrold(j+1)+ycrold(j))*0.5d0
      enddo
      ymrold(n2or) = 2.d0*ymrold(n2omr)-ymrold(n2omr-1)

!  z (axial) direction, with different method

      if (istr3o.eq.0) then
        do k=1,n3o
          x3=dble(k-1)/dble(n3om)
          etaz(k)=alx3*x3
          zcold(k)=etaz(k)
        enddo
        do k=1,n3om
          do l=1,mref3o
            zcrold(mref3o*(k-1)+l) = zcold(k) &
                  +(zcold(k+1)-zcold(k))*dble(l-1)/dble(mref3o)
          end do
        end do
        zcrold(n3or)=zcold(n3o)
      endif

      if (istr3o.eq.4) then
        tstr3 = dtanh(str3o)
        zcold(1)=0.0d0
        do k=2,n3o
          z2dp=float(2*k-n3o-1)/float(n3om)
          zcold(k)=(1.0+dtanh(str3o*z2dp)/tstr3)*0.5*alx3
          if(zcold(k).lt.0.0 .or. zcold(k).gt.alx3)then
            write(*,*)'Forza la griglia: ','zc(',k,')=',zcold(k)
            stop
          endif
        end do
      end if

      if(istr3o.eq.6) then
        nclipr = int(str3o)*mref3o
        n3mallr = n3or+nclipr+nclipr
        do k=1,n3mallr
          etazm(k) = dcos(pi*(dble(k)-0.5d0)/dble(n3mallr))
        end do
        do k=1,n3or
          etaz(k)=etazm(k+nclipr)
        end do
        delet = etaz(1)-etaz(n3or)
        do k=1,n3or
          etaz(k)=etaz(k)/(0.5d0*delet)
        end do
        zcrold(1) = 0.d0
        do k=2,n3omr
          zcrold(k) = alx3*(1.d0-etaz(k))*0.5d0
        end do
        zcrold(n3or) = alx3

        do k=1,n3o
          zcold(k) = zcrold((k-1)*mref3o+1)
        enddo
      endif

      if(istr3o.eq.7) then
        nclipr = int(str3o)*mref3o
        n3mallr = n3omr+nclipr+nclipr
        wdth = 2.d0
        stfn = 4.d0
        dxi = wdth / dble(n3mallr)
        do k=1,n3or
          delet = -1.d0+dxi*dble(k+nclipr-1)
          zcrold(k) = erf(stfn*delet/2.d0)/erf(stfn/2.d0)
        enddo

        zcrold = zcrold - zcrold(1)
        delet = alx3/zcrold(n3or)
        do k=1,n3or
          zcrold(k) = zcrold(k)*delet
        enddo

        do k=1,n3o
          zcold(k) = zcrold((k-1)*mref3o+1)
        end do
      endif

      if(istr3o.ge.8)then
        nclipr = int(str3o)*mref3o
        n3mallr = n3omr+nclipr+nclipr
        wdth = dble(istr3o)
        stfn = 1.d0
        dxi = wdth/dble(n3mallr)
        do k=1,n3mallr/2+1
          delet = dxi*dble(k)
          dzr(k) = dtanh(stfn*delet)-dtanh(stfn*(delet-wdth))
        enddo
        dzr = dzr - dtanh(0.d0)+dtanh(stfn*(-wdth))

        zcrold(1) = 0.d0
        do k=2,n3omr/2+1
          zcrold(k) = zcrold(k-1) + dzr(k+nclipr-1)
        enddo

        delet = 0.5d0*alx3/zcrold(n3omr/2+1)
        do k=1,n3omr/2+1
          zcrold(k) = zcrold(k)*delet
        enddo
        do k=1,n3omr/2
          zcrold(n3or-k+1) = alx3 - zcrold(k)
        enddo

        do k=1,n3o
          zcold(k) = zcrold((k-1)*mref3o+1)
        end do
      endif

      do k=1,n3om
        zmold(k)=(zcold(k)+zcold(k+1))*0.5d0
      enddo
      zmold(n3o) = alx3
      do k=1,n3omr
       zmrold(k) = (zcrold(k+1)+zcrold(k))*0.5d0
      end do
      zmrold(n3or) = alx3

      deallocate(etaz,etazm,dzr)

      end subroutine iniitp_grid

!================================================
!   index for interpolation 
      subroutine iniitp_indc
      use param
      use mpih
      use init_itp
      implicit none
      integer :: i,j,k,l

      real(DP),allocatable,dimension(:) :: xmoa,ymoa,zmoa
      real(DP),allocatable,dimension(:) :: xmroa,ymroa,zmroa

!  for de

      allocate(xmoa(0:n1o),ymoa(0:n2o),zmoa(0:n3o))

      xmoa(1:n1o)=xmold(1:n1o)
      xmoa(0) = 2.d0*xmoa(1) - xmoa(2)
      do i=1,n1m
        do l=0,n1om
          if(xm(i).ge.xmoa(l) .and. xm(i).lt.xmoa(l+1))then
            iidec1(i) = l
          endif
        enddo
      enddo
      
      ymoa(1:n2o)=ymold(1:n2o)
      ymoa(0) = 2.d0*ymoa(1) - ymoa(2)
      do j=1,n2m
        do l=0,n2om
          if(ym(j).ge.ymoa(l) .and. ym(j).lt.ymoa(l+1))then
            iidec2(j) = l
          endif
        enddo
      enddo

      zmoa(1:n3om)=zmold(1:n3om)
      zmoa(0) = 0.d0
      zmoa(n3o) = alx3
      do k=1,n3m
        do l=0,n3om
          if(zm(k).ge.zmoa(l) .and. zm(k).lt.zmoa(l+1))then
            iidec3(k) = l
          endif
        enddo
      enddo

      deallocate(xmoa,ymoa,zmoa)

!  for q1

      iiq1c1 = 1
      do i=1,n1m
        do l=1,n1om
          if(xc(i).ge.xcold(l) .and. xc(i).lt.xcold(l+1))then
            iiq1c1(i) = l
          endif
        enddo
      enddo

      iiq1c2 = iidec2
      iiq1c3 = iidec3

!  for q2

      iiq2c1 = iidec1

      iiq2c2 = 1
      do j=1,n2m
        do l=1,n2om
          if(yc(j).ge.ycold(l) .and. yc(j).lt.ycold(l+1))then
            iiq2c2(j) = l
          endif
        enddo
      enddo

      iiq2c3 = iidec3

!  for q3

      iiq3c1 = iidec1
      iiq3c2 = iidec2

      iiq3c3 = 1
      do k=1,n3m
        do l=1,n3om
          if(zc(k).ge.zcold(l) .and. zc(k).lt.zcold(l+1))then
            iiq3c3(k) = l
          endif
        enddo
      enddo

!  for sa 

      allocate(xmroa(0:n1or),ymroa(0:n2or),zmroa(0:n3or))

      xmroa(1:n1or)=xmrold(1:n1or)
      xmroa(0) = 2.d0*xmroa(1) - xmroa(2)
      do i=1,n1mr
        do l=0,n1omr
          if(xmr(i).ge.xmroa(l) .and. xmr(i).lt.xmroa(l+1))then
            iisac1(i) = l
          endif
        enddo
      enddo

      ymroa(1:n2or)=ymrold(1:n2or)
      ymroa(0) = 2.d0*ymroa(1) - ymroa(2)
      do j=1,n2mr
        do l=0,n2omr
          if(ymr(j).ge.ymroa(l) .and. ymr(j).lt.ymroa(l+1))then
            iisac2(j) = l
          endif
        enddo
      enddo

      zmroa(1:n3omr) = zmrold(1:n3omr)
      zmroa(0) = 0.d0
      zmroa(n3or) = alx3
      do k=1,n3mr
        do l=0,n3omr
          if(zmr(k).ge.zmroa(l) .and. zmr(k).lt.zmroa(l+1))then
            iisac3(k) = l
          endif
        enddo
      enddo

      deallocate(xmroa,ymroa,zmroa)

      return
      end subroutine iniitp_indc
