!
!     This routine finds the indices in the computational grid
!     close to the physical location of the body.
!
      subroutine topogr
      use param
      use ibm_param
      use mpi_param, only: kstart,kend

      implicit none
      integer :: i,j,k,l,n,kstartp
      integer :: km,kp,jm,jp,im,ip

      real(DP)    :: xe, xem, xep
      real(DP)    :: ye, yem, yep
      real(DP)    :: ze, zem, zep
      real(DP)    :: delta1x, delta2x
      real(DP)    :: radsph

! constants 
      radsph = 0.1

      npunt=0
      npunr=0
      npunz=0
!
!     IDENTIFICATION OF THE GRID POINTS IN THE BODY
!     (BOUNDARY + INNER PART)
!
!     X_3 vertical direction
!     X_2 radial direction
!     X_1 azimuthal direction
!
!      allocate(forclo(1:n1,1:n2,kstart:kend))

      do l = 1,3 !{ start do over the 3 velocity components
      n=0
!      forclo = 0.0d0

!     l = 1   Q_1 vel. component
!     l = 2   Q_2 vel. component
!     l = 3   Q_3 vel. component
!

      do k=kstart,kend
       do j=1,n2m
        do i=1,n1m
         km=kmv(k)
         kp=kpv(k)
         xe=tm(i)
         ye=rm(j)
         ze=zm(k)
         zem=zm(km)
         zep=zm(kp)
         if(l.eq.3) then
           ze=zz(k)
           zem=zz(km)
           zep=zz(kp)
         end if
!
!    SOLID PART
!

         if(((ze-0.5)**2+(ye-0.5)**2+(xe-0.5)**2).lt.radsph**2) then 

          n=n+1
          indgeo(l,n,1)=i
          indgeo(l,n,2)=j
          indgeo(l,n,3)=k
          indgeoe(l,n,1)=i
          indgeoe(l,n,2)=j
          indgeoe(l,n,3)=k
          distb(l,n)= 0.
!          forclo(i,j,k) = 1.

         end if

          end do
        end do
      end do

      if(l.eq.1) then
        if(n.gt.mpun) write(*,*) 'Max dimension of indgeot has exceeded n=',n
        npunt= n
        npunfx(1) = npunt
        npunifx(1) = npunt
        write(6,332)npunt
 332  format(5x,'For Q_1 N ='i7)
      end if
      if(l.eq.2) then
        if(n.gt.mpun) write(*,*) 'Max dimension of indgeor has exceeded n=',n
        npunr= n
        npunfx(2) = npunr
        npunifx(2) = npunr
        write(6,331)npunr
 331  format(5x,'For Q_2 N ='i7)
      end if
      if(l.eq.3) then
        if(n.gt.mpun) write(*,*) 'Max dimension of indgeoz has exceeded n=',n
        npunz= n
        npunfx(3) = npunz
        npunifx(3) = npunz
        write(6,330)npunz
 330  format(5x,'For Q_3 N ='i7)
      end if
      end do   !} end do over the 3 velocity components 

  
!      if(allocated(forclo)) deallocate(forclo)


      return
      end

