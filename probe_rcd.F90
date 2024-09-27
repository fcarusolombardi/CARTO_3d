      subroutine probe_rcd
      use param
      use local_arrays, only: q1,q2,q3,dens,dsal
      use mpih
      use mpi_param
      use slab_param
      implicit none

      integer is,ic,jc,kc,ir,jr,kr
      integer i0,j0,i1,j1
      real(DP),dimension(4) :: uxp,uyp,uzp,tep,sap

      i0 = 5
      j0 = 5
      i1 = i0 + 5
      j1 = j0 + 5

      IF(myid.eq.idslab(4))then
       kc = kslab(4)
       kr = kslabr(4)
       write(804,'(1x,f10.4,20(1x,ES20.8))') time,               &    
       q1(i0,j0,kc),q2(i0,j0,kc),q3(i0,j0,kc),                   &
       dens(i0,j0,kc),dsal((i0-1)*mref1+1,(j0-1)*mref2+1,kr),    &
       q1(i0,j1,kc),q2(i0,j1,kc),q3(i0,j1,kc),                   &
       dens(i0,j1,kc),dsal((i0-1)*mref1+1,(j1-1)*mref2+1,kr),    &
       q1(i1,j1,kc),q2(i1,j1,kc),q3(i1,j1,kc),                   &
       dens(i1,j1,kc),dsal((i1-1)*mref1+1,(j1-1)*mref2+1,kr),    &
       q1(i1,j0,kc),q2(i1,j0,kc),q3(i1,j0,kc),                   &
       dens(i1,j0,kc),dsal((i1-1)*mref1+1,(j0-1)*mref2+1,kr)
      ENDIF

      IF(myid.eq.idslab(5))then
       kc = kslab(5)
       kr = kslabr(5)
       write(805,'(1x,f10.4,20(1x,ES20.8))') time,               &
       q1(i0,j0,kc),q2(i0,j0,kc),q3(i0,j0,kc),                   &
       dens(i0,j0,kc),dsal((i0-1)*mref1+1,(j0-1)*mref2+1,kr),    &
       q1(i0,j1,kc),q2(i0,j1,kc),q3(i0,j1,kc),                   &
       dens(i0,j1,kc),dsal((i0-1)*mref1+1,(j1-1)*mref2+1,kr),    &
       q1(i1,j1,kc),q2(i1,j1,kc),q3(i1,j1,kc),                   &
       dens(i1,j1,kc),dsal((i1-1)*mref1+1,(j1-1)*mref2+1,kr),    &
       q1(i1,j0,kc),q2(i1,j0,kc),q3(i1,j0,kc),                   &
       dens(i1,j0,kc),dsal((i1-1)*mref1+1,(j0-1)*mref2+1,kr)
      ENDIF

      IF(myid.eq.idslab(8))then
       kc = kslab(8)
       kr = kslabr(8)
       write(808,'(1x,f10.4,20(1x,ES20.8))') time,               &
       q1(i0,j0,kc),q2(i0,j0,kc),q3(i0,j0,kc),                   &
       dens(i0,j0,kc),dsal((i0-1)*mref1+1,(j0-1)*mref2+1,kr),    &
       q1(i0,j1,kc),q2(i0,j1,kc),q3(i0,j1,kc),                   &
       dens(i0,j1,kc),dsal((i0-1)*mref1+1,(j1-1)*mref2+1,kr),    &
       q1(i1,j1,kc),q2(i1,j1,kc),q3(i1,j1,kc),                   &
       dens(i1,j1,kc),dsal((i1-1)*mref1+1,(j1-1)*mref2+1,kr),    &
       q1(i1,j0,kc),q2(i1,j0,kc),q3(i1,j0,kc),                   &
       dens(i1,j0,kc),dsal((i1-1)*mref1+1,(j0-1)*mref2+1,kr)
      ENDIF

      IF(myid.eq.idslab(11))then
       kc = kslab(11)
       kr = kslabr(11)
       write(811,'(1x,f10.4,20(1x,ES20.8))') time,               &
       q1(i0,j0,kc),q2(i0,j0,kc),q3(i0,j0,kc),                   &
       dens(i0,j0,kc),dsal((i0-1)*mref1+1,(j0-1)*mref2+1,kr),    &
       q1(i0,j1,kc),q2(i0,j1,kc),q3(i0,j1,kc),                   &
       dens(i0,j1,kc),dsal((i0-1)*mref1+1,(j1-1)*mref2+1,kr),    &
       q1(i1,j1,kc),q2(i1,j1,kc),q3(i1,j1,kc),                   &
       dens(i1,j1,kc),dsal((i1-1)*mref1+1,(j1-1)*mref2+1,kr),    &
       q1(i1,j0,kc),q2(i1,j0,kc),q3(i1,j0,kc),                   &
       dens(i1,j0,kc),dsal((i1-1)*mref1+1,(j0-1)*mref2+1,kr)
      ENDIF

      IF(myid.eq.idslab(12))then
       kc = kslab(12)
       kr = kslabr(12)
       write(812,'(1x,f10.4,20(1x,ES20.8))') time,              &
       q1(i0,j0,kc),q2(i0,j0,kc),q3(i0,j0,kc),                  &
       dens(i0,j0,kc),dsal((i0-1)*mref1+1,(j0-1)*mref2+1,kr),   &
       q1(i0,j1,kc),q2(i0,j1,kc),q3(i0,j1,kc),                  &
       dens(i0,j1,kc),dsal((i0-1)*mref1+1,(j1-1)*mref2+1,kr),   &
       q1(i1,j1,kc),q2(i1,j1,kc),q3(i1,j1,kc),                  &
       dens(i1,j1,kc),dsal((i1-1)*mref1+1,(j1-1)*mref2+1,kr),   &
       q1(i1,j0,kc),q2(i1,j0,kc),q3(i1,j0,kc),                  &
       dens(i1,j0,kc),dsal((i1-1)*mref1+1,(j0-1)*mref2+1,kr)
      ENDIF

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end subroutine probe_rcd

      subroutine probe_open
      use param
      use mpih
      use slab_param
      implicit none

      IF(myid.eq.idslab(4))then
        open(804, file='data/probe_z04.out',           &
              status='unknown',                        &
              access='sequential',                     &
              position='append')
        if(ireset.eq.1 .or. nread.eq.0)rewind(804)
      ENDIF

      IF(myid.eq.idslab(5))then
        open(805, file='data/probe_z05.out',          &
              status='unknown',                       &
              access='sequential',                    &
              position='append')
        if(ireset.eq.1 .or. nread.eq.0)rewind(805)
      ENDIF

      IF(myid.eq.idslab(8))then
        open(808, file='data/probe_z08.out',         &
              status='unknown',                      &
              access='sequential',                   &
              position='append')
        if(ireset.eq.1 .or. nread.eq.0)rewind(808)
      ENDIF

      IF(myid.eq.idslab(11))then
        open(811, file='data/probe_z11.out',        &
              status='unknown',                     &
              access='sequential',                  &
              position='append')
        if(ireset.eq.1 .or. nread.eq.0)rewind(811)
      ENDIF

      IF(myid.eq.idslab(12))then
        open(812, file='data/probe_z12.out',        &
              status='unknown',                     &
              access='sequential',                  &
              position='append')
        if(ireset.eq.1 .or. nread.eq.0)rewind(812)
      ENDIF

      return
      end subroutine probe_open

      subroutine probe_close
      use mpih
      use slab_param 
      implicit none

      IF(myid.eq.idslab(4))close(804)
      IF(myid.eq.idslab(5))close(805)
      IF(myid.eq.idslab(8))close(808)
      IF(myid.eq.idslab(11))close(811)
      IF(myid.eq.idslab(12))close(812)

      return
      end subroutine probe_close

