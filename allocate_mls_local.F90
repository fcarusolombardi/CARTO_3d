      subroutine allocate_mls_local
      use param
      use mls_local
      use mpih
      use outflow_vars
      use mpi_param, only: kstart, kend, kstartr, kendr
      use mls_param!, only: Hboxx,celvol,Hboxxr,celvolr
      implicit none 
      integer :: merr
      !-------------------------------------------------
      allocate(for_xc(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo-1),stat=merr)
      
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory: forxc"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      !-------------------------------------------------
      allocate(for_yc(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo-1),stat=merr)

                                                                 
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory: foryc"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      !-------------------------------------------------
      allocate(for_zc(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo-1),stat=merr)
      
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory: forzc"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      !-------------------------------------------------
      allocate(for_te(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo-1),stat=merr)
                                                                 
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory: forte"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 

      !-------------------------------------------------
      allocate(coll(1:n1,1:n2,kstartr-1:kendr),stat=merr)
      
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory: coll"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
     endif
      !-------------------------------------------------
#ifdef CONTACTV3
      allocate(coll3(1:n1,1:n2,kstartr-1:kendr,3),stat=merr)
      
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory: coll3"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
     endif
#endif     
      !-------------------------------------------------
      allocate(for_sc(1:n1r,1:n2r,kstartr-lvlhalo:kendr+lvlhalo-1),stat=merr)
      
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory: forsc"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 

!     Allocate outlfow variables
      allocate(qb1s(1:n1,1:n2),qb1n(1:n1,1:n2))
      allocate(qb2s(1:n1,1:n2),qb2n(1:n1,1:n2))
      allocate(qb3s(1:n1,1:n2),qb3n(1:n1,1:n2))
#if defined(STRONG) || defined(STRONGRV)
      allocate(qb1sg(1:n1,1:n2),qb1ng(1:n1,1:n2))
      allocate(qb2sg(1:n1,1:n2),qb2ng(1:n1,1:n2))
      allocate(qb3sg(1:n1,1:n2),qb3ng(1:n1,1:n2))
#endif
      
      allocate(dqb1s(1:n1,1:n2),dqb1n(1:n1,1:n2))
      allocate(dqb2s(1:n1,1:n2),dqb2n(1:n1,1:n2))
      allocate(dqb3s(1:n1,1:n2),dqb3n(1:n1,1:n2))

      allocate(dq1x2o(1:m1,1:m2),dq2x2o(1:m1,1:m2),dq3x2o(1:m1,1:m2))
#if defined(STRONG) || defined(STRONGRV)
      allocate(dq1x2og(1:n1,1:n2),dq2x2og(1:n1,1:n2),dq3x2og(1:n1,1:n2))
#endif      
      !     spostate qui da allocate_trigeo
      allocate(Hboxx(0:n3),celvol(0:n3))
      allocate(Hboxxr(0:n3r),celvolr(0:n3r))

      !moving probe
      allocate(tpcnt(1))      
      allocate(mypr(Nmlsprobe))
      allocate(xyzMLSprobe(3,Nmlsprobe), outMLSprobe(4,Nmlsprobe))
      allocate(pindMLSprobe(6,Nmlsprobe))
      
      return
      end 
