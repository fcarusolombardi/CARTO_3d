        subroutine SetTiling
        use mpih
        use tile_arrays
        use param, only : n1m,n2m,n3m,alx3,rext1,rext2,tilemax
        use mls_param

        implicit none

        integer :: ii,l,v1,v2,v3, Nreq, Nreqmax, Nava, qualeriga, rigamax
        integer :: ist, iet,inp,tilemaxeff
        real R0
        integer,dimension(tilemax,4) :: TABELLA
        real,    dimension(:,:), allocatable :: alphaBarco, betaBarco, gammaBarco


        Lgrid = ((rext1/dble(n1m))*(rext2/dble(n2m))*(alx3/dble(n3m)))**(1./3.)
        if(myid.eq.0) write(*,*) "Lgrid = ", Lgrid
!Read TABELLA
        open(unit=17,file='TABELLA.in',status='old')
        do ii=1,tilemax
        read(17,*)TABELLA(ii,1:4)
        enddo
        close(17)

!Find baricentre coordinate coefficients
        rigamax = TABELLA(tilemax,3)
        tilemaxeff = TABELLA(tilemax,2) !new
        allocate( alphaBarco(rigamax,tilemaxeff),betaBarco(rigamax,tilemaxeff),gammaBarco(rigamax,tilemaxeff))
        open(unit=17,file='DBbarco.in',status='old')
        do ii=1,rigamax
        read(17,*)alphaBarco(ii,1:tilemaxeff)
        enddo
        do ii=1,rigamax
        read(17,*)betaBarco(ii,1:tilemaxeff)
        enddo
        do ii=1,rigamax
        read(17,*)gammaBarco(ii,1:tilemaxeff)
        enddo
        close(17)


!Check what is needed
        allocate(tri_tiling(nftot,2))
        Nreqmax = -1
        do l=1,nftot
           R0   = sqrt(sur0(l))
           Nreq = ceiling((R0/(0.7*Lgrid))**2)

#ifdef SCALINGTEST
           Nava=Nreq
           qualeriga= nint(sqrt(real(Nava)))
#else
           if (Nreq.GT.tilemax) then
              write(*,*) "R0, Nreq",R0, Nreq
              write(*,*) "ERROR: Nreq > tilemax",l,face_to_part(l)
              stop
           endif
           Nava      = TABELLA(Nreq,2)
           qualeriga = TABELLA(Nreq,3)
#endif
           tri_tiling(l,1) = Nava  
           tri_tiling(l,2) = qualeriga

           Nreqmax = max(Nreqmax,Nreq)
        enddo
#ifdef SCALINGTEST
        Navamax = Nreqmax
        rigamax = nint(sqrt(real(Nreqmax)))
#else
        Navamax = TABELLA(Nreqmax,2)
        rigamax = TABELLA(Nreqmax,3)
#endif
        if(myid.eq.0) write(*,*) tilemax, Nreq, rigamax
        
        if(myid.eq.0) write(*,*) "Navamax1", Navamax

!allocate      
        allocate( albegaBar(rigamax,Navamax,3) )        
#ifdef SCALINGTEST
        albegaBar(1:rigamax,1:Navamax,1) = 0.333D0!alphaBarco(1:rigamax,1:Navamax)
        albegaBar(1:rigamax,1:Navamax,2) = 0.333D0!betaBarco(1:rigamax,1:Navamax)
        albegaBar(1:rigamax,1:Navamax,3) = 0.333D0!gammaBarco(1:rigamax,1:Navamax) !PJ
#else
        albegaBar(1:rigamax,1:Navamax,1) = alphaBarco(1:rigamax,1:Navamax)
        albegaBar(1:rigamax,1:Navamax,2) = betaBarco(1:rigamax,1:Navamax)
        albegaBar(1:rigamax,1:Navamax,3) = gammaBarco(1:rigamax,1:Navamax)
#endif

        nttot = 0

        allocate(ntilei(Nparticle),tilestart(Nparticle),tileend(Nparticle))
        ntilei(:)=0.d0
        do l=1,nftot
            nttot = nttot + tri_tiling(l,1)
            inp = face_to_part(l)
            ntilei(inp) = ntilei(inp) + tri_tiling(l,1)
        end do
        tilestart(1) = 1
        tileend(1) = ntilei(1)
        do inp=2,Nparticle
            tilestart(inp) = tileend(inp-1)+1
            tileend(inp)   = tilestart(inp) + ntilei(inp)-1
        enddo

        allocate(faceid_t(nttot)) ! face id of a specific tile      
        allocate(pindt(6,nttot))
        allocate(tstart(nftot))
        allocate(tend(nftot))

        allocate(tlcnt(1))
        allocate(mytl(nttot))

        iet = 0
        do l=1,nftot
            ist = iet + 1
            iet = iet + tri_tiling(l,1)
            faceid_t(ist:iet) = l
            tstart(l)=ist
            tend(l)=iet
        end do

#ifdef OFFSETBODY
      dOff=2.D0*Lgrid
      inpstartOff = 1
      inpendOff = 1
#ifdef SOLOUNO
      inpendOff = 1
#endif
      ntOffstart = tilestart(inpstartOff)
      ntOffend   = tileend(inpendOff)
      ntOfftot   = ntOffend-ntOffstart+1
      allocate(pindtOff(6,nttot))
#endif
      deallocate(alphaBarco,betaBarco,gammaBarco)

        return
        end subroutine SetTiling


