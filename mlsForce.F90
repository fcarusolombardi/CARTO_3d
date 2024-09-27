!-----------------------------------------------------------------------
!     read in position of a individual trial marker for MLS and compute 
!     compute support domain, shape function and interpolate
!------------------------------------------------------------------------

      subroutine mlsForce
      USE mpih
      USE param
      USE mls_param
      USE local_arrays, only: q1,q2,q3
      USE mpi_param, only: kstart, kend
      USE mls_local, only: for_xc, for_yc, for_zc, coll
      USE ieee_arithmetic
      IMPLICIT NONE

      integer :: inp,seed,merr,ntr
      real(DP) :: pos_MLS(3)
      integer :: pind_i(3),pind_o(3)

      integer :: inw,i,j,k,ii,jj,kk,cmp,ind_pal
      real(DP) :: norp(3),el_cx(3)
      real(DP) :: elmag,norpd,normd,epsw
      real(DP) :: pinvA(4,4),invA(4,4),B(4,nel),pxk(4,1)
      real(DP) :: ui(nel,3),Wt(3,nel),Wtx(nel),Hbox(3)
      real(DP) :: ptx(1,4),Bu(4,1),ABu(4,1),pABu(1,1)
      real(DP) :: ptxA(1,4),ptxAB(1,nel),ptxABu(1,1)
!     --------------------------------------------------
      real(DP) :: Um1,Um2,Um3
      real(DP) :: tcelvol,cfac,cel_fx,cel_fy,cel_fz
      real(DP) :: Forx,Fory,Forz
      real(DP) :: blend

      integer :: coll_dum,ci,cj,ck
!     --------------------------------------------------------
      mvol(:)=0.0
      dum_for(:) = -1
!     Initialise collision array and allocate coll values for walls
!     -1 to activate ; 0 for free; inp for single body occupancy

      coll(:,:,:) = 0
      !if(myid.eq.0)coll(:,:,4) = 0
      !if(myid.eq.numtasks-1)coll(:,:,n3m-4) = 0
      coll_dum = -1
      
      
      do ntr=1,nftot
         inp=face_to_part(ntr)
      if (body_nsl(inp) .eq. 1) then
         blend = 1.0D0
      elseif (body_nsl(inp) .eq. 0) then
        blend = 0.0D0
      else
        blend = 0.0D0
      endif

      ind_pal = pind(6,ntr)

      ci=pind(1,ntr); cj=pind(2,ntr); ck=pind(3,ntr)
      coll_dum = inp

      pind_i(:)=-1 ; pind_o(:)=-1
      
      if(ind_pal.ge.kstart-1.and.ind_pal.le.kend-1)then
      dum_for(ntr)=0
      
        pos_MLS(1)=tri_bar(1,ntr)
        pos_MLS(2)=tri_bar(2,ntr)
        pos_MLS(3)=tri_bar(3,ntr)



!     check for collision from brother face
      inp = face_to_part(ntr)

      if(inp.ge.3.and.inp.le.7) then !only for bodies between 3 and 7
      if(coll(ci,cj,ck).eq.0) coll(ci,cj,ck)=inp
      end if

!     initialise pre-factor matrix

      ptx(1,1)=1.d0;ptx(1,2)=pos_MLS(1);ptx(1,3)=pos_MLS(2)
      ptx(1,4)=pos_MLS(3)

      pind_i(1)=pind(1,ntr)-1;pind_o(1)=pind(1,ntr)+1
      pind_i(2)=pind(2,ntr)-1;pind_o(2)=pind(2,ntr)+1
      pind_i(3)=pind(3,ntr)-1;pind_o(3)=pind(3,ntr)+1

!     --------------------------------------------------------

      Hbox(1)=dx1
      Hbox(2)=dx2
      Hbox(3)=dx3

!     --------------------------------------------------------

      inw=1
      epsw = 1e-8

!     spline weights for all three components
!     and initialise ui matrix
      ui(:,:) = 0.0d0
      !dismax(ntr,inp) = dismax(ntr,inp)+epsw 
      
      do k=pind_i(3),pind_o(3)
       do j=pind_i(2),pind_o(2)
        do i=pind_i(1),pind_o(1)
        
          el_cx(1)=xm(i)-pos_MLS(1)+epsw
          el_cx(2)=ym(j)-pos_MLS(2)+epsw
          el_cx(3)=zm(k)-pos_MLS(3)+epsw

          elmag=sqrt(el_cx(1)**2+el_cx(2)**2+el_cx(3)**2)

          do cmp=1,3
            norp(cmp) = abs(el_cx(cmp))*Hbox(cmp)/wscl
!     ----------------EXPONENTIAL SPLINES------------------
            if(wexp.eq.1)then 
              if(norp(cmp).le.1.0d0)then
                Wt(cmp,inw)=exp(-(norp(cmp)/wcon)**2)
              else
                Wt(cmp,inw)=0.d0
                write(*,*)'Something wrong in support domain-mlsForce'
            write(*,*)cmp,norp(cmp),el_cx(:)
            write(*,*)'posmls', pos_MLS(:)
            write(*,*)'indx ',ntr
            stop
              end if
            end if
!     ----------------CUBIC SPLINES---------------------
            if(wcub.eq.1)then
              if(norp(cmp).le.0.5d0)then
                Wt(cmp,inw)=(2.d0/3.d0)-4.d0*(norp(cmp)**2)+ &
                                   4.d0*(norp(cmp)**3)
              elseif(norp(cmp).le.1.d0)then
                Wt(cmp,inw)=(4.d0/3.d0)*(1.d0-norp(cmp)**3)- &
                             4.d0*(norp(cmp)-norp(cmp)**2)
              else
                Wt(cmp,inw)=0.d0
              end if
            end if
!     ------------------------------------------------
          end do

          Wtx(inw) = Wt(1,inw)*Wt(2,inw)*Wt(3,inw)

!     -------------------------------------------------

          ui(inw,1) = 0.5*(q1(i,j,k)+q1(i+1,j,k))
          ui(inw,2) = 0.5*(q2(i,j,k)+q2(i,j+1,k))
          ui(inw,3) = 0.5*(q3(i,j,k)+q3(i,j,k+1))
          
          inw = inw + 1
         end do
        end do
       end do

!     --------------------------------------------------------
      
      B(:,:)=0.0d0; pinvA(:,:)=0.0d0; invA(:,:)=0.0d0
      inw = 1

!     pre-inverse matrix A(x) and B(x)
      do k=pind_i(3),pind_o(3)      
       do j=pind_i(2),pind_o(2)
        do i=pind_i(1),pind_o(1)

        pxk(1,1)=1.0d0;pxk(2,1)=xm(i)
        pxk(3,1)=ym(j);pxk(4,1)=zm(k)

        call DGEMM('N','T',4,4,1,Wtx(inw),pxk,4,pxk,4, &
                                1.0d0,pinvA,4)

!     -------------------------------------------------------------------

        B(:,inw)=Wtx(inw)*pxk(:,1)

!     -----------------------------------------------------------------
        inw = inw + 1
  
        end do
       end do
      end do

!     calling routine to compute inverse
      call inverseLU(pinvA,invA)
            
!     ------------------------------------------------------
!     matrix multiplications for final interpolation
!     DGEMM(transA,transB,m,n,k,alpha,A,LDA,b,LDB,beta,c,LDC)
!     C = alpha * A * B + beta * C

!     ---------------Shape function calculation---------------
      call DGEMM('N','N',1,4,4,1.0d0,ptx,1,invA,4,0.0d0,ptxA,1) 
      call DGEMM('N','N',1,nel,4,1.0d0,ptxA,1,B,4,0.0d0,ptxAB,1) 

!     --------------Velocity Interpolation from \phi---------
      do cmp=1,3
      
      call DGEMM('N','N',1,1,nel,1.0d0,ptxAB,1,ui(:,cmp),nel,0.0d0, &
                                                     ptxABu,1) 

      if(wcheck.eq.1)then
        if(sum(ptxAB).le.0.95.or.sum(ptxAB).gt.1.05)then
        write(*,*)'Shape function sum mlsForce',sum(ptxAB),pos_MLS(1:3)
        call MPI_ABORT(MPI_COMM_WORLD,ierr)
        end if
        if(ieee_is_nan(sum(ptxAB))) then
        write(*,*)'NaN Shape function mlsForce',pos_MLS(1:3)
        call MPI_ABORT(MPI_COMM_WORLD,ierr)
        end if
      end if

      select CASE(cmp)
            CASE(1)
              Um1 = ptxABu(1,1)
            CASE(2)
              Um2 = ptxABu(1,1)
            CASE(3)
              Um3 = ptxABu(1,1)
      end select


      end do
!     ------------------------------------------------------

!     For moving particles
      
      Forx = blend * (vel_tri(1,ntr)-Um1)
      Fory = blend * (vel_tri(2,ntr)-Um2)
      Forz = blend * (vel_tri(3,ntr)-Um3)

!     For stationary particles
!      Forx = -Um1
!      Fory = -Um2
!      Forz = -Um3

!!     For zero forcing
!      Forx = 0.0
!      Fory = 0.0
!      Forz = 0.0

!     -------------FORCING FUNCTION------------------------
!     volume of a face with a specific marker - thickness taken as average
!     of grid spacing

      inw = 1 ; tcelvol = 0.0
      do k=pind_i(3),pind_o(3)
       do j=pind_i(2),pind_o(2)
        do i=pind_i(1),pind_o(1)

        mvol(ntr)=mvol(ntr)+ptxAB(1,inw)*Hboxx(k)
        tcelvol = tcelvol+ptxAB(1,inw)*celvol(k)
        inw = inw+1

        end do
       end do
      end do

      cfac = (1.0-0.75*float(imlsref))*(mvol(ntr)* &
                                   sur(ntr))/tcelvol
!     ------------------------------------------------------------------
!     Resolutions need to be such that cfac is close to 1
!     cfac is ratio of volume of face to the volume of cell it resides in 
!     For test case - can be artificially set to 1

      cel_fx = cfac*Forx
      cel_fy = cfac*Fory
      cel_fz = cfac*Forz

!     ------------------------------------------------------

      inw = 1
      do k=pind_i(3),pind_o(3)
       do j=pind_i(2),pind_o(2)
        do i=pind_i(1),pind_o(1)

            for_xc(i,j,k)=for_xc(i,j,k)+ptxAB(1,inw)*cel_fx
            for_yc(i,j,k)=for_yc(i,j,k)+ptxAB(1,inw)*cel_fy
            for_zc(i,j,k)=for_zc(i,j,k)+ptxAB(1,inw)*cel_fz

            inw = inw+1
        end do
       end do
      end do 

!     -------------------------------------------------------


      end if

       end do
!     --------------------------------------------------------
 
      return
      end
