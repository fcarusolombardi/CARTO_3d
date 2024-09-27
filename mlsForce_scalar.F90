!-----------------------------------------------------------------------
!     read in position of a individual trial marker for MLS and compute 
!     compute support domain, shape function and interpolate
!------------------------------------------------------------------------

      subroutine mlsForce_scalar
      USE mpih
      USE param
      USE mls_param
      USE local_arrays, only: dsal
      USE mpi_param, only: kstart,kend,kstartr,kendr
      USE mls_local, only: for_xc, for_yc, for_zc, coll, for_sc
      USE ieee_arithmetic
      IMPLICIT NONE

      integer :: inp,seed,merr,ntr
      real(DP) :: pos_MLS(3)
      integer :: pind_i(3),pind_o(3)

      integer :: inw,i,j,k,ii,jj,kk,cmp,ind_pal
      real(DP) :: norp(3),el_cx(3)
      real(DP) :: elmag,norpd,normd,epsw
      real(DP) :: pinvA(4,4),invA(4,4),B(4,nel),pxk(4,1)
      real(DP) :: sci(nel),Wt(3,nel),Wtx(nel),Hbox(3)
      real(DP) :: ptx(1,4),Bu(4,1),ABu(4,1),pABu(1,1)
      real(DP) :: ptxA(1,4),ptxAB(1,nel),ptxABu(1,1)
!     --------------------------------------------------
      real(DP) :: tcelvol,cfac,cel_f
      real(DP) :: for_sc_local,sc_int
!     --------------------------------------------------------
      mvol(:,:)=0.0
      
      do inp=1,Nparticle
      
      do ntr=1,maxnf

      ind_pal = pindr(6,ntr,inp)

      pind_i(:)=-1 ; pind_o(:)=-1
      
      
      if(ind_pal.ge.kstartr-1.and.ind_pal.le.kendr-1)then

        pos_MLS(1)=tri_bar(1,ntr,inp)
        pos_MLS(2)=tri_bar(2,ntr,inp)
        pos_MLS(3)=tri_bar(3,ntr,inp)

!     correct for boundaries
      if(pos_MLS(1).ge.xc(n1))pos_MLS(1)=pos_MLS(1)-xc(n1)
      if(pos_MLS(1).lt.xc(1)) pos_MLS(1)=pos_MLS(1)+xc(n1)

!     initialise pre-factor matrix

      ptx(1,1)=1.d0;ptx(1,2)=pos_MLS(1);ptx(1,3)=pos_MLS(2)
      ptx(1,4)=pos_MLS(3)

      if(nel.eq.27)then

      pind_i(1)=pindr(1,ntr,inp)-1;pind_o(1)=pindr(1,ntr,inp)+1
      pind_i(2)=pindr(2,ntr,inp)-1;pind_o(2)=pindr(2,ntr,inp)+1
      pind_i(3)=pindr(3,ntr,inp)-1;pind_o(3)=pindr(3,ntr,inp)+1

      elseif(nel.eq.125)then

      pind_i(1)=pindr(1,ntr,inp)-2;pind_o(1)=pindr(1,ntr,inp)+2
      pind_i(2)=pindr(2,ntr,inp)-2;pind_o(2)=pindr(2,ntr,inp)+2
      pind_i(3)=pindr(3,ntr,inp)-2;pind_o(3)=pindr(3,ntr,inp)+2

      end if
!     --------------------------------------------------------

      Hbox(1)=dx1
      Hbox(2)=dx2
      if(istr3.eq.6)then
      Hbox(3)=udx3m(pind(inp,ntr,3))
      elseif(istr3.eq.0)then
      Hbox(3)=dx3
      end if

!     --------------------------------------------------------

      inw=1
      epsw = 1e-10

!     spline weights for all three components
!     and initialise ui matrix
      sci(:) = 0.0d0
      dismax(ntr,inp) = dismax(ntr,inp)+epsw 
      
      do k=pind_i(3),pind_o(3)
       do j=pind_i(2),pind_o(2)
        do i=pind_i(1),pind_o(1)
        
          el_cx(1)=xmr(i)-pos_MLS(1)
          el_cx(2)=ymr(j)-pos_MLS(2)
          el_cx(3)=zmr(k)-pos_MLS(3)

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

      if(i.gt.n1m) then
            ii = i-n1m
      else
            ii = i
      end if

      if(ii.lt.1)then
            ii = n1m
      else
            ii = ii
      end if 

          sci(inw) = dsal(ii,j,k)
          
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

        pxk(1,1)=1.0d0;pxk(2,1)=xmr(i)
        pxk(3,1)=ymr(j);pxk(4,1)=zmr(k)

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
      call invert4(pinvA,invA)
            
!     ------------------------------------------------------
!     matrix multiplications for final interpolation
!     DGEMM(transA,transB,m,n,k,alpha,A,LDA,b,LDB,beta,c,LDC)
!     C = alpha * A * B + beta * C

!     ---------------Shape function calculation---------------
      call DGEMM('N','N',1,4,4,1.0d0,ptx,1,invA,4,0.0d0,ptxA,1) 
      call DGEMM('N','N',1,nel,4,1.0d0,ptxA,1,B,4,0.0d0,ptxAB,1) 

!     --------------Scalar Interpolation from \phi---------
      
      call DGEMM('N','N',1,1,nel,1.0d0,ptxAB,1,sci(:),nel,0.0d0, &
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

              sc_int = ptxABu(1,1)
!     ------------------------------------------------------

!     For moving particles
      
      for_sc_local = 1.0-sc_int

!     For stationary particles
!      Forx = -Um1
!      Fory = -Um2
!      Forz = -Um3

!     For zero forcing
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

        mvol(ntr,inp)=mvol(ntr,inp)+ptxAB(1,inw)*Hboxxr(k)
        tcelvol = tcelvol+ptxAB(1,inw)*celvolr(k)
        inw = inw+1

        end do
       end do
      end do

      cfac = (1.0-0.75*float(imlsref))*(mvol(ntr,inp)* &
                                   sur(ntr,inp))/tcelvol
!     ------------------------------------------------------------------
!     Resolutions need to be such that cfac is close to 1
!     cfac is ratio of volume of face to the volume of cell it resides in 
!     For test case - can be artificially set to 1

      cel_f = cfac*for_sc_local

!     ------------------------------------------------------

      inw = 1
      do k=pind_i(3),pind_o(3)
       do j=pind_i(2),pind_o(2)
        do i=pind_i(1),pind_o(1)

            if(i.gt.n1m) then
                  ii = i-n1m
            else
                  ii = i
            end if

            if(ii.lt.1)then
                  ii = n1m
            else
                  ii = ii
            end if 


            for_sc(ii,j,k)=for_sc(ii,j,k)+ptxAB(1,inw)*cel_f

            inw = inw+1
        end do
       end do
      end do 

!     -------------------------------------------------------


      end if

       end do
      end do
!     --------------------------------------------------------
 
      return
      end
