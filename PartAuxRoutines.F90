!-----------------------------------------------------
!     Compute mean position for all particles
!-----------------------------------------------------

      subroutine meanpos
      USE mls_param

      IMPLICIT NONE

      integer inp


      do inp = 1,Nparticle
       
        meanxyz(1,inp) = sum(xyz(1,:,inp))/float(maxnv)
        meanxyz(2,inp) = sum(xyz(2,:,inp))/float(maxnv)
        meanxyz(3,inp) = sum(xyz(3,:,inp))/float(maxnv)

        meanvel(1,inp) = sum(xyzv(1,:,inp))/float(maxnv)
        meanvel(2,inp) = sum(xyzv(2,:,inp))/float(maxnv)
        meanvel(3,inp) = sum(xyzv(3,:,inp))/float(maxnv)

      end do

      return
      end

!-----------------------------------------------------
!     Reset the particle shape to its initial shape
!-----------------------------------------------------

      subroutine partshape_reset
      USE mls_param

      integer dummy,i,inp,inp_reset
      character*50 filename
      real(DP) init_pos(3,maxnv)

      filename = gtsfx

!     read the position of vertices from sphere gts
      open(109,file=filename)
      read(109,*)dummy,dummy,dummy
      do i=1,maxnv
        read(109,*)init_pos(1,i),init_pos(2,i),init_pos(3,i)
      end do
      close(109)

!     Check which particle has exploded
      do inp = 1,Nparticle
      inp_reset = 0
       do ntr = 1,maxnv
         if(abs(maxval(xyzv(:,ntr,inp))).gt.3.0)then
           inp_reset = inp
           if(myid.eq.0)write(*,*)'Reseting particle ',inp
         end if
       end do
      end do

!     translate the sphere to its reset position
      if(inp_reset.gt.0)then      
       do i=1,maxnv
        xyz(1,i,inp_reset)=init_pos(1,i)*sclf+meanxyz(1,inp_reset)
        xyz(2,i,inp_reset)=init_pos(2,i)*sclf+meanxyz(2,inp_reset)
        xyz(3,i,inp_reset)=init_pos(3,i)*sclf+meanxyz(3,inp_reset)

        xyzv(:,i,inp_reset) = 0.0
        xyzv0(:,i,inp_reset) = 0.0
        xyza(:,i,inp_reset) = 0.0
        xyza0(:,i,inp_reset) = 0.0
       end do
      end if

      return
      end


