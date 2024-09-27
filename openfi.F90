      subroutine openfi
      use param
      use mpih
      implicit none

      IF(myid.eq.0)then

!   volume integration of heat flux
        open(95, file='pressure.out', &
              status='unknown',           &
              access='sequential',        &
              position='append')

#ifdef ELECTRO
!   volume integration of heat flux
        open(195, file='electroEF.out', &
              status='unknown',           &
              access='sequential',        &
              position='append')
!   ECG
        open(196, file='ECG.out', &
              status='unknown',           &
              access='sequential',        &
              position='append')

!   ECG_Atr
        open(197, file='ECG_Atr.out', &
              status='unknown',           &
              access='sequential',        &
              position='append')

!   ECG_Ventr
        open(198, file='ECG_Ventr.out', &
              status='unknown',           &
              access='sequential',        &
              position='append')

!   ECG_2D
        open(199, file='ECG_2D.out', &
              status='unknown',           &
              access='sequential',        &
              position='append')
! EGMs CARTO
        open(200, file='EGM_CARTO.out', &
             !  !status='unknown',           &
             !  access='stream',        &
             !  position='append',      &
             ! ! status='unknown',          &
             !  form='unformatted')
             ! ! access='sequential',           &
             ! ! position='append')
              status='unknown',           &
              access='sequential',        &
              position='append')

        
#endif

!   average over top and bottom wall
      open(97, file="nusse_walls.out", &
            status='unknown',         &
            access='sequential',      &
            position='append')

!   rms_vel.out in stst.F
      open(94, file='rms_vel.out', &
            status='unknown',     &
            access='sequential',  &
            position='append')

!   exact relations of GL 
      open(92, file='nusse_balance.out', &
            status='unknown',           &
            access='sequential',        &
            position='append')

!   balance of base and refined mesh
      open(98, file='kediss_bothmesh.out', &
            status='unknown',             &
            access='sequential',          &
            position='append')

!   total kinetic energy
        open(96, file='total_KE.out', &
              status='unknown',      &
              access='sequential',   &
              position='append')

! reset the time history
        if(ireset.eq.1 .or. nread.eq.0)then
          rewind(95)
          rewind(97)
          rewind(94)
          rewind(92)
          rewind(98)
          rewind(96)
        endif


      ENDIF

      return
      end   
      
!==============================================

      subroutine closefi
      use mpih
      implicit none
      
      if(myid.eq.0)then

!    nusse_volume.out in vmaxv.F
      close(95)
!    nusse_wall.out  in densmc.F
      close(97)
!    rms_vel.out in stst.F
      close(94)
!    nusse_balance.out in balance.F
      close(92)
!    kediss_bothmesh.out
      close(98)
!    total kinetic energy
        close(96)


      endif
     
      return      
      end

