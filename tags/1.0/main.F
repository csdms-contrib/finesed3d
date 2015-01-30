      PROGRAM main
#include "common.inc"
c
      INTEGER ch_fin,m_fin,iters,itfirst,iend,i,j,k
      INTEGER tt_yes,twc_yes,otau_yes,nbsq_yes,md_rhs_yes,ngbal,n1d,code
      REAL*8  dt0,dt_p,KE,Ep,diss,ddtEp
      REAL*8  espz(nz) 
      CHARACTER*32 identifier
c
      REAL*8 espx(nxh,nz),dspx(nxh,nz),espy(nyh,nz),dspy(nyh,nz)
      REAL*8 fmean(nz),xrms(nz),yrms(nz),pmean(nz),prms(nz)
      REAL*8 tfmean,txrms,tyrms,tpmean,tprms,igmod
      REAL*8 auxdepo(nx0,ny0),auxeros(nx0,ny0)
      REAL*8 ushearx(nx0,ny0),usheary(nx0,ny0)
      logical :: exist
c
      print *, '*************************** FineSed3D **************
     &****************'
      print *, 'A turbulence-resolving numerical model'
      print *, 'for fine sediment transport in bottom boundary layer'
      print *, ''
      print *, 'Developed by (1)Zhen Cheng, (1)Xiao Yu, (2)Celalettin 
     &E. Ozdemir,'
      print *, '(1)Tian-Jian Hsu, (3)S. Balachandar'
      print *, ''
      print *, '(1)Center for Applied Coastal Research, University 
     &of Delaware'
      print *, '(2)Civil and Environmental Engineering, Louisiana 
     &State University'
      print *, '(3)Mechanical and Aerospace engineering, University 
     &of Florida'
      print *, 'For questions, please contact: zcheng@udel.edu (Zhen 
     &Cheng)'
      print *, '*****************************************************
     &****************'
      print *, 'Program running'
c
c     read in the parameters for the time-stepping and the filenames
c     ch_fin = 0 generates an initial solution
c     ch_fin = 1 reads the initial data from ch_file.(istart)
c
c     outfirst = whether to output data (other than c and u) at
c     first time-step
c
      read '(a)', identifier
      if (identifier(1:5).ne.'#MAIN') then
        print *, 'Expecting identifier #MAIN in input file'
        stop
      endif
      read '(a)', ch_file
      read(*,*) istart,iters,ch_fin,m_fin
      read(*,*) oi_chan,out_press,limit
      read(*,*) oi_timer,oi_cfl,oi_spec,oi_mean,oi_gbal,oi_1d
      read(*,*) re,gradpbar,dt,cflmax
      read(*,*) tt_yes,twc_yes, otau_yes
c      read(*,*) prop

c      re=redel*30.0
c
c      dt = dt_p/re
     
      dt = 4.d0*atan(1.d0)/2880.d0
cz make it more stable
cz      dt = dt/4.d0
      FLTIMER = oi_timer.gt.0
      FLAVER  = oi_mean.gt.0
      FLTHRM  = tt_yes.ne.0
      FLTWC   = twc_yes.eq.1
      FLOTAU  = otau_yes.eq.1
      print *,'1 chnl input read'
c
      if (FLTIMER) then
        call init_timers
      else
        oi_timer = 1
      endif
c
      read '(a)', identifier
      print *, identifier
      if (identifier(1:9).ne.'#END_MAIN') then
        print *, 'Expecting identifier #END_MAIN in input file'
        stop
      endif
c
      if (FLTHRM) then
         read '(a)', identifier
         print *, identifier
         if (identifier(1:5).ne.'#THRM') then
            print *, 'Expecting identifier #THRM in input file'
            stop
         endif
         read '(a)', tt_file
         read(*,*) sc,ri,Tmax,RR,nbsq_yes,md_rhs_yes
         read(*,*) ttab(1),ttat(1) 
         read(*,*) ttbb(1),ttbt(1)
         read(*,*) ttgb(1),ttgt(1)
         read(*,*) w_settling,tau_p,beta_p,rey_p,phiave,Ero,ushearc
         read '(a)', identifier
         print *, identifier
         if (identifier(1:9).ne.'#END_THRM') then
            print *, 'Expecting identifier #END_THRM in input file'
            stop
         endif
      endif
c
c      print *,'1 thermal input read'
      FLNBSQ    = nbsq_yes.eq.1
      FLSETTLE  = w_settling.ne.0.0
      FLINERTIA = tau_p.ne.0.0
      FLPART    = w_settling.ne.0.0.or.tau_p.ne.0.0
      FLMDRHS   = md_rhs_yes.eq.1
      FLHINDER = phiave.ne.0.0
      FLEro = Ero.ne.0.0
      phind = 0.d0
      if(FLHINDER)then
        print *, 'Hindered settling is included'
c  power in hindered settling
c  Richardson & Zaki (1954)
        if(rey_p.le.2d-1)then
        phind = 4.6d0
        else if(rey_p.gt.2.d-1.and.rey_p.le.1.d0)then
        phind = 4.4d0*rey_p**(-3.d-2)
        else if(rey_p.gt.1.d0.and.rey_p.le.5.d2)then
        phind = 4.4d0*rey_p**(-1.d-1)
        else
        phind = 2.39d0
        endif
      endif

      print *,'1 logical variables are assigned'
c
      read '(a)', identifier
      print *, identifier
      if (identifier(1:8).ne.'#SP_VISC') then
         print *, 'Expecting identifier #SP_VISC in input file'
         stop
      endif
      read(*,*) epssvx,epssvy
      print *, 'epssvx apssvy are read'
      read(*,*) msvx,msvy 
      print *, 'msvx msvy are read'
      read '(a)', identifier
      if (identifier(1:12).ne.'#END_SP_VISC') then
         print *, 'Expecting identifier #END_SP_VISC in input file'
         stop
      endif
      print *,'2 viscous input read'
c
      FLSV_NO=.false.
      FLSV_YES=epssvx.gt.0.or.epssvy.gt.0
c
      print *, 'Input read'
c
c make sure the bc coeff is right
      ttat(1)=re*sc*abs(w_settling*gz(1,1)/sqrt(gy(1,1)**2
     &     +gz(1,1)**2))
      ttbt(1)=1.d0
      ttgt(1)=0.d0
      ttab(1)=ttat(1)
      ttbb(1)=ttbt(1)
      ttgb(1)=ttgt(1)
c  done with the re-def of bc coeff

      itfirst = istart+1
      iend = istart + iters
      itime = istart
c
      print *,              'Flow file root:           ', ch_file
      if (FLTHRM)  print *, 'Temperature file root:    ', tt_file
c
      open(26,file='logfile')
      open(31,file='gbal')
      open(60,file='ucl')
cz to write u_shear
      inquire(file='ushear.dat',exist=exist)
      if (exist) then
      open(80,file='ushear.dat',Access = 'append',Status='old')
      else
      open(80,file='ushear.dat',action = 'write',Status='new')
      end if
c      print *,'3 before rewind'
      rewind(31)
c      print *,'4 after rewind'
      ngbal = 0
      open(30,file='1dout')
      rewind(30)
      n1d = 0
c
      call initial(dt0,ch_fin,m_fin)
c      print *,'5 initial subroutine is called'
c produce complete output for initial condition if ...
      if(ch_fin.eq.0)then
         call make_gbalrhs(KE,Ep,diss,ddtEp)
c         print *,'6 make_gbalrhs subroutine is called'
         write(31,200) rtime,KE,Ep,diss,ddtEp
         ngbal = ngbal + 1
c     
         call output
c         print *,'7 output subroutine is called'
c     
         if(FLAVER)then
         endif
      endif !(ch_fin.eq.0)
c     
      call divg
c      print *,'8 divg subroutine is called'
c
      write(26,*)' rlenx,rleny'
      write(26,*)  rlenx,rleny
      write(26,*)' deltax,deltay,dt0'
      write(26,*)  deltax,deltay,dt0
      write(26,*)' istart,iters,itfirst,iend'
      write(26,*)  istart,iters,itfirst,iend
      write(26,*)' oi_chan,out_press'
      write(26,*)  oi_chan,out_press
      write(26,*)' oi_timer,oi_cfl,oi_spec,oi_mean,oi_gbal,oi_1d'
      write(26,*)  oi_timer,oi_cfl,oi_spec,oi_mean,oi_gbal,oi_1d
      write(26,*)' nx, ny, nz, nxl, nyl'
      write(26,*)  nx, ny, nz, nxl, nyl
      write(26,*)' rlenx, rleny, dt'
      write(26,*)  rlenx, rleny, dt
      write(26,*)' Re,gradpbar,dt_p,cflmax'
      write(26,*)  re,gradpbar,dt_p,cflmax
      write(26,*)' tt_yes'
      write(26,*)  tt_yes
      if (FLTHRM) then
        write(26,*)' Sc, Ri, Tmax, RR, nbsq_yes, md_rhs_yes'
        write(26,*)  sc, ri, Tmax, RR, nbsq_yes, md_rhs_yes
        write(26,*)' BC: ttat(1), ttab(1), ttbt(1), ttbb(1)
     &, ttgt(1), ttgb(1)'
        write(26,*)      ttat(1), ttab(1), ttbt(1), ttbb(1)
     &, ttgt(1), ttgb(1)
        write(26,*)' w_settling, tau_p, beta_p, rey_p'
        write(26,*)  w_settling, tau_p, beta_p, rey_p 
      endif
      write(26,*)
c
c     transform variables to Fourier space
      call fft_r2f_2d_new(u)    !call fft_r2f_2d(u,ru)
      call fft_r2f_2d_new(v)    !call fft_r2f_2d(v,rv)
      call fft_r2f_2d_new(w)    !call fft_r2f_2d(w,rw)
      if (FLTHRM) call fft_r2f_2d_new(tt) !call fft_r2f_2d(tt,rtt)
c      print *,'9 fft_r2f_2d_new subroutine is called'
c
c     transform bc's (that were computed in initial)
cNOTE: to do time varying bc see subroutine advance
      call fft_r2f_2d_new_slice(ubct)
      call fft_r2f_2d_new_slice(ubcb)
      call fft_r2f_2d_new_slice(vbct)
      call fft_r2f_2d_new_slice(vbcb)
      call fft_r2f_2d_new_slice(wbct)
      call fft_r2f_2d_new_slice(wbcb)
      if (FLTHRM) call fft_r2f_2d_new_slice(ttbct)
      if (FLTHRM) call fft_r2f_2d_new_slice(ttbcb)
      if (FLINERTIA.or.FLHINDER) then
         call dcopy(nx0y0,tt(1,1,1),1,ttotn(1,1),1)
         call dcopy(nx0y0,tt(1,1,nz),1,ttobn(1,1),1)
         call dcopy(nx0y0,ttotn(1,1),1,ttoto(1,1),1)
         call dcopy(nx0y0,ttobn(1,1),1,ttobo(1,1),1)
      endif
c      print *,'10 fft_r2f_2d_new_slice subroutine is called'
c
cNOTE: DDts and old vel are not fft, they are always in Fourier space
c     
c     loop over the time steps
      do itime = itfirst,iend
c    
c     write the time step number to logfile
        if (mod((itime-istart),oi_cfl).eq.0) then
          write(26,*)' '
          write(26,*)' itime =  ',itime
          write(26,*)' '
        endif
c
c     perform one complete time step
        call advance(0)
        print *,'11 advance subroutine is called'
        rtime  = rtime+dt
c
c     compute mean values
        call make_means
c
c     output means
        if(FLAVER)then
c     keep these two "if" like this because oi_mean=0 -> FLAVER=.false.
           if(mod((itime-istart),oi_mean).eq.0)call output_means
        endif
c
c     compute deposit, must be done every time step        
        if(FLEro.and.w_settling.ne.0)then
           call comp_bott_shear(ushearx,usheary)
cz returned in fourier space, just output the mean
           write(80,104) rtime,ushearx(1,1),usheary(1,1)
c         if(.false..and.w_settling.ne.0)then
cz           call dcopy(nx0y0,tt(1,1,nz-1),1,auxdepo,1) !copy bottom concentration
c           print *,'12 dcopy subroutine is called'
cz           call fft_f2r_2d_new_slice(auxdepo) !bring it to real space
c           print *,'13 fft_f2r_2d_new_slice subroutine is called'
cz           call comp_erosion(auxeros) !erosion returns in real space
ccz   need to make sure of w_s, should include inertia effect and
ccz   hindered settling
cz           do j=1,ny
cz             do i=1,nx
cz                erosion(i,j)= erosion(i,j)
ccz note we are using different strategy, comment the following 1 line
ccz     &                       - abs(w_settling)*auxeros(i,j)*dt
cz     &                        - auxeros(i,j)*dt
c
cz                deposit(i,j)= deposit(i,j)
cz     &                       + abs(w_settling)*auxdepo(i,j)*dt
cz             enddo
cz           enddo
c           print *, erosion(nxh,nyh)
c           print *, deposit(nxh,nyh)
        endif
104   format(e15.8,2f12.5)
c
c     produce outputs in Fourier space
c     output spectra
        if((mod(itime-istart,oi_spec).eq.0.or.
     &      mod(itime-istart,oi_mean).eq.0)) then
c u 
           call make_spec_new(u,espx,dspx,espy,dspy,xrms,
     &          txrms,yrms,tyrms)
c     if(mod(itime,oi_mean).eq.0)call make_meanf(u,fmean,tfmean)
           if(mod(itime-istart,oi_spec).eq.0)
     &          call spec_out(espx,dspx,espy,dspy,'dspu',itime)
c v
           call make_spec_new(v,espx,dspx,espy,dspy,xrms,
     &          txrms,yrms,tyrms)
c     if (mod(itime,oi_mean).eq.0)call make_meanf(v,fmean,tfmean)
           if (mod(itime-istart,oi_spec).eq.0)
     &          call spec_out(espx,dspx,espy,dspy,'dspv',itime)
c w
           call make_spec_new(w,espx,dspx,espy,dspy,xrms,
     &          txrms,yrms,tyrms)
c     if (mod(itime,oi_mean).eq.0)call make_meanf(w,fmean,tfmean)
           if (mod(itime-istart,oi_spec).eq.0)
     &          call spec_out(espx,dspx,espy,dspy,'dspw',itime)
c tt
           if (FLTHRM) then
              call make_spec_new(tt,espx,dspx,espy,dspy,xrms,
     &             txrms,yrms,tyrms)
c     if (mod(itime,oi_mean).eq.0)call make_meanf(tt,fmean,tfmean)
              if (mod(itime-istart,oi_spec).eq.0)
     &             call spec_out(espx,dspx,espy,dspy,'dspT',itime)

           endif
        endif
c         print *,'14 make_spec_new subroutine is called'
c
c     produce outputs in real space
        if( mod((itime-istart),oi_chan).eq.0.or.
     &      mod((itime-istart),oi_cfl ).eq.0.or.
     &      mod((itime-istart),oi_mean).eq.0.or.
     &      mod((itime-istart),oi_spec).eq.0.or.
     &      mod((itime-istart),oi_gbal).eq.0.or.
     &      mod((itime-istart),oi_1d  ).eq.0)then
c
c     transform to real space
           call fft_f2r_2d_new(u)
           call fft_f2r_2d_new(v)
           call fft_f2r_2d_new(w)
           if (FLTHRM)call fft_f2r_2d_new(tt)
c NOTE: old vel are not fft, they are always in Fourier space
c
c     output a restart file
           if(itime.gt.limit.and.mod((itime-istart),oi_chan).eq.0) 
     &          call output
c
c     output pressure (press is in Fourier space in var p allways)
           if(oi_chan.gt.limit) then
           if(out_press.gt.0.and.mod((itime-istart),oi_chan).eq.0)
     &          call output_press(out_press)
           endif
c
c     compute the divergence and the cfl condition 
c     output the vars at (2,2,3)
           if (mod((itime-istart),oi_cfl).eq.0) then
              call divg
              print *,u(1,1,1),v(1,1,1),w(1,1,1)
              call courant(dt0,u,v,w)
c              print *,'15 courant subroutine is called'
              write(26,*)' '
c
              if(FLSETTLE)then
C     recall w_settling has its own sign (<0 for heavier than fluid)
                 igmod=sqrt(gy(1,1)**2+gz(1,1)**2)
                 igmod=1.d0/igmod
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
                 do k=1,nz
                    do j=1,ny
                       do i=1,nx
                          rv(i,j,k)=v(i,j,k)-w_settling*igmod*gy(1,1)
                          rw(i,j,k)=w(i,j,k)-w_settling*igmod*gz(1,1)
                       enddo
                    enddo
                 enddo                 
C$OMP END PARALLEL DO
                 write(26,*)' tt(nxh,nyh,nzm):',tt(nxh,nyh,nzm)
                 write(*,*)' tt(nxh,nyh,nzm):',tt(nxh,nyh,nzm)
                 write(26,*)'  Disperse phase info:'
                 write( *,*)'  Disperse phase info:'
                 call courant(dt0,u,rv,rw)
                 write(26,*)' '
              endif
           endif
c
c     output deposition and erosion
           if(FLEro.and.mod((itime-istart),oi_mean).eq.0)
     &       call output_2d_section('deposit',deposit,erosion,deposit,2)
c         print *, '16 output_2d_section subroutine is called'
c
c     output Chebyshev spectrum
           if (mod((itime-istart),oi_spec).eq.0)then
              call make_chb_spec(u,espz)
              call chb_spec_out(espz,'chspu',itime)
              call make_chb_spec(v,espz)
              call chb_spec_out(espz,'chspv',itime)
              call make_chb_spec(w,espz)
              call chb_spec_out(espz,'chspw',itime)
              if(FLTHRM)then
                 call make_chb_spec(tt,espz)
                 call chb_spec_out(espz,'chspT',itime)
              endif
           endif
c           print *,'17 chebyshev subroutines are called'
c
c     output 1d stats
           if (mod((itime-istart),oi_1d).eq.0) then
              call make_output_1d(30)
              n1d = n1d + 1
           endif
c     
c     output global balances
           if (mod((itime-istart),oi_gbal).eq.0) then
              call make_gbalrhs(KE,Ep,diss,ddtEp)
              write(31,200) rtime,KE,Ep,diss,ddtEp
              ngbal = ngbal + 1
           endif
c
c     transform back to Fourier space
           call fft_r2f_2d_new(u)
           call fft_r2f_2d_new(v)
           call fft_r2f_2d_new(w)
           if (FLTHRM) call fft_r2f_2d_new(tt)
c NOTE: old vel are not fft, they are always in Fourier space
c     
        endif!outputs in real space
c
c     close and reopen gbal and 1dout when main data written, 
c     so no data is lost during a crash and restart
        if (mod((itime-istart),oi_chan).eq.0) then
           close(30)
           open(30,file='1dout')
           rewind(30)
           do j=1,n1d
c check in make_output_1d # of variables to read
              read(30,*) diss,diss,diss,diss,diss
           enddo
           close(31)
           open(31,file='gbal')
           rewind(31)
           do j=1,ngbal
              read(31,*) diss,diss,diss,diss,diss
           enddo
        endif
c
c     if the time step is divisible by oi_timer, print timers
        if (FLTIMER.and.mod((itime-istart),oi_timer).eq.0) then
           print *, 'Intermediate timings, time = ', itime
           call print_timers
        endif
c     
c     all done with one complete time step
      enddo!itime
c
c     fill in and output remaining data for the final time-step
      call advance(1)
c
      close(26)
      close(60)
      close(30)
      close(31)
      close(80)
      write(40,110) (zets(k),u(nxh,nyh,k),k=1,nz)
      if (FLTIMER) then
        print *, 'Final timings'
        call print_timers
      endif
c
 101  FORMAT(7i5,3f10.5)
 102  FORMAT(3f10.5,2i5)
 103  FORMAT(4i5) 
 109  FORMAT(i5,2f10.5)
 110  FORMAT(f12.8,5x,f12.8)
 200  FORMAT(1p5e18.10)
c 200  FORMAT(i9,1p6e18.10)
      END
