      SUBROUTINE input
#include "common.inc"
c
      INTEGER itold
      REAL*8  rtimet,dt_oldt
      INTEGER idum
      REAL*8  dum
      CHARACTER*32 anum

      open(11,file=anum(ch_file,istart),form=format_mode)
      rewind(11)
      read(11) rtime,itold,dt_old,idum,idum,idum,dum,dum
      read(11) u,v,w
      read(11) pbxnt,pcxnt,pbynt,pcynt,pbxot,pcxot,pbyot,pcyot
      read(11) pbxnb,pcxnb,pbynb,pcynb,pbxob,pcxob,pbyob,pcyob
      close(11)
      write(26,*)'input: rtime,itold,dt_old ', rtime,itold,dt_old
      if (itold.ne.istart) then
        write(26,*)'Error in input: itold is not istart'
        stop
      endif
c
      if (FLTHRM) then
        open(11,file=anum(tt_file,istart),form=format_mode)
        rewind(11)
        read(11) rtimet,itold,dt_oldt,idum,idum,idum,dum,dum
        read(11) tt
        close(11)
        if (itold.ne.istart) then
          write(26,*)'Error in input: itold is not istart'
          stop
        endif
        if ((rtime.ne.rtimet).or.(dt_old.ne.dt_oldt)) then
          write(26,*)'Error in input: time mismatch with thermal file',
     &      rtime, rtimet, dt_old, dt_oldt
          stop
        endif
      endif
c
c if non-boussinesq or particle phase read material derivative
      if(FLNBSQ.or.FLINERTIA)then
       if(FLMDRHS)then
         open(11,file=anum('DDt',istart),form=format_mode)
         rewind(11)
         read(11) rtimet,itold,dt_oldt,idum,idum,idum,dum,dum
         read(11) DuDt,DvDt,DwDt
         close(11)
         if(itold.ne.istart) then
            write(26,*)'Error in input: itold is not istart'
            stop
         endif
         if((rtime.ne.rtimet).or.(dt_old.ne.dt_oldt)) then
           write(26,*)'Error in input: time mismatch with DDt file',
     &           rtime, rtimet, dt_old, dt_oldt
            stop
         endif
c
       else!if(FLMDLHS)then
c read old_vel1
         open(11,file=anum('old_vel1',istart),form=format_mode)
         rewind(11)
         read(11) rtimet,itold,dt_oldt,idum,idum,idum,dum,dum
         read(11) ou1,ov1,ow1
         close(11)
         if(itold.ne.istart) then
            write(26,*)'Error in input: itold is not istart'
            stop
         endif
         if((rtime.ne.rtimet).or.(dt_old.ne.dt_oldt)) then
           write(26,*)'Error in input: time mismatch with old_vel file',
     &           rtime, rtimet, dt_old, dt_oldt
            stop
         endif
c read old_vel2
         open(11,file=anum('old_vel2',istart),form=format_mode)
         rewind(11)
         read(11) rtimet,itold,dt_oldt,idum,idum,idum,dum,dum
         read(11) ou2,ov2,ow2
         close(11)
         if(itold.ne.istart) then
            write(26,*)'Error in input: itold is not istart'
            stop
         endif
         if((rtime.ne.rtimet).or.(dt_old.ne.dt_oldt)) then
           write(26,*)'Error in input: time mismatch with old_vel file',
     &           rtime, rtimet, dt_old, dt_oldt
            stop
         endif
c read old_vel3
         open(11,file=anum('old_vel3',istart),form=format_mode)
         rewind(11)
         read(11) rtimet,itold,dt_oldt,idum,idum,idum,dum,dum
         read(11) ou3,ov3,ow3
         close(11)
         if (itold.ne.istart) then
            write(26,*)'Error in input: itold is not istart'
            stop
         endif
         if((rtime.ne.rtimet).or.(dt_old.ne.dt_oldt)) then
           write(26,*)'Error in input: time mismatch with old_vel file',
     &           rtime, rtimet, dt_old, dt_oldt
            stop
         endif
       endif!end FLMDXXX
      endif!end (FLNBSQ.or.FLINERTIA)
c
c if w_settling.ne.0 input deposit
c      if(FLEro.and.w_settling.ne.0)then 
       if(.false..and.w_settling.ne.0)then
         open(11,file=anum('deposit',istart),form=format_mode)
         rewind(11)
         read(11) rtimet,itold,dt_oldt,idum,idum,idum,dum,dum,idum
         read(11) deposit,erosion
         close(11)
      endif
c
      return
      END
c
      SUBROUTINE output
#include "common.inc"
c
      CHARACTER*32 anum
c      integer dum1,dum2,dum3
c      real*8 dum4,dum5
c
      if (FLTIMER) call begin_timer(30)
      open(21,file=anum(ch_file,itime),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay
      write(21) u,v,w
      write(21) pbxnt,pcxnt,pbynt,pcynt,pbxot,pcxot,pbyot,pcyot
      write(21) pbxnb,pcxnb,pbynb,pcynb,pbxob,pcxob,pbyob,pcyob
      close(21)
c     
      if(FLTHRM)then         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC the line below is used for tricking the initial sediment   CCCC
CCC concentration. for the case in which the fully             CCCC
CCC developed clear fluid turbulence is used as the            CCCC
CCC initial condition we have to give an initial               CCCC
CCC concentration profile. this will be accomplished           CCCC
CCC by outputting a concentration profile after an iteration   CCCC
CCC at t=0.0 which will be used as an initial condition for    CCCC
CCC turbulent flow with sediment at t=to (the latest output    CCCC
CCC time for clear fluid turbulent flow). Therefore the time   CCCC
CCC in the output file should be changed.                      CCCC
CCC                              to=itime(turb)*dt             CCCC
CCC                                                            CCCC
CCC                                       EMRE                 CCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
         
c         open(45,file='conc.2741001',form=format_mode)
c         read(45) rtime,itime,dt,dum1,dum2,dum3,dum4,dum5
c         read(45) tt        
c         close(45)
c
c         open(44,file='vel.1481001',form=format_mode)
c         read(44) rtime,itime,dt,dum1,dum2,dum3,dum4,dum5        
c         close(44)
  
         open(21,file=anum(tt_file,itime),form=format_mode)
         write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay
 
C         open(21,file=anum(tt_file,itime),form=format_mode)
C         write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay
c         write(21) tt*0.01
         write(21) tt
         close(21)
      endif
c     
c     if non-boussinesq or particle phase
c     write material derivative or old ve in Fourier space
      if((FLNBSQ.or.FLINERTIA))then
       if(FLMDRHS)then
         open(21,file=anum('DDt',itime),form=format_mode)
         write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay
         write(21) DuDt,DvDt,DwDt
         close(21)
c
       else!if(FLMDLHS)then
         call output_old
       endif!end FLMDXXX
      endif
c     
      if (FLTIMER) call end_timer(30)
      return
      END
c
      SUBROUTINE output_old
#include "common.inc"
c
      CHARACTER*32 anum
c
      if (FLTIMER) call begin_timer(35)
c write old_x1
      open(21,file=anum('old_vel1',itime),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay
      write(21) ou1,ov1,ow1
      close(21)
c write old_x2
      open(21,file=anum('old_vel2',itime),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay
      write(21) ou2,ov2,ow2
      close(21)
c write old_x3
      open(21,file=anum('old_vel3',itime),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay
      write(21) ou3,ov3,ow3
      close(21)
c
      if (FLTIMER) call end_timer(35)
      return
      END
c
      SUBROUTINE output_up(aup,avp,awp)
#include "common.inc"
      REAL*8 aup(nx0,ny0,nz),avp(nx0,ny0,nz),awp(nx0,ny0,nz)
c
      CHARACTER*32 anum
c
c      if (FLTIMER) call begin_timer()
      open(21,file=anum('vel_p',itime),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay
      write(21) aup,avp,awp
      close(21)
c     
c      if (FLTIMER) call end_timer()
      return
      END
c
      SUBROUTINE output_2d_section(name,val2d_1,val2d_2,val2d_3,dim)
#include "common.inc"
      REAL*8 val2d_1(nx0,ny0),val2d_2(nx0,ny0),val2d_3(nx0,ny0)
      INTEGER dim
      CHARACTER*(*) name
c
      CHARACTER*32 anum
c
      if (FLTIMER) call begin_timer(36)
c write averaged value
      open(21,file=anum(name,itime),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,0
      if(dim.eq.1) write(21) val2d_1
      if(dim.eq.2) write(21) val2d_1,val2d_2
      if(dim.eq.3) write(21) val2d_1,val2d_2,val2d_3
      close(21)
c
      if (FLTIMER) call end_timer(36)
      return
      END
c
      SUBROUTINE input_means
#include "common.inc"
c
      CHARACTER*32 anum
      INTEGER nz_old,idum
      REAL*8 dum
c
c      if (FLTIMER) call begin_timer()
      open(21,file=anum('means',itime),form=format_mode)
      read(21) NT,dum,idum,dum,idum,idum,nz_old,dum,dum
      if(nz.ne.nz_old)then
         write(26,*)'Error in input of mean: nz_old ne nz'
         stop
      endif
      read(21)ust2,usb2,vst2,vsb2,dum
      read(21)um,vm,wm,ttm,u2m,v2m,w2m,tt2m,u3m,v3m,w3m,tt3m,u4m,v4m
     &,w4m,tt4m,uvm,uwm,uttm,vwm,vttm,wttm,uxm,uym,uzm,ux2m,uy2m,uz2m
     &,ux3m,uy3m,uz3m,ux4m,uy4m,uz4m,vxm,vym,vzm,vx2m,vy2m,vz2m,vx3m
     &,vy3m,vz3m,vx4m,vy4m,vz4m,wxm,wym,wzm,wx2m,wy2m,wz2m,wx3m,wy3m
     &,wz3m,wx4m,wy4m,wz4m,ttxm,ttym,ttzm,ttx2m,tty2m,ttz2m,ttx3m,tty3m
     &,ttz3m,ttx4m,tty4m,ttz4m
      close(21)
c
c      if (FLTIMER) call end_timer()
      return
      END
c
      SUBROUTINE output_means
#include "common.inc"
c
      CHARACTER*32 anum
c
c      if (FLTIMER) call begin_timer()
      open(21,file=anum('means',itime),form=format_mode)
      write(21)NT,rtime,itime,dt,nx0,ny0,nz,deltax,deltay
      write(21)ust2,usb2,vst2,vsb2,re
      write(21)um,vm,wm,ttm,u2m,v2m,w2m,tt2m,u3m,v3m,w3m,tt3m,u4m,v4m
     &,w4m,tt4m,uvm,uwm,uttm,vwm,vttm,wttm,uxm,uym,uzm,ux2m,uy2m,uz2m
     &,ux3m,uy3m,uz3m,ux4m,uy4m,uz4m,vxm,vym,vzm,vx2m,vy2m,vz2m,vx3m
     &,vy3m,vz3m,vx4m,vy4m,vz4m,wxm,wym,wzm,wx2m,wy2m,wz2m,wx3m,wy3m
     &,wz3m,wx4m,wy4m,wz4m,ttxm,ttym,ttzm,ttx2m,tty2m,ttz2m,ttx3m,tty3m
     &,ttz3m,ttx4m,tty4m,ttz4m
      close(21)
c
c      if (FLTIMER) call end_timer()
      return
      END
c
      SUBROUTINE mean5_out(cl1,cl2,cl3,cl4,cl5,m1,m2,m3,m4,m5,outfile,t)
#include "common.inc"
      INTEGER t
      REAL*8 cl1(nz),cl2(nz),cl3(nz),cl4(nz),cl5(nz),m1,m2,m3,m4,m5
      CHARACTER*(*) outfile
c
      INTEGER k
      CHARACTER*32 anum
c
      open(21,file=anum(outfile,t))
      write(21,1000) m1, m2, m3, m4, m5
      do k=1,nz
        write(21,2000) zets(k), cl1(k), cl2(k), cl3(k), cl4(k), cl5(k)
      enddo
      close(21)
c
      return
 1000 FORMAT('#   means:        ', 1p9e18.10)
 2000 FORMAT(1p9e18.10)
      END
c
      SUBROUTINE spec_out(espx,dspx,espy,dspy,outfile,t)
#include "common.inc"
      INTEGER t
      REAL*8 espx(nxh,nz),dspx(nxh,nz),espy(nyh,nz),dspy(nyh,nz)
      CHARACTER*(*) outfile
c
      CHARACTER*32 anum
c
      open(21,file=anum(outfile,t),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay
      write(21) espx
      write(21) dspx
      write(21) espy
      write(21) dspy
      close(21)
c
      return
      END
c
      SUBROUTINE chb_spec_out(espz,outfile,t)
#include "common.inc"
      INTEGER t
      REAL*8 espz(nz)
      CHARACTER*(*) outfile
c
      CHARACTER*32 anum
c
      open(21,file=anum(outfile,t),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay
      write(21) espz
      close(21)
c
      return
      END
c
      SUBROUTINE output_press(code)
#include "common.inc"
      INTEGER code
c
c     This subroutine postprocesses and outputs pressure for output 
c     according to code:
c     code=1 -> output dyn press + hydro press due to strat by variable tt
c     code=2 -> output dyn press
c     Needs tt in physical space
c
      real*8 tmeanp,hydp,s
      integer i,j,k,l
      CHARACTER*32 anum
c
c      if (FLTIMER) call begin_timer(36)
c
c     copy p to rw
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(k)
      do k=1,nz
        call dcopy(nx0y,p(1,1,k),1,rw(1,1,k),1)
      enddo

C$OMP END PARALLEL DO
c
c     IFFT rw -> rw is pressure
      call fft_f2r_2d_new(rw)
c
c     substract stratification pressure 
c     and compute total mean dynamic pressure
      tmeanp = 0.d0
      if(code.eq.2)then
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,l,hydp,s)
C$OMP&                           ,REDUCTION(+:tmeanp)
         do k=1,nz
            s=0.d0
            do j=1,ny
               do i=1,nx
                  hydp=0.d0
                  do l=1,k
                     hydp=hydp+tt(i,j,l)*(2.d0*wzet(l))
                  enddo
                  rw(i,j,k)=rw(i,j,k)-hydp
                  s=s+rw(i,j,k)
               enddo
            enddo
            s=s/dble(nx*ny)
            tmeanp=tmeanp+s*wzet(k)
         enddo
C$OMP END PARALLEL DO
      elseif(code.eq.1)then
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,s)
C$OMP&                           ,REDUCTION(+:tmeanp)
         do k=1,nz
            s=0.d0
            do j=1,ny
               do i=1,nx
                  s=s+rw(i,j,k)
               enddo
            enddo
            s=s/dble(nx*ny)
            tmeanp=tmeanp+s*wzet(k)
         enddo
C$OMP END PARALLEL DO
      endif
c
c     substract mean pressure to normalize output
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  rw(i,j,k)=rw(i,j,k)-tmeanp
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO
c
      print *, 'before writing' 
      open(21,file=anum('press',itime),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay
      write(21) rw    
      print *,' the output is written to the file'
      close(21)
      print *,'file is closed'
c
c      if (FLTIMER) call end_timer(36)
      return
      END
c
      FUNCTION anum(string,num)
      IMPLICIT none
      CHARACTER*32 anum
      CHARACTER*(*) string
      INTEGER num
c
      INTEGER i,n,dot,fin
c
      anum = string
      dot = index(anum,' ')
      if (dot.gt.0.and.num.ge.0) then
        if (dot.le.32) anum(dot:dot) = '.'
        fin = dot
        n = num
 100    fin = fin + 1
        n = n/10
        if (n.gt.0) go to 100
        n = num
        do i=fin,dot+1,-1
          if (i.le.32) anum(i:i) = char(ichar('0')+mod(n,10))
          n = n/10
        enddo
      endif
c
      return
      END
