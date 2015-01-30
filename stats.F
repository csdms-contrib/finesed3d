c     This computes the x- and y- spectra for energy and dissipation.
c     Contributions from the mean profile aren't included.
c     Also outputs the rms energy as a function of z, and the rms
c     average of this over the box
c
c     xrms and yrms should be the same, as should txrms and tyrms.
c     (Both are included as a check.)
c
      SUBROUTINE make_spec_new(uf,espx,dspx,espy,dspy,xrms,txrms,
     &  yrms,tyrms)
#include "common.inc"
      REAL*8 uf(nx0,ny0,nz)
      REAL*8 espx(nxh,nz),dspx(nxh,nz),espy(nyh,nz),dspy(nyh,nz)
      REAL*8 xrms(nz),txrms,yrms(nz),tyrms
c
      INTEGER i,ii,j,jimag,k
      REAL*8 ensum,dssum,temp
c
      if (FLTIMER) call begin_timer(110)
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP&            PRIVATE(i,ii,j,jimag,k,ensum,dssum,temp)
      do k=1,nz
c
c     compute x-spectra
c
        ensum = 0.d0
        dssum = 0.d0
        do j=2,nyh
          jimag=nyh+j
          temp = uf(1,j,k)**2 + uf(1,jimag,k)**2
          ensum = ensum + temp
          dssum = dssum + waveys(j)*temp
        enddo
        espx(1,k) = 2.d0*ensum
        dspx(1,k) = 2.d0*dssum
        xrms(k) = espx(1,k)
        do i=2,nxh
          ii = nx+2-i
          ensum = 0.5d0*(uf(i, 1,k)**2 + uf(i, nyhp,k)**2 +
     &                   uf(ii,1,k)**2 + uf(ii,nyhp,k)**2)
          dssum = 0.d0
          do j=2,nyh
            jimag=nyh+j
            temp = uf(i, j,k)**2 + uf(i, jimag,k)**2 +
     &             uf(ii,j,k)**2 + uf(ii,jimag,k)**2
            ensum = ensum + temp
            dssum = dssum + waveys(j)*temp
          enddo
          espx(i,k) = 2.d0*ensum
          dspx(i,k) = 2.d0*(wavexs(i)*ensum + dssum)
          xrms(k) = xrms(k) + espx(i,k)
        enddo
c
c     compute y-spectra
c
        ensum = 0.d0
        dssum = 0.d0
        do i=2,nxh
          ii = nx+2-i
          temp = uf(i, 1,k)**2 + uf(i, nyhp,k)**2 +
     &           uf(ii,1,k)**2 + uf(ii,nyhp,k)**2
          ensum = ensum + temp
          dssum = dssum + wavexs(i)*temp
        enddo
        espy(1,k) = ensum
        dspy(1,k) = dssum
        yrms(k) = espy(1,k)
        do j=2,nyh
          jimag=nyh+j
          ensum = uf(1,j,k)**2 + uf(1,jimag,k)**2
          dssum = 0.d0
          do i=2,nxh
            ii = nx+2-i
            temp = uf(i, j,k)**2 + uf(i, jimag,k)**2 +
     &             uf(ii,j,k)**2 + uf(ii,jimag,k)**2
            ensum = ensum + temp
            dssum = dssum + wavexs(i)*temp
          enddo
          espy(j,k) = 2.d0*ensum
          dspy(j,k) = 2.d0*(waveys(j)*ensum + dssum)
          yrms(k) = yrms(k) + espy(j,k)
        enddo
      enddo
C$OMP END PARALLEL DO
c
      txrms = 0.d0
      tyrms = 0.d0
      do k=1,nz
        txrms = txrms + xrms(k)*wzet(k)
        tyrms = tyrms + yrms(k)*wzet(k)
        xrms(k) = sqrt(xrms(k))
        yrms(k) = sqrt(yrms(k))
      enddo
      txrms = sqrt(txrms)
      tyrms = sqrt(tyrms)
c
      if (FLTIMER) call end_timer(110)
      return
      END
c
c     This computes Chebyshev energy spectra 
c
      SUBROUTINE make_chb_spec(val,espz)
#include "common.inc"
      REAL*8 val(nx0,ny0,nz)
      REAL*8 espz(nz)
c
      INTEGER i,j,jcomp,k
c
c     First half
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j)
      do j = 1,nyh
c     transpose data
        do k = 1,nz
          call dcopy(nx,val(1,j,k),1,auxt1(k,1,j),nz0)
        enddo
c     compute
        call dgemm('N','N',nz,nx,nz,1.d0,chb,nz0,auxt1(1,1,j),nz0,
     &                                      0.d0,auxt2(1,1,j),nz0)
      enddo
C$OMP END PARALLEL DO
c     
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
      do k=1,nz
       espz(k)=0.d0
       do j=1,nyh 
        do i=1,nx
         espz(k)=espz(k)+auxt2(k,i,j)**2
        enddo
       enddo
       espz(k)=espz(k)/dble(nx*ny)
      enddo
C$OMP END PARALLEL DO
c     Second half
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j,jcomp)
      do j = 1,nyh
        jcomp=nyh+j
c     transpose data
        do k = 1,nz
          call dcopy(nx,val(1,jcomp,k),1,auxt1(k,1,j),nz0)
        enddo
c     compute
        call dgemm('N','N',nz,nx,nz,1.d0,chb,nz0,auxt1(1,1,j),nz0,
     &                                      0.d0,auxt2(1,1,j),nz0)
      enddo
C$OMP END PARALLEL DO
c     
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
      do k=1,nz
       do j=1,nyh
        do i=1,nx
         espz(k)=espz(k)+auxt2(k,i,j)**2
        enddo
       enddo
       espz(k)=espz(k)/dble(nx*ny)
      enddo
C$OMP END PARALLEL DO
c
      return
      END 
c
c     This computes the mean as a function of z
c     and the total mean over the box.
c     Uses Fourier space information.
c
      SUBROUTINE make_meanf(ur,mean,tmean)
#include "common.inc"
      REAL*8 ur(nx0,ny0,nz),mean(nz),tmean
c
      INTEGER k
c
      tmean = 0.d0
      do k=1,nz
        mean(k) = ur(1,1,k)
        tmean = tmean + mean(k)*wzet(k)
      enddo
c
      return
      END
c
c     This computes the mean and rms as a function of z
c     and the total mean and rms over the box.
c     Uses physical space information.
c
      SUBROUTINE make_meanp(val,mean,tmean,rms,trms)
#include "common.inc"
      REAL*8 val(nx0,ny0,nz),mean(nz),tmean,rms(nz),trms
c
      INTEGER i,j,k
      REAL*8 s,su
c
      if (FLTIMER) call begin_timer(111)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,s,su)
      do k=1,nz
        su = 0.d0
        do j=1,ny
          s = 0.d0
          do i=1,nx
            s = s + val(i,j,k)
          enddo
          su = su + s
        enddo
        mean(k) = su/dble(nx*ny)
c
        su = 0.d0
        do j=1,ny
          s = 0.d0
          do i=1,nx
            s = s + (val(i,j,k)-mean(k))**2
          enddo
          su = su + s
        enddo
        rms(k) = su/dble(nx*ny)
      enddo
C$OMP END PARALLEL DO
c
      tmean = 0.d0
      trms = 0.d0
      do k=1,nz
        tmean = tmean + mean(k)*wzet(k)
        trms = trms + rms(k)*wzet(k)
        rms(k) = sqrt(rms(k))
      enddo
      trms = sqrt(trms)
c
      if (FLTIMER) call end_timer(111)
      return
      END
c
c     This computes the mean square of a quantity
c
      SUBROUTINE make_msq(val,msq)
#include "common.inc"
      REAL*8 val(nx0,ny0,nz),msq
c
      INTEGER i,j,k
      REAL*8 s,su,msqz(nz)
c
      if (FLTIMER) call begin_timer(112)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,s,su)
      do k=1,nz
        su = 0.d0
        do j=1,ny
          s = 0.d0
          do i=1,nx
            s = s + val(i,j,k)**2
          enddo
          su = su + s
        enddo
        msqz(k) = su/dble(nx*ny)
      enddo
C$OMP END PARALLEL DO
c
      msq = 0.d0
      do k=1,nz
        msq = msq + msqz(k)*wzet(k)
      enddo
c
      if (FLTIMER) call end_timer(112)
      return
      END
c
c
      SUBROUTINE make_meant(val,valm,tvalm)
#include "common.inc"
      REAL*8 val(nz),valm(nz)
      REAL*8 tvalm
c
      INTEGER k
c
c      if (FLTIMER) call begin_timer()
      tvalm=0.d0
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(k),REDUCTION(+:tvalm)
      do k=1,nz
         valm(k)=(dble(NT-1)*valm(k)+val(k))/dble(NT)
         tvalm=tvalm+valm(k)*wzet(k)
      enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer()
      return
      END
c
c
c     This computes the quantities on the RHS of global balances
c     Uses physical space quantities
c
      SUBROUTINE make_gbalrhs(KE,Ep,diss,ddtEp)
#include "common.inc"
      REAL*8 KE,Ep,diss,ddtEp
c
      INTEGER i,j,k
      REAL*8 mean(nz),tmean,s,su
c
      if (FLTIMER) call begin_timer(115)
c
c     compute KE and 
c
      KE = 0.d0
      call make_msq(u,tmean)
      KE = KE + tmean
      call make_msq(v,tmean)
      KE = KE + tmean
      call make_msq(w,tmean)
      KE = KE + tmean
c 
      KE = KE*rlenx*rleny ! recall 2*rlenx*rleny is the vol of domain
c
c
c     compute KE dissipation rate
c
      diss = 0.d0
c
      call ddx(u,ru)
      call make_msq(ru,tmean)
      diss = diss + tmean
      call ddy(v,ru)
      call make_msq(ru,tmean)
      diss = diss + tmean
      call ddz(w,ru)
      call make_msq(ru,tmean)
      diss = diss + tmean
c
      call ddy(u,ru)
      call ddx(v,rv)
      call ddz(u,rw)
      call ddx(w,rtt)
      call ddy(w,rp)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
       do k=1,nz
        do j=1,ny
         do i=1,nx
           ru(i,j,k)=(ru(i,j,k)+ rv(i,j,k))/2.d0
           rw(i,j,k)=(rw(i,j,k)+rtt(i,j,k))/2.d0
         enddo
        enddo
       enddo
C$OMP END PARALLEL DO
      call ddz(v,rv)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
       do k=1,nz
        do j=1,ny
         do i=1,nx
            rv(i,j,k)=( rv(i,j,k)+  rp(i,j,k))/2.d0
         enddo
        enddo
       enddo
C$OMP END PARALLEL DO
      call make_msq(ru,tmean)
      diss = diss + 2.d0*tmean
      call make_msq(rw,tmean)
      diss = diss + 2.d0*tmean
      call make_msq(rv,tmean)
      diss = diss + 2.d0*tmean
c
      diss = 2.d0*diss/re*(2.d0*rlenx*rleny)! recall 2*rlenx*rleny is the vol of domain
c
c
c     compute Ep 
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
      do k=1,nz
        mean(k)=0.d0
        do j=1,ny
          do i=1,nx
            mean(k) = mean(k) + (zets(k)+1.d0)*tt(i,j,k)
          enddo
        enddo
        mean(k) = mean(k)*deltax*deltay
      enddo
C$OMP END PARALLEL DO
c
      Ep = 0.d0
      do k=1,nz
        Ep = Ep + mean(k)*wzet(k)
      enddo
      Ep = Ep * 2.d0      
c
c
c      compute ddtEp
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
      do k=1,nz
        mean(k)=0.d0
        do j=1,ny
          do i=1,nx
            mean(k) = mean(k) + w(i,j,k)*tt(i,j,k)
          enddo
        enddo
        mean(k) = mean(k)*deltax*deltay
      enddo
C$OMP END PARALLEL DO
c
      ddtEp = 0.d0
      do k=1,nz
        ddtEp = ddtEp + mean(k)*wzet(k)
      enddo
      ddtEp = ddtEp * 2.d0
c
      if(FLSETTLE)then
       call ddz(tt,rtt)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
       do k=1,nz
         mean(k)=0.d0
         do j=1,ny
           do i=1,nx
             mean(k) = mean(k) + (zets(k)+1.d0)*(-w_settling)*rtt(i,j,k)
           enddo
         enddo
         mean(k) = mean(k)*deltax*deltay
       enddo
C$OMP END PARALLEL DO        
c
       do k=1,nz
         ddtEp = ddtEp + mean(k)*wzet(k)*2.d0
       enddo
      endif
c
      if (FLTIMER) call end_timer(115)
      return
      END
c
c
c     This computes the accumulative means
c     Assumes variable in Fourier space and all reusable variables free
c
      SUBROUTINE make_means
#include "common.inc"
c     
      INTEGER i,j,k
      REAL*8 sum,svm,swm,sttm,su2m,sv2m,sw2m,stt2m,su3m,sv3m,sw3m,stt3m,
     &       su4m,sv4m,sw4m,stt4m,suvm,suwm,suttm,svwm,svttm,swttm
c
      REAL*8 ddot
      EXTERNAL ddot
c     
      NT=NT+1
c     
c     variables
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(k)
      do k=1,nz
         call dcopy(nx0y0,u( 1,1,k),1,ru( 1,1,k),1)
         call dcopy(nx0y0,v( 1,1,k),1,rv( 1,1,k),1)
         call dcopy(nx0y0,w( 1,1,k),1,rw( 1,1,k),1)
         call dcopy(nx0y0,tt(1,1,k),1,rtt(1,1,k),1)
      enddo
C$OMP END PARALLEL DO
      call fft_f2r_2d_new(ru)
      call fft_f2r_2d_new(rv)
      call fft_f2r_2d_new(rw)
      call fft_f2r_2d_new(rtt)
c
c     v* u* at top and bottom
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do j = 1,ny
        do i = 1,nx
           vsb2=vsb2+ddot(nz,chbd1(nz,1),nz0,rv(i,j,1),nx0y0)/dble(nxny)
           vst2=vst2+ddot(nz,chbd1(1,1),nz0,rv(i,j,1),nx0y0)/dble(nxny)
c
           usb2=usb2+ddot(nz,chbd1(nz,1),nz0,ru(i,j,1),nx0y0)/dble(nxny)
           ust2=ust2+ddot(nz,chbd1(1,1),nz0,ru(i,j,1),nx0y0)/dble(nxny)
        enddo
      enddo
C$OMP END PARALLEL DO
cccccccccccccccc
c      vsb2=vsb2+ddot(nz,chbd1(nz,1),nz0,v(1,1,1),nx0y0)
c      vst2=vst2+ddot(nz,chbd1(1,1),nz0,v(1,1,1),nx0y0)
cc     
c      usb2=usb2+ddot(nz,chbd1(nz,1),nz0,u(1,1,1),nx0y0)
c      ust2=ust2+ddot(nz,chbd1(1,1),nz0,u(1,1,1),nx0y0)
ccccccccccccccccc
      write(11,11) rtime-dt,vsb2/dble(NT)/re,vst2/dble(NT)/re,usb2/dble(
     &NT)/re,ust2/dble(NT)/re
 11   format(g25.15,g25.15,g25.15,g25.15,g25.15)
c
c     means
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP& PRIVATE(i,j,k,sum,svm,swm,sttm,su2m,sv2m,sw2m,stt2m,
C$OMP&               su3m,sv3m,sw3m,stt3m,su4m,sv4m,sw4m,stt4m,
C$OMP&               suvm,suwm,suttm,svwm,svttm,swttm)
      do k=1,nz
         sum   = 0.d0 
         svm   = 0.d0
         swm   = 0.d0
         sttm  = 0.d0
         su2m  = 0.d0
         sv2m  = 0.d0
         sw2m  = 0.d0
         stt2m = 0.d0
         su3m  = 0.d0
         sv3m  = 0.d0
         sw3m  = 0.d0
         stt3m = 0.d0
         su4m  = 0.d0
         sv4m  = 0.d0
         sw4m  = 0.d0
         stt4m = 0.d0
         suvm  = 0.d0
         suwm  = 0.d0
         suttm = 0.d0
         svwm  = 0.d0
         svttm = 0.d0
         swttm = 0.d0 
         do j=1,ny
            do i=1,nx
               sum   = sum   + ru( i,j,k) 
               svm   = svm   + rv( i,j,k)
               swm   = swm   + rw( i,j,k)
               sttm  = sttm  + rtt(i,j,k)
               su2m  = su2m  + ru( i,j,k)**2
               sv2m  = sv2m  + rv( i,j,k)**2
               sw2m  = sw2m  + rw( i,j,k)**2
               stt2m = stt2m + rtt(i,j,k)**2
               su3m  = su3m  + ru( i,j,k)**3
               sv3m  = sv3m  + rv( i,j,k)**3
               sw3m  = sw3m  + rw( i,j,k)**3
               stt3m = stt3m + rtt(i,j,k)**3
               su4m  = su4m  + ru( i,j,k)**4
               sv4m  = sv4m  + rv( i,j,k)**4
               sw4m  = sw4m  + rw( i,j,k)**4
               stt4m = stt4m + rtt(i,j,k)**4
               suvm  = suvm  + ru( i,j,k)*rv( i,j,k)
               suwm  = suwm  + ru( i,j,k)*rw( i,j,k)
               suttm = suttm + ru( i,j,k)*rtt(i,j,k)
               svwm  = svwm  + rv( i,j,k)*rw( i,j,k)
               svttm = svttm + rv( i,j,k)*rtt(i,j,k)
               swttm = swttm + rw( i,j,k)*rtt(i,j,k)
            enddo
         enddo
         um(  k) = um(  k) + sum   / dble(nxny) 
         vm(  k) = vm(  k) + svm   / dble(nxny) 
         wm(  k) = wm(  k) + swm   / dble(nxny) 
         ttm( k) = ttm( k) + sttm  / dble(nxny) 
         u2m( k) = u2m( k) + su2m  / dble(nxny)
         v2m( k) = v2m( k) + sv2m  / dble(nxny)
         w2m( k) = w2m( k) + sw2m  / dble(nxny)
         tt2m(k) = tt2m(k) + stt2m / dble(nxny)
         u3m( k) = u3m( k) + su3m  / dble(nxny)
         v3m( k) = v3m( k) + sv3m  / dble(nxny)
         w3m( k) = w3m( k) + sw3m  / dble(nxny)
         tt3m(k) = tt3m(k) + stt3m / dble(nxny)
         u4m( k) = u4m( k) + su4m  / dble(nxny)
         v4m( k) = v4m( k) + sv4m  / dble(nxny)
         w4m( k) = w4m( k) + sw4m  / dble(nxny)
         tt4m(k) = tt4m(k) + stt4m / dble(nxny)
         uvm( k) = uvm( k) + suvm  / dble(nxny)
         uwm( k) = uwm( k) + suwm  / dble(nxny)
         uttm(k) = uttm(k) + suttm / dble(nxny)
         vwm( k) = vwm( k) + svwm  / dble(nxny)
         vttm(k) = vttm(k) + svttm / dble(nxny)
         wttm(k) = wttm(k) + swttm / dble(nxny)
      enddo
C$OMP END PARALLEL DO
      write(10,10) rtime-dt,vm(3*nz/4)/dble(NT),ttm(3*nz/4)/dble(NT)
 10   format(g25.15,g25.15,g25.15)
c
c     z-component of gradients
      call ddz(ru,rp)
      call ddz(rv,ru)
      call ddz(rw,rv)
      call ddz(rtt,rw)
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP& PRIVATE(i,j,k,sum,svm,swm,sttm,su2m,sv2m,sw2m,stt2m,
C$OMP&               su3m,sv3m,sw3m,stt3m,su4m,sv4m,sw4m,stt4m)
      do k=1,nz
         sum   = 0.d0 
         su2m  = 0.d0
         su3m  = 0.d0
         su4m  = 0.d0
         svm   = 0.d0
         sv2m  = 0.d0
         sv3m  = 0.d0
         sv4m  = 0.d0
         swm   = 0.d0
         sw2m  = 0.d0
         sw3m  = 0.d0
         sw4m  = 0.d0
         sttm  = 0.d0
         stt2m = 0.d0
         stt3m = 0.d0
         stt4m = 0.d0
         do j=1,ny
            do i=1,nx
               sum   = sum   + rp( i,j,k) 
               su2m  = su2m  + rp( i,j,k)**2
               su3m  = su3m  + rp( i,j,k)**3
               su4m  = su4m  + rp( i,j,k)**4
c
               svm   = svm   + ru( i,j,k)
               sv2m  = sv2m  + ru( i,j,k)**2
               sv3m  = sv3m  + ru( i,j,k)**3
               sv4m  = sv4m  + ru( i,j,k)**4
c     
               swm   = swm   + rv( i,j,k)
               sw2m  = sw2m  + rv( i,j,k)**2
               sw3m  = sw3m  + rv( i,j,k)**3
               sw4m  = sw4m  + rv( i,j,k)**4
c     
               sttm  = sttm  + rw(i,j,k)
               stt2m = stt2m + rw(i,j,k)**2
               stt3m = stt3m + rw(i,j,k)**3
               stt4m = stt4m + rw(i,j,k)**4
            enddo
         enddo
         uzm(  k) = uzm(  k) + sum   / dble(nxny) 
         uz2m( k) = uz2m( k) + su2m  / dble(nxny)
         uz3m( k) = uz3m( k) + su3m  / dble(nxny)
         uz4m( k) = uz4m( k) + su4m  / dble(nxny)
c     
         vzm(  k) = vzm(  k) + svm   / dble(nxny) 
         vz2m( k) = vz2m( k) + sv2m  / dble(nxny)
         vz3m( k) = vz3m( k) + sv3m  / dble(nxny)
         vz4m( k) = vz4m( k) + sv4m  / dble(nxny)
c     
         wzm(  k) = wzm(  k) + swm   / dble(nxny) 
         wz2m( k) = wz2m( k) + sw2m  / dble(nxny)
         wz3m( k) = wz3m( k) + sw3m  / dble(nxny)
         wz4m( k) = wz4m( k) + sw4m  / dble(nxny)
c     
         ttzm( k) = ttzm( k) + sttm  / dble(nxny) 
         ttz2m(k) = ttz2m(k) + stt2m / dble(nxny)
         ttz3m(k) = ttz3m(k) + stt3m / dble(nxny)
         ttz4m(k) = ttz4m(k) + stt4m / dble(nxny)
      enddo
C$OMP END PARALLEL DO
c
c     x-component of gradients
      call ddx_f_new(u,ru)
      call ddx_f_new(v,rv)
      call ddx_f_new(w,rw)
      call ddx_f_new(tt,rtt)
c     
      call fft_f2r_2d_new(ru)
      call fft_f2r_2d_new(rv)
      call fft_f2r_2d_new(rw)
      call fft_f2r_2d_new(rtt)
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP& PRIVATE(i,j,k,sum,svm,swm,sttm,su2m,sv2m,sw2m,stt2m,
C$OMP&               su3m,sv3m,sw3m,stt3m,su4m,sv4m,sw4m,stt4m)
      do k=1,nz
         sum   = 0.d0 
         su2m  = 0.d0
         su3m  = 0.d0
         su4m  = 0.d0
         svm   = 0.d0
         sv2m  = 0.d0
         sv3m  = 0.d0
         sv4m  = 0.d0
         swm   = 0.d0
         sw2m  = 0.d0
         sw3m  = 0.d0
         sw4m  = 0.d0
         sttm  = 0.d0
         stt2m = 0.d0
         stt3m = 0.d0
         stt4m = 0.d0
         do j=1,ny
            do i=1,nx
               sum   = sum   + ru( i,j,k) 
               su2m  = su2m  + ru( i,j,k)**2
               su3m  = su3m  + ru( i,j,k)**3
               su4m  = su4m  + ru( i,j,k)**4
c     
               svm   = svm   + rv( i,j,k)
               sv2m  = sv2m  + rv( i,j,k)**2
               sv3m  = sv3m  + rv( i,j,k)**3
               sv4m  = sv4m  + rv( i,j,k)**4
c     
               swm   = swm   + rw( i,j,k)
               sw2m  = sw2m  + rw( i,j,k)**2
               sw3m  = sw3m  + rw( i,j,k)**3
               sw4m  = sw4m  + rw( i,j,k)**4
c     
               sttm  = sttm  + rtt(i,j,k)
               stt2m = stt2m + rtt(i,j,k)**2
               stt3m = stt3m + rtt(i,j,k)**3
               stt4m = stt4m + rtt(i,j,k)**4
            enddo
         enddo
         uxm(  k) = uxm(  k) + sum   / dble(nxny) 
         ux2m( k) = ux2m( k) + su2m  / dble(nxny)
         ux3m( k) = ux3m( k) + su3m  / dble(nxny)
         ux4m( k) = ux4m( k) + su4m  / dble(nxny)
c     
         vxm(  k) = vxm(  k) + svm   / dble(nxny) 
         vx2m( k) = vx2m( k) + sv2m  / dble(nxny)
         vx3m( k) = vx3m( k) + sv3m  / dble(nxny)
         vx4m( k) = vx4m( k) + sv4m  / dble(nxny)
c     
         wxm(  k) = wxm(  k) + swm   / dble(nxny) 
         wx2m( k) = wx2m( k) + sw2m  / dble(nxny)
         wx3m( k) = wx3m( k) + sw3m  / dble(nxny)
         wx4m( k) = wx4m( k) + sw4m  / dble(nxny)
c     
         ttxm( k) = ttxm( k) + sttm  / dble(nxny) 
         ttx2m(k) = ttx2m(k) + stt2m / dble(nxny)
         ttx3m(k) = ttx3m(k) + stt3m / dble(nxny)
         ttx4m(k) = ttx4m(k) + stt4m / dble(nxny)
      enddo
C$OMP END PARALLEL DO
c
c     y-component of gradients
      call ddy_f_new(u,ru)
      call ddy_f_new(v,rv)
      call ddy_f_new(w,rw)
      call ddy_f_new(tt,rtt)
c     
      call fft_f2r_2d_new(ru)
      call fft_f2r_2d_new(rv)
      call fft_f2r_2d_new(rw)
      call fft_f2r_2d_new(rtt)
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP& PRIVATE(i,j,k,sum,svm,swm,sttm,su2m,sv2m,sw2m,stt2m,
C$OMP&               su3m,sv3m,sw3m,stt3m,su4m,sv4m,sw4m,stt4m)
      do k=1,nz
         sum   = 0.d0 
         su2m  = 0.d0
         su3m  = 0.d0
         su4m  = 0.d0
         svm   = 0.d0
         sv2m  = 0.d0
         sv3m  = 0.d0
         sv4m  = 0.d0
         swm   = 0.d0
         sw2m  = 0.d0
         sw3m  = 0.d0
         sw4m  = 0.d0
         sttm  = 0.d0
         stt2m = 0.d0
         stt3m = 0.d0
         stt4m = 0.d0
         do j=1,ny
            do i=1,nx
               sum   = sum   + ru( i,j,k) 
               su2m  = su2m  + ru( i,j,k)**2
               su3m  = su3m  + ru( i,j,k)**3
               su4m  = su4m  + ru( i,j,k)**4
c     
               svm   = svm   + rv( i,j,k)
               sv2m  = sv2m  + rv( i,j,k)**2
               sv3m  = sv3m  + rv( i,j,k)**3
               sv4m  = sv4m  + rv( i,j,k)**4
c     
               swm   = swm   + rw( i,j,k)
               sw2m  = sw2m  + rw( i,j,k)**2
               sw3m  = sw3m  + rw( i,j,k)**3
               sw4m  = sw4m  + rw( i,j,k)**4
c     
               sttm  = sttm  + rtt(i,j,k)
               stt2m = stt2m + rtt(i,j,k)**2
               stt3m = stt3m + rtt(i,j,k)**3
               stt4m = stt4m + rtt(i,j,k)**4
            enddo
         enddo
         uym(  k) = uym(  k) + sum   / dble(nxny) 
         uy2m( k) = uy2m( k) + su2m  / dble(nxny)
         uy3m( k) = uy3m( k) + su3m  / dble(nxny)
         uy4m( k) = uy4m( k) + su4m  / dble(nxny)
c     
         vym(  k) = vym(  k) + svm   / dble(nxny) 
         vy2m( k) = vy2m( k) + sv2m  / dble(nxny)
         vy3m( k) = vy3m( k) + sv3m  / dble(nxny)
         vy4m( k) = vy4m( k) + sv4m  / dble(nxny)
c     
         wym(  k) = wym(  k) + swm   / dble(nxny) 
         wy2m( k) = wy2m( k) + sw2m  / dble(nxny)
         wy3m( k) = wy3m( k) + sw3m  / dble(nxny)
         wy4m( k) = wy4m( k) + sw4m  / dble(nxny)
c     
         ttym( k) = ttym( k) + sttm  / dble(nxny) 
         tty2m(k) = tty2m(k) + stt2m / dble(nxny)
         tty3m(k) = tty3m(k) + stt3m / dble(nxny)
         tty4m(k) = tty4m(k) + stt4m / dble(nxny)
      enddo
C$OMP END PARALLEL DO
c
      return
      END





