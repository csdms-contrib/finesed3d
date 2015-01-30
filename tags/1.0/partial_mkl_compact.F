      SUBROUTINE ddx(val,dval)
#include "common.inc"
      include 'fftw3.f'
      REAL*8 val(nx0,ny0,nz),dval(nx0,ny0,nz)
c
      INTEGER i,j,k
      REAL*8     rtmp(tnxhp0)
      COMPLEX*16 ctmp( nxhp0)
#ifdef COPPER
c for ESSL
c fft_work must be 8-byte real array
      REAL*8 fft_work(nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
      REAL*8 fft_work(tnxhp0)
      INTEGER*8 fftw_plan
#endif
c
c      if (FLTIMER) call begin_timer(80)
c
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,rtmp,ctmp)
C$OMP+ PRIVATE(fft_work,fftw_plan)
      do k=1,nz
        do j = 1,ny
          call dcopy(nx,val(1,j,k),1,rtmp,1)
#ifdef COPPER
c for ESSL
          call drcft(0,rtmp,tnxhp0,ctmp,nxhp0,nx,1,1,1.d0/dble(nx),
     &    rx_tb_f,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
          call dfftw_plan_dft_r2c_1d(fftw_plan,nx,rtmp,ctmp,
     &         FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call zdscal(nxhp,1.d0/dble(nx),ctmp,1)
          call dfftw_destroy_plan(fftw_plan)
#endif
c
          do i=1,nxh
c          do i=1,nxhp
c           call zscal(1,(0.d0,1.d0)*wavex(i),ctmp(i),0)
c           ctmp(i)=dcmplx(((0.d0,1.d0)*wavex(i))*ctmp(i))
             ctmp(i)=dcmplx(0.d0,wavex(i))*ctmp(i)
          enddo
          ctmp(nxhp)=(0.d0,0.d0)
c
#ifdef COPPER
c for ESSL
          call dcrft(0,ctmp,nxhp0,rtmp,tnxhp0,nx,1,-1,1.d0,
     &    rx_tb_b,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
          call dfftw_plan_dft_c2r_1d(fftw_plan,nx,ctmp,rtmp,
     &         FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call dfftw_destroy_plan(fftw_plan)
c          call dscal(nx,1.d0/dble(nx),rtmp,1)
#endif
c
          call dcopy(nx,rtmp,1,dval(1,j,k),1)
        enddo
      enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer(80)
      return
      END
c
c
      SUBROUTINE ddy(val,dval)
#include "common.inc"
      include 'fftw3.f'
      REAL*8 val(nx0,ny0,nz),dval(nx0,ny0,nz)
c
      INTEGER i,j,k
      REAL*8     rtmp(tnyhp0)
      COMPLEX*16 ctmp( nyhp0)
#ifdef COPPER
c for ESSL
c fft_work must be 8-byte real array
      REAL*8 fft_work(nfft_aux)
#endif
#ifdef COBALT
c for SCSL
      INTEGER*8 fftw_plan
      REAL*8 fft_work(tnyhp0)
#endif
c
c      if (FLTIMER) call begin_timer(82)
c
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,rtmp,ctmp)
C$OMP+ PRIVATE(fft_work,fftw_plan)
      do k=1,nz
        do i=1,nx
          call dcopy(ny,val(i,1,k),nx0,rtmp,1)
c
#ifdef COPPER
c for ESSL
          call drcft(0,rtmp,tnyhp0,ctmp,nyhp0,ny,1,1,1.d0/dble(ny),
     &    ry_tb_f,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
          call dfftw_plan_dft_r2c_1d(fftw_plan,ny,rtmp,ctmp,
     &         FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call zdscal(nyhp,1.d0/dble(ny),ctmp,1)
          call dfftw_destroy_plan(fftw_plan)
#endif
c
          do j=1,nyh
c          do j=1,nyhp
c           call zscal(1,(0.d0,1.d0)*wavey(j),ctmp(j),0)
c           ctmp(j)=dcmplx(((0.d0,1.d0)*wavey(j))*ctmp(j))
             ctmp(j)=dcmplx(0.d0,wavey(j))*ctmp(j)
          enddo
          ctmp(nyhp)=(0.d0,0.d0)
c
#ifdef COPPER
c for ESSL
          call dcrft(0,ctmp,nyhp0,rtmp,tnyhp0,ny,1,-1,1.d0,
     &    ry_tb_b,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
          call dfftw_plan_dft_c2r_1d(fftw_plan,ny,ctmp,rtmp,
     &         FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call dfftw_destroy_plan(fftw_plan)
c          call dscal(ny,1.d0/dble(ny),rtmp,1)
#endif
c
          call dcopy(ny,rtmp,1,dval(i,1,k),nx0)
        enddo
      enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer(82)
      return
      END
c
c
c$$$      SUBROUTINE ddz(val,dval)!ddz_f_new(val,dval)
c$$$#include "common.inc"
c$$$      REAL*8 val(nx0,ny0,nz),dval(nx0,ny0,nz)
c$$$c
c$$$c     This subroutine computes ddz in real and Fourier spaces
c$$$c     There is not difference in its implementation
c$$$c
c$$$      INTEGER i,j,jcomp,k
c$$$c
c$$$c      if (FLTIMER) call begin_timer()
c$$$c
c$$$C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,jcomp,k)
c$$$      do j = 1,nyh
c$$$c     First part. This is done this way to reuse auxt1/auxt2(nz0,nx0,nyh)
c$$$c     traspose data
c$$$        do k = 1,nz
c$$$          call dcopy(nx,val(1,j,k),1,auxt1(k,1,j),nz0)
c$$$        enddo
c$$$c     cumpute
c$$$       call dgemm('N','N',nz,nx,nz,1.d0,chbd1,nz0,auxt1(1,1,j),nz0,
c$$$     &                                       0.d0,auxt2(1,1,j),nz0)
c$$$c     transpose back data
c$$$        do i=1,nx
c$$$          call dcopy(nz,auxt2(1,i,j),1,dval(i,j,1),nx0y0)
c$$$        enddo
c$$$c     Second part.
c$$$        jcomp=nyh+j
c$$$c     traspose data
c$$$        do k = 1,nz
c$$$          call dcopy(nx,val(1,jcomp,k),1,auxt1(k,1,j),nz0)
c$$$        enddo
c$$$c     compute
c$$$       call dgemm('N','N',nz,nx,nz,1.d0,chbd1,nz0,auxt1(1,1,j),nz0,
c$$$     &                                       0.d0,auxt2(1,1,j),nz0)
c$$$c     transpose back data
c$$$        do i=1,nx
c$$$          call dcopy(nz,auxt2(1,i,j),1,dval(i,jcomp,1),nx0y0)
c$$$        enddo
c$$$      enddo
c$$$C$OMP END PARALLEL DO
c$$$c
c$$$c      if (FLTIMER) call end_timer()
c$$$      return
c$$$      END
c$$$c
c

      subroutine ddz(val,dval)
#include "common.inc"
      real*8 val(nx0,ny0,nz),dval(nx0,ny0,nz)
      real*8 tmp(nz),tmp2(nz)
      real*8 diagu(nzm),diagl(nzm),diag(nz)
C     This subroutine computes ddz in both real and Fourire spaces
      integer i,j,k

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,tmp,tmp2,diag,diagu,diagl)
      do j=1,ny
         do i=1,nx
            call dcopy(nz,CA1(2,1),3,diag,1)
            call dcopy(nzm,CA1(1,2),3,diagu,1)
            call dcopy(nzm,CA1(3,1),3,diagl,1)
c     Copy data into tmp
            call dcopy(nz,val(i,j,1),nx0y0,tmp,1)
c     B*tmp
            call dgbmv('n',nz,nz,5,5,1.d0,CB1,11,tmp,1,0.d0,
     &           tmp2,1)
c     solve A*tmp2=tmp2
            call dgtsv(nz,1,diagl,diag,diagu,tmp2,nz,k)
c     copy back to dval
            call dcopy(nz,tmp2,1,dval(i,j,1),nx0y0)
         enddo 
      enddo
C$OMP END PARALLEL DO    
      return
      end

      SUBROUTINE ddx_f_new(var,dvar)
#include "common.inc"
      REAL*8 var(nx0,ny0,nz)
      REAL*8 dvar(nx0,ny0,nz)
c     The input var is preserved
      INTEGER i,j,k,jimag
c
c      if (FLTIMER) call begin_timer()
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,k)
        do k=1,nz
         do j=1,ny
           call dcopy(nx,var(1,j,k),1,dvar(1,j,k),1)
         enddo
         do i=1,nx
           call dscal(nyh, wavex(i),dvar(i,1   ,k),nx0)
           call dscal(nyh,-wavex(i),dvar(i,nyhp,k),nx0)
         enddo
         do j=1,nyh
           jimag=nyh+j
           call dswap(nx,dvar(1,j,k),1,dvar(1,jimag,k),1)
         enddo
        enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer()
      return
      END
c
c
      SUBROUTINE ddy_f_new(var,dvar)
#include "common.inc"
      REAL*8 var(nx0,ny0,nz)
      REAL*8 dvar(nx0,ny0,nz)
c     The input var is preserved
      INTEGER j,k,jimag
c
c      if (FLTIMER) call begin_timer()
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(j,jimag,k)
        do k=1,nz
         do j=1,nyh
           jimag=nyh+j
           call dcopy(nx,var(1,j    ,k),1,dvar(1,j    ,k),1)
           call dcopy(nx,var(1,jimag,k),1,dvar(1,jimag,k),1)
c
           call dscal(nx, wavey(j),dvar(1,j    ,k),1)
           call dscal(nx,-wavey(j),dvar(1,jimag,k),1)
c
           call dswap(nx,dvar(1,j,k),1,dvar(1,jimag,k),1)
         enddo
        enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer()
      return
      END
c
c
      SUBROUTINE ddt(irk,val0,val1,val2,val3,dval)
#include "common.inc"
      INTEGER irk
      REAL*8 val0(nx0,ny0,nz),val1(nx0,ny0,nz),val2(nx0,ny0,nz)
     &       ,val3(nx0,ny0,nz),dval(nx0,ny0,nz)
c
      REAL*8 hh0,hh1,hh2,hh3,a0,a1,a2,a3
      INTEGER i,j,k
c
c      if (FLTIMER) call begin_timer(86)
c initialize values 
      if(irk.eq.1)then
        hh0=0.d0
      elseif(irk.eq.2)then
        hh0=dt/3.d0
      elseif(irk.eq.3)then
        hh0=3.d0*dt/4.d0
      endif
      hh1=dt+hh0
      hh2=dt+hh1
      hh3=dt+hh2
c for restart with different dt
      if (istart.gt.0.and.itime.eq.istart+1)then
        hh1=dt_old+hh0
        hh2=dt_old+hh1
        hh3=dt_old+hh2
      elseif(istart.gt.0.and.itime.eq.istart+2)then
        hh2=dt_old+hh1
        hh3=dt_old+hh2
      elseif(istart.gt.0.and.itime.eq.istart+3)then
        hh3=dt_old+hh2
      endif     
c 3rd order
      a0=-(hh2*hh3+hh1*(hh2+hh3))/(hh0-hh1)/(hh0-hh2)/(hh0-hh3)
      a1=-(hh2*hh3+hh0*(hh2+hh3))/(hh1-hh0)/(hh1-hh2)/(hh1-hh3)
      a2=-(hh1*hh3+hh0*(hh1+hh3))/(hh2-hh0)/(hh2-hh1)/(hh2-hh3)
      a3=-(hh1*hh2+hh0*(hh1+hh2))/(hh3-hh0)/(hh3-hh1)/(hh3-hh2)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
       do k=1,nz
        do j=1,ny
         do i=1,nx
c 3rd order 
          dval(i,j,k)=a0*val0(i,j,k)+a1*val1(i,j,k)
     &               +a2*val2(i,j,k)+a3*val3(i,j,k)
         enddo
        enddo
       enddo
C$OMP END PARALLEL DO
c      if (FLTIMER) call end_timer(86)
      return
      END
