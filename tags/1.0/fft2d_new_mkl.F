c     For mkl, use fftw3x library, initfft is not need and can be deleted
c
c
      subroutine fft_r2f_2d_new(ur)
#include "common.inc"
      include 'fftw3.f'
      REAL*8 ur(nx0,ny0,nz)
c
c     2-d real to fourier transform in x and y directions
c     On input:  ur contains data on 1:nx,1:ny,1:nz
c     On output: ur contains data on (1:nx),1:nyh,1:nz   (real component)
c     On output: ur contains data on (1:nx),nyhp:ny,1:nz (imag component)
c
      INTEGER i,j,k,jimag
      REAL*8 tp(tnyhp0)
      COMPLEX*16 tfw(nyhp)
      COMPLEX*16 tf(nx0)
c      COMPLEX*16 tmp(nyhp0,nx0,nz)!this is a common block now ->tmpf
#ifdef COPPER
c for ESSL
c fft_work must be 8-byte real array
      REAL*8 fft_work(nfft_aux)
#endif
#ifdef COBALT
c for SCSL
      INTEGER*8 fftw_plan
      COMPLEX*16 tx(nx0)
#endif
c
c      if (FLTIMER) call begin_timer(70)
c
#ifdef COPPER
c for ESSL
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,k,tp,tf,fft_work)
#endif
#ifdef COBALT
c for MKL
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP&            PRIVATE(i,j,jimag,k,tp,tf,tfw,fftw_plan,tx)
#endif
      do k=1,nz
c compute fft in y-direction
        do i=1,nx
          call dcopy(ny,ur(i,1,k),nx0,tp,1)
#ifdef COPPER
c ESSL
          call drcft(0,tp,tnyhp0,tmpf(1,i,k),nyhp0,ny,1,1,1.d0/dble(ny),
     &    ry_tb_f,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c MKL
          call dfftw_plan_dft_r2c_1d(fftw_plan,ny,tp,tfw,FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call zdscal(nyhp,1.d0/dble(ny),tfw,1)
          call zcopy(nyhp,tfw,1,tmpf(1,i,k),1)
          call dfftw_destroy_plan(fftw_plan)
c          write(35,*) tp
c          write(36,*) tfw
#endif
        enddo
c compute fft in x-direction
        do j=1,nyh 
#ifdef COPPER
c for ESSL
          call dcft(0,tmpf(j,1,k),nyhp0,1,tf,1,nx0,nx,1,1,1.d0/dble(nx),
     &    cx_tb_f,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c for MKL
          call zcopy(nx,tmpf(j,1,k),nyhp0,tx,1)
          call dfftw_plan_dft_1d(fftw_plan,nx,tx,tf,FFTW_FORWARD,
     &         FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call zdscal(nx,1.d0/dble(nx),tf,1)
          call dfftw_destroy_plan(fftw_plan)
#endif
c
          jimag=nyh+j
          do i=1,nx
            ur(i,j    ,k) = dreal(tf(i))
            ur(i,jimag,k) = dimag(tf(i))
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer(70)
      return
      end
c
c
      subroutine fft_f2r_2d_new(ur)
#include "common.inc"
      include 'fftw3.f'
      REAL*8 ur(nx0,ny0,nz)
c
c     2-d fourier to real transform in x and y directions
c     On input:  ur contains data on 1:nx,   1:nyh,1:nz (real component)
c     On input:  ur contains data on 1:nx,nyhp:ny ,1:nz (imag component)
c     On output: ur contains data on 1:nx,1:ny,1:nz
c
      INTEGER i,j,k,jimag
      REAL*8 tp(tnyhp0)
      COMPLEX*16 tf(nx0)
      COMPLEX*16 tfw(nyhp)
c      COMPLEX*16 tmp(nyhp0,nx0,nz)!this is a common block now ->tmpf
#ifdef COPPER
c for ESSL
c fft_work must be 8-byte real array
      REAL*8 fft_work(nfft_aux)
#endif
#ifdef COBALT
c for SCSL
      COMPLEX*16 tx(nx0)
      INTEGER*8 fftw_plan
#endif
c
c      if (FLTIMER) call begin_timer(75)
c
#ifdef COPPER
c for ESSL
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,k,tf,tp,fft_work)
#endif
#ifdef COBALT
c for SCSL
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP&            PRIVATE(i,j,jimag,k,tf,tp,tx,tfw,fftw_plan)
#endif
      do k=1,nz
c cumpute fft in x-direction
        do j=1,nyh
          jimag=nyh+j
          do i=1,nx
            tf(i) = dcmplx(ur(i,j,k),ur(i,jimag,k))
          enddo
          tf(nxhp)=(0.d0,0.d0)
c
#ifdef COPPER
c for ESSL
          call dcft(0,tf,1,nx0,tmpf(j,1,k),nyhp0,1,nx,1,-1,1.d0,
     &    cx_tb_b,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
          call dfftw_plan_dft_1d(fftw_plan,nx,tf,tx,FFTW_BACKWARD,
     &         FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call dfftw_destroy_plan(fftw_plan)
          call zcopy(nx,tx,1,tmpf(j,1,k),nyhp0)
#endif
        enddo
c cumpute fft in y-direction
        do i=1,nx
          tmpf(nyhp,i,k)=(0.d0,0.d0)
c
#ifdef COPPER
c for ESSL
          call dcrft(0,tmpf(1,i,k),nyhp0,tp,tnyhp0,ny,1,-1,1.d0,
     &    ry_tb_b,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
          call zcopy(nyhp,tmpf(1,i,k),1,tfw,1)
          call dfftw_plan_dft_c2r_1d(fftw_plan,ny,tfw,tp,FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call dfftw_destroy_plan(fftw_plan)
c          call dscal(ny,1.d0/dble(nxny),tp,1)
#endif
          call dcopy(ny,tp,1,ur(i,1,k),nx0)
        enddo
      enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer(75)
      return
      end
c
c
      subroutine fftl_r2f_2d_new(ur)
#include "common.inc"
      include 'fftw3.f'
      REAL*8 ur(nxl0,nyl0,nz)
c
c     2-d real to fourier transform in x and y directions
c     On input:  ur contains data on 1:nxl,1:nyl,1:nz       (real component)
c     On output: ur contains data on (1:nxl),nylhp:nyl,1:nz (imag component)
c
      INTEGER i,j,k,jimag
      REAL*8 tp(tnylhp0)
      COMPLEX*16 tf(nxl0)
      COMPLEX*16 tfw(nylhp)
c      COMPLEX*16 tmp(nylhp0,nxl0,nz)!this is a common block now ->tmplf
#ifdef COPPER
c for ESSL
c fft_work must be 8-byte real array
      REAL*8 fft_work(nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
      COMPLEX*16 tx(nxl0)
      INTEGER*8 fftw_plan
#endif
c
c      if (FLTIMER) call begin_timer(71)
c
#ifdef COPPER
c for ESSL
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,k,tp,tf,fft_work)
#endif
#ifdef COBALT
c for SCSL
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP&            PRIVATE(i,j,jimag,k,tp,tf,tx,tfw,fftw_plan)
#endif
      do k=1,nz
c compute fft in y-direction
        do i=1,nxl
          call dcopy(nyl,ur(i,1,k),nxl0,tp,1)
#ifdef COPPER
c ESSL
          call drcft(0,tp,tnylhp0,tmplf(1,i,k),nylhp0,nyl,1,1,
     &    1.d0/dble(nyl),ry_ltb_f,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c MKL FFTW3XF
          call dfftw_plan_dft_r2c_1d(fftw_plan,nyl,tp,tfw,FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call zdscal(nylhp,1.d0/dble(nyl),tfw,1)
          call zcopy(nylhp,tfw,1,tmplf(1,i,k),1)
          call dfftw_destroy_plan(fftw_plan)
#endif
        enddo
c compute fft in x-direction
        do j=1,nylh 
#ifdef COPPER
c for ESSL
          call dcft(0,tmplf(j,1,k),nylhp0,1,tf,1,nxl0,nxl,1,1,
     &    1.d0/dble(nxl),cx_ltb_f,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
          call zcopy(nxl,tmplf(j,1,k),nylhp0,tx,1)
          call dfftw_plan_dft_1d(fftw_plan,nxl,tx,tf,FFTW_FORWARD,
     &         FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call zdscal(nxl,1.d0/dble(nxl),tf,1)
          call dfftw_destroy_plan(fftw_plan)
#endif
c
          jimag=nylh+j
          do i=1,nxl
            ur(i,j    ,k) = dreal(tf(i))
            ur(i,jimag,k) = dimag(tf(i))   
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer(71)
      return
      end
c
c
      subroutine fftl_f2r_2d_new(ur)
#include "common.inc"
      include 'fftw3.f'
      REAL*8 ur(nxl0,nyl0,nz)
c
c     2-d fourier to real transform in x and y directions
c     On input:  ur contains data on 1:nxl,    1:nylh,1:nz (real component)
c     On input:  ur contains data on 1:nxl,nylhp:nyl ,1:nz (imag component)
c     On output: ur contains data on 1:nxl,1:nyl ,1:nz
c
      INTEGER i,j,k,jimag
      REAL*8 tp(tnylhp0)
      COMPLEX*16 tfw(nylhp)
      COMPLEX*16 tf(nxl0)
c      COMPLEX*16 tmp(nylhp0,nxl0,nz)!this is a common block now ->tmplf
#ifdef COPPER
c for ESSL
c fft_work must be 8-byte real array
      REAL*8 fft_work(nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
      COMPLEX*16 tx(nxl0)
      INTEGER*8 fftw_plan
#endif
c
c      if (FLTIMER) call begin_timer(76)
c
#ifdef COPPER
c for ESSL
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,k,tf,tp,fft_work)
#endif
#ifdef COBALT
c for SCSL
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP&            PRIVATE(i,j,jimag,k,tf,tp,tx,tfw,fftw_plan)
#endif
      do k=1,nz
c cumpute fft in x-direction
        do j=1,nylh
          jimag=nylh+j
          do i=1,nxl
            tf(i) = dcmplx(ur(i,j,k),ur(i,jimag,k))
          enddo
          tf(nxlhp)=(0.d0,0.d0)
c
#ifdef COPPER
c for ESSL
          call dcft(0,tf,1,nxl0,tmplf(j,1,k),nylhp0,1,nxl,1,-1,1.d0,
     &    cx_ltb_b,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
          call dfftw_plan_dft_1d(fftw_plan,nxl,tf,tx,FFTW_BACKWARD,
     &         FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call dfftw_destroy_plan(fftw_plan)
c          call zzfftm(1,nxl,1,1.d0,tf,nxl0,tx,nxl0,cx_ltb,cx_wk,1)
          call zcopy(nxl,tx,1,tmplf(j,1,k),nylhp0)
#endif
        enddo
c cumpute fft in y-direction
        do i=1,nxl
          tmplf(nylhp,i,k)=(0.d0,0.d0)
c
#ifdef COPPER
c for ESSL
          call dcrft(0,tmplf(1,i,k),nylhp0,tp,tnylhp0,nyl,1,-1,1.d0,
     &    ry_ltb_b,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
          call zcopy(nylhp,tmplf(1,i,k),1,tfw,1)
          call dfftw_plan_dft_c2r_1d(fftw_plan,nyl,tfw,tp,
     &         FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call dfftw_destroy_plan(fftw_plan)
c          call dscal(nyl,1.d0/dble(nxl*nyl),tp,1)
#endif
          call dcopy(nyl,tp,1,ur(i,1,k),nxl0)
        enddo
      enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer(76)
      return
      end
c
c
      subroutine fft_r2f_2d_new_slice(ur)
#include "common.inc"
      include 'fftw3.f'
      REAL*8 ur(nx0,ny0)
c
c     2-d real to fourier transform in x and y directions
c     On input:  ur contains data on 1:nx,1:ny
c     On output: ur contains data on (1:nx),1:nyh   (real component)
c     On output: ur contains data on (1:nx),nyhp:ny (imag component)
c
      INTEGER i,j,jimag
      REAL*8 tp(tnyhp0)
      COMPLEX*16 tfw(nyhp)
      COMPLEX*16 tf(nx0)
c      COMPLEX*16 tmp(nyhp0,nx0,nz)!this is a common block now ->tmpf
#ifdef COPPER
c for ESSL
c fft_work must be 8-byte real array
      REAL*8 fft_work(nfft_aux)
#endif
#ifdef COBALT
c for SCSL
      COMPLEX*16 tx(nx0)
      INTEGER*8 fftw_plan
#endif
c
c      if (FLTIMER) call begin_timer()
c
CC#ifdef COPPER
CCc for ESSL
CCC$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,k,tp,tf,fft_work)
CC#endif
CC#ifdef COBALT
CCc for SCSL
CCC$OMP PARALLEL DO DEFAULT(SHARED),
CCC$OMP&            PRIVATE(i,j,jimag,k,tp,tf,tx,tfw,fftw_plan)
CC#endif
CC      do k=1,nz
c compute fft in y-direction
        do i=1,nx
          call dcopy(ny,ur(i,1),nx0,tp,1)
#ifdef COPPER
c ESSL
          call drcft(0,tp,tnyhp0,tmpf(1,i,1),nyhp0,ny,1,1,1.d0/dble(ny),
     &    ry_tb_f,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c MKL FFTW3XF
          call dfftw_plan_dft_r2c_1d(fftw_plan,ny,tp,tfw,FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call zdscal(nyhp,1.d0/dble(ny),tfw,1)
          call zcopy(nyhp,tfw,1,tmpf(1,i,1),1)
          call dfftw_destroy_plan(fftw_plan)
#endif
        enddo
c compute fft in x-direction
        do j=1,nyh 
#ifdef COPPER
c for ESSL
          call dcft(0,tmpf(j,1,1),nyhp0,1,tf,1,nx0,nx,1,1,1.d0/dble(nx),
     &    cx_tb_f,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c for SCSL
          call zcopy(nx,tmpf(j,1,1),nyhp0,tx,1)
          call dfftw_plan_dft_1d(fftw_plan,nx,tx,tf,FFTW_FORWARD,
     &         FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call zdscal(nx,1.d0/dble(nx),tf,1)
          call dfftw_destroy_plan(fftw_plan)
#endif
c
          jimag=nyh+j
          do i=1,nx
            ur(i,j    ) = dreal(tf(i))
            ur(i,jimag) = dimag(tf(i))
          enddo
        enddo
CCC      enddo
CCCC$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer()
      return
      end
c
c
      subroutine fft_f2r_2d_new_slice(ur)
#include "common.inc"
      include 'fftw3.f'
      REAL*8 ur(nx0,ny0)
c
c     2-d fourier to real transform in x and y directions
c     On input:  ur contains data on 1:nx,   1:nyh,1:nz (real component)
c     On input:  ur contains data on 1:nx,nyhp:ny ,1:nz (imag component)
c     On output: ur contains data on 1:nx,1:ny,1:nz
c
      INTEGER i,j,jimag
      REAL*8 tp(tnyhp0)
      COMPLEX*16 tfw(nyhp)
      COMPLEX*16 tf(nx0)
c      COMPLEX*16 tmp(nyhp0,nx0,nz)!this is a common block now ->tmpf
#ifdef COPPER
c for ESSL
c fft_work must be 8-byte real array
      REAL*8 fft_work(nfft_aux)
#endif
#ifdef COBALT
c for SCSL
      COMPLEX*16 tx(nx0)
      INTEGER*8 fftw_plan
#endif
c
c      if (FLTIMER) call begin_timer()
c
CCC#ifdef COPPER
CCCc for ESSL
CCCC$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jimag,k,tf,tp,fft_work)
CCC#endif
CCC#ifdef COBALT
CCCc for SCSL
CCCC$OMP PARALLEL DO DEFAULT(SHARED),
CCCC$OMP&            PRIVATE(i,j,jimag,k,tf,tp,tx,tfw,fftw_plan)
CCC#endif
CCC      do k=1,nz
c cumpute fft in x-direction
        do j=1,nyh
          jimag=nyh+j
          do i=1,nx
            tf(i) = dcmplx(ur(i,j),ur(i,jimag))
          enddo
          tf(nxhp)=(0.d0,0.d0)
c
#ifdef COPPER
c for ESSL
          call dcft(0,tf,1,nx0,tmpf(j,1,1),nyhp0,1,nx,1,-1,1.d0,
     &    cx_tb_b,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
          call dfftw_plan_dft_1d(fftw_plan,nx,tf,tx,FFTW_BACKWARD,
     &         FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call dfftw_destroy_plan(fftw_plan)
          call zcopy(nx,tx,1,tmpf(j,1,1),nyhp0)
#endif
        enddo
c cumpute fft in y-direction
        do i=1,nx
          tmpf(nyhp,i,1)=(0.d0,0.d0)
c
#ifdef COPPER
c for ESSL
          call dcrft(0,tmpf(1,i,1),nyhp0,tp,tnyhp0,ny,1,-1,1.d0,
     &    ry_tb_b,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
          call zcopy(nyhp,tmpf(1,i,1),1,tfw,1)
          call dfftw_plan_dft_c2r_1d(fftw_plan,ny,tfw,tp,FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call dfftw_destroy_plan(fftw_plan)
c          call dscal(ny,1.d0/dble(nxny),tp,1)
#endif
          call dcopy(ny,tp,1,ur(i,1),nx0)
        enddo
CCC      enddo
CCCC$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer()
      return
      end

      subroutine fftl_r2f_2d_new_slice(ur)
#include "common.inc"
      include 'fftw3.f'
      REAL*8 ur(nxl0,nyl0)
c     2-d real to fourier transform in x-y directions
c     On input: ur contains data on 1:nxl, 1:nyl
c     On output: ur contains data on 1:nxl, 1:nylh (real component)
c     On output: ur contains data on 1:nxl, nylhp:nyl (imag component)
      INTEGER i,j,jimag
      REAL*8 tp(tnylhp0)
      COMPLEX*16 tfw(nylhp)
      COMPLEX*16 tf(nxl0)
      COMPLEX*16 tx(nxl0)
      INTEGER*8 fftw_plan
c     compute fft in y-direction
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,tp,tfw,fftw_plan)
      do i=1,nxl
         call dcopy(nyl,ur(i,1),nxl0,tp,1)
         call dfftw_plan_dft_r2c_1d(fftw_plan,nyl,tp,tfw,FFTW_ESTIMATE)
         call dfftw_execute(fftw_plan)
         call zdscal(nylhp,1.d0/dble(nyl),tfw,1)
         call zcopy(nylhp,tfw,1,tmplf(1,i,1),1)
         call dfftw_destroy_plan(fftw_plan)
      enddo
C$OMP END PARALLEL DO
c     compute fft in x-direction
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(j,i,jimag,tx,tf,fftw_plan)
      do j=1,nylh
         call zcopy(nxl,tmplf(j,1,1),nylhp0,tx,1)
         call dfftw_plan_dft_1d(fftw_plan,nxl,tx,tf,FFTW_FORWARD,
     &        FFTW_ESTIMATE)
         call dfftw_execute(fftw_plan)
         call zdscal(nxl,1.d0/dble(nxl),tf,1)
         call dfftw_destroy_plan(fftw_plan)
         jimag = nylh+j
         do i=1,nxl
            ur(i,j)     = dreal(tf(i))
            ur(i,jimag) = dimag(tf(i))
         enddo
      enddo
C$OMP END PARALLEL DO
      return
      end

      subroutine fftl_f2r_2d_new_slice(ur)
#include "common.inc"
      include 'fftw3.f'
      REAL*8 ur(nxl0,nyl0)
      
      INTEGER i,j,jimag
      REAL*8 tp(tnylhp0)
      COMPLEX*16 tfw(nylhp)
      COMPLEX*16 tf(nxl0)
      COMPLEX*16 tx(nxl0)
      INTEGER*8 fftw_plan
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(j,i,jimag,tf,tx,fftw_plan)
      do j=1,nylh
         jimag = nylh+j
         do i=1,nxl
            tf(i)=dcmplx(ur(i,j),ur(i,jimag))
         enddo
         tf(nxlhp) = (0.d0,0.d0)
         call dfftw_plan_dft_1d(fftw_plan,nxl,tf,tx,FFTW_BACKWARD,
     &        FFTW_ESTIMATE)
         call dfftw_execute(fftw_plan)
         call dfftw_destroy_plan(fftw_plan)
         call zcopy(nxl,tx,1,tmplf(j,1,1),nylhp0)
      enddo
C$OMP END PARALLEL DO
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,tfw,tp,fftw_plan)
      do i=1,nxl
         call zcopy(nylhp,tmplf(1,i,1),1,tfw,1)
         call dfftw_plan_dft_c2r_1d(fftw_plan,nyl,tfw,tp,FFTW_ESTIMATE)
         call dfftw_execute(fftw_plan)
         call dfftw_destroy_plan(fftw_plan)
         call dcopy(nyl,tp,1,ur(i,1),nxl0)
      enddo
C$OMP END PARALLEL DO
      return
      end
