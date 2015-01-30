c  Common blocks for temporary storage in ffts
      COMPLEX*16 tmpf(nyhp0,nx0,nz),tmplf(nylhp0,nxl0,nz)
      COMMON /tmpfft/ tmpf
      COMMON /tmplfft/ tmplf

c  Support for FFT routines from ESSL library
#ifdef COPPER
      INTEGER nfft_aux
      PARAMETER (nfft_aux = 20000)
c
c  all of these arrays must be 8-byte reals
      REAL*8 cx_tb_f(nfft_aux),cx_tb_b(nfft_aux)
      REAL*8 rx_tb_f(nfft_aux),rx_tb_b(nfft_aux)
      REAL*8 ry_tb_f(nfft_aux),ry_tb_b(nfft_aux)
c
      REAL*8 cx_ltb_f(nfft_aux),cx_ltb_b(nfft_aux)
      REAL*8 ry_ltb_f(nfft_aux),ry_ltb_b(nfft_aux)
c
      COMMON /fft/ cx_tb_f,cx_tb_b,rx_tb_f,rx_tb_b,ry_tb_f,ry_tb_b,
     &	           cx_ltb_f,cx_ltb_b,ry_ltb_f,ry_ltb_b
#endif
c
c  Support for FFT routines from SCSL library
#ifdef COBALT
      INTEGER nxpnfr,nypnfr,tnxpfr
      PARAMETER (nxpnfr=nx+256, nypnfr=ny+256)
      PARAMETER (tnxpfr=2*nx+256)
c
      INTEGER nylpnfr,tnxlpfr
      PARAMETER (nylpnfr=nyl+256)
      PARAMETER (tnxlpfr=2*nxl+256)
c
      INTEGER tnx,tnxl
      PARAMETER(tnx=2*nx,tnxl=2*nxl)
c
c  all of these arrays must be 8-byte reals
      REAL*8 ry_tb_f(nypnfr),ry_tb_b(nypnfr)
      REAL*8 cx_tb(tnxpfr)
      REAL*8 rx_tb_f(nxpnfr),rx_tb_b(nxpnfr)
c
      REAL*8 ry_ltb_f(nylpnfr),ry_ltb_b(nylpnfr)
      REAL*8 cx_ltb(tnxlpfr)
c
      COMMON /fft/  ry_tb_b,ry_tb_f,cx_tb,rx_tb_b,rx_tb_f,
     &                              ry_ltb_b,ry_ltb_f,cx_ltb
#endif
