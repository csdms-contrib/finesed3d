c
cjpf  Collapsed all subroutines to fix error in compiling with -O2
c
      SUBROUTINE viscxyz(val,dval,reyn)
#include "common.inc"
      include 'fftw3.f'
      REAL*8 val(nx0,ny0,nz),dval(nx0,ny0,nz),reyn
c
c     compute dval = laplacian(val)/reyn
c
      INTEGER i,j,jcomp,k
      REAL*8 visc
      REAL*8 tmp(nz),tmp2(nz)
      REAL*8 diagu(nzm),diagl(nzm),diag(nz)
      REAL*8     rxtmp(tnxhp0),rytmp(tnyhp0)
      COMPLEX*16 cxtmp( nxhp0),cytmp( nyhp0)
#ifdef COPPER
c for ESSL
c  fft_work must be an 8-byte real array
      REAL*8 fft_work(nfft_aux)
#endif
#ifdef COBALT
c for SCSL
      REAL*8 rx_wk(nx+2),ry_wk(ny+2)
      INTEGER*8 fftw_plan
#endif
c
c      if (FLTIMER) call begin_timer()
c
      visc = 1.d0/reyn
c
C$OMP PARALLEL DO DEFAULT(SHARED),SHARED(visc)
C$OMP+ PRIVATE(i,j,tmp,tmp2,diag,diagu,diagl)
      do j=1,ny
         do i=1,nx
C     Get diagonals
            call dcopy(nz,CA2(2,1),3,diag,1)
            call dcopy(nzm,CA2(1,2),3,diagu,1)
            call dcopy(nzm,CA2(3,1),3,diagl,1)
C     Copy data into tmp
            call dcopy(nz,val(i,j,1),nx0y0,tmp,1)
c     Get visc*B*tmp
            call dgbmv('n',nz,nz,5,5,visc,CB2,11,tmp,1,0.d0,
     &           tmp2,1)
C     Solve A*tmp2 = tmp2
            call dgtsv(nz,1,diagl,diag,diagu,tmp2,nz,k)
c     Copy back to dval
            call dcopy(nz,tmp2,1,dval(i,j,1),nx0y0)
         enddo
      enddo
C$OMP END PARALLEL DO
c
C$OMP PARALLEL DEFAULT(SHARED),SHARED(visc)
#ifdef COPPER
c for ESSL
C$OMP DO PRIVATE(i,j,k,rxtmp,cxtmp,rytmp,cytmp,fft_work)
#endif
#ifdef COBALT
c for MKL FFTW3XF
C$OMP DO PRIVATE(i,j,k,rxtmp,cxtmp,rytmp,cytmp,fftw_plan)
#endif
      do k=1,nz
c     old d2dx
        do j = 1,ny
          call dcopy(nx,val(1,j,k),1,rxtmp,1)
c
#ifdef COPPER
c for ESSL
          call drcft(0,rxtmp,tnxhp0,cxtmp,nxhp0,nx,1,1,1.d0/dble(nx),
     &    rx_tb_f,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
          call dfftw_plan_dft_r2c_1d(fftw_plan,nx,rxtmp,cxtmp,
     &         FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call zdscal(nxhp,1.d0/dble(nx),cxtmp,1)
          call dfftw_destroy_plan(fftw_plan)
#endif
c
          do i=1,nxh
c          do i=1,nxhp
c            call zdscal(1,-wavexs(i),cxtmp(i),0)
            cxtmp(i)=dcmplx(-wavexs(i)*cxtmp(i))
          enddo
          cxtmp(nxhp)=(0.d0,0.d0)
c
#ifdef COPPER
c for ESSL
          call dcrft(0,cxtmp,nxhp0,rxtmp,tnxhp0,nx,1,-1,1.d0,
     &    rx_tb_b,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
          call dfftw_plan_dft_c2r_1d(fftw_plan,nx,cxtmp,rxtmp,
     &         FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call dfftw_destroy_plan(fftw_plan)
c          call dscal(nx,1.d0/dble(nx),rxtmp,1)
#endif
c
c     add to laplacian
          call daxpy(nx,visc,rxtmp,1,dval(1,j,k),1)
        enddo
c     end d2dx
c
c     old d2dy
        do i = 1,nx
          call dcopy(ny,val(i,1,k),nx0,rytmp,1)
c
#ifdef COPPER
c for ESSL
          call drcft(0,rytmp,tnyhp0,cytmp,nyhp0,ny,1,1,1.d0/dble(ny),
     &    ry_tb_f,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
          call dfftw_plan_dft_r2c_1d(fftw_plan,ny,rytmp,cytmp,
     &         FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call zdscal(nyhp,1.d0/dble(ny),cytmp,1)
          call dfftw_destroy_plan(fftw_plan)
#endif
c
          do j=1,nyh
c          do j=1,nyhp
c            call zdscal(1,-waveys(j),cytmp(j),0)
            cytmp(j)=dcmplx(-waveys(j)*cytmp(j))
          enddo
          cytmp(nyhp)=(0.d0,0.d0)
c
#ifdef COPPER
c for ESSL
          call dcrft(0,cytmp,nyhp0,rytmp,tnyhp0,ny,1,-1,1.d0,
     &    ry_tb_b,nfft_aux,fft_work,nfft_aux)
#endif
#ifdef COBALT
c for MKL FFTW3XF
          call dfftw_plan_dft_c2r_1d(fftw_plan,ny,cytmp,rytmp,
     &         FFTW_ESTIMATE)
          call dfftw_execute(fftw_plan)
          call dfftw_destroy_plan(fftw_plan)
c          call dscal(ny,1.d0/dble(ny),rytmp,1)
#endif
c
c     add to laplacian
          call daxpy(ny,visc,rytmp,1,dval(i,1,k),nx0)
        enddo
c     end d2dy
c     rezero above nx
        do i=nx+1,nx0
          do j=1,ny
            dval(i,j,k) = 0.d0
          enddo
        enddo
      enddo
C$OMP END DO
C$OMP END PARALLEL
c
c      if (FLTIMER) call end_timer()
      return
      END
c
c
      SUBROUTINE viscxyz_f(val,dval,reyn,flgsv)
#include "common.inc"
      REAL*8 val(nx0,ny0,nz),dval(nx0,ny0,nz),reyn
      REAL*8 tmp(nz),tmp2(nz)
      REAL*8 diag(nz),diagu(nzm),diagl(nzm)
      LOGICAL flgsv !this allows to use the same subroutine for vel and tt
c
c     compute dval = laplacian(val)/reyn in Fourier space
c
      INTEGER i,j,jimag,k
      REAL*8 visc,mltx,mlty
c
c      if (FLTIMER) call begin_timer(90)
c
      visc = 1.d0/reyn
C$OMP PARALLEL DO DEFAULT(SHARED),SHARED(visc)
C$OMP+ PRIVATE(i,j,tmp,tmp2,diag,diagu,diagl)
      do j=1,ny
         do i=1,nx
C     Get diagonals
            call dcopy(nz,CA2(2,1),3,diag,1)
            call dcopy(nzm,CA2(1,2),3,diagu,1)
            call dcopy(nzm,CA2(3,1),3,diagl,1)
C     Copy data into tmp
            call dcopy(nz,val(i,j,1),nx0y0,tmp,1)
c     Get visc*B*tmp
            call dgbmv('n',nz,nz,5,5,visc,CB2,11,tmp,1,0.d0,
     &           tmp2,1)
C     Solve A*tmp2 = tmp2
            call dgtsv(nz,1,diagl,diag,diagu,tmp2,nz,k)
c     Copy back to dval
            call dcopy(nz,tmp2,1,dval(i,j,1),nx0y0)
         enddo
      enddo
C$OMP END PARALLEL DO
c
C$OMP PARALLEL DEFAULT(SHARED),SHARED(visc)
      if(flgsv)then !use anisotropic spectral viscosity
C$OMP DO PRIVATE(i,j,jimag,k,mltx,mlty)
         do k=1,nz
            do j=1,nyh
               jimag=nyh+j
               mlty=1.d0+svy(j)
               do i=1,nx
                  mltx=1.d0+svx(i)
                  dval(i,j    ,k)=dval(i,j    ,k)
     &                            -visc*(mltx*wavexs(i)+mlty*waveys(j))
     &                                  *val(i,j    ,k)
                  dval(i,jimag,k)=dval(i,jimag,k)
     &                            -visc*(mltx*wavexs(i)+mlty*waveys(j))
     &                                  *val(i,jimag,k)
               enddo
            enddo
         enddo
C$OMP END DO
      else
C$OMP DO PRIVATE(i,j,jimag,k)
         do k=1,nz
            do j=1,nyh
               jimag=nyh+j
               do i=1,nx
                  dval(i,j    ,k)=dval(i,j    ,k)
     &                            -visc*(wavexs(i)+waveys(j))
     &                                  *val(i,j    ,k)
                  dval(i,jimag,k)=dval(i,jimag,k)
     &                            -visc*(wavexs(i)+waveys(j))
     &                                  *val(i,jimag,k)
               enddo
            enddo
         enddo
C$OMP END DO
      endif
C$OMP END PARALLEL
c
c      if (FLTIMER) call end_timer(90)
      return
      END
