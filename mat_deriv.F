      SUBROUTINE mat_deriv_f_rhs(irk)
#include "common.inc"
      INTEGER irk
c
c     This subroutine computes FFT[ Dui/Dt ] i=x, y, z and
c     return values in Fourier space stored in 
c     DuDt, DvDt and DwDt, respectively
c     This subroutine assumes state variables u, v, w, p and tt 
c     are in Fourier space   
      INTEGER i,j,k
      REAL*8 thrmfac,alpha,saat
c
c      if (FLTIMER) call begin_timer()
c
      thrmfac = 0.d0
      if (FLTHRM) thrmfac = -1.d0*ri   !ra/(16.d0*re*re*pr)!ra < 0 
      alpha = 1.8d3/re
c
c     add -grad(p)+lap(u)/re+buoyancy
c
c     x-component
      call viscxyz_f(u,ru,re,FLSV_NO)
      call ddx_f_new(p,DuDt)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            ru(i,j,k)=-DuDt(i,j,k)+ru(i,j,k)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
c
c     y-component
c NOTE: here we also add gradpbar in Fourier space
      call viscxyz_f(v,rv,re,FLSV_NO)
      call ddy_f_new(p,DvDt)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rv(i,j,k)=-DvDt(i,j,k)+rv(i,j,k)
          enddo
        enddo
c       add gradpbar to zero mode only
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC For the unsteady forcing the pressure gradient should
CC be modified @ each three time step. This is not necessary
CC as the pressure gradient is constant in steady channel flow 
CC however it is necessary in the unsteady flow
CC                                      EMRE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        if(irk.eq.1) saat=rtime+dble(1.d0/3.d0*dt)
        if(irk.eq.2) saat=rtime+dble(0.75d0*dt)
        if(irk.eq.3) saat=rtime+dble(dt)
        rv(1,1,k)=rv(1,1,k)-alpha*sin(alpha*saat) !gradpbar!0.06*sin(0.06*saat)          !20.0*sin(20.0*rtime)   !0.06*cos(0.06*rtime)   !gradpbar
      enddo
C$OMP END PARALLEL DO
c
c     z-component
      call viscxyz_f(w,rw,re,FLSV_NO)
      call ddz(p,DwDt)!ddz computes in real and Fourier spaces
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rw(i,j,k)=-DwDt(i,j,k)+rw(i,j,k)+thrmfac*tt(i,j,k)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
c
c     multiply by 1/(1+Tmax*RR*tt) (in de-aliased way)
c
c     enlarge tt --> ddxul
      call enlarge_f_new(tt,ddxul)
      call fftl_f2r_2d_new(ddxul)
c
c     x-component
      call enlarge_f_new(ru,ddyul)
      call fftl_f2r_2d_new(ddyul)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
        do j=1,nyl
          do i=1,nxl
            ddyul(i,j,k)=ddyul(i,j,k)/(1.d0+Tmax*RR*ddxul(i,j,k))
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
      call fftl_r2f_2d_new(ddyul)
      call reduce_f_new(ddyul,DuDt)
c
c     y-component
      call enlarge_f_new(rv,ddyul)
      call fftl_f2r_2d_new(ddyul)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
        do j=1,nyl
          do i=1,nxl
            ddyul(i,j,k)=ddyul(i,j,k)/(1.d0+Tmax*RR*ddxul(i,j,k))
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
      call fftl_r2f_2d_new(ddyul)
      call reduce_f_new(ddyul,DvDt)
c
c     z-component
      call enlarge_f_new(rw,ddyul)
      call fftl_f2r_2d_new(ddyul)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
        do j=1,nyl
          do i=1,nxl
            ddyul(i,j,k)=ddyul(i,j,k)/(1.d0+Tmax*RR*ddxul(i,j,k))
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
      call fftl_r2f_2d_new(ddyul)
      call reduce_f_new(ddyul,DwDt)
c
c      if (FLTIMER) call end_timer()
      return
      END
c
c
      SUBROUTINE mat_deriv_f_lhs(irk)
#include "common.inc"
      INTEGER irk
c
c     This subroutine computes FFT[ Dui/Dt ] i=x, y, z and
c     return values in Fourier space stored in 
c     DuDt, DvDt and DwDt, respectively
c     This subroutine assumes variables in Fourier space and
c     non-linear terms in ru, rv and rw, respectively.
      INTEGER i,j,k
c
      call ddt(irk,u,ou1,ou2,ou3,rtt)
      call ddt(irk,v,ov1,ov2,ov3,rp)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            DuDt(i,j,k)=rtt(i,j,k)+ru(i,j,k)
            DvDt(i,j,k)=rp( i,j,k)+rv(i,j,k)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
      call ddt(irk,w,ow1,ow2,ow3,rtt)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            DwDt(i,j,k)=rtt(i,j,k)+rw(i,j,k)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
c
      return
      END
