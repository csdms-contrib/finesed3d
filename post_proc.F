      SUBROUTINE comp_vort(varx,vary,varz)
#include "common.inc"
      REAL*8 varx(nx0,ny0,nz),vary(nx0,ny0,nz),varz(nx0,ny0,nz)
c
      INTEGER i,j,k
      if (FLTIMER) call begin_timer(87)
c
      call ddy(w,varx)
      call ddz(v,ru)
c
      call ddz(u,vary)
      call ddx(w,rv)
c
      call ddx(v,varz)
      call ddy(u,rw)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
       do k=1,nz
        do j=1,ny
          do i=1,nx
          varx(i,j,k)=varx(i,j,k)-ru(i,j,k)
          vary(i,j,k)=vary(i,j,k)-rv(i,j,k)
          varz(i,j,k)=varz(i,j,k)-rw(i,j,k)
          enddo
        enddo
       enddo
C$OMP END PARALLEL DO
c
      if (FLTIMER) call end_timer(87)
      return
      END
c
c
      SUBROUTINE comp_bott_shear(bott_shear_zx,bott_shear_zy)
cz note that k=nz correspond to bottom!!!
cz so use chbd1(nz,1) to begin, and u(i,j,1) to start
#include "common.inc"
      REAL*8 bott_shear_zx(nx0,ny0),bott_shear_zy(nx0,ny0)
c
      INTEGER i,j
c
      REAL*8 ddot
      EXTERNAL ddot
c
      if (FLTIMER) call begin_timer(93)
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do j = 1,ny
        do i = 1,nx
         bott_shear_zx(i,j)= ddot(nz,chbd1(nz,1),nz0,u(i,j,1),nx0*ny0)
c
         bott_shear_zy(i,j)= ddot(nz,chbd1(nz,1),nz0,v(i,j,1),nx0*ny0)
        enddo
      enddo
C$OMP END PARALLEL DO
c
      if (FLTIMER) call end_timer(93)
      return
      END
c
c     
      SUBROUTINE comp_erosion(auxeros)
#include "common.inc"
      REAL*8 auxeros(nx0,ny0)
c     
c     This subroutine computes local erosion assuming variables 
c     in Fourier space
c     Local erosion is returned in physical space 
c     The entrainment function used is Garcia and Parker (1993), JGR-Oceans
c     To get to use this subroutine (w_settling.ne.0)
c
      REAL*8 bott_shear_zx(nx0,ny0),bott_shear_zy(nx0,ny0)
      REAL*8 ushear,ttvm,auxzu5,ushearcc,Ndtt,ddelta,ufree,dmass
      REAL*8 alpha,epsilon,sum,pin,hpin,Eros,rrhof
      INTEGER i,j,k
c
c      if (FLTIMER) call begin_timer()
c
c     compute bottom shear stress
c     here it is assumed that u,v,w,tt are in fourier space
      call comp_bott_shear(bott_shear_zx,bott_shear_zy)
c     transform bottom shear to physical space
      call fft_f2r_2d_new_slice(bott_shear_zx)
      call fft_f2r_2d_new_slice(bott_shear_zy)
cz critical shear stress=0.4 N/m^2, E0 = 5*10^(-5)kg/(m^2*s)
cz see the physical value in Mulder and Udink (1991)
cz      Ndtt = 1.088d-3
cz      ddelta = 1.784d-3
      ufree = 5.6d-1
cz      Eros = Ero/(2.65d3*Ndtt*ufree)
      Eros = Ero
      rrhof = 1.d-3
      ushearcc = sqrt(ushearc*rrhof)/ufree
cz in reality, critical shear stress my be function of eroded mass
cz see Sanford and Maa, 2001
cz      ushearc = ushearc + sqrt(6.7d-4*(dmass**0.5d0))/ufree
cz      print *, 'updated critical shear vel is', ushearc
cz  Use mean shear rate for stability
cz first transfer the shear rate to Fourier space to get the mean
cz      call fft_r2f_2d_new_slice(bott_shear_zy) !transform to Forier space

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,ushear)
      do j = 1,ny
        do i = 1,nx
           ushear = sqrt((1.d0/re)*sqrt(bott_shear_zy(i,j)**2.d0))
           if(ushear.le.ushearcc) then
                auxeros(i,j) = 0.d0
           else
                auxeros(i,j) =Eros*((ushear/ushearcc)**2.d0-1.d0)
           end if
        enddo
      enddo
C$OMP END PARALLEL DO
c
c      if (FLTIMER) call end_timer()
      return
      END
c
      SUBROUTINE make_output_1d(LU)
#include "common.inc"
      INTEGER LU
c
      REAL*8 bott_shear_zx,bott_shear_zy
      INTEGER iout,jout,kout
c
      REAL*8 ddot
      EXTERNAL ddot
c
      if (FLTIMER) call begin_timer(96)
c
      iout=nxh+1
      jout=nyh+1+nyh/2
      kout=nzmm
c compute bottom shear stress
      bott_shear_zx=ddot(nz,chbd1(nz,1),nz0,u(iout,jout,1),nx0*ny0)
      bott_shear_zy=ddot(nz,chbd1(nz,1),nz0,v(iout,jout,1),nx0*ny0)
c
      write(LU,300)rtime,u(iout,jout,kout),v(iout,jout,kout),
     &  bott_shear_zx,bott_shear_zy
c
      if (FLTIMER) call end_timer(96)
 300  FORMAT(1p5e18.10)
      return
      END
