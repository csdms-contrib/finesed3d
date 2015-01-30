      SUBROUTINE comp_h(cfnl1,irk)
#include "common.inc"
      REAL*8 cfnl1
      INTEGER irk
c
c     This subroutine computes h for RK-3 for the three components 
c     of velocity in real space
c     Uses the non-linear terms computed by nltrms
c     with convective and divergence forms alternately,
c     these are stored in ru, rv and rw
c     Use a variable mean pressure gradient in the streamwise direction (y), 
c     set gradpbar = 1.d0 in main (or shell script) to keep force balance
c
c
      REAL*8 thrmfac,saat,alpha
      INTEGER i,j,k
      INTEGER ntaux
c
      if (FLTIMER) call begin_timer(50)
c
      thrmfac = 0.d0
      alpha   = 1.8d3/re
c      alpha = 30.d0
      if (FLTHRM) thrmfac = -1.0*ri    !ra/(16.d0*re*re*pr)!ra < 0 
c
c     add non-linear terms and mean pressure gradient
c     compute de-aliased non-linear terms in Fourier space, 
c     results in ru, rv and rw
      call nltrms_f
c      print *, pi
c
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
        do j=1,ny
          do i=1,nx
           h1(i,j,k)=cfnl1*h1(i,j,k)-dt*ru(i,j,k)
           h2(i,j,k)=cfnl1*h2(i,j,k)-dt*rv(i,j,k)
           h3(i,j,k)=cfnl1*h3(i,j,k)-dt*rw(i,j,k)
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
c        h2(1,1,k)=h2(1,1,k)+dt*gradpbar!-dt*0.06*sin(0.06*saat)       !20.0*sin(20.0*rtime)    !0.06*cos(0.06*rtime)   !+dt*gradpbar
        h2(1,1,k) = h2(1,1,k)-dt*alpha*sin(alpha*saat)
c         h2(1,1,k) = h2(1,1,k) + dt*gradpbar
      enddo
C$OMP END PARALLEL DO
c
c     do not use ru, rv and rw until after computing 
c     the material derivatives because it is used in
c     mat_deriv_f_lhs
c
c     compute material derivative in Fourier space
c     results stored in DuDt, DvDt and DwDt
      if(FLNBSQ.or.FLINERTIA)then
         if(FLMDRHS)then 
C            write(26,*)'Check implementation of DDt (rhs.F)'
C            write(26,*)'which is not suited for sloping channels'
C            write(*,*)'Check implementation of DDt (rhs.F)'
C            write(*,*)'which is not suited for sloping channels'
C            stop
            call mat_deriv_f_rhs(irk)
c            write(*,*)'using RHS'
         else
C            write(26,*)'Check implementation of DDt (rhs.F)'
C            write(26,*)'which is not suited for sloping channels'
C            write(*,*)'Check implementation of DDt (rhs.F)'
C            write(*,*)'which is not suited for sloping channels'
C            stop
            call mat_deriv_f_lhs(irk)
c            write(*,*)'using LHS'
         endif
      endif
c     
c     ru, rv, rw, rtt, rp can be re-used now
c     ul, vl, and wl cannot be reused until tt_nltrms has been called
c     DuDt, DvDt and DwDt have the values of Dui/Dt
c
c     add non-Boussinesq terms
      if(FLNBSQ)then
c     compute de-aliased non-Boussinesq terms (tt*Dui/Dt) in Fourier space 
c     using material derivative (DuDt, DvDt and DwDt), 
c     results in ru, rv and rw
c
c     enlarge tt --> ddxul
         call enlarge_f_new(tt,ddxul)
         call fftl_f2r_2d_new(ddxul)
c     x-component
         call enlarge_f_new(DuDt,ddyul)
         call fftl_f2r_2d_new(ddyul)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
         do k=1,nz
            do j=1,nyl
               do i=1,nxl
                  ddyul(i,j,k)=ddxul(i,j,k)*ddyul(i,j,k)
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO
         call fftl_r2f_2d_new(ddyul)
         call reduce_f_new(ddyul,ru)
c     y-component
         call enlarge_f_new(DvDt,ddyul)
         call fftl_f2r_2d_new(ddyul)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
         do k=1,nz
            do j=1,nyl
               do i=1,nxl
                  ddyul(i,j,k)=ddxul(i,j,k)*ddyul(i,j,k)
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO
         call fftl_r2f_2d_new(ddyul)
         call reduce_f_new(ddyul,rv)
c     z-component
         call enlarge_f_new(DwDt,ddyul)
         call fftl_f2r_2d_new(ddyul)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
         do k=1,nz
            do j=1,nyl
               do i=1,nxl
                  ddyul(i,j,k)=ddxul(i,j,k)*ddyul(i,j,k)
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO
         call fftl_r2f_2d_new(ddyul)
         call reduce_f_new(ddyul,rw)
c
c     add to h
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  h1(i,j,k)=h1(i,j,k)-dt*Tmax*RR*ru(i,j,k)
                  h2(i,j,k)=h2(i,j,k)-dt*Tmax*RR*rv(i,j,k)
                  h3(i,j,k)=h3(i,j,k)-dt*Tmax*RR*rw(i,j,k)
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO
      endif!FLNBSQ
c
c     ru, rv, rw, rtt, rp can be re-used now
c     ul, vl, and wl cannot be reused until tt_nltrms has been called
c
c     add body force terms computed in de-aliased way if gravity varies
      if(FLGVAR)then
c     
         write(26,*)'Check implementation of varying gravity'
         write(26,*)'which is not completely suited' 
         write(26,*)'for sloping channels'
         stop
c     enlarge tt --> ddxul
         call enlarge_f_new(tt,ddxul)
         call fftl_f2r_2d_new(ddxul)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
         do k=1,nz
            do j=1,nyl
               do i=1,nxl
                  ddyul(i,j,k)=ddxul(i,j,k)*gz(j,k)
                  ddxul(i,j,k)=ddxul(i,j,k)*gy(j,k)
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO
         call fftl_r2f_2d_new(ddxul)
         call reduce_f_new(ddxul,rv)
         call fftl_r2f_2d_new(ddyul)
         call reduce_f_new(ddyul,rw)
c
c     add to h
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  h2(i,j,k)=h2(i,j,k)-dt*thrmfac*rv(i,j,k)
                  h3(i,j,k)=h3(i,j,k)-dt*thrmfac*rw(i,j,k)
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO
      else!if gravity does not vay with y or z not need de-aliasing
!         ntaux=NT
!         if(ntaux.eq.0.d0)ntaux=1
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  h2(i,j,k)=h2(i,j,k)-dt*thrmfac*tt(i,j,k)*gy(1,1)!grav fixed
                  h3(i,j,k)=h3(i,j,k)-dt*thrmfac*tt(i,j,k)*gz(1,1)!grav fixed
               enddo
            enddo
C NOTE by mic: here i substract the mean of the tt to avoid press BC problem 
!            h3(1,1,k)=h3(1,1,k)-dt*thrmfac*(-ttm(k)/dble(ntaux))*gz(1,1) 
             h3(1,1,k)=h3(1,1,k)-dt*thrmfac*(-tt(1,1,k))*gz(1,1)
         enddo
         
C$OMP END PARALLEL DO         
      endif
c
c     ru, rv, rw, rtt, rp can be re-used now
c     ul, vl, and wl cannot be reused until tt_nltrms has been called
c
c     add two-way coupling terms
      if(FLTWC)then
c     compute de-aliased terms in Fourier space, 
c     results in ru, rv and rw
c     here call (and programe) 
c        call twctrms_f
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
         do k=1,nz
            do j=1,ny
               do i=1,nx
c                  h1(i,j,k)=h1(i,j,k)-dt*ru(i,j,k)
c                  h2(i,j,k)=h2(i,j,k)-dt*rv(i,j,k)
c                  h3(i,j,k)=h3(i,j,k)-dt*rw(i,j,k)
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO
      endif!FLTWC
c
c     ru, rv, rw, rtt, rp can be re-used now
c     ul, vl, and wl cannot be reused until tt_nltrms has been called
c
      if (FLTIMER) call end_timer(50)
      return
      END
c
c
      SUBROUTINE rhs(cfnl2,cfdif)
#include "common.inc"
      REAL*8 cfnl2,cfdif
c
c     find right hand side of helmholtz eqn for the three components of 
c     velocity, in real space, storing results in u, v and w
c     uses h (advection + forces) computed before in comp_h
c     h has not aliasing arrors here
c
      REAL*8 rhsfac
      INTEGER i,j,k
c
      rhsfac  = -re/cfdif     ! check this one
c
      if (FLTIMER) call begin_timer(53)
c
c     save old velocity values for time derivative if using LHS for DDt
      if((FLNBSQ.or.FLINERTIA).and.(.not.FLMDRHS))then
c        write(*,*)'using LHS, saving old vels'
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  ou3(i,j,k)= ou2(i,j,k)
                  ou2(i,j,k)= ou1(i,j,k)
                  ou1(i,j,k)= u(  i,j,k)
                  ov3(i,j,k)= ov2(i,j,k)
                  ov2(i,j,k)= ov1(i,j,k)
                  ov1(i,j,k)= v(  i,j,k)
                  ow3(i,j,k)= ow2(i,j,k)
                  ow2(i,j,k)= ow1(i,j,k)
                  ow1(i,j,k)= w(  i,j,k)
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO
      endif
c
c     compute 1/re*div grad u for the i-momentum equation and
c     sum the right hand side terms for the i-momentum equation, i=x,y,z
c
c      call viscxyz(u,ru,re)
c      call viscxyz(v,rv,re)
c      call viscxyz(w,rw,re)
      call viscxyz_f(u,ru,re,FLSV_NO)
      call viscxyz_f(v,rv,re,FLSV_NO)
      call viscxyz_f(w,rw,re,FLSV_NO)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            u(i,j,k)=(u(i,j,k)+cfnl2*h1(i,j,k)+cfdif*ru(i,j,k))*rhsfac
            v(i,j,k)=(v(i,j,k)+cfnl2*h2(i,j,k)+cfdif*rv(i,j,k))*rhsfac
            w(i,j,k)=(w(i,j,k)+cfnl2*h3(i,j,k)+cfdif*rw(i,j,k))*rhsfac
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
c
      if (FLTIMER) call end_timer(53)
      return
      END
