      SUBROUTINE advance(code)
#include "common.inc"
      INTEGER code
c
      REAL*8 cfnl1,cfnl2,cfdif,cfdifpv
      REAL*8 dtcr,idtcr,dtpv,idtpv,prfc1,prfc2,temp,divfac,dt0
      INTEGER i,j,jimag,k,irk,nthreads
      REAL*8 chbd1t(nz),chbd1b(nz),tmp(nz)
      REAL*8 auxeros(nx0,ny0)
      REAL*8 ddot
      EXTERNAL ddot
      INTEGER  idmin,idmax
      EXTERNAL idmin,idmax
c
      INTEGER  OMP_GET_NUM_THREADS
      EXTERNAL OMP_GET_NUM_THREADS
c
      if (FLTIMER.and.code.eq.0) call begin_timer(40)
c
      divfac=0.d0
      if(FLSETTLE.and.(.not.FLOTAU))divfac=Tmax*w_settling
c
      do irk = 1,3
        if (FLTIMER) call begin_timer(42)
c
        cfnl1   = cfnl1v(irk)
        cfnl2   = cfnl2v(irk)
        cfdif   = cfdifv(irk)
        dtcr    = cfdif*2.d0
        idtcr   = 1.d0/dtcr
c
        cfdifpv = cfdifv(mod(irk+1,3)+1)
c     for restarting with different time step
        if (istart.gt.0.and.itime.eq.istart+1.and.irk.eq.1)
     &   cfdifpv= dt_old/8.d0
        dtpv    = cfdifpv*2.d0
        idtpv   = 1.d0/dtpv !=0.5d0/cfdifpv
c
        prfc1  = 1.d0+(dtcr/dtpv)
        prfc2  = dtcr/dtpv
c
c     compute erosion at bottom now that we have velocity
cNOTE:in this way, erosion bc is put explicitely as afunction of the 
c     previous time shear stress 
c     we may need to compute first the velocities and the compute erosion
c     to do this, there is a problem with the volumetric coupling
c     since div(u)=f(tt) 
c     see below before pressure step
c     another idea is to extrapolate the erosion bc
        if (FLSETTLE.and.ttgb(1).eq.-1.d0) then !compute erosion at bottom
c           call comp_erosion(auxeros) !erosion returns in real space
c           call fft_r2f_2d_new_slice(auxeros) !transform to Forier space
        endif
c     compute h of eqn for tt
        if(FLTHRM) call tt_comp_h(cfnl1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC THE ORDER IN THIS PORTION HAS BEEN CHANGED. THE REASON .FOR THAT 
CCC IS TO USE THE UPDATED CONCENTRATION OR TEMPERATURE ON THE 
CCC CALCULATION OF THE RHS OF THE MOMENTUM EQUATION. THEREFORE 
CCC RK3 SCHEME IS FULLY COMPLETED AND THE OBTAINED CONCENTRATION
CCC VALUE IS THEN USED ON THE RHS OF THE MOMENTUM EQUATION
CCC                                   EMRE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
c     compute the right hand side for tt, result stored in tt
        if (FLTHRM) call tt_rhs(cfnl2,cfdif)
        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC  FOLLOWING SOLVES THE HELMHOLZ EQUATION FOR CONCENTRATION OR
CCC  TEMPERATURE.           EMRE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     solve the helmholtz equation for tt
        if (FLTHRM) call gen_helmholz(re*sc/4.0/cfdif,tt,ttbct,ttbcb,
     &       tte_mimi,ttci_mimi,ttcn_mimi,
     &       ttprbc1,ttmbct1,ttmbcb1,ttprbcn,ttmbctn,ttmbcbn)
c
c     compute h of eqns for u, v and w 
        call comp_h(cfnl1,irk)
c

       
c
c     compute the right hand sides for u, v and w, results stored in u, v and w
        call rhs(cfnl2,cfdif)
c     compute the right hand side for tt, result stored in tt
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC BELOW COMMAND IS BLOCKED AND MOVED TO LINE 69
CC                                   EMRE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C        if (FLTHRM) call tt_rhs(cfnl2,cfdif)
c
c     update bc's in Fourier space
c     velocity bc's
        nthreads = 1
C$OMP PARALLEL DEFAULT(SHARED) SHARED(nthreads) 
        nthreads = OMP_GET_NUM_THREADS()
C$OMP DO PRIVATE(i,j,jimag) 
        do j = 1,nyh
          jimag=nyh+j
          do i = 1,nx
c u top bc
           ubctax(i,j    ) = (prfc1*pbxnt(i,j)-prfc2*pbxot(i,j))*dtcr
     &                        + ubct(i,j    )
           ubctax(i,jimag) = (prfc1*pcxnt(i,j)-prfc2*pcxot(i,j))*dtcr
     &                        + ubct(i,jimag)
c u bottom bc
           ubcbax(i,j    ) = (prfc1*pbxnb(i,j)-prfc2*pbxob(i,j))*dtcr
     &                        + ubcb(i,j    )
           ubcbax(i,jimag) = (prfc1*pcxnb(i,j)-prfc2*pcxob(i,j))*dtcr
     &                        + ubcb(i,jimag)
c v top bc
           vbctax(i,j    ) = (prfc1*pbynt(i,j)-prfc2*pbyot(i,j))*dtcr
     &                        + vbct(i,j    )
           vbctax(i,jimag) = (prfc1*pcynt(i,j)-prfc2*pcyot(i,j))*dtcr
     &                        + vbct(i,jimag)
c v bottom bc
           vbcbax(i,j    ) = (prfc1*pbynb(i,j)-prfc2*pbyob(i,j))*dtcr
     &                        + vbcb(i,j    )
           vbcbax(i,jimag) = (prfc1*pcynb(i,j)-prfc2*pcyob(i,j))*dtcr
     &                        + vbcb(i,jimag)
          enddo
        enddo
C$OMP END DO
C$OMP END PARALLEL
c
c     tt bc's
cNOTE:tt bc's - Implemented only erosion at bottom
c     ttab(1)=0.d0, ttbb(1)=1.d0 and ttgb(1)=1.d0
        if (FLSETTLE.and.ttgb(1).eq.-1.d0) then !compute erosion at bottom
c           call comp_erosion(auxeros) !erosion returns in real space
c           call fft_r2f_2d_new_slice(auxeros) !transform to Forier space
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,jimag) 
           do j = 1,nyh
              jimag=nyh+j
              do i = 1,nx
c top bc
c                ttbct(i,j    )=
c                ttbct(i,jimag)=
c bottom bc
c                ttbcb(i,j    )=-re*pr*abs(w_settling)*auxeros(i,j    )
c                ttbcb(i,jimag)=-re*pr*abs(w_settling)*auxeros(i,jimag)
              enddo
           enddo
C$OMP END PARALLEL DO
        endif !FLTHRM.and.ttgb(1).eq.1.d0
c
c     solve the helmholtz equations for u, v, w and tt
        if (mod(itime-istart,oi_cfl).eq.0.and.irk.eq.3) !prnthreads = 1
     &       print *, 'helmholz: time, nthreads = ',itime, nthreads
c     solve the helmholtz equation for u
c        call gen_helmholz(re/cfdif,u,ubctax,ubcbax,
c     &       uve_mimi,uvci_mimi,uvcn_mimi,
c     &       uvprbc1,uvmbct1,uvmbcb1,uvprbcn,uvmbctn,uvmbcbn)
         call tustep(re/cfdif,u)
c     solve the helmholtz equation for v
c        call gen_helmholz(re/cfdif,v,vbctax,vbcbax,
c     &       uve_mimi,uvci_mimi,uvcn_mimi,
c     &       uvprbc1,uvmbct1,uvmbcb1,uvprbcn,uvmbctn,uvmbcbn)
         call tustep(re/cfdif,v)
c     solve the helmholtz equation for w
        call gen_helmholz(re/cfdif,w,wbct,wbcb,
     &       we_mimi,wci_mimi,wcn_mimi,
     &       wprbc1,wmbct1,wmbcb1,wprbcn,wmbctn,wmbcbn)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC BELOW COMMAND IS BLOCKED AND MOVED TO LINE 76
CCC                                      EMRE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
c     solve the helmholtz equation for tt
C        if (FLTHRM) call gen_helmholz(re*sc/4.0/cfdif,tt,ttbct,ttbcb,
C     &       tte_mimi,ttci_mimi,ttcn_mimi,
C     &       ttprbc1,ttmbct1,ttmbcb1,ttprbcn,ttmbctn,ttmbcbn)
c
c     compute the divergence for the pressure step
c     the right hand side of eqn for p is stored in rw
        call ddz(w,rw)!ddz computes in real and Fourier spaces
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,jimag,k,temp),
C$OMP$                             SHARED(idtcr)
        do k=1,nz
          do j=1,nyh
            jimag=nyh+j
c add dw/dz
            call dscal(nx,idtcr,rw(1,j    ,k),1)
            call dscal(nx,idtcr,rw(1,jimag,k),1)
c add dv/dy
            temp = idtcr*wavey(j)
            call daxpy(nx,-temp,v(1,jimag,k),1,rw(1,j,    k),1)
            call daxpy(nx, temp,v(1,j    ,k),1,rw(1,jimag,k),1)
          enddo
c add du/dx
          do i=1,nx
            temp = idtcr*wavex(i)
            call daxpy(nyh,-temp,u(i,nyhp,k),nx0,rw(i,1   ,k),nx0)
            call daxpy(nyh, temp,u(i,1   ,k),nx0,rw(i,nyhp,k),nx0)
          enddo
        enddo
C$OMP END PARALLEL DO
c
c add two-phase flow compressibility effects
        if(divfac.ne.0.d0)then
           call ddz(tt,rtt)
           temp = idtcr*divfac
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j,jimag,k), SHARED(temp)
           do k=1,nz
              do j=1,nyh
                 jimag=nyh+j
c add w_s * dtt/dz
                 call daxpy(nx,temp,rtt(1,j    ,k),1,rw(1,j    ,k),1)
                 call daxpy(nx,temp,rtt(1,jimag,k),1,rw(1,jimag,k),1)
              enddo
           enddo
C$OMP END PARALLEL DO
        endif!divfac
c
c     solve the helmholtz equation for pressure p
c     rhs for pressure equation is in rw
c     pressure returned in Fourie space in variable p
        call pstep
c     compute u(i) from u(*i) and dp/dx
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,k,temp)
        do k=1,nz
          do i=1,nx
            temp = dtcr*wavex(i)
            call daxpy(nyh, temp,p(i,nyhp,k),nx0, u(i,1,k),nx0)
            call daxpy(nyh,-temp,p(i,1,k),nx0,u(i,nyhp,k),nx0)
          enddo
        enddo
C$OMP END PARALLEL DO
c
c     compute v(i) from v(*i) and dp/dy
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j,jimag,k,temp)
        do k=1,nz
          do j=1,nyh
            jimag=nyh+j
c
            temp = dtcr*wavey(j)
            call daxpy(nx, temp,p(1,jimag,k),1, v(1,j,k),1)
            call daxpy(nx,-temp,p(1,j,k),1,v(1,jimag,k),1)
          enddo
        enddo
C$OMP END PARALLEL DO
c
c     compute w(i) from w(*i) and dp/dz
        call ddz(p,ru)!ddz computes in real and Fourier spaces
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j,jimag,k)
        do k = 1,nz
          do j=1,nyh
            jimag=nyh+j
            call daxpy(nx,-dtcr,ru(1,j,k),1, w(1,j,k),1)
            call daxpy(nx,-dtcr,ru(1,jimag,k),1,w(1,jimag,k),1)
          enddo
cjpf  zero out the mean of w (which can grow very slowly over time)
          w( 1,1,k) = 0.d0
          w(1,nyhp,k) = 0.d0
        enddo
C$OMP END PARALLEL DO
c
        if (FLTIMER) call end_timer(42)
c
        if (FLTIMER) call begin_timer(43)
c     update pressure gradient for bc's for tangencial velocity
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,jimag)
        do j = 1,nyh
           jimag=nyh+j
           do i = 1,nx
c top 
              pbxot(i,j) = pbxnt(i,j)
              pcxot(i,j) = pcxnt(i,j)
              pbyot(i,j) = pbynt(i,j)
              pcyot(i,j) = pcynt(i,j)
c bottom
              pbxob(i,j) = pbxnb(i,j)
              pcxob(i,j) = pcxnb(i,j)
              pbyob(i,j) = pbynb(i,j)
              pcyob(i,j) = pcynb(i,j)
           enddo
        enddo
C$OMP END PARALLEL DO
c
c update values of aux pressure at boundaries for aux vel bc's
c if needed ddz is computed, but ddx and ddy are computed below
c top boundary 
        if(uvat.eq.1.d0)then
c no-slip
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j,jimag)
           do j = 1,nyh
              jimag=nyh+j
              call dcopy(nx,p(1,jimag,1),1,pbxnt(1,j),1)
              call dcopy(nx,p(1,jimag,1),1,pbynt(1,j),1)
              call dcopy(nx,p(1,j    ,1),1,pcxnt(1,j),1)
              call dcopy(nx,p(1,j    ,1),1,pcynt(1,j),1)
           enddo
C$OMP END PARALLEL DO
        elseif(uvbt.eq.1.d0)then
c slip
           call dcopy(nz,chbd1( 1,1),nz0,chbd1t,1)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,jimag,tmp)
           do j = 1,nyh
              jimag=nyh+j
              do i=1,nx
c     ddz and place into imag part
                 call dcopy(nz,p(i,j,1),nx0y0,tmp,1)
                 pcynt(i,j) = ddot(nz,chbd1t,1,tmp,1)
                 pcxnt(i,j) = pcynt(i,j)
c     ddz and place into real part
                 call dcopy(nz,p(i,jimag,1),nx0y0,tmp,1)
                 pbynt(i,j) = ddot(nz,chbd1t,1,tmp,1)
                 pbxnt(i,j) = pbynt(i,j)
              enddo
           enddo
C$OMP END PARALLEL DO
        else
           write(26,*)'Wrong BC for velocity for top boundary'
           write( *,*)'Wrong BC for velocity for top boundary'
           stop
        endif
c
c bottom boundary
        if(uvab.eq.1.d0)then
c no-slip
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j,jimag)
           do j = 1,nyh
              jimag=nyh+j
              call dcopy(nx,p(1,jimag,nz),1,pbxnb(1,j),1)
              call dcopy(nx,p(1,jimag,nz),1,pbynb(1,j),1)
              call dcopy(nx,p(1,j    ,nz),1,pcxnb(1,j),1)
              call dcopy(nx,p(1,j    ,nz),1,pcynb(1,j),1)
           enddo
C$OMP END PARALLEL DO
        elseif(uvbb.eq.1.d0)then
c slip
           call dcopy(nz,chbd1(nz,1),nz0,chbd1b,1)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,jimag,tmp)
           do j = 1,nyh
              jimag=nyh+j
              do i=1,nx
c     ddz and place into imag part
                 call dcopy(nz,p(i,j,1),nx0y0,tmp,1)
                 pcynb(i,j) = ddot(nz,chbd1b,1,tmp,1)
                 pcxnb(i,j) = pcynb(i,j)
c     ddz and place into real part
                 call dcopy(nz,p(i,jimag,1),nx0y0,tmp,1)
                 pbynb(i,j) = ddot(nz,chbd1b,1,tmp,1)   
                 pbxnb(i,j) = pbynb(i,j)
              enddo
           enddo
C$OMP END PARALLEL DO
        else
           write(26,*)'Wrong BC for velocity for bottom boundary'
           write( *,*)'Wrong BC for velocity for bottom boundary'
           stop
        endif
c
c     Compute ddx and ddy of aux pressure for aux vel bc's
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,jimag)
        do j = 1,nyh
           jimag=nyh+j
           do i = 1,nx
c top
              pbxnt(i,j) = -wavex(i)*pbxnt(i,j)
              pcxnt(i,j) =  wavex(i)*pcxnt(i,j)
              pbynt(i,j) = -wavey(j)*pbynt(i,j)
              pcynt(i,j) =  wavey(j)*pcynt(i,j)
c bottom
              pbxnb(i,j) = -wavex(i)*pbxnb(i,j)
              pcxnb(i,j) =  wavex(i)*pcxnb(i,j)
              pbynb(i,j) = -wavey(j)*pbynb(i,j)
              pcynb(i,j) =  wavey(j)*pcynb(i,j)
           enddo
        enddo
C$OMP END PARALLEL DO
        if (FLTIMER) call end_timer(43)
c     
c     compute pressure for output
        if(out_press.ne.0.and.
     &       mod((itime-istart),oi_chan).eq.0.and.irk.eq.3)then
           call viscxyz_f(p,rw,re,FLSV_NO)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
           do k=1,nz
              do j=1,ny
                 do i=1,nx
                    p(i,j,k)=p(i,j,k)-cfdif*rw(i,j,k)
                 enddo
              enddo
           enddo
C$OMP END PARALLEL DO          
        endif
c
c     output particles velocity (stored in up, vp and wp)
        if(FLINERTIA.and.irk.eq.3.and.
     &       (mod((itime-istart),oi_chan).eq.0.or.
     &        mod((itime-istart),oi_cfl ).eq.0))then
           call fft_f2r_2d_new(up)
           call fft_f2r_2d_new(vp)
           call fft_f2r_2d_new(wp)
           if(mod((itime-istart),oi_chan).eq.0)
     &          call output_up(up,vp,wp)
           if(mod((itime-istart),oi_cfl).eq.0)then
              write(26,*)'  Disperse phase info:'
              write( *,*)'  Disperse phase info:'
              call courant(dt0,up,vp,wp)
           endif
        endif
c     
      enddo !irk
c     
      if (FLTIMER) call end_timer(40)
c     
      return
      END
