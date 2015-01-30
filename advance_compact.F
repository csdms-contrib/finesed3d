      SUBROUTINE advance(code)
#include "common.inc"
      INTEGER code
c
      REAL*8 cfnl1,cfnl2,cfdif,cfdifpv
      REAL*8 dtcr,idtcr,dtpv,idtpv,prfc1,prfc2,temp,divfac,dt0
      INTEGER i,j,jimag,k,irk,nthreads
      REAL*8 chbd1t(nz),chbd1b(nz),tmp(nz)
      REAL*8 wb,wt
      REAL*8 auxeros(nx0,ny0),auxdep2(nx0,ny0)
      REAL*8 ddot
      REAL*8 ttvm
      REAL*8 igmod
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

cNOTE:tt bc's - Implemented only erosion at bottom
c     ttab(1)=0.d0, ttbb(1)=1.d0 and ttgb(1)=1.d0
ccz to include the erosion, uncommented
        if (FLSETTLE.and.FLEro) then
c     Calculate the settling vel based on the given hindered settling
c     model
c     This is done in the real space (can be a nonlinear function)
      if(FLINERTIA.or.FLHINDER)then
C use w and wp to calculate the effective settling velocity 
cz         igmod=sqrt(gy(1,1)**2+gz(1,1)**2)
cz         igmod=1.d0/igmod
cz in Fourier Space
c First calculate the deposition rate
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
        do j = 1,ny0
           do i = 1,nx0
              wsh(i,j) = 0.d0 
           enddo
        enddo
C$OMP END PARALLEL DO
      wsh(1,1) = -abs(w_settling)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
        do j = 1,ny0
           do i = 1,nx0
              wsh(i,j) = wsh(i,j) + wp(i,j,nz-1)-w(i,j,nz-1)
              ttwb(i,j) = tt(i,j,nz-1)
           enddo
        enddo
C$OMP END PARALLEL DO
        call enlarge_f_slice_new(ttwb,ttwbl)
        call enlarge_f_slice_new(wsh,wshl)
        call fftl_f2r_2d_new_slice(ttwbl)
        call fftl_f2r_2d_new_slice(wshl)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
        do j=1,nyl
           do i=1,nxl
              auxdepl(i,j)=wshl(i,j)*ttwbl(i,j)
           enddo
        enddo
C$OMP END PARALLEL DO
        call fftl_r2f_2d_new_slice(auxdepl)
        call reduce_f_slice_new(auxdepl,auxdep)
      else
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
        do j=1,ny
           do i=1,nx
              auxdep(i,j)=-abs(w_settling)*tt(i,j,nz-1)
           enddo
        enddo
C$OMP END PARALLEL DO
c  End of deposition calculation
       endif
       print *, 'average depos flux is', auxdep(1,1)
C Then calculate the erosion rate
ccz     bottom
           call comp_erosion(auxeros) !erosion returns in real space
cz add the deposition rate at the bottom
           call fft_r2f_2d_new_slice(auxeros) !transform to Forier space
        print *, 'average Erosion flux is', auxeros(1,1)
c End of erosion rate 
c combine dep and eros to get R.H.S. of B.C.
cC$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
c           do j=1,ny
c              do i=1,nx
cc Now we decided to eliminate deposition flux
cc                 auxeros(i,j) = auxeros(i,j) + auxdep(i,j)
c              enddo
c           enddo
cC$OMP END PARALLEL DO
c        print *, 'average net flux is', auxeros(1,1)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,jimag) 
           do j = 1,nyh
              jimag=nyh+j
              do i = 1,nx
c bottom bc
                ttbcber(i,j    ) = -re*sc*auxeros(i,j    )
                ttbcber(i,jimag) = -re*sc*auxeros(i,jimag)
              enddo
           enddo
C$OMP END PARALLEL DO
        endif !FLTHRM.and.ttgb(1).eq.1.d0
cz end of erosion computation

c     compute h of eqns for u, v and w 
        call comp_h(cfnl1,irk)
c
c     compute h of eqn for tt
        if(FLTHRM) call tt_comp_h(cfnl1)
c
c     compute the right hand sides for u, v and w, results stored in u, v and w
        call rhs(cfnl2,cfdif)
c     compute the right hand side for tt, result stored in tt
        if (FLTHRM) call tt_rhs(cfnl2,cfdif)
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
c     solve the helmholtz equations for u, v, w and tt
        if (mod(itime-istart,oi_cfl).eq.0.and.irk.eq.3) !prnthreads = 1
     &       print *, 'helmholz: time, nthreads = ',itime, nthreads
c     solve the helmholtz equation for u
c     use Neumann-Neumann boundary condition for mixing layer simulation
c     use Dirichlet-Dirichlet boundary condition for Channel flow simulation

c$$$        call gen_helmholz(re/cfdif,u,ubctax,ubcbax,
c$$$     &       uve_mimi,uvci_mimi,uvcn_mimi,
c$$$     &       uvprbc1,uvmbct1,uvmbcb1,uvprbcn,uvmbctn,uvmbcbn)
        call helmholtz_compact(re/cfdif,u,ubctax,ubcbax,
     &       CA2_drdr,3,1,1,CB2_drdr,11,5,5)
c     solve the helmholtz equation for v
c$$$        call gen_helmholz(re/cfdif,v,vbctax,vbcbax,
c$$$     &       uve_mimi,uvci_mimi,uvcn_mimi,
c$$$     &       uvprbc1,uvmbct1,uvmbcb1,uvprbcn,uvmbctn,uvmbcbn)
        call helmholtz_compact(re/cfdif,v,vbctax,vbcbax,
     &       CA2_drdr,3,1,1,CB2_drdr,11,5,5)
c     solve the helmholtz equation for w, use Dirichlet-Dirichlet boundary condition
c$$$        call gen_helmholz(re/cfdif,w,wbct,wbcb,
c$$$     &       we_mimi,wci_mimi,wcn_mimi,
c$$$     &       wprbc1,wmbct1,wmbcb1,wprbcn,wmbctn,wmbcbn)
        call helmholtz_compact(re/cfdif,w,wbct,wbcb,
     &       CA2_drdr,3,1,1,CB2_drdr,11,5,5)
c     solve the helmholtz equation for tt, use no-flux boundary condition first
c$$$        if (FLTHRM) call gen_helmholz(re*sc/cfdif,tt,ttbct,ttbcb,
c$$$     &       tte_mimi,ttci_mimi,ttcn_mimi,
c$$$     &       ttprbc1,ttmbct1,ttmbcb1,ttprbcn,ttmbctn,ttmbcbn)

C     Sediment concentration
        if(FLINERTIA.or.FLHINDER) then
C for the boundary condition
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
           do j=1,ny
              do i=1,nx
                 ttot(i,j)=prfc1*ttotn(i,j)-prfc2*ttoto(i,j)
                 ttob(i,j)=prfc1*ttobn(i,j)-prfc2*ttobo(i,j)
              enddo
           enddo
C$OMP END PARALLEL DO
           call enlarge_f_slice_new(ttot,ttotl)
           call enlarge_f_slice_new(ttob,ttobl)
           wb = dwbcb(1,1)
           wt = dwbct(1,1)
c           print *, 'Wb and wt are ',wb,wt
           dwbct(1,1) = 0.d0
           dwbcb(1,1) = 0.d0
           call enlarge_f_slice_new(dwbct,dwbctl)
           call enlarge_f_slice_new(dwbcb,dwbcbl)

           call fftl_f2r_2d_new_slice(ttotl)
           call fftl_f2r_2d_new_slice(ttobl)
           call fftl_f2r_2d_new_slice(dwbctl)
           call fftl_f2r_2d_new_slice(dwbcbl)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
           do j=1,nyl
              do i=1,nxl
                 dwbctl(i,j)=re*sc*dwbctl(i,j)*ttotl(i,j)
                 dwbcbl(i,j)=re*sc*dwbcbl(i,j)*ttobl(i,j)
              enddo
           enddo
C$OMP END PARALLEL DO
           call fftl_r2f_2d_new_slice(dwbctl)
           call fftl_r2f_2d_new_slice(dwbcbl)
           call reduce_f_slice_new(dwbctl,ttbct)
           call reduce_f_slice_new(dwbcbl,ttbcb)
           ttat(1) = -re*sc*(w_settling+wt)
           ttab(1) = -re*sc*(w_settling+wb)
           if(FLEro) then
           ttab(1) = 0.d0
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
           do j=1,ny0
              do i=1,nx0
                 ttbcb(i,j)= 0.d0
              enddo
           enddo
C$OMP END PARALLEL DO
           endif
           ttbt(1) = 1.d0
           ttbb(1) = 1.d0
           CB2_mimi(6,1)  = CB2_nono(6,1)*ttbt(1)+ttat(1)
           CB2_mimi(6,nz) = CB2_nono(6,nz)*ttbb(1)+ttab(1)
           else
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
           do j=1,ny0
              do i=1,nx0
                 ttbct(i,j)= 0.d0
                 ttbcb(i,j)= 0.d0
              enddo
           enddo
C$OMP END PARALLEL DO
           ttat(1) = -re*sc*w_settling
           ttab(1) = -re*sc*w_settling
           if(FLEro) then
           ttab(1) = 0.d0
           endif
           ttbt(1) = 1.d0
           ttbb(1) = 1.d0
           CB2_mimi(6,1)  = CB2_nono(6,1)*ttbt(1)+ttat(1)
           CB2_mimi(6,nz) = CB2_nono(6,nz)*ttbb(1)+ttab(1)
        endif !FLINERTIA.or.FLHINDER
ccz update ttbcb to include erosinal flux
ccz to include the erosion, uncommented
c        if (FLSETTLE.and.ttgb(1).eq.-1.d0) then !compute erosion at
        if (FLSETTLE.and.FLEro) then
c     bottom
c to compensate the settling flux with the deposition flux
cz        ttab(1) = 0.d0
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,jimag) 
           do j = 1,nyh
              jimag=nyh+j
              do i = 1,nx
c bottom bc
                ttbcb(i,j    ) = ttbcb(i,j    ) 
     &                          + ttbcber(i,j    )
                ttbcb(i,jimag) = ttbcb(i,jimag) 
     &                          + ttbcber(i,jimag)
              enddo
           enddo
C$OMP END PARALLEL DO
        endif !FLTHRM.and.ttgb(1).eq.1.d0
ccz end of including erosional flux

        if(FLTHRM) call helmholtz_compact(re*sc/cfdif,tt,ttbct,ttbcb,
     &       CA2_mimi,3,1,1,CB2_mimi,11,5,5)

C After update sediment concentration
        if(FLINERTIA.or.FLHINDER.or.FLEro) then
cC calculate the volume integrated concentration
cz after solving helmholtz equation, tt is in Fourier space
           ttvm = 0.d0
           do k=1,nz
              ttvm = ttvm+wzet(k)*tt(1,1,k)
           enddo        
           print *, 'Volume integrated sediment concentration',ttvm
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
           do j=1,ny
              do i=1,nx
                 ttoto(i,j)=ttotn(i,j)
                 ttobo(i,j)=ttobn(i,j)
cz right, k=1 correspond to top, and k=nz at bottom
                 ttotn(i,j)=tt(i,j,1)
                 ttobn(i,j)=tt(i,j,nz)
              enddo
           enddo
C$OMP END PARALLEL DO
        endif ! FLINERTIA or FLHINDER or FEro
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
            call dscal(nx,idtcr,rw(1,j ,k),1)
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
c           call dcopy(nz,chbd1( 1,1),nz0,chbd1t,1)
c     Here, use explicit finite difference for dp/dz at the top
           do k=1,nz
              chbd1t(k) = 0.d0
           enddo
           chbd1t(1) = CB2_nono(6,1)
           chbd1t(2) = CB2_nono(5,2)
           chbd1t(3) = CB2_nono(4,3)
           chbd1t(4) = CB2_nono(3,4)
           chbd1t(5) = CB2_nono(2,5)
           chbd1t(6) = CB2_nono(1,6)
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
c           call dcopy(nz,chbd1(nz,1),nz0,chbd1b,1)
c     Here, use explicit finite difference for ddz
           do k=1,nz
              chbd1b(k) = 0.d0
           enddo
           chbd1b(nz)   = CB2_nono(6,nz)
           chbd1b(nzm)  = CB2_nono(7,nzm)
           chbd1b(nz-2) = CB2_nono(8,nz-2)
           chbd1b(nz-3) = CB2_nono(9,nz-3)
           chbd1b(nz-4) = CB2_nono(10,nz-4)
           chbd1b(nz-5) = CB2_nono(11,nz-5)
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
        if((FLINERTIA.or.FLHINDER).and.irk.eq.3.and.
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
