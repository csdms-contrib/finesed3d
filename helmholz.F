c  Helmholz solver for velocity
c
      SUBROUTINE helmholz(cfdif,cfdifpv,prnthreads)
#include "common.inc"
      REAL*8 cfdif,cfdifpv
      INTEGER prnthreads
c
      REAL*8 dtcr,dtpv,prfc1,prfc2,fac,fac2,fac3
      INTEGER i,j,jimag,k,nthreads
c
      INTEGER  OMP_GET_NUM_THREADS
      EXTERNAL OMP_GET_NUM_THREADS
c
      if (FLTIMER) call begin_timer(60)
c
      dtcr   = cfdif*2.d0
      dtpv   = cfdifpv*2.d0
      prfc1  = 1.d0+(dtcr/dtpv)
      prfc2  = dtcr/dtpv
c
C$OMP PARALLEL DEFAULT(SHARED) SHARED(nthreads)
      nthreads = OMP_GET_NUM_THREADS()
C$OMP SECTIONS PRIVATE(i,j,jimag,fac,fac2)
      do j = 1,nyh
        jimag=nyh+j
c
c        fac2 = 1.d0 - waveys(j) - re/cfdif
        fac2 =  e_drdr(1) - waveys(j) - re/cfdif
        do i = 1,nx
          fac = (fac2 - wavexs(i))*dtcr
          u(i,j    ,1) = (prfc1*pbxnt(i,j)-prfc2*pbxot(i,j))*fac
          u(i,jimag,1) = (prfc1*pcxnt(i,j)-prfc2*pcxot(i,j))*fac
c
          v(i,j    ,1) = (prfc1*pbynb(i,j)-prfc2*pbyob(i,j))*fac
          v(i,jimag,1) = (prfc1*pcynb(i,j)-prfc2*pcyob(i,j))*fac
c
          w(i,j    ,1) = 0.d0
          w(i,jimag,1) = 0.d0
        enddo
      enddo
C$OMP SECTION
      do j = 1,nyh
        jimag=nyh+j
c
c        fac2 = 1.d0 - waveys(j) - re/cfdif
        fac2 =  e_drdr(nz) - waveys(j) - re/cfdif
        do i = 1,nx
          fac = (fac2 - wavexs(i))*dtcr
          u(i,j    ,nz) = (prfc1*pbxnb(i,j)-prfc2*pbxob(i,j))*fac
          u(i,jimag,nz) = (prfc1*pcxnb(i,j)-prfc2*pcxob(i,j))*fac
c
          v(i,j    ,nz) = (prfc1*pbynt(i,j)-prfc2*pbyot(i,j))*fac
          v(i,jimag,nz) = (prfc1*pcynt(i,j)-prfc2*pcyot(i,j))*fac
c
          w(i,j    ,nz) = 0.d0
          w(i,jimag,nz) = 0.d0
        enddo
      enddo
C$OMP END SECTIONS
c
C$OMP DO PRIVATE(i,j,jimag,k,fac,fac2,fac3)
      do j=1,nyh
        jimag=nyh+j
c
        fac3 = waveys(j) + re/cfdif
        do k=1,nz
          fac = 1.d0/(e_drdr(k) - fac3)
          mlt(k,1,j) = fac
        enddo
        do i=2,nxh
          fac2 = wavexs(i) + fac3
          do k=1,nz
            fac = 1.d0/(e_drdr(k) - fac2)
            mlt(k,i,j)      = fac
            mlt(k,nx+2-i,j) = fac
          enddo
        enddo
cc
cc
c     Solving for u
c     Real component
c     transpose data
        do k = 1,nz
          call dcopy(nx,u(1,j,k),1,auxt1(k,1,j),nz0)
        enddo
c     solve
        call helm_aux_2d(auxt1,auxt2,mlt,j)
c     transpose back data
        do i=1,nx
          call dcopy(nz,auxt1(1,i,j),1,u(i,j,1),nx0y0)
        enddo
c     Imaginary component
c     transpose data
        do k = 1,nz
          call dcopy(nx,u(1,jimag,k),1,auxt1(k,1,j),nz0)
        enddo
c     solve
        call helm_aux_2d(auxt1,auxt2,mlt,j)
c     transpose back data
        do i=1,nx
          call dcopy(nz,auxt1(1,i,j),1,u(i,jimag,1),nx0y0)
        enddo
c
c     Solving for v
c     Real component
c     transpose data
         do k = 1,nz
          call dcopy(nx,v(1,j,k),1,auxt1(k,1,j),nz0)
        enddo
c     solve
        call helm_aux_2d(auxt1,auxt2,mlt,j)
c     transpose back data
        do i=1,nx
          call dcopy(nz,auxt1(1,i,j),1,v(i,j,1),nx0y0)
        enddo
c     Imaginary component
c     transpose data
         do k = 1,nz
          call dcopy(nx,v(1,jimag,k),1,auxt1(k,1,j),nz0)
        enddo
c     solve
        call helm_aux_2d(auxt1,auxt2,mlt,j)
c     transpose back data
        do i=1,nx
          call dcopy(nz,auxt1(1,i,j),1,v(i,jimag,1),nx0y0)
        enddo
c
c     Solving for w
c     Real component
c     transpose data
         do k = 1,nz
          call dcopy(nx,w(1,j,k),1,auxt1(k,1,j),nz0)
        enddo
c     solve
        call helm_aux_2d(auxt1,auxt2,mlt,j)
c     transpose back data
        do i=1,nx
          call dcopy(nz,auxt1(1,i,j),1,w(i,j,1),nx0y0)
        enddo
c     Imaginary component
c     transpose data
         do k = 1,nz
          call dcopy(nx,w(1,jimag,k),1,auxt1(k,1,j),nz0)
        enddo
c     solve
        call helm_aux_2d(auxt1,auxt2,mlt,j)
c     transpose back data
        do i=1,nx
          call dcopy(nz,auxt1(1,i,j),1,w(i,jimag,1),nx0y0)
        enddo
c
      enddo
C$OMP END DO
C$OMP END PARALLEL
c
      if (prnthreads.ne.0) print *, 'helmholz: time, nthreads = ',
     &  itime, nthreads
c
      if (FLTIMER) call end_timer(60)
      return
      END
c
c
      SUBROUTINE helm_aux_2d(tmp1,tmp2,mlt_drdr,j)
#include "common.inc"
      REAL*8 mlt_drdr(nz0,nx0,nyh),tmp1(nz0,nx0,nyh),tmp2(nz0,nx0,nyh)
c
      INTEGER i,j,k
c
        call dgemm('N','N',nz,nx,nz,1.d0,ci_drdr,nz0,tmp1(1,1,j),
     &    nz0,0.d0,tmp2(1,1,j),nz0)
        do i = 1,nx
          do k = 1,nz
            tmp2(k,i,j) = tmp2(k,i,j)*mlt_drdr(k,i,j)
          enddo
        enddo
        call dgemm('N','N',nz,nx,nz,1.d0,cn_drdr,nz0,tmp2(1,1,j),
     &    nz0,0.d0,tmp1(1,1,j),nz0)
c
      return
      END
