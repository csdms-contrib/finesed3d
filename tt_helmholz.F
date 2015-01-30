c  Helmholz solver for mixed-mixed boundary conditions on temperature
c
      SUBROUTINE tt_helmholz(cfdif)
#include "common.inc"
      REAL*8 cfdif
c
      REAL*8 fac,fac2,fac3
      INTEGER i,j,k,jimag
c
      REAL*8 ddot
      EXTERNAL ddot
c
      if (FLTIMER) call begin_timer(61)
c
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,jimag,k,fac,fac2,fac3)
      do j=1,nyh
        jimag=nyh+j
c
        fac3 = waveys(j) + re*sc/4.0/cfdif
        do k=1,nzmm
          fac = 1.d0/(tte_mimi(k) - fac3)
          mlt(k,1,j) = fac
        enddo
        do i=2,nxh
          fac2 = wavexs(i) + fac3
          do k=1,nzmm
            fac = 1.d0/(tte_mimi(k) - fac2)
            mlt(k,i,j)      = fac
            mlt(k,nx+2-i,j) = fac
          enddo
        enddo
c
c     Real component
c     transpose data
        do k = 1,nzmm
          call dcopy(nx,tt(1,j,k+1),1,auxt1(k,1,j),nz0)
        enddo
c     solve
        call tt_helm_aux_2d(auxt1,auxt2,mlt,j)
c     transpose back data
        do i=1,nx
          call dcopy(nzmm,auxt1(1,i,j),1,tt(i,j,2),nx0y0)
        enddo
c     compute top and bottom (Newmann condition)
        do i = 1,nx
          tt(i,j,1)  = ddot(nzmm,ttprbc1,1,auxt1(1,i,j),1)
          tt(i,j,nz) = ddot(nzmm,ttprbcn,1,auxt1(1,i,j),1)
        enddo
c     Imaginary component
c     transpose data
        do k = 1,nzmm
          call dcopy(nx,tt(1,jimag,k+1),1,auxt1(k,1,j),nz0)
        enddo
c     solve
          call tt_helm_aux_2d(auxt1,auxt2,mlt,j)
c     transpose back data
        do i=1,nx
          call dcopy(nzmm,auxt1(1,i,j),1,tt(i,jimag,2),nx0y0)
        enddo
c     compute top and bottom (Newmann condition)
        do i = 1,nx
          tt(i,jimag,1)  = ddot(nzmm,ttprbc1,1,auxt1(1,i,j),1)
          tt(i,jimag,nz) = ddot(nzmm,ttprbcn,1,auxt1(1,i,j),1)
        enddo
c
      enddo
C$OMP END PARALLEL DO
c
      if (FLTIMER) call end_timer(61)
      return
      END
c
c
      SUBROUTINE tt_helm_aux_2d(tmp1,tmp2,mlt_mimi,j)
#include "common.inc"
      REAL*8 mlt_mimi(nz0,nx0,nyh),tmp1(nz0,nx0,nyh),tmp2(nz0,nx0,nyh)
c
      INTEGER i,j,k
c
        call dgemm('N','N',nzmm,nx,nzmm,1.d0,ttci_mimi,nzmm0,
     &    tmp1(1,1,j),nz0,0.d0,tmp2(1,1,j),nz0)
        do i = 1,nx
          do k = 1,nzmm
            tmp2(k,i,j) = tmp2(k,i,j)*mlt_mimi(k,i,j)
          enddo
        enddo
        call dgemm('N','N',nzmm,nx,nzmm,1.d0,ttcn_mimi,nzmm0,
     &    tmp2(1,1,j),nz0,0.d0,tmp1(1,1,j),nz0)
c
      return
      END
