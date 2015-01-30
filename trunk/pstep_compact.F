      SUBROUTINE pstep_compact
#include "common.inc"

c     Solve Poisson equation for pressure in Fourier space
C     for i!=1 and j!=1 the governing equation is Helmholtz
C     equation, which has a unique solution using Neumann-Neumann
c     boundary condition. For i=1 and j=1, the governing equation
c     is Poisson equation, and Neumann-Neumann boundary condition
c     is not a applausible boundary condition.

      real*8 tmp1(nz),tmp2(nz),tmp(nz),tmp3(nz)
      real*8 p1(nz),p2(nz),dp1dz,dp2dz
      real*8 diag(nz),diagl(nz),diagu(nz)
      real*8 coef(11,nz)
      integer i,j,jimag,k,l,info
      real*8 fac
      real*8 alpha,beta
C     External functions
      real*8 ddot
      external ddot

      if(FLTIMER) call begin_timer(65)
c     setting BC for pressure at wall = 0 
      rw(1,1,1)     = 0.d0
      rw(1,nyhp,1)  = 0.d0
      rw(1,1,nz)    = 0.d0
      rw(1,nyhp,nz) = 0.d0
C     
      do k=2,nzm
         tmp1(k) = rw(1,1,k)
         tmp2(k) = rw(1,1,k) - 1.d0
      enddo
      tmp1(1)  = 0.d0
      tmp2(1)  = 0.d0
      tmp1(nz) = 0.d0
      tmp2(nz) = 0.d0
C     Solve equation using Gauss Elimination method
c     p1
      call dcopy(nz,CA2_nodr(2,1),3,diag,1)
      call dcopy(nzm,CA2_nodr(1,2),3,diagu,1)
      call dcopy(nzm,CA2_nodr(3,1),3,diagl,1)
      call dgbmv('N',nz,nz,1,1,1.d0,CA2_nodr,3,tmp1,1,0.d0,tmp,1)
      do k=1,nz
         do l=1,11
            coef(l,k) = CB2_nodr(l,k)
         enddo
      enddo
      call daxpy(nzmm,0.d0,diag(2),1,coef(6,2),11)
      call daxpy(nzmm,0.d0,diagu(2),1,coef(5,3),11)
      call daxpy(nzmm,0.d0,diagl(1),1,coef(7,1),11)
      call GaussEliminate(nz,coef,5,5,11,tmp,p1,info)
      if(info.ne.0) then
         print *, 'Error Solving Pressure Poisson Equation'
         stop
      endif
c     p2
      call dcopy(nz,CA2_nodr(2,1),3,diag,1)
      call dcopy(nzm,CA2_nodr(1,2),3,diagu,1)
      call dcopy(nzm,CA2_nodr(3,1),3,diagl,1)
      call dgbmv('N',nz,nz,1,1,1.d0,CA2_nodr,3,tmp2,1,0.d0,tmp,1)
      do k=1,nz
         do l=1,11
            coef(l,k) = CB2_nodr(l,k)
         enddo
      enddo
      call daxpy(nzmm,0.d0,diag(2),1,coef(6,2),11)
      call daxpy(nzmm,0.d0,diagu(2),1,coef(5,3),11)
      call daxpy(nzmm,0.d0,diagl(1),1,coef(7,1),11)
      call GaussEliminate(nz,coef,5,5,11,tmp,p2,info)
      if(info.ne.0) then
         print *, 'Error Solving Pressure Poisson Equation'
         stop
      endif
C     Compute dp1dz at z = 1 and dp2dz at z = 1
c     use one sided explicit finite difference scheme
c     
      do k=1,nz
         tmp3(k) = 0.d0
      enddo
      tmp3(nz)   = CB2_nono(6,nz)
      tmp3(nzm)  = CB2_nono(7,nzm)
      tmp3(nz-2) = CB2_nono(8,nz-2)
      tmp3(nz-3) = CB2_nono(9,nz-3)
      tmp3(nz-4) = CB2_nono(10,nz-4)
      tmp3(nz-5) = CB2_nono(11,nz-5)
C     
      dp1dz = ddot(nz,tmp3,1,p1,1)
      dp2dz = ddot(nz,tmp3,1,p2,1)
      print *, 'dp1dz,dp2dz'
      print *, dp1dz,dp2dz
C     
      alpha = dp2dz/(dp2dz-dp1dz)
      beta  = dp1dz/(dp1dz-dp2dz)
c      if(mod(itime-istart,oi_cfl).eq.0) write(26,900) alpha,beta
      
c     
      do k=1,nz
         p(1,1,k)    = alpha*p1(k)+beta*p2(k)
         p(1,nyhp,k) = 0.d0
      enddo

      call dcopy(nz,CA2_nono(2,1),3,diag,1)
      call dcopy(nzm,CA2_nono(1,2),3,diagu,1)
      call dcopy(nzm,CA2_nono(3,1),3,diagl,1)

C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP+ PRIVATE(i,j,jimag,k,l,coef,tmp,tmp2,fac,info)
      do j = 1,nyh
         jimag = nyh + j
         do i = 1,nx
            if(i.ne.1.or.j.ne.1) then
               call dcopy(nz,rw(i,j,1),nx0y0,tmp,1)
C     Set boundary condition
               tmp(1)  = 0.d0
               tmp(nz) = 0.d0
               call dgbmv('N',nz,nz,1,1,1.d0,CA2_nono,3,tmp,1,0.d0,
     &              tmp2,1)
               do k=1,nz
                  do l=1,11
                     coef(l,k) = CB2_nono(l,k)
                  enddo
               enddo
               fac = -wavexs(i)-waveys(j)
               call daxpy(nzmm,fac,diag(2),1,coef(6,2),11)
               call daxpy(nzmm,fac,diagu(2),1,coef(5,3),11)
               call daxpy(nzmm,fac,diagl(1),1,coef(7,1),11)
               call GaussEliminate(nz,coef,5,5,11,tmp2,tmp,info)
               if(info.ne.0) then
                  print *, 'Error solving pressure poisson equation'
                  stop
               endif
               call dcopy(nz,tmp,1,p(i,j,1),nx0y0)
C     Imaginary part
               call dcopy(nz,rw(i,jimag,1),nx0y0,tmp,1)
               tmp(1)  = 0.d0
               tmp(nz) = 0.d0
               call dgbmv('N',nz,nz,1,1,1.d0,CA2_nono,3,tmp,1,0.d0,
     &              tmp2,1)
               do k=1,nz
                  do l=1,11
                     coef(l,k) = CB2_nono(l,k)
                  enddo
               enddo
               call daxpy(nzmm,fac,diag(2),1,coef(6,2),11)
               call daxpy(nzmm,fac,diagu(2),1,coef(5,3),11)
               call daxpy(nzmm,fac,diagl(1),1,coef(7,1),11)
               call GaussEliminate(nz,coef,5,5,11,tmp2,tmp,info)
               if(info.ne.0) then
                  print *, 'Error solving pressure poisson equation'
                  stop
               endif
               call dcopy(nz,tmp,1,p(i,jimag,1),nx0y0)
            endif
         enddo
      enddo
C$OMP END PARALLEL DO

 900  FORMAT('pstep: alpha = ',f12.8,' beta = ',f12.8)

      if(FLTIMER) call end_timer(65)
      return
      end
