      subroutine helmholtz_compact(cfre,var,bct,bcb,A,lda,kla,kua,
     &     B,ldb,klb,kub)
#include "common.inc"
      integer lda,kua,kla,ldb,kub,klb
      real*8 cfre,var(nx0,ny0,nz),bct(nx0,ny0),bcb(nx0,ny0)
      real*8 A(lda,*),B(ldb,*)
      real*8 coef(ldb,nz)
      real*8 tmp(nz),tmp2(nz)
      real*8 diag(nz),diagu(nzm),diagl(nzm)
      real*8 fac
      integer i,j,k,l,info,jimag

      
      call dcopy(nz,A(2,1),lda,diag,1)
      call dcopy(nzm,A(1,2),lda,diagu,1)
      call dcopy(nzm,A(3,1),lda,diagl,1)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,l,coef,tmp)
C$OMP+ PRIVATE(tmp2,fac,jimag,info)
      do i=1,nx
         do j=1,nyh
C     Get the right hand side
            call dcopy(nz,var(i,j,1),nx0y0,tmp,1)
C     Set Boundary condition for tmp
            tmp(1)  = bct(i,j)
            tmp(nz) = bcb(i,j)
c     tmp now is the right hand side of linear system A^{-1}B = b
C     multiply both hand side with A
            call dgbmv('N',nz,nz,kla,kua,1.d0,A,lda,tmp,1,0.d0,tmp2,1)
C     Assemble the coefficient matrix
            do k = 1,nz
               do l = 1,ldb
                  coef(l,k) = B(l,k)
               enddo
            enddo
C     Compute multiplicative factor for constant fac
            fac = -wavexs(i)-waveys(j)-cfre
C     For constant fac, use the following code 
C     Main Diagonal
            call daxpy(nzmm,fac,diag(2),1,coef(6,2),ldb)
C     Upper diagonal
            call daxpy(nzmm,fac,diagu(2),1,coef(5,3),ldb)
C     Lower diagonal
            call daxpy(nzmm,fac,diagl(1),1,coef(7,1),ldb)

c$$$  c     For fac as a function of z, use following code
c$$$            fac2 = -wavexs(i)+wavexs(j)
c$$$            do k = 1,nz
c$$$               fac(k) = fac2 - cfre(k)
c$$$            enddo
c$$$  ind1 = kub+1
c$$$  ind2 = kua+1
c$$$  do k = 2,nzm
c$$$  C     Main diagonal
c$$$  coef(ind1,k) = coef(ind1,k)-A(ind2,k)*fac(k)
c$$$  C     Upper diagonal
c$$$  coef(ind1-1,k+1) = coef(ind1-1,k+1)-A(ind2-1,k+1)*fac(k+1)
c$$$  C     Lower diagonal
c$$$  coef(ind1+1,k-1) = coef(ind1+1,k-1)-A(ind2+1,k-1)*fac(k-1)
c$$$  enddo

C     Call GuassEliminate subroutine to solve the linear system
            call GaussEliminate(nz,coef,klb,kub,ldb,tmp2,tmp,info)
            if(info.ne.0) then
               print *, 'Error solving linear equation'
               stop
            endif
            call dcopy(nz,tmp,1,var(i,j,1),nx0y0)

C     Imaginary part
            jimag = j+nyh
            call dcopy(nz,var(i,jimag,1),nx0y0,tmp,1)
            tmp(1)  = bct(i,jimag)
            tmp(nz) = bcb(i,jimag)
            call dgbmv('N',nz,nz,kla,kua,1.d0,A,lda,tmp,1,0.d0,tmp2,1)
            do k = 1,nz
               do l = 1,ldb
                  coef(l,k) = B(l,k)
               enddo
            enddo
            fac = -wavexs(i)-waveys(j)-cfre
C     For constant fac, use the following code 
C     Main Diagonal
            call daxpy(nzmm,fac,diag(2),1,coef(6,2),ldb)
C     Upper diagonal
            call daxpy(nzmm,fac,diagu(2),1,coef(5,3),ldb)
C     Lower diagonal
            call daxpy(nzmm,fac,diagl(1),1,coef(7,1),ldb)
            call GaussEliminate(nz,coef,klb,kub,ldb,tmp2,tmp,info)
            if(info.ne.0) then
               print *, 'Error solving linear equation'
               stop
            endif
            call dcopy(nz,tmp,1,var(i,jimag,1),nx0y0)
         enddo
      enddo
C$OMP END PARALLEL DO
      return
      end
