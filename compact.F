C     This subroutine get the first order and second order differentiation
C     matrix based on polynomial interpolation formula for non-uniform grid

      SUBROUTINE compact
#include "common.inc"
c     alpha is the stretching coefficient used for the non-uniform grid
c     in the form of z=asin(alpha*cos(pi*i/N))/asin(alpha)
c     For alpha = 1, we have uniform grid
c     For alpha = 0, we have Chebyshev collocation points
      INTEGER i,j,k
      INTEGER ku,kl,ind,ind1,ind2
c      REAL*8 alpha,beta
c     Thses are temporary variables used to calculate the coefficient
      REAL*8 a1,a2,a3,a4,a5,b1,b2,b,b1m,b1p,b2m,b2p,D
      REAL*8 a1p,a2p,a1m,a2m,c1m,c2m
      REAL*8 t1,t2,t3,t4,t5,t6,t7
      REAL*8 a11,a12,a21,a22
      REAL*8 c1,c2,c

C     Initialize CB1,CB2,CA1 and CA2
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(j,k)
      do k = 1,nz
         do j = 1,11
            CB1(j,k) = 0.d0
            CB2(j,k) = 0.d0
            CB2_nono(j,k) = 0.d0
            CB2_drdr(j,k) = 0.d0
            CB2_nodr(j,k) = 0.d0
            CB2_mimi(j,k) = 0.d0
         enddo
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(j,k)
      do k = 1,nz
         do j =1,3
            CA1(j,k) = 0.d0
            CA2(j,k) = 0.d0
            CA2_nono(j,k) = 0.d0
            CA2_drdr(j,k) = 0.d0
            CA2_nodr(j,k) = 0.d0
            CA2_mimi(j,k) = 0.d0
         enddo
      enddo
C$OMP END PARALLEL DO

C     Get the first order differentiation matrix
C     such that A du/dz = B u
C     The tridiagonal 6th order scheme is used here
C     For the boundary, 6-th order biased scheme is still in used
C     The non-uniform grid should cure the instability introduced by the 
C     higher order boundary condition
C     A, B are of the general band matrix storage scheme
      ku  = 1
      kl  = 1
C     For interior grid points of matrix A
      ind = ku+1
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(a1,a2,a3,i)
      do i=3,nzm-1
         a1  = (zets(i)-zets(i-2))/(zets(i-1)-zets(i-2))
         a2  = (zets(i)-zets(i+2))/(zets(i-1)-zets(i+2))
         a3  = (zets(i)-zets(i+1))/(zets(i-1)-zets(i+1))
         CA1(ind+1,i-1) = a1*a2*a3**2
         CA1(ind,i)   = 1.d0
         a1  = (zets(i)-zets(i-2))/(zets(i+1)-zets(i-2))
         a2  = (zets(i)-zets(i+2))/(zets(i+1)-zets(i+2))
         a3  = (zets(i)-zets(i-1))/(zets(i+1)-zets(i-1))
         CA1(ind-1,i+1) = a1*a2*a3**2
      enddo
C$OMP END PARALLEL DO

C     For interior grid points of matrix B
      kl  = 5
      ku  = 5
      ind = ku+1
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,a1,a2,a3,a4)
      do i = 3,nzm-1
         a1  = (zets(i)-zets(i-1))/(zets(i-2)-zets(i-1))
         a2  = (zets(i)-zets(i+1))/(zets(i-2)-zets(i+1))
         a3  = (zets(i)-zets(i+2))/((zets(i-2)-zets(i))
     &        *(zets(i-2)-zets(i+2)))
         CB1(ind+2,i-2) = (a1*a2)**2*a3
         a1  = (zets(i)-zets(i+1))/(zets(i-1)-zets(i+1))
         a2  = (zets(i)-zets(i-2))/(zets(i-1)-zets(i-2))
         a3  = (zets(i)-zets(i+2))/(zets(i-1)-zets(i+2))
         a4  = 2.d0/(zets(i-1)-zets(i)) + 2.d0/(zets(i-1)-zets(i+1))
     &        +1.d0/(zets(i-1)-zets(i-2))+1.d0/(zets(i-1)-zets(i+2))
         CB1(ind+1,i-1) = a1**2*a2*a3*a4
         CB1(ind,i)   =1.d0/(zets(i)-zets(i-2))+1.d0/(zets(i)-zets(i+2))
     &               + 2.d0/(zets(i)-zets(i-1))+2.d0/(zets(i)-zets(i+1))
         a1  = (zets(i)-zets(i-1))/(zets(i+1)-zets(i-1))
         a2  = (zets(i)-zets(i+2))/(zets(i+1)-zets(i+2))
         a3  = (zets(i)-zets(i-2))/(zets(i+1)-zets(i-2))
         a4  = 2.d0/(zets(i+1)-zets(i))+2.d0/(zets(i+1)-zets(i-1))
     &        +1.d0/(zets(i+1)-zets(i+2))+1.d0/(zets(i+1)-zets(i-2))
         CB1(ind-1,i+1) = a1**2*a2*a3*a4
         a1  = (zets(i)-zets(i-1))/(zets(i+2)-zets(i-1))
         a2  = (zets(i)-zets(i+1))/(zets(i+2)-zets(i+1))
         a3  = (zets(i)-zets(i-2))/((zets(i+2)-zets(i))
     &        *(zets(i+2)-zets(i-2)))
         CB1(ind-2,i+2) = (a1*a2)**2*a3
      enddo
C$OMP END PARALLEL DO

C     For i=N
C$OMP PARALLEL SECTIONS DEFAULT(SHARED)
C$OMP+ PRIVATE(a1,a2,a3,a4,ku,ind)
      ku  = 1
      ind = 1+ku
      CA1(ind,nzm)   = 1.d0
      a1 = (zets(nzm)-zets(nz))/(zets(nzm-1)-zets(nz))
      a2 = (zets(nzm)-zets(nz-3))/(zets(nzm-1)-zets(nz-3))
      a3 = (zets(nzm)-zets(nz-4))/(zets(nzm-1)-zets(nz-4))
      CA1(ind+1,nz-2) = a1**2*a2*a3
      a1 = (zets(nzm)-zets(nzm-1))/(zets(nz)-zets(nzm-1))
      a2 = (zets(nzm)-zets(nz-3))/(zets(nz)-zets(nz-3))
      a3 = (zets(nzm)-zets(nz-4))/(zets(nz)-zets(nz-4))
      CA1(ind-1,nz) = a1**2*a2*a3
C$OMP SECTION
      ku  = 5
      ind = ku+1
      a1  = (zets(nzm)-zets(nzm-1))/(zets(nzm-3)-zets(nzm-1)) 
      a2  = (zets(nzm)-zets(nz))/(zets(nzm-3)-zets(nz))
      a3  = (zets(nzm)-zets(nz-3))/((zets(nz-4)-zets(nzm))
     &     *(zets(nz-4)-zets(nz-3)))
      CB1(ind+3,nz-4) = (a1*a2)**2*a3
C$OMP SECTION
      ku  = 5
      ind = ku+1
      a1  = (zets(nzm)-zets(nzm-1))/(zets(nzm-2)-zets(nzm-1))
      a2  = (zets(nzm)-zets(nz))/(zets(nzm-2)-zets(nz))
      a3  = (zets(nzm)-zets(nzm-3))/((zets(nzm-2)-zets(nzm))
     &     *(zets(nzm-2)-zets(nzm-3)))
      CB1(ind+2,nz-3) = (a1*a2)**2*a3
C$OMP SECTION
      ku  = 5
      ind = ku+1
      a1  = (zets(nzm)-zets(nzm-2))/(zets(nzm-1)-zets(nz-3))
      a2  = (zets(nzm)-zets(nzm-3))/(zets(nzm-1)-zets(nz-4))
      a3  = (zets(nzm)-zets(nz))/(zets(nzm-1)-zets(nz))
      a4  = 2.d0/(zets(nzm-1)-zets(nzm))+2.d0/(zets(nzm-1)-zets(nz))
     &     +1.d0/(zets(nzm-1)-zets(nz-3))+1.d0/(zets(nzm-1)-zets(nz-4))
      CB1(ind+1,nz-2) = a1*a2*a3**2*a4
C$OMP SECTION
      ku  = 5
      ind = ku+1
      a1  = (zets(nzm)-zets(nz-3))/(zets(nz)-zets(nz-3))
      a2  = (zets(nzm)-zets(nz-4))/(zets(nz)-zets(nz-4))
      a3  = (zets(nzm)-zets(nzm-1))/(zets(nz)-zets(nzm-1))
      a4  = 2.d0/(zets(nz)-zets(nzm))+2.d0/(zets(nz)-zets(nzm-1))
     &     +1.d0/(zets(nz)-zets(nzm-2))+1.d0/(zets(nz)-zets(nzm-3))
      CB1(ind-1,nz) = a1*a2*a3**2*a4
      CB1(ind,nzm)   = 1.d0/(zets(nzm)-zets(nzm-2))
     &     +1.d0/(zets(nzm)-zets(nzm-3)) + 2.d0/(zets(nzm)-zets(nzm-1))
     &     +2.d0/(zets(nzm)-zets(nz))
C$OMP SECTION
C     For i = N+1
      ku  = 1
      ind = ku+1
      CA1(ind,nz)   = 1.d0
      a1  = (zets(nz)-zets(nz-2))/(zets(nzm)-zets(nz-2)) 
      a2  = (zets(nz)-zets(nz-3))/(zets(nzm)-zets(nz-3))
      a3  = (zets(nz)-zets(nz-4))/(zets(nzm)-zets(nz-4))
      a4  = (zets(nz)-zets(nz-5))/(zets(nzm)-zets(nz-5))
      CA1(ind+1,nzm) = a1*a2*a3*a4
C$OMP SECTION
      ku  = 5
      ind = ku+1
      a1  = (zets(nz)-zets(nzm))/(zets(nz-5)-zets(nzm))
      a2  = (zets(nz)-zets(nz-2))/(zets(nz-5)-zets(nz-2))
      a3  = (zets(nz)-zets(nz-3))/(zets(nz-5)-zets(nz-3))
      a4  = (zets(nz)-zets(nz-4))/((zets(nz-5)-zets(nz-4))
     &     *(zets(nz-5)-zets(nz)))
      CB1(ind+5,nz-5) = a1**2*a2*a3*a4
C$OMP SECTION
      ku  = 5
      ind = ku+1
      a1  = (zets(nz)-zets(nzm))/(zets(nz-4)-zets(nzm))
      a2  = (zets(nz)-zets(nz-2))/(zets(nz-4)-zets(nz-2))
      a3  = (zets(nz)-zets(nz-3))/(zets(nz-4)-zets(nz-3))
      a4  = (zets(nz)-zets(nz-5))/((zets(nz-4)-zets(nz-5))
     &     *(zets(nz-4)-zets(nz)))
      CB1(ind+4,nz-4) = a1**2*a2*a3*a4
C$OMP SECTION
      ku  = 5
      ind = ku+1
      a1  = (zets(nz)-zets(nzm))/(zets(nz-3)-zets(nzm))
      a2  = (zets(nz)-zets(nz-2))/(zets(nz-3)-zets(nz-2))
      a3  = (zets(nz)-zets(nz-4))/(zets(nz-3)-zets(nz-4))
      a4  = (zets(nz)-zets(nz-5))/((zets(nz-3)-zets(nz-5))
     &     *(zets(nz-3)-zets(nz)))
      CB1(ind+3,nz-3) = a1**2*a2*a3*a4
C$OMP SECTION
      ku  = 5
      ind = ku+1
      a1  = (zets(nz)-zets(nzm))/(zets(nz-2)-zets(nzm))
      a2  = (zets(nz)-zets(nz-3))/(zets(nz-2)-zets(nz-3))
      a3  = (zets(nz)-zets(nz-4))/(zets(nz-2)-zets(nz-4))
      a4  = (zets(nz)-zets(nz-5))/((zets(nz-2)-zets(nz-5))
     &     *(zets(nz-2)-zets(nz)))
      CB1(ind+2,nz-2) = a1**2*a2*a3*a4
C$OMP SECTION
      ku  = 5
      ind = ku+1
      a1  = (zets(nz)-zets(nz-2))/(zets(nzm)-zets(nz-2))
      a2  = (zets(nz)-zets(nz-3))/(zets(nzm)-zets(nz-3))
      a3  = (zets(nz)-zets(nz-4))/(zets(nzm)-zets(nz-4))
      a4  = (zets(nz)-zets(nz-5))/(zets(nzm)-zets(nz-5))
      CB1(ind+1,nzm) = a1*a2*a3*a4*(2.d0/(zets(nzm)-zets(nz))
     &     +1.d0/(zets(nzm)-zets(nz-2))+1.d0/(zets(nzm)-zets(nz-3))
     &     +1.d0/(zets(nzm)-zets(nz-4))+1.d0/(zets(nzm)-zets(nz-5)))
      CB1(ind,nz)   = 2.d0/(zets(nz)-zets(nzm))
     &     +1.d0/(zets(nz)-zets(nz-2))+1.d0/(zets(nz)-zets(nz-3))
     &     +1.d0/(zets(nz)-zets(nz-4))+1.d0/(zets(nz)-zets(nz-5))
C$OMP END PARALLEL SECTIONS

C     No boundary condition involved, use biased scheme at those boundary
C     points
C     For i = 2
      ku  = 1
      ind = ku+1
      CA1(ind,2)   = 1.d0
      a1 = (zets(2)-zets(4))*(zets(2)-zets(5))
      a2 = (zets(1)-zets(4))*(zets(1)-zets(5))
      a3 = (zets(3)-zets(4))*(zets(3)-zets(5))
      a4 = (zets(2)-zets(3))/(zets(1)-zets(3))
      CA1(ind+1,1) = a1/a2*a4**2
      a4 = (zets(2)-zets(1))/(zets(3)-zets(1))
      CA1(ind-1,3) = a1/a3*a4**2
     
      ku  = 5
      ind = ku+1
      CB1(ind+1,1) = CA1(3,1)*(2.d0/(zets(1)-zets(3))
     &     +2.d0/(zets(1)-zets(2))+1.d0/(zets(1)-zets(4))
     &     +1.d0/(zets(1)-zets(5)))
      CB1(ind-1,3)   = CA1(1,3)*(2.d0/(zets(3)-zets(1))
     &     +2.d0/(zets(3)-zets(2))+1.d0/(zets(3)-zets(4))
     &     +1.d0/(zets(3)-zets(5)))
      CB1(ind,2) = 2.d0/(zets(2)-zets(1))+2.d0/(zets(2)-zets(3))
     &     +1.d0/(zets(2)-zets(4))+1.d0/(zets(2)-zets(5))
      a1 = (zets(2)-zets(1))*(zets(2)-zets(3))
      a1 = a1**2
      a2 = (zets(4)-zets(1))*(zets(4)-zets(3))
      a2 = a2**2
      a3 = (zets(5)-zets(1))*(zets(5)-zets(3))
      a3 = a3**2
      CB1(ind-2,4) = a1/a2*(zets(2)-zets(5))/((zets(4)-zets(2))
     &     *(zets(4)-zets(5))) 
      CB1(ind-3,5) = a1/a3*(zets(2)-zets(4))/((zets(5)-zets(2))
     &     *(zets(5)-zets(4)))

C     For i = 1
      ku  = 1
      ind = ku+1
      CA1(ind,1)   = 1.d0
      a1 = (zets(1)-zets(3))/(zets(2)-zets(3))
      a2 = (zets(1)-zets(4))/(zets(2)-zets(4))
      a3 = (zets(1)-zets(5))/(zets(2)-zets(5))
      a4 = (zets(1)-zets(6))/(zets(2)-zets(6))
      CA1(ind-1,2) = a1*a2*a3*a4

      ku  = 5
      ind = ku+1
C$OMP PARALLEL SECTIONS DEFAULT(SHARED)
C$OMP+ PRIVATE(a1,a2,a3)      
      CB1(ind-1,2) = CA1(1,2)*(1.d0/(zets(2)-zets(3))
     &     +1.d0/(zets(2)-zets(4))+1.d0/(zets(2)-zets(5))
     &     +1.d0/(zets(2)-zets(6))+2.d0/(zets(2)-zets(1)))
      CB1(ind,1)   = 1.d0/(zets(1)-zets(3))+1.d0/(zets(1)-zets(4))
     &     +1.d0/(zets(1)-zets(5))+1.d0/(zets(1)-zets(6))
     &     +2.d0/(zets(1)-zets(2))
C$OMP SECTION
      a1 = (zets(1)-zets(2))/(zets(3)-zets(2))
      a2 = (zets(1)-zets(4))*(zets(1)-zets(5))*(zets(1)-zets(6))
      a3 = (zets(3)-zets(4))*(zets(3)-zets(5))*(zets(3)-zets(6))
      CB1(ind-2,3) = a1**2*a2/(a3*(zets(3)-zets(1)))
C$OMP SECTION
      a1 = (zets(1)-zets(2))/(zets(4)-zets(2))
      a2 = (zets(1)-zets(3))*(zets(1)-zets(5))*(zets(1)-zets(6))
      a3 = (zets(4)-zets(3))*(zets(4)-zets(5))*(zets(4)-zets(6))
      CB1(ind-3,4) = a1**2*a2/(a3*(zets(4)-zets(1)))
C$OMP SECTION
      a1 = (zets(1)-zets(2))/(zets(5)-zets(2))
      a2 = (zets(1)-zets(3))*(zets(1)-zets(4))*(zets(1)-zets(6)) 
      a3 = (zets(5)-zets(3))*(zets(5)-zets(4))*(zets(5)-zets(6))
      CB1(ind-4,5) = a1**2*a2/(a3*(zets(5)-zets(1)))
C$OMP SECTION
      a1 = (zets(1)-zets(2))/(zets(6)-zets(2))
      a2 = (zets(1)-zets(3))*(zets(1)-zets(4))*(zets(1)-zets(5))
      a3 = (zets(6)-zets(3))*(zets(6)-zets(4))*(zets(6)-zets(5))
      CB1(ind-5,6) = a1**2*a2/(a3*(zets(6)-zets(1)))
C$OMP END PARALLEL SECTIONS

C     Get the diagonal form of matrix CA1, where d1 is the main diagonal
c     d1u is the upper diagonal and d1l is the lower diagonal
      call dcopy(nz,CA1(2,1),3,d1,1)
      call dcopy(nzm,CA1(1,2),3,d1u,1)
      call dcopy(nzm,CA1(3,1),3,d1l,1)


C     For second order compact sheme
C     A d^2u/dz^2 = B u
C     For interior points

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,a1,a2,a3,a4,b1,b2,b,D)
C$OMP+ PRIVATE(b1m,b1p,b2m,b2p,c1,c2,c,a1p,a2p,a1m,a2m,t1,t2,t3)
C$OMP+ PRIVATE(t4,t5,t6,a11,a12,a21,a22,t7,c1m,c2m,ku,ind)
      do i = 3,nzm-1
         ku  = 1
         ind = ku+1
         a1  = 1.d0/(zets(i+1)-zets(i-2))+1.d0/(zets(i+1)-zets(i))
     &        +1.d0/(zets(i+1)-zets(i+2))
         a2  = 1.d0/(zets(i-1)-zets(i-2))+1.d0/(zets(i-1)-zets(i))
     &        +1.d0/(zets(i-1)-zets(i+2))
         b1  = (zets(i+1)-zets(i+2))*(zets(i+1)-zets(i-2))
     &        *(zets(i+1)-zets(i))
         b2  = (zets(i-1)-zets(i+2))*(zets(i-1)-zets(i-2))
     &        *(zets(i-1)-zets(i))
         a3  = 2.d0*(3*zets(i+1)-zets(i-2)-zets(i+2)-zets(i))/b1
         a4  = 2.d0*(3*zets(i-1)-zets(i-2)-zets(i+2)-zets(i))/b2
         b   = zets(i+1)-zets(i-1)
         D   = 6.d0+4.d0*b*(a1-a2)-2.d0*b**2*a1*a2
         b1m = -(2.d0/b+a1)*b**2/D
         b1p = -(-2.d0/b+a2)*b**2/D
         b2m = (1.d0+b*a1)/D
         b2p = (1.d0-b*a2)/D
         c1  = (zets(i)-zets(i+2))*(zets(i)-zets(i-2))
         c2  = 2.d0*(2.d0*zets(i)-zets(i+2)-zets(i-2))
         c   = (zets(i)-zets(i+1))*(zets(i)-zets(i-1))
         CA2(ind,i)   = 1.d0
         CA2(ind-1,i+1) =-(2.d0*c1*(b1p
     &        *(2.d0*zets(i)-zets(i+1)-zets(i-1))
     &        +b2p*(zets(i)-zets(i+1))*(3.d0*zets(i)-2.d0*zets(i-1)
     &        -zets(i+1)))+c2*c*(b1p+b2p*(zets(i)-zets(i+1))))/b/b1
         CA2(ind+1,i-1) = (2.d0*c1*(b1m
     &        *(2.d0*zets(i)-zets(i+1)-zets(i-1))
     &        +b2m*(zets(i)-zets(i-1))*(3.d0*zets(i)-2.d0*zets(i+1)
     &        -zets(i-1)))+c2*c*(b1m+b2m*(zets(i)-zets(i-1))))/b/b2
         ku  = 5
         ind = ku+1
         a1p = (-4.d0*a1-2.d0*a2+2.d0*b*(a1*a2-a3)+b**2*a3*a2)/D
         a2p = (4.d0*a1*a2-a3-2.d0*(a1-a2)/b+b*a2*a3)/D
         a1m = (-4.d0*a2-2.d0*a1-2.d0*b*(a1*a2-a4)+b**2*a4*a1)/D
         a2m = (4.d0*a1*a2-a4+2.d0*(a2-a1)/b-b*a4*a1)/D
         CB2(ind-1,i+1) = (2.d0*c1*(1.d0+a1p*(2.d0*zets(i)-zets(i+1)
     &        -zets(i-1))+a2p*(zets(i)-zets(i+1))*(3.d0*zets(i)
     &        -2.d0*zets(i-1)-zets(i+1)))+c2*c*
     &        (1.d0/(zets(i)-zets(i+1))+a1p+a2p*(zets(i)-zets(i+1))))
     &        /b/b1
         CB2(ind+1,i-1) = -(2.d0*c1*(1.d0+a1m*(2.d0*zets(i)-zets(i+1)
     &        -zets(i-1))+a2m*(zets(i)-zets(i-1))*(3.d0*zets(i)
     &        -2.d0*zets(i+1)-zets(i-1)))+c2*c*
     &        (1.d0/(zets(i)-zets(i-1))+a1m+a2m*(zets(i)-zets(i-1))))
     &        /b/b2
         
         t1 = zets(i+1)-zets(i)
         t2 = zets(i-1)-zets(i)
         t3 = (zets(i-1)-zets(i-2))*(zets(i-1)-zets(i+1))
     &        *(zets(i-1)-zets(i+2))
         t4 = (zets(i+1)-zets(i-2))*(zets(i+1)-zets(i-1))
     &        *(zets(i+1)-zets(i+2))
         t5 = 2.d0*((zets(i-1)-zets(i-2))*(zets(i-1)-zets(i+1))
     &        +(zets(i-1)-zets(i+2))*(zets(i-1)-zets(i+1))
     &        +(zets(i-1)-zets(i-2))*(zets(i-1)-zets(i+2)))
         t6 = 2.d0*((zets(i+1)-zets(i-2))*(zets(i+1)-zets(i+2))
     &        +(zets(i+1)-zets(i-2))*(zets(i+1)-zets(i-1))
     &        +(zets(i+1)-zets(i+2))*(zets(i+1)-zets(i-1)))
         a11 = (zets(i-1)-zets(i))*t5+2*t3
         a12 = (zets(i-1)-zets(i))**2*t5+4.d0*(zets(i-1)-zets(i))*t3
         a21 = (zets(i+1)-zets(i))*t6+2.d0*t4
         a22 = (zets(i+1)-zets(i))**2*t6+4.d0*(zets(i+1)-zets(i))*t4
         t7  = a11*a22-a12*a21
         c1m = (a12*t6-a22*t5)/t7
         c2m = (a21*t5-a11*t6)/t7
         t3  = 1.d0/(zets(i)-zets(i-2))+1.d0/(zets(i)-zets(i+2))
         t4  = 2.d0/((zets(i)-zets(i-2))*(zets(i)-zets(i+2)))
         CB2(ind,i) = 2.d0*c2m+2.d0*c1m*(t3-(1.d0/t1+1.d0/t2))+t4
     &        +(2.d0-2.d0*t3*(t1+t2))/t1/t2
         
         t3  = (zets(i-1)-zets(i-2))*(zets(i-1)-zets(i+1))
     &        *(zets(i-1)-zets(i))
         t4  = (zets(i+1)-zets(i-2))*(zets(i+1)-zets(i-1))
     &        *(zets(i+1)-zets(i))
         t5  = 2.d0*((zets(i-1)-zets(i-2))*(zets(i-1)-zets(i+1))
     &        +(zets(i-1)-zets(i))*(zets(i-1)-zets(i+1))
     &        +(zets(i-1)-zets(i-2))*(zets(i-1)-zets(i)))
         t6  = 2.d0*((zets(i+1)-zets(i-2))*(zets(i+1)-zets(i-1))
     &        +(zets(i+1)-zets(i))*(zets(i+1)-zets(i-1))
     &        +(zets(i+1)-zets(i-2))*(zets(i+1)-zets(i)))
         a11 = (zets(i-1)-zets(i+2))*t5+2.d0*t3
         a12 = (zets(i-1)-zets(i+2))**2*t5+4.d0*(zets(i-1)-zets(i+2))*t3
         a21 = (zets(i+1)-zets(i+2))*t6+2.d0*t4
         a22 = (zets(i+1)-zets(i+2))**2*t6+4.d0*(zets(i+1)-zets(i+2))*t4
         t7  = a11*a22-a12*a21
         c1m = (a12*t6-a22*t5)/t7
         c2m = (a21*t5-a11*t6)/t7
         t3  = zets(i+1)-zets(i)
         t4  = zets(i-1)-zets(i)
         t5  = (zets(i)-zets(i-2))/((zets(i+2)-zets(i-2))
     &        *(zets(i+2)-zets(i)))
         t6  = t3*t4/((zets(i+2)-zets(i+1))*(zets(i+2)-zets(i-1)))
         t1  = zets(i+2)-zets(i)
         t7  = 2.d0/((zets(i+2)-zets(i))*(zets(i+2)-zets(i-2)))
         CB2(ind-2,i+2) =(c1m*(1.d0+t1/t3+t1/t4)
     &        -c2m*(2.d0+t1/t3+t1/t4)*t1
     &        -1.d0/t3-1.d0/t4)*2.d0*t5*t6+t6*t7*(1.d0-c1m*t1+c2m*t1**2)
         
         t3  = (zets(i-1)-zets(i+2))*(zets(i-1)-zets(i+1))
     &        *(zets(i-1)-zets(i))
         t4  = (zets(i+1)-zets(i+2))*(zets(i+1)-zets(i-1))
     &        *(zets(i+1)-zets(i))
         t5  = 2.d0*((zets(i-1)-zets(i+2))*(zets(i-1)-zets(i+1))
     &        +(zets(i-1)-zets(i))*(zets(i-1)-zets(i+1))
     &        +(zets(i-1)-zets(i+2))*(zets(i-1)-zets(i)))
         t6  = 2.d0*((zets(i+1)-zets(i+2))*(zets(i+1)-zets(i-1))
     &        +(zets(i+1)-zets(i-1))*(zets(i+1)-zets(i))
     &        +(zets(i+1)-zets(i+2))*(zets(i+1)-zets(i)))
         a11 = (zets(i-1)-zets(i-2))*t5+2.d0*t3
         a12 = (zets(i-1)-zets(i-2))**2*t5+4.d0*(zets(i-1)-zets(i-2))*t3
         a21 = (zets(i+1)-zets(i-2))*t6+2.d0*t4
         a22 = (zets(i+1)-zets(i-2))**2*t6+4.d0*(zets(i+1)-zets(i-2))*t4
         t7  = a11*a22-a12*a21
         c1m = (a12*t6-a22*t5)/t7 
         c2m = (a21*t5-a11*t6)/t7
         t3  = zets(i+1)-zets(i)
         t4  = zets(i-1)-zets(i)
         t1  = zets(i-2)-zets(i)
         t5  = (zets(i)-zets(i+2))/((zets(i-2)-zets(i+2))
     &        *(zets(i-2)-zets(i)))
         t6  = t3*t4/((zets(i-2)-zets(i+1))*(zets(i-2)-zets(i-1)))
         t7  = 2.d0/((zets(i-2)-zets(i))*(zets(i-2)-zets(i+2)))
         CB2(ind+2,i-2) =(c1m*(1.d0+t1/t3+t1/t4)
     &        -c2m*(2.d0+t1/t3+t1/t4)*t1
     &        -1.d0/t3-1.d0/t4)*2.d0*t5*t6+t6*t7*(1.d0-c1m*t1+c2m*t1**2)
      enddo
C$OMP END PARALLEL DO

C     Get CA2,CB2 without boundary treatment, i.e. i=1,2,N-1,N
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(j,k)
      do k=1,nz
         do j=1,3
            CA2_drdr(j,k) = CA2(j,k)
         enddo
      enddo
C$OMP END PARALLEL DO
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(j,k)
      do k=1,nz
         do j=1,11
            CB2_drdr(j,k) = CB2(j,k)
         enddo
      enddo
C$OMP END PARALLEL DO

c     For i=1
      b  = zets(2)-zets(1)
CC$OMP PARALLEL SECTIONS DEFAULT(SHARED),PRIVATE(ku,ind,b2,a1,a2,a3)
      ku  = 1
      ind = ku+1
      CA2(ind,1) = 1.d0
      a1 = (zets(1)-zets(3))*(zets(1)-zets(4))*(zets(1)-zets(5))
     &     *(zets(1)-zets(6))
      a2 = (zets(2)-zets(1))*(zets(2)-zets(3))*(zets(2)-zets(4))
     &     *(zets(2)-zets(5))*(zets(2)-zets(6))*(1.d0/(zets(2)-zets(1))
     &     +1.d0/(zets(2)-zets(3))+1.d0/(zets(2)-zets(4))
     &     +1.d0/(zets(2)-zets(5))+1.d0/(zets(2)-zets(6)))
      a3 = (zets(1)-zets(3))*(zets(1)-zets(4))*(zets(1)-zets(5))
     &     +(zets(1)-zets(3))*(zets(1)-zets(4))*(zets(1)-zets(6))
     &     +(zets(1)-zets(3))*(zets(1)-zets(5))*(zets(1)-zets(6))
     &     +(zets(1)-zets(4))*(zets(1)-zets(5))*(zets(1)-zets(6))
      a3 = 2.d0*a3
      CA2(ind-1,2) = 0.5d0*b*a3/a2-a1/a2

      a4 = (zets(2)-zets(1))*(zets(2)-zets(3))*(zets(2)-zets(4))
     &     +(zets(2)-zets(1))*(zets(2)-zets(3))*(zets(2)-zets(5))
     &     +(zets(2)-zets(1))*(zets(2)-zets(3))*(zets(2)-zets(6))
     &     +(zets(2)-zets(1))*(zets(2)-zets(4))*(zets(2)-zets(5))
     &     +(zets(2)-zets(1))*(zets(2)-zets(4))*(zets(2)-zets(6))
     &     +(zets(2)-zets(1))*(zets(2)-zets(5))*(zets(2)-zets(6))
     &     +(zets(2)-zets(3))*(zets(2)-zets(4))*(zets(2)-zets(5))
     &     +(zets(2)-zets(3))*(zets(2)-zets(4))*(zets(2)-zets(6))
     &     +(zets(2)-zets(3))*(zets(2)-zets(5))*(zets(2)-zets(6))
     &     +(zets(2)-zets(4))*(zets(2)-zets(5))*(zets(2)-zets(6))
      a4 = 2.d0*a4
      a5 = (zets(2)-zets(1))*(zets(2)-zets(3))*(zets(2)-zets(4))
     &     *(zets(2)-zets(5))*(zets(2)-zets(6))
      ku = 5
      ind= ku+1
      CB2(ind-1,2) = a3/a5+0.5d0*b*a3*a4/(a2*a5)-a1*a4/(a2*a5)
CC$OMP SECTION
      ku  = 5
      ind = ku+1
      a1  = 1.d0/(zets(1)-zets(3))+1.d0/(zets(1)-zets(4))
     &     +1.d0/(zets(1)-zets(5))+1.d0/(zets(1)-zets(6))
      a2  = (zets(2)-zets(3))*(zets(2)-zets(4))*(zets(2)-zets(5))
     &     *(zets(2)-zets(6))/((zets(1)-zets(3))*(zets(1)-zets(4))
     &     *(zets(1)-zets(5))*(zets(1)-zets(6)))
      a3  = a2*(1.d0/(zets(2)-zets(3))+1.d0/(zets(2)-zets(4))
     &     +1.d0/(zets(2)-zets(5))+1.d0/(zets(2)-zets(6)))
      a4  = 1.d0/((zets(1)-zets(3))*(zets(1)-zets(4)))
     &     +1.d0/((zets(1)-zets(3))*(zets(1)-zets(5)))
     &     +1.d0/((zets(1)-zets(3))*(zets(1)-zets(6)))
     &     +1.d0/((zets(1)-zets(4))*(zets(1)-zets(5)))
     &     +1.d0/((zets(1)-zets(4))*(zets(1)-zets(6)))
     &     +1.d0/((zets(1)-zets(5))*(zets(1)-zets(6)))
      a4  = 2.d0*a4
      CB2(ind,1) = a4-2.d0*a1*(a2+2.d0*b*a3)/(b*(a2+b*a3))
     &     +2.d0*a3/(b*a2+b**2*a3)
CC$OMP SECTION
      ku  = 5
      ind = ku+1
      b2  = zets(3)-zets(2)
      a1  = (zets(3)-zets(1))*(zets(3)-zets(4))*(zets(3)-zets(5))
     &     *(zets(3)-zets(6))
      a1  = 1.d0/a1
      a2  = a1*(zets(2)-zets(1))*(zets(2)-zets(4))*(zets(2)-zets(5))
     &     *(zets(2)-zets(6))
      a3  = a2*(1.d0/(zets(2)-zets(1))+1.d0/(zets(2)-zets(4))
     &     +1.d0/(zets(2)-zets(5))+1.d0/(zets(2)-zets(6)))
      a4  = a1*(zets(1)-zets(4))*(zets(1)-zets(5))*(zets(1)-zets(6))
      a5  = 2.d0*((zets(1)-zets(4))*(zets(1)-zets(5))+(zets(1)-zets(4))
     &     *(zets(1)-zets(6))+(zets(1)-zets(5))*(zets(1)-zets(6)))*a1
      CB2(ind-2,3) = a5*(-b*a2-b**2*a3)/(b2*a2-b2**2*a3)
     &     +2.d0*a4*(a2+2.d0*b*a3)/(a2-b2*a3)/b2
CC$OMP SECTION   
      ku  = 5
      ind = ku+1
      b2 = zets(2)-zets(4)
      a1  = (zets(4)-zets(1))*(zets(4)-zets(3))*(zets(4)-zets(5))
     &     *(zets(4)-zets(6))
      a1  = 1.d0/a1
      a2  = a1*(zets(2)-zets(1))*(zets(2)-zets(3))*(zets(2)-zets(5))
     &     *(zets(2)-zets(6))
      a3  = a2*(1.d0/(zets(2)-zets(1))+1.d0/(zets(2)-zets(3))
     &     +1.d0/(zets(2)-zets(5))+1.d0/(zets(2)-zets(6)))
      a4  = a1*(zets(1)-zets(3))*(zets(1)-zets(5))*(zets(1)-zets(6))
      a5  = 2.d0*((zets(1)-zets(3))*(zets(1)-zets(5))+(zets(1)-zets(3))
     &     *(zets(1)-zets(6))+(zets(1)-zets(5))*(zets(1)-zets(6)))*a1
      CB2(ind-3,4) = a5*(-b*a2-b**2*a3)/(-b2*a2-b2**2*a3)
     &     -2.d0*a4*(a2+2.d0*b*a3)/(a2+b2*a3)/b2
CC$OMP SECTION
      ku  = 5
      ind = ku+1
      b2  = zets(2)-zets(5)
      a1  = (zets(5)-zets(1))*(zets(5)-zets(3))*(zets(5)-zets(4))
     &     *(zets(5)-zets(6))
      a1  = 1.d0/a1
      a2  = a1*(zets(2)-zets(1))*(zets(2)-zets(3))*(zets(2)-zets(4))
     &     *(zets(2)-zets(6))
      a3  = a2*(1.d0/(zets(2)-zets(1))+1.d0/(zets(2)-zets(3))
     &     +1.d0/(zets(2)-zets(4))+1.d0/(zets(2)-zets(6)))
      a4  = a1*(zets(1)-zets(3))*(zets(1)-zets(4))*(zets(1)-zets(6))
      a5  = 2.d0*((zets(1)-zets(3))*(zets(1)-zets(4))+(zets(1)-zets(3))
     &     *(zets(1)-zets(6))+(zets(1)-zets(4))*(zets(1)-zets(6)))*a1
      CB2(ind-4,5) = a5*(-b*a2-b**2*a3)/(-b2*a2-b2**2*a3)
     &     -2.d0*a4*(a2+2.d0*b*a3)/(a2+b2*a3)/b2
CC$OMP SECTION
      ku  = 5
      ind = ku+1
      b2  = zets(2)-zets(6)
      a1  = (zets(6)-zets(1))*(zets(6)-zets(3))*(zets(6)-zets(4))
     &     *(zets(6)-zets(5))
      a1  = 1.d0/a1
      a2  = a1*(zets(2)-zets(1))*(zets(2)-zets(3))*(zets(2)-zets(4))
     &     *(zets(2)-zets(5))
      a3  = a2*(1.d0/(zets(2)-zets(1))+1.d0/(zets(2)-zets(3))
     &     +1.d0/(zets(2)-zets(4))+1.d0/(zets(2)-zets(5)))
      a4  = a1*(zets(1)-zets(3))*(zets(1)-zets(4))*(zets(1)-zets(5))
      a5  = 2.d0*((zets(1)-zets(3))*(zets(1)-zets(4))+(zets(1)-zets(3))
     &     *(zets(1)-zets(5))+(zets(1)-zets(4))*(zets(1)-zets(5)))*a1
      CB2(ind-5,6) = a5*(-b*a2-b**2*a3)/(-b2*a2-b2**2*a3)
     &     -2.d0*a4*(a2+2.d0*b*a3)/(a2+b2*a3)/b2
CC$OMP END PARALLEL SECTIONS
      
C     For i=N
      b = zets(nzm)-zets(nz)
CC$OMP PARALLEL SECTIONS DEFAULT(SHARED),
CC$OMP+ PRIVATE(ku,ind,a1,a2,a3,a4,a5,b2)
      ku  = 1
      ind = ku+1
      CA2(ind,nz) = 1.d0
      a1 = (zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-3))
     &     *(zets(nz)-zets(nz-4))*(zets(nz)-zets(nz-5))
      a2 = (zets(nzm)-zets(nz-2))*(zets(nzm)-zets(nz-3))
     &     *(zets(nzm)-zets(nz-4))*(zets(nzm)-zets(nz-5))
     &     *(zets(nzm)-zets(nz))*(1.d0/(zets(nzm)-zets(nz))
     &     +1.d0/(zets(nzm)-zets(nz-2))+1.d0/(zets(nzm)-zets(nz-3))
     &     +1.d0/(zets(nzm)-zets(nz-4))+1.d0/(zets(nzm)-zets(nz-5)))
      a3 = (zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-3))
     &     *(zets(nz)-zets(nz-4))+(zets(nz)-zets(nz-2))
     &     *(zets(nz)-zets(nz-3))*(zets(nz)-zets(nz-5))
     &     +(zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-4))
     &     *(zets(nz)-zets(nz-5))+(zets(nz)-zets(nz-3))
     &     *(zets(nz)-zets(nz-4))*(zets(nz)-zets(nz-5))
      a3 = 2.d0*a3
      CA2(ind+1,nzm) = 0.5d0*b*a3/a2-a1/a2

      ku  = 5
      ind = ku+1
      a4  = (zets(nzm)-zets(nz))*(zets(nzm)-zets(nz-2))
     &     *(zets(nzm)-zets(nz-3))+(zets(nzm)-zets(nz))
     &     *(zets(nzm)-zets(nz-2))*(zets(nzm)-zets(nz-4))
     &     +(zets(nzm)-zets(nz))*(zets(nzm)-zets(nz-2))
     &     *(zets(nzm)-zets(nz-5))+(zets(nzm)-zets(nz))
     &     *(zets(nzm)-zets(nz-3))*(zets(nzm)-zets(nz-4))
     &     +(zets(nzm)-zets(nz))*(zets(nzm)-zets(nz-3))
     &     *(zets(nzm)-zets(nz-5))+(zets(nzm)-zets(nz))
     &     *(zets(nzm)-zets(nz-4))*(zets(nzm)-zets(nz-5))
     &     +(zets(nzm)-zets(nz-2))*(zets(nzm)-zets(nz-3))
     &     *(zets(nzm)-zets(nz-4))+(zets(nzm)-zets(nz-2))
     &     *(zets(nzm)-zets(nz-3))*(zets(nzm)-zets(nz-5))
     &     +(zets(nzm)-zets(nz-2))*(zets(nzm)-zets(nz-4))
     &     *(zets(nzm)-zets(nz-5))+(zets(nzm)-zets(nz-3))
     &     *(zets(nzm)-zets(nz-4))*(zets(nzm)-zets(nz-5))
      a4  = 2.d0*a4
      a5  = (zets(nzm)-zets(nz))*(zets(nzm)-zets(nz-2))
     &     *(zets(nzm)-zets(nz-3))*(zets(nzm)-zets(nz-4))
     &     *(zets(nzm)-zets(nz-5))
      CB2(ind+1,nzm) = a3/a5+0.5d0*b*a3*a4/(a2*a5)-a1*a4/(a2*a5)
CC$OMP SECTION
      ku  = 5
      ind = ku+1
      a1  = (zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-3))
     &     *(zets(nz)-zets(nz-4))*(zets(nz)-zets(nz-5))
      a1  = 1.d0/a1
      a2  = a1*(zets(nzm)-zets(nz-2))*(zets(nzm)-zets(nz-3))
     &     *(zets(nzm)-zets(nz-4))*(zets(nzm)-zets(nz-5))
      a3  = a2*(1.d0/(zets(nzm)-zets(nz-2))+1.d0/(zets(nzm)-zets(nz-3))
     &     +1.d0/(zets(nzm)-zets(nz-4))+1.d0/(zets(nzm)-zets(nz-5)))
      a4  = (zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-3))
     &     +(zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-4))
     &     +(zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-5))
     &     +(zets(nz)-zets(nz-3))*(zets(nz)-zets(nz-4))
     &     +(zets(nz)-zets(nz-3))*(zets(nz)-zets(nz-5))
     &     +(zets(nz)-zets(nz-4))*(zets(nz)-zets(nz-5))
      a4  = 2.d0*a4*a1
      a1  = 1.d0/(zets(nz)-zets(nz-2))+1.d0/(zets(nz)-zets(nz-3))
     &     +1.d0/(zets(nz)-zets(nz-4))+1.d0/(zets(nz)-zets(nz-5))
      CB2(ind,nz) = a4-2.d0*a1*(a2+2.d0*b*a3)/(b*(a2+b*a3))
     &     +2.d0*a3/(b*a2+b**2*a3)
CC$OMP SECTION
      ku  = 5
      ind = ku+1
      b2  = zets(nzm)-zets(nz-2)
      a1  = (zets(nz-2)-zets(nz))*(zets(nz-2)-zets(nz-3))
     &     *(zets(nz-2)-zets(nz-4))*(zets(nz-2)-zets(nz-5))
      a1  = 1.d0/a1
      a2  = a1*(zets(nzm)-zets(nz))*(zets(nzm)-zets(nz-3))
     &     *(zets(nzm)-zets(nz-4))*(zets(nzm)-zets(nz-5))
      a3  = a2*(1.d0/(zets(nzm)-zets(nz))+1.d0/(zets(nzm)-zets(nz-3))
     &     +1.d0/(zets(nzm)-zets(nz-4))+1.d0/(zets(nzm)-zets(nz-5)))
      a4  = a1*(zets(nz)-zets(nz-3))*(zets(nz)-zets(nz-4))
     &     *(zets(nz)-zets(nz-5))
      a5  = (zets(nz)-zets(nz-3))*(zets(nz)-zets(nz-4))
     &     +(zets(nz)-zets(nz-3))*(zets(nz)-zets(nz-5))
     &     +(zets(nz)-zets(nz-4))*(zets(nz)-zets(nz-5))
      a5  = 2.d0*a5*a1
      CB2(ind+2,nz-2)=a5*(-b*a2-b**2*a3)/(-b2*a2-b2**2*a3)
     &     -2.d0*a4*(a2+2.d0*b*a3)/(a2+b2*a3)/b2
CC$OMP SECTION
      ku  = 5
      ind = ku+1
      b2  = zets(nzm)-zets(nz-3)
      a1  = (zets(nz-3)-zets(nz))*(zets(nz-3)-zets(nz-2))
     &     *(zets(nz-3)-zets(nz-4))*(zets(nz-3)-zets(nz-5))
      a1  = 1.d0/a1
      a2  = a1*(zets(nzm)-zets(nz))*(zets(nzm)-zets(nz-2))
     &     *(zets(nzm)-zets(nz-4))*(zets(nzm)-zets(nz-5))
      a3  = a2*(1.d0/(zets(nzm)-zets(nz))+1.d0/(zets(nzm)-zets(nz-2))
     &     +1.d0/(zets(nzm)-zets(nz-4))+1.d0/(zets(nzm)-zets(nz-5)))
      a4  = a1*(zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-4))
     &     *(zets(nz)-zets(nz-5))
      a5  = (zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-4))
     &     +(zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-5))
     &     +(zets(nz)-zets(nz-4))*(zets(nz)-zets(nz-5))
      a5  = 2.d0*a5*a1
      CB2(ind+3,nz-3)=a5*(-b*a2-b**2*a3)/(-b2*a2-b2**2*a3)
     &     -2.d0*a4*(a2+2.d0*b*a3)/(a2+b2*a3)/b2
CC$OMP SECTION
      ku  = 5
      ind = ku+1
      b2  = zets(nzm)-zets(nz-4)
      a1  = (zets(nz-4)-zets(nz))*(zets(nz-4)-zets(nz-2))
     &     *(zets(nz-4)-zets(nz-3))*(zets(nz-4)-zets(nz-5))
      a1  = 1.d0/a1
      a2  = a1*(zets(nzm)-zets(nz))*(zets(nzm)-zets(nz-2))
     &     *(zets(nzm)-zets(nz-3))*(zets(nzm)-zets(nz-5))
      a3  = a2*(1.d0/(zets(nzm)-zets(nz))+1.d0/(zets(nzm)-zets(nz-2))
     &     +1.d0/(zets(nzm)-zets(nz-3))+1.d0/(zets(nzm)-zets(nz-5)))
      a4  = a1*(zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-3))
     &     *(zets(nz)-zets(nz-5))
      a5  = (zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-3))
     &     +(zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-5))
     &     +(zets(nz)-zets(nz-3))*(zets(nz)-zets(nz-5))
      a5  = 2.d0*a5*a1
      CB2(ind+4,nz-4)=a5*(-b*a2-b**2*a3)/(-b2*a2-b2**2*a3)
     &     -2.d0*a4*(a2+2.d0*b*a3)/(a2+b2*a3)/b2
CC$OMP SECTION
      ku  = 5
      ind = ku+1
      b2  = zets(nzm)-zets(nz-5)
      a1  = (zets(nz-5)-zets(nz))*(zets(nz-5)-zets(nz-2))
     &     *(zets(nz-5)-zets(nz-3))*(zets(nz-5)-zets(nz-4))
      a1  = 1.d0/a1
      a2  = a1*(zets(nzm)-zets(nz))*(zets(nzm)-zets(nz-2))
     &     *(zets(nzm)-zets(nz-3))*(zets(nzm)-zets(nz-4))
      a3  = a2*(1.d0/(zets(nzm)-zets(nz))+1.d0/(zets(nzm)-zets(nz-2))
     &     +1.d0/(zets(nzm)-zets(nz-3))+1.d0/(zets(nzm)-zets(nz-4)))
      a4  = a1*(zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-3))
     &     *(zets(nz)-zets(nz-4))
      a5  = (zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-3))
     &     +(zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-4))
     &     +(zets(nz)-zets(nz-3))*(zets(nz)-zets(nz-4))
      a5  = 2.d0*a5*a1
      CB2(ind+5,nz-5)=a5*(-b*a2-b**2*a3)/(-b2*a2-b2**2*a3)
     &     -2.d0*a4*(a2+2.d0*b*a3)/(a2+b2*a3)/b2
CC$OMP END PARALLEL SECTIONS

C     For i=2
C     Get a1m,a2m,b1m,b2m
CC$OMP PARALLEL SECTIONS DEFAULT(SHARED),PRIVATE(a1m,a2m,a1p,a2p)
CC$OMP+ PRIVATE(b,t1,a1,a2,a3,a4,b1,b2,a11,a12,a21,a22,ind,ku)
      b   = zets(3)-zets(1)
      t1  = (zets(1)-zets(2))*(zets(1)-zets(3))*(zets(1)-zets(4))
     &     *(zets(1)-zets(5))
      t1  = 1.d0/t1
      a1  =  (zets(1)-zets(2))*(zets(1)-zets(3))*(zets(1)-zets(4))
     &     +(zets(1)-zets(2))*(zets(1)-zets(3))*(zets(1)-zets(5))
     &     +(zets(1)-zets(2))*(zets(1)-zets(4))*(zets(1)-zets(5))
     &     +(zets(1)-zets(3))*(zets(1)-zets(4))*(zets(1)-zets(5))
      a1  = a1*t1
      a2  = (zets(1)-zets(2))*(zets(1)-zets(3))
     &     +(zets(1)-zets(2))*(zets(1)-zets(4))
     &     +(zets(1)-zets(2))*(zets(1)-zets(5))
     &     +(zets(1)-zets(3))*(zets(1)-zets(4))
     &     +(zets(1)-zets(3))*(zets(1)-zets(5))
     &     +(zets(1)-zets(4))*(zets(1)-zets(5))
      a2  = 2.d0*a2*t1
      a3  = (zets(3)-zets(2))*(zets(3)-zets(4))*(zets(3)-zets(5))*t1
      a4  = (zets(3)-zets(2))*(zets(3)-zets(4))
     &     +(zets(3)-zets(2))*(zets(3)-zets(5))
     &     +(zets(3)-zets(4))*(zets(3)-zets(5))
      a4  = 2.d0*a4*t1
C     Solve for a1m,a2m
      a11 = 2.d0*a1
      a12 = 2.d0
      a21 = b*a4+2.d0*a3
      a22 = b**2*a4+4.d0*b*a3
      b1  = -a2
      b2  = -a4
      a1m = (b1*a22-b2*a12)/(a11*a22-a12*a21)
      a2m = (a11*b2-a21*b1)/(a11*a22-a12*a21)
c     Solve for b1m,b2m
      a11 = 2.d0*a1
      a12 = 2.d0
      b1  = 1.d0
      b2  = 0.d0
      b1m = (b1*a22-b2*a12)/(a11*a22-a12*a21)
      b2m = (a11*b2-a21*b1)/(a11*a22-a12*a21) 
      b1  = (zets(1)-zets(2))*(zets(1)-zets(4))*(zets(1)-zets(5))
      c1  = (zets(2)-zets(4))*(zets(2)-zets(5))
      c2  = 2.d0*((zets(2)-zets(4))+(zets(2)-zets(5)))
      c   = (zets(2)-zets(1))*(zets(2)-zets(3))
c     Get coefficient for i=1
      ku  = 1
      ind = ku+1
      CA2(ind+1,1) = (2.d0*c1*(b1m
     &     *(2.d0*zets(2)-zets(3)-zets(1))
     &     +b2m*(zets(2)-zets(1))*(3.d0*zets(2)-2.d0*zets(3)
     &     -zets(1)))+c2*c*(b1m+b2m*(zets(2)-zets(1))))/b/b1
      ku  = 5
      ind = ku+1
      CB2(ind+1,1) = -(2.d0*c1*(1.d0+a1m*(2.d0*zets(2)-zets(3)
     &     -zets(1))+a2m*(zets(2)-zets(1))*(3.d0*zets(2)
     &     -2.d0*zets(3)-zets(1)))+c2*c*
     &     (1.d0/(zets(2)-zets(1))+a1m+a2m*(zets(2)-zets(1))))
     &     /b/b1
CC$OMP SECTION 
C     Get a1p,a2p,b1p,b2p
      b   = zets(1)-zets(3)
      t1  = (zets(3)-zets(1))*(zets(3)-zets(2))*(zets(3)-zets(4))
     &     *(zets(3)-zets(5))
      t1  = 1.d0/t1
      a1  =  (zets(3)-zets(1))*(zets(3)-zets(2))*(zets(3)-zets(4))
     &     +(zets(3)-zets(1))*(zets(3)-zets(2))*(zets(3)-zets(5))
     &     +(zets(3)-zets(1))*(zets(3)-zets(4))*(zets(3)-zets(5))
     &     +(zets(3)-zets(2))*(zets(3)-zets(4))*(zets(3)-zets(5))
      a1  = a1*t1
      a2  = (zets(3)-zets(1))*(zets(3)-zets(2))
     &     +(zets(3)-zets(1))*(zets(3)-zets(4))
     &     +(zets(3)-zets(1))*(zets(3)-zets(5))
     &     +(zets(3)-zets(2))*(zets(3)-zets(4))
     &     +(zets(3)-zets(2))*(zets(3)-zets(5))
     &     +(zets(3)-zets(4))*(zets(3)-zets(5))
      a2  = 2.d0*a2*t1
      a3  = (zets(1)-zets(2))*(zets(1)-zets(4))
     &     *(zets(1)-zets(5))*t1
      a4  = (zets(1)-zets(2))*(zets(1)-zets(4))
     &     +(zets(1)-zets(2))*(zets(1)-zets(5))
     &     +(zets(1)-zets(4))*(zets(1)-zets(5))
      a4  = 2.d0*a4*t1
C     Solve for a1p,a2p
      a11 = 2.d0*a1
      a12 = 2.d0
      a21 = b*a4+2.d0*a3
      a22 = b**2*a4+4.d0*b*a3
      b1  = -a2
      b2  = -a4
      a1p = (b1*a22-b2*a12)/(a11*a22-a12*a21)
      a2p = (a11*b2-a21*b1)/(a11*a22-a12*a21)
c     Solve for b1p,b2p
      a11 = 2.d0*a1
      a12 = 2.d0
      b1  = 1.d0
      b2  = 0.d0
      b1p = (b1*a22-b2*a12)/(a11*a22-a12*a21)
      b2p = (a11*b2-a21*b1)/(a11*a22-a12*a21)
      b   = zets(3)-zets(1)
      b1  = (zets(3)-zets(2))*(zets(3)-zets(4))*(zets(3)-zets(5))
      c1  = (zets(2)-zets(4))*(zets(2)-zets(5))
      c2  = 2.d0*(2.d0*zets(2)-zets(4)-zets(5))
      c   = (zets(2)-zets(1))*(zets(2)-zets(3))
C     Get coefficient for i=3
      ku  = 1
      ind = ku+1
      CA2(ind,2)   = 1.d0
      CA2(ind-1,3) =-(2.d0*c1*(b1p
     &     *(2.d0*zets(2)-zets(1)-zets(3))
     &     +b2p*(zets(2)-zets(3))*(3.d0*zets(2)-2.d0*zets(1)
     &     -zets(3)))+c2*c*(b1p+b2p*(zets(2)-zets(3))))/b/b1
      ku  = 5
      ind = ku+1
      CB2(ind-1,3) = (2.d0*c1*(1.d0+a1p*(2.d0*zets(2)-zets(3)
     &     -zets(1))+a2p*(zets(2)-zets(3))*(3.d0*zets(2)
     &     -2.d0*zets(1)-zets(3)))+c2*c*
     &     (1.d0/(zets(2)-zets(3))+a1p+a2p*(zets(2)-zets(3))))
     &     /b/b1
CC$OMP END PARALLEL SECTIONS

      ku  = 5
      ind = ku+1
CC$OMP PARALLEL SECTIONS DEFAULT(SHARED),PRIVATE(t1,t2,t3,t4,t5)
CC$OMP+ PRIVATE(t6,t7,a11,a12,a21,a22,c1m,c2m)
      t1 = zets(3)-zets(2)
      t2 = zets(1)-zets(2)
      t3 = (zets(1)-zets(3))*(zets(1)-zets(4))
     &     *(zets(1)-zets(5))
      t4 = (zets(3)-zets(1))*(zets(3)-zets(4))
     &     *(zets(3)-zets(5))
      t5 = 2.d0*((zets(1)-zets(3))*(zets(1)-zets(4))
     &     +(zets(1)-zets(3))*(zets(1)-zets(5))
     &     +(zets(1)-zets(4))*(zets(1)-zets(5)))
      t6 = 2.d0*((zets(3)-zets(1))*(zets(3)-zets(4))
     &     +(zets(3)-zets(1))*(zets(3)-zets(5))
     &     +(zets(3)-zets(4))*(zets(3)-zets(5)))
      a11 = (zets(1)-zets(2))*t5+2*t3
      a12 = (zets(1)-zets(2))**2*t5+4.d0*(zets(1)-zets(2))*t3
      a21 = (zets(3)-zets(2))*t6+2.d0*t4
      a22 = (zets(3)-zets(2))**2*t6+4.d0*(zets(3)-zets(2))*t4
      t7  = a11*a22-a12*a21
      c1m = (a12*t6-a22*t5)/t7
      c2m = (a21*t5-a11*t6)/t7
      t3  = 1.d0/(zets(2)-zets(4))+1.d0/(zets(2)-zets(5))
      t4  = 2.d0/((zets(2)-zets(4))*(zets(2)-zets(5)))
      CB2(ind,2) = 2.d0*c2m+2.d0*c1m*(t3-(1.d0/t1+1.d0/t2))+t4
     &     +(2.d0-2.d0*t3*(t1+t2))/t1/t2
CC$OMP SECTION      
      t3  = (zets(1)-zets(2))*(zets(1)-zets(3))
     &     *(zets(1)-zets(5))
      t4  = (zets(3)-zets(1))*(zets(3)-zets(2))
     &     *(zets(3)-zets(5))
      t5  = 2.d0*((zets(1)-zets(2))*(zets(1)-zets(3))
     &     +(zets(1)-zets(2))*(zets(1)-zets(5))
     &     +(zets(1)-zets(3))*(zets(1)-zets(5)))
      t6  = 2.d0*((zets(3)-zets(2))*(zets(3)-zets(1))
     &     +(zets(3)-zets(2))*(zets(3)-zets(5))
     &     +(zets(3)-zets(1))*(zets(3)-zets(5)))
      a11 = (zets(1)-zets(4))*t5+2.d0*t3
      a12 = (zets(1)-zets(4))**2*t5+4.d0*(zets(1)-zets(4))*t3
      a21 = (zets(3)-zets(4))*t6+2.d0*t4
      a22 = (zets(3)-zets(4))**2*t6+4.d0*(zets(3)-zets(4))*t4
      t7  = a11*a22-a12*a21
      c1m = (a12*t6-a22*t5)/t7
      c2m = (a21*t5-a11*t6)/t7
      t3  = zets(3)-zets(2)
      t4  = zets(1)-zets(2)
      t1  = zets(4)-zets(2)
      t5  = (zets(2)-zets(5))/((zets(4)-zets(5))
     &     *(zets(4)-zets(2)))
      t6  = t3*t4/((zets(4)-zets(1))*(zets(4)-zets(3)))
      t7  = 2.d0/((zets(4)-zets(5))*(zets(4)-zets(2)))
      CB2(ind-2,4) =(c1m*(1.d0+t1/t3+t1/t4)
     &     -c2m*(2.d0+t1/t3+t1/t4)*t1
     &     -1.d0/t3-1.d0/t4)*2.d0*t5*t6+t6*t7*(1.d0-c1m*t1+c2m*t1**2)
CC$OMP SECTION      
      t3  = (zets(1)-zets(2))*(zets(1)-zets(3))
     &     *(zets(1)-zets(4))
      t4  = (zets(3)-zets(1))*(zets(3)-zets(2))
     &     *(zets(3)-zets(4))
      t5  = 2.d0*((zets(1)-zets(2))*(zets(1)-zets(3))
     &     +(zets(1)-zets(2))*(zets(1)-zets(4))
     &     +(zets(1)-zets(3))*(zets(1)-zets(4)))
      t6  = 2.d0*((zets(3)-zets(1))*(zets(3)-zets(2))
     &     +(zets(3)-zets(1))*(zets(3)-zets(4))
     &     +(zets(3)-zets(2))*(zets(3)-zets(4)))
      a11 = (zets(1)-zets(5))*t5+2.d0*t3
      a12 = (zets(1)-zets(5))**2*t5+4.d0*(zets(1)-zets(5))*t3
      a21 = (zets(3)-zets(5))*t6+2.d0*t4
      a22 = (zets(3)-zets(5))**2*t6+4.d0*(zets(3)-zets(5))*t4
      t7  = a11*a22-a12*a21
      c1m = (a12*t6-a22*t5)/t7 
      c2m = (a21*t5-a11*t6)/t7
      t3  = zets(3)-zets(2)
      t4  = zets(1)-zets(2)
      t1  = zets(5)-zets(2)
      t5  = (zets(2)-zets(4))/((zets(5)-zets(4))
     &     *(zets(5)-zets(2)))
      t6  = t3*t4/((zets(5)-zets(3))*(zets(5)-zets(1)))
      t7  = 2.d0/((zets(5)-zets(2))*(zets(5)-zets(4)))
      CB2(ind-3,5) =(c1m*(1.d0+t1/t3+t1/t4)
     &     -c2m*(2.d0+t1/t3+t1/t4)*t1
     &     -1.d0/t3-1.d0/t4)*2.d0*t5*t6+t6*t7*(1.d0-c1m*t1+c2m*t1**2)
CC$OMP END PARALLEL SECTIONS

c     For i = N-1
C     Get a1m,a2m,b1m,b2m
CC$OMP PARALLEL SECTIONS DEFAULT(SHARED),PRIVATE(a1m,a2m,a1p,a2p)
CC$OMP+ PRIVATE(b,t1,a1,a2,a3,a4,b1,b2,a11,a12,a21,a22,ind,ku)
      b   = zets(nz)-zets(nz-2)
      t1  = (zets(nz-2)-zets(nz))*(zets(nz-2)-zets(nzm))
     &     *(zets(nz-2)-zets(nz-3))*(zets(nz-2)-zets(nz-4))
      t1  = 1.d0/t1
      a1  =  (zets(nz-2)-zets(nz))*(zets(nz-2)-zets(nzm))
     &     *(zets(nz-2)-zets(nz-3))
     &     +(zets(nz-2)-zets(nz))*(zets(nz-2)-zets(nzm))
     &     *(zets(nz-2)-zets(nz-4))
     &     +(zets(nz-2)-zets(nz))*(zets(nz-2)-zets(nz-3))
     &     *(zets(nz-2)-zets(nz-4))
     &     +(zets(nz-2)-zets(nzm))*(zets(nz-2)-zets(nz-3))
     &     *(zets(nz-2)-zets(nz-4))
      a1  = a1*t1
      a2  = (zets(nz-2)-zets(nz))*(zets(nz-2)-zets(nzm))
     &     +(zets(nz-2)-zets(nz))*(zets(nz-2)-zets(nz-3))
     &     +(zets(nz-2)-zets(nz))*(zets(nz-2)-zets(nz-4))
     &     +(zets(nz-2)-zets(nzm))*(zets(nz-2)-zets(nz-3))
     &     +(zets(nz-2)-zets(nzm))*(zets(nz-2)-zets(nz-4))
     &     +(zets(nz-2)-zets(nz-3))*(zets(nz-2)-zets(nz-4))
      a2  = 2.d0*a2*t1
      a3  = (zets(nz)-zets(nzm))*(zets(nz)-zets(nz-3))
     &     *(zets(nz)-zets(nz-4))*t1
      a4  = (zets(nz)-zets(nzm))*(zets(nz)-zets(nz-3))
     &     +(zets(nz)-zets(nzm))*(zets(nz)-zets(nz-4))
     &     +(zets(nz)-zets(nz-3))*(zets(nz)-zets(nz-4))
      a4  = 2.d0*a4*t1
C     Solve for a1m,a2m
      a11 = 2.d0*a1
      a12 = 2.d0
      a21 = b*a4+2.d0*a3
      a22 = b**2*a4+4.d0*b*a3
      b1  = -a2
      b2  = -a4
      a1m = (b1*a22-b2*a12)/(a11*a22-a12*a21)
      a2m = (a11*b2-a21*b1)/(a11*a22-a12*a21)
c     Solve for b1m,b2m
      a11 = 2.d0*a1
      a12 = 2.d0
      b1  = 1.d0
      b2  = 0.d0
      b1m = (b1*a22-b2*a12)/(a11*a22-a12*a21)
      b2m = (a11*b2-a21*b1)/(a11*a22-a12*a21) 
      b   = zets(nz-2)-zets(nz)
      b1  = (zets(nz-2)-zets(nzm))*(zets(nz-2)-zets(nz-3))
     &     *(zets(nz-2)-zets(nz-4))
      c1  = (zets(nzm)-zets(nz-3))*(zets(nzm)-zets(nz-4))
      c2  = 2.d0*(2.d0*zets(nzm)-zets(nz-3)-zets(nz-4))
      c   = (zets(nzm)-zets(nz-2))*(zets(nzm)-zets(nz))
c     Get coefficient for i=nz-2
      ku  = 1
      ind = ku+1
      CA2(ind+1,nz-2) = -(2.d0*c1*(b1m
     &     *(2.d0*zets(nzm)-zets(nz-2)-zets(nz))
     &     +b2m*(zets(nzm)-zets(nz-2))*(3.d0*zets(nzm)-2.d0*zets(nz)
     &     -zets(nz-2)))+c2*c*(b1m+b2m*(zets(nzm)-zets(nz-2))))/b/b1
      ku  = 5
      ind = ku+1
      CB2(ind+1,nz-2) = (2.d0*c1*(1.d0+a1m*(2.d0*zets(nzm)-zets(nz)
     &     -zets(nz-2))+a2m*(zets(nzm)-zets(nz-2))*(3.d0*zets(nzm)
     &     -2.d0*zets(nz)-zets(nz-2)))+c2*c*
     &     (1.d0/(zets(nzm)-zets(nz-2))+a1m+a2m*(zets(nzm)-zets(nz-2))))
     &     /b/b1
CC$OMP SECTION 
C     Get a1p,a2p,b1p,b2p
      b   = zets(nz-2)-zets(nz)
      t1  = (zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-3))
     &     *(zets(nz)-zets(nz-4))*(zets(nz)-zets(nzm))
      t1  = 1.d0/t1
      a1  =  (zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-3))
     &     *(zets(nz)-zets(nz-4))
     &     +(zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-3))
     &     *(zets(nz)-zets(nzm))
     &     +(zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-4))
     &     *(zets(nz)-zets(nzm))
     &     +(zets(nz)-zets(nz-3))*(zets(nz)-zets(nz-4))
     &     *(zets(nz)-zets(nzm))
      a1  = a1*t1
      a2  = (zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-3))
     &     +(zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-4))
     &     +(zets(nz)-zets(nz-2))*(zets(nz)-zets(nzm))
     &     +(zets(nz)-zets(nz-3))*(zets(nz)-zets(nz-4))
     &     +(zets(nz)-zets(nz-3))*(zets(nz)-zets(nzm))
     &     +(zets(nz)-zets(nz-4))*(zets(nz)-zets(nzm))
      a2  = 2.d0*a2*t1
      a3  = (zets(nz-2)-zets(nzm))*(zets(nz-2)-zets(nz-3))
     &     *(zets(nz-2)-zets(nz-4))*t1
      a4  = (zets(nz-2)-zets(nzm))*(zets(nz-2)-zets(nz-3))
     &     +(zets(nz-2)-zets(nzm))*(zets(nz-2)-zets(nz-4))
     &     +(zets(nz-2)-zets(nz-3))*(zets(nz-2)-zets(nz-4))
      a4  = 2.d0*a4*t1
C     Solve for a1p,a2p
      a11 = 2.d0*a1
      a12 = 2.d0
      a21 = b*a4+2.d0*a3
      a22 = b**2*a4+4.d0*b*a3
      b1  = -a2
      b2  = -a4
      a1p = (b1*a22-b2*a12)/(a11*a22-a12*a21)
      a2p = (a11*b2-a21*b1)/(a11*a22-a12*a21)
c     Solve for b1p,b2p
      a11 = 2.d0*a1
      a12 = 2.d0
      b1  = 1.d0
      b2  = 0.d0
      b1p = (b1*a22-b2*a12)/(a11*a22-a12*a21)
      b2p = (a11*b2-a21*b1)/(a11*a22-a12*a21)
      b   = zets(nz)-zets(nz-2)
      b1  = (zets(nz)-zets(nzm))*(zets(nz)-zets(nz-3))
     &     *(zets(nz)-zets(nz-4))
      c1  = (zets(nzm)-zets(nz-3))*(zets(nzm)-zets(nz-4))
      c2  = 2.d0*(2.d0*zets(nzm)-zets(nz-3)-zets(nz-4))
      c   = (zets(nzm)-zets(nz))*(zets(nzm)-zets(nz-2))
C     Get coefficient for i=nz
      ku  = 1
      ind = ku+1
      CA2(ind,nzm)  = 1.d0
      CA2(ind-1,nz) = -(2.d0*c1*(b1p
     &     *(2.d0*zets(nzm)-zets(nz)-zets(nz-2))
     &     +b2p*(zets(nzm)-zets(nz))*(3.d0*zets(nzm)-2.d0*zets(nz-2)
     &     -zets(nz)))+c2*c*(b1p+b2p*(zets(nzm)-zets(nz))))/b/b1
      ku  = 5
      ind = ku+1
      CB2(ind-1,nz) = (2.d0*c1*(1.d0+a1p*(2.d0*zets(nzm)-zets(nz)
     &     -zets(nz-2))+a2p*(zets(nzm)-zets(nz))*(3.d0*zets(nzm)
     &     -2.d0*zets(nz-2)-zets(nz)))+c2*c*
     &     (1.d0/(zets(nzm)-zets(nz))+a1p+a2p*(zets(nzm)-zets(nz))))
     &     /b/b1
CC$OMP END PARALLEL SECTIONS

      ku  = 5
      ind = ku+1
CC$OMP PARALLEL SECTIONS DEFAULT(SHARED),PRIVATE(t1,t2,t3,t4,t5)
CC$OMP+ PRIVATE(t6,t7,a11,a12,a21,a22,c1m,c2m)
      t1 = zets(nz)-zets(nzm)
      t2 = zets(nz-2)-zets(nzm)
      t3 = (zets(nz-2)-zets(nz))*(zets(nz-2)-zets(nz-3))
     &     *(zets(nz-2)-zets(nz-4))
      t4 = (zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-3))
     &     *(zets(nz)-zets(nz-4))
      t5 = 2.d0*((zets(nz-2)-zets(nz))*(zets(nz-2)-zets(nz-3))
     &     +(zets(nz-2)-zets(nz))*(zets(nz-2)-zets(nz-4))
     &     +(zets(nz-2)-zets(nz-3))*(zets(nz-2)-zets(nz-4)))
      t6 = 2.d0*((zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-3))
     &     +(zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-4))
     &     +(zets(nz)-zets(nz-3))*(zets(nz)-zets(nz-4)))
      a11 = (zets(nz-2)-zets(nzm))*t5+2*t3
      a12 = (zets(nz-2)-zets(nzm))**2*t5+4.d0*(zets(nz-2)-zets(nzm))*t3
      a21 = (zets(nz)-zets(nzm))*t6+2.d0*t4
      a22 = (zets(nz)-zets(nzm))**2*t6+4.d0*(zets(nz)-zets(nzm))*t4
      t7  = a11*a22-a12*a21
      c1m = (a12*t6-a22*t5)/t7
      c2m = (a21*t5-a11*t6)/t7
      t3  = 1.d0/(zets(nzm)-zets(nz-3))+1.d0/(zets(nzm)-zets(nz-4))
      t4  = 2.d0/((zets(nzm)-zets(nz-3))*(zets(nzm)-zets(nz-4)))
      CB2(ind,nzm) = 2.d0*c2m+2.d0*c1m*(t3-(1.d0/t1+1.d0/t2))+t4
     &     +(2.d0-2.d0*t3*(t1+t2))/t1/t2
CC$OMP SECTION      
      t3  = (zets(nz-2)-zets(nzm))*(zets(nz-2)-zets(nz))
     &     *(zets(nz-2)-zets(nz-4))
      t4  = (zets(nz)-zets(nzm))*(zets(nz)-zets(nz-2))
     &     *(zets(nz)-zets(nz-4))
      t5  = 2.d0*((zets(nz-2)-zets(nzm))*(zets(nz-2)-zets(nz))
     &     +(zets(nz-2)-zets(nzm))*(zets(nz-2)-zets(nz-4))
     &     +(zets(nz-2)-zets(nz))*(zets(nz-2)-zets(nz-4)))
      t6  = 2.d0*((zets(nz)-zets(nzm))*(zets(nz)-zets(nz-2))
     &     +(zets(nz)-zets(nzm))*(zets(nz)-zets(nz-4))
     &     +(zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-4)))
      a11 = (zets(nz-2)-zets(nz-3))*t5+2.d0*t3
      a12 = (zets(nz-2)-zets(nz-3))**2*t5
     &     +4.d0*(zets(nz-2)-zets(nz-3))*t3
      a21 = (zets(nz)-zets(nz-3))*t6+2.d0*t4
      a22 = (zets(nz)-zets(nz-3))**2*t6+4.d0*(zets(nz)-zets(nz-3))*t4
      t7  = a11*a22-a12*a21
      c1m = (a12*t6-a22*t5)/t7
      c2m = (a21*t5-a11*t6)/t7
      t3  = zets(nz)-zets(nzm)
      t4  = zets(nz-2)-zets(nzm)
      t5  = (zets(nzm)-zets(nz-4))/((zets(nz-3)-zets(nzm))
     &     *(zets(nz-3)-zets(nz-4)))
      t6  = t3*t4/((zets(nz-3)-zets(nz))*(zets(nz-3)-zets(nz-2)))
      t1  = zets(nz-3)-zets(nzm)
      t7  = 2.d0/((zets(nz-3)-zets(nzm))*(zets(nz-3)-zets(nz-4)))
      CB2(ind+2,nz-3) =(c1m*(1.d0+t1/t3+t1/t4)
     &     -c2m*(2.d0+t1/t3+t1/t4)*t1
     &     -1.d0/t3-1.d0/t4)*2.d0*t5*t6+t6*t7*(1.d0-c1m*t1+c2m*t1**2)
CC$OMP SECTION      
      t3  = (zets(nz-2)-zets(nz))*(zets(nz-2)-zets(nzm))
     &     *(zets(nz-2)-zets(nz-3))
      t4  = (zets(nz)-zets(nz-2))*(zets(nz)-zets(nzm))
     &     *(zets(nz)-zets(nz-3))
      t5  = 2.d0*((zets(nz-2)-zets(nz))*(zets(nz-2)-zets(nzm))
     &     +(zets(nz-2)-zets(nz))*(zets(nz-2)-zets(nz-3))
     &     +(zets(nz-2)-zets(nzm))*(zets(nz-2)-zets(nz-3)))
      t6  = 2.d0*((zets(nz)-zets(nz-2))*(zets(nz)-zets(nzm))
     &     +(zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-3))
     &     +(zets(nz)-zets(nzm))*(zets(nz)-zets(nz-3)))
      a11 = (zets(nz-2)-zets(nz-4))*t5+2.d0*t3
      a12 = (zets(nz-2)-zets(nz-4))**2*t5
     &     +4.d0*(zets(nz-2)-zets(nz-4))*t3
      a21 = (zets(nz)-zets(nz-4))*t6+2.d0*t4
      a22 = (zets(nz)-zets(nz-4))**2*t6+4.d0*(zets(nz)-zets(nz-4))*t4
      t7  = a11*a22-a12*a21
      c1m = (a12*t6-a22*t5)/t7 
      c2m = (a21*t5-a11*t6)/t7
      t3  = zets(nz)-zets(nzm)
      t4  = zets(nz-2)-zets(nzm)
      t1  = zets(nz-4)-zets(nzm)
      t5  = (zets(nzm)-zets(nz-3))/((zets(nz-4)-zets(nzm))
     &     *(zets(nz-4)-zets(nz-3)))
      t6  = t3*t4/((zets(nz-4)-zets(nz))*(zets(nz-4)-zets(nz-2)))
      t7  = 2.d0/((zets(nz-4)-zets(nzm))*(zets(nz-4)-zets(nz-3)))
      CB2(ind+3,nz-4) =(c1m*(1.d0+t1/t3+t1/t4)
     &     -c2m*(2.d0+t1/t3+t1/t4)*t1
     &     -1.d0/t3-1.d0/t4)*2.d0*t5*t6+t6*t7*(1.d0-c1m*t1+c2m*t1**2)
CC$OMP END PARALLEL SECTIONS

C     Get the diagonal form of matrix CA2, where d2 is the main diagonal
c     d2u is the upper diagonal and d2l is the lower diagonal
      call dcopy(nz,CA2(2,1),3,d2,1)
      call dcopy(nzm,CA2(1,2),3,d2u,1)
      call dcopy(nzm,CA2(3,1),3,d2l,1)
 
C     For Dirichlet-Dirichlet Boundary condition
C     For i = 2
      b  = zets(3)-zets(2)
CC$OMP PARALLEL SECTIONS DEFAULT(SHARED),PRIVATE(ku,ind,b2,a1,a2,a3)
      ku  = 1
      ind = ku+1
      CA2_drdr(ind,2) = 1.d0
      a1 = (zets(2)-zets(1))*(zets(2)-zets(4))*(zets(2)-zets(5))
     &     *(zets(2)-zets(6))
      a2 = (zets(3)-zets(1))*(zets(3)-zets(2))*(zets(3)-zets(4))
     &     *(zets(3)-zets(5))*(zets(3)-zets(6))*(1.d0/(zets(3)-zets(1))
     &     +1.d0/(zets(3)-zets(2))+1.d0/(zets(3)-zets(4))
     &     +1.d0/(zets(3)-zets(5))+1.d0/(zets(3)-zets(6)))
      a3 = (zets(2)-zets(1))*(zets(2)-zets(4))*(zets(2)-zets(5))
     &     +(zets(2)-zets(1))*(zets(2)-zets(4))*(zets(2)-zets(6))
     &     +(zets(2)-zets(1))*(zets(2)-zets(5))*(zets(2)-zets(6))
     &     +(zets(2)-zets(4))*(zets(2)-zets(5))*(zets(2)-zets(6))
      a3 = 2.d0*a3
      CA2_drdr(ind-1,3) = 0.5d0*b*a3/a2-a1/a2
      a4 = (zets(3)-zets(1))*(zets(3)-zets(2))*(zets(3)-zets(4))
     &     +(zets(3)-zets(1))*(zets(3)-zets(2))*(zets(3)-zets(5))
     &     +(zets(3)-zets(1))*(zets(3)-zets(2))*(zets(3)-zets(6))
     &     +(zets(3)-zets(1))*(zets(3)-zets(4))*(zets(3)-zets(5))
     &     +(zets(3)-zets(1))*(zets(3)-zets(4))*(zets(3)-zets(6))
     &     +(zets(3)-zets(1))*(zets(3)-zets(5))*(zets(3)-zets(6))
     &     +(zets(3)-zets(2))*(zets(3)-zets(4))*(zets(3)-zets(5))
     &     +(zets(3)-zets(2))*(zets(3)-zets(4))*(zets(3)-zets(6))
     &     +(zets(3)-zets(2))*(zets(3)-zets(5))*(zets(3)-zets(6))
     &     +(zets(3)-zets(4))*(zets(3)-zets(5))*(zets(3)-zets(6))
      a4 = 2.d0*a4
      a5 = (zets(3)-zets(1))*(zets(3)-zets(2))*(zets(3)-zets(4))
     &     *(zets(3)-zets(5))*(zets(3)-zets(6))
      ku  = 5
      ind = ku+1
      CB2_drdr(ind-1,3) = a3/a5+0.5d0*b*a3*a4/(a2*a5)-a1*a4/(a2*a5)
CC$OMP SECTION
      ku  = 5
      ind = ku+1
      a1 = 1.d0/(zets(2)-zets(1))+1.d0/(zets(2)-zets(4))
     &     +1.d0/(zets(2)-zets(5))+1.d0/(zets(2)-zets(6))
      a2 = (zets(3)-zets(1))*(zets(3)-zets(4))*(zets(3)-zets(5))
     &     *(zets(3)-zets(6))/((zets(2)-zets(1))*(zets(2)-zets(4))
     &     *(zets(2)-zets(5))*(zets(2)-zets(6)))
      a3 = a2*(1.d0/(zets(3)-zets(1))+1.d0/(zets(3)-zets(4))
     &     +1.d0/(zets(3)-zets(5))+1.d0/(zets(3)-zets(6)))
      a4 = (zets(2)-zets(1))*(zets(2)-zets(4))+(zets(2)-zets(1))
     &     *(zets(2)-zets(5))+(zets(2)-zets(1))*(zets(2)-zets(6))
     &     +(zets(2)-zets(4))*(zets(2)-zets(5))+(zets(2)-zets(4))
     &     *(zets(2)-zets(6))+(zets(2)-zets(5))*(zets(2)-zets(6))
      a4 = 2.d0*a4/((zets(2)-zets(1))*(zets(2)-zets(4))
     &     *(zets(2)-zets(5))*(zets(2)-zets(6)))
      CB2_drdr(ind,2) = a4-2.d0*a1*(a2+2.d0*b*a3)/(b*(a2+b*a3))
     &     +2.d0*a3/(b*a2+b**2*a3)
CC$OMP SECTION
      ku  = 5
      ind = ku+1
      b2 = zets(3)-zets(1)
      a1 = (zets(1)-zets(2))*(zets(1)-zets(4))
     &     *(zets(1)-zets(5))*(zets(1)-zets(6))
      a2 = (zets(3)-zets(2))*(zets(3)-zets(4))*(zets(3)-zets(5))
     &     *(zets(3)-zets(6))/a1
      a3 = a2*(1.d0/(zets(3)-zets(2))+1.d0/(zets(3)-zets(4))
     &     +1.d0/(zets(3)-zets(5))+1.d0/(zets(3)-zets(6)))
      a4 = (zets(2)-zets(4))*(zets(2)-zets(5))*(zets(2)-zets(6))/a1
      a5 = 2.d0*((zets(2)-zets(4))*(zets(2)-zets(5))+(zets(2)-zets(4))
     &     *(zets(2)-zets(6))+(zets(2)-zets(5))*(zets(2)-zets(6)))/a1
      CB2_drdr(ind+1,1) = a5*(-b*a2-b**2*a3)/(-b2*a2-b2**2*a3)
     &     -2.d0*a4*(a2+2.d0*b*a3)/(a2+b2*a3)/b2
CC$OMP SECTION
      ku  = 5
      ind = ku+1
      b2 = zets(3)-zets(4)
      a1 = (zets(4)-zets(1))*(zets(4)-zets(2))*(zets(4)-zets(5))
     &     *(zets(4)-zets(6))
      a2 = (zets(3)-zets(1))*(zets(3)-zets(2))*(zets(3)-zets(5))
     &     *(zets(3)-zets(6))/a1
      a3 = a2*(1.d0/(zets(3)-zets(1))+1.d0/(zets(3)-zets(2))
     &     +1.d0/(zets(3)-zets(5))+1.d0/(zets(3)-zets(6)))
      a4 = (zets(2)-zets(1))*(zets(2)-zets(5))*(zets(2)-zets(6))/a1
      a5 = 2.d0*((zets(2)-zets(1))*(zets(2)-zets(5))+(zets(2)-zets(1))
     &     *(zets(2)-zets(6))+(zets(2)-zets(5))*(zets(2)-zets(6)))/a1
      CB2_drdr(ind-2,4) = a5*(-b*a2-b**2*a3)/(-b2*a2-b2**2*a3)
     &     -2.d0*a4*(a2+2.d0*b*a3)/(a2+b2*a3)/b2
CC$OMP SECTION
      ku  = 5
      ind = ku+1
      b2 = zets(3)-zets(5)
      a1 = (zets(5)-zets(1))*(zets(5)-zets(2))*(zets(5)-zets(4))
     &     *(zets(5)-zets(6))
      a2 = (zets(3)-zets(1))*(zets(3)-zets(2))*(zets(3)-zets(4))
     &     *(zets(3)-zets(6))/a1
      a3 = a2*(1.d0/(zets(3)-zets(1))+1.d0/(zets(3)-zets(2))
     &     +1.d0/(zets(3)-zets(4))+1.d0/(zets(3)-zets(6)))
      a4 = (zets(2)-zets(1))*(zets(2)-zets(4))*(zets(2)-zets(6))/a1
      a5 = 2.d0*((zets(2)-zets(1))*(zets(2)-zets(4))+(zets(2)-zets(1))
     &     *(zets(2)-zets(6))+(zets(2)-zets(4))*(zets(2)-zets(6)))/a1
      CB2_drdr(ind-3,5) = a5*(-b*a2-b**2*a3)/(-b2*a2-b2**2*a3)
     &     -2.d0*a4*(a2+2.d0*b*a3)/(a2+b2*a3)/b2
CC$OMP SECTION
      ku  = 5
      ind = ku+1
      b2 = zets(3)-zets(6)
      a1 = (zets(6)-zets(1))*(zets(6)-zets(2))*(zets(6)-zets(4))
     &     *(zets(6)-zets(5))
      a2 = (zets(3)-zets(1))*(zets(3)-zets(2))*(zets(3)-zets(4))
     &     *(zets(3)-zets(5))/a1
      a3 = a2*(1.d0/(zets(3)-zets(1))+1.d0/(zets(3)-zets(2))
     &     +1.d0/(zets(3)-zets(4))+1.d0/(zets(3)-zets(5)))
      a4 = (zets(2)-zets(1))*(zets(2)-zets(4))*(zets(2)-zets(5))/a1
      a5 = 2.d0*((zets(2)-zets(1))*(zets(2)-zets(4))+(zets(2)-zets(1))
     &     *(zets(2)-zets(5))+(zets(2)-zets(4))*(zets(2)-zets(5)))/a1
      CB2_drdr(ind-4,6) = a5*(-b*a2-b**2*a3)/(-b2*a2-b2**2*a3)
     &     -2.d0*a4*(a2+2.d0*b*a3)/(a2+b2*a3)/b2
CC$OMP END PARALLEL SECTIONS

C     For i = N-1
CC$OMP PARALLEL SECTIONS DEFAULT(SHARED), 
CC$OMP+ PRIVATE(ku,ind,a1,a2,a3,a4,a5,b2)
      ku  = 1
      ind = ku+1
      CA2_drdr(ind,nzm) = 1.d0
      b  = zets(nz-2)-zets(nzm)
      a1 = (zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz-4))
     &     *(zets(nzm)-zets(nz-3))*(zets(nzm)-zets(nz))
      a2 = (zets(nz-2)-zets(nz-5))*(zets(nz-2)-zets(nz-4))
     &     *(zets(nz-2)-zets(nz-3))*(zets(nz-2)-zets(nzm))*
     &     (zets(nz-2)-zets(nz))*(1.d0/(zets(nz-2)-zets(nz-5))
     &     +1.d0/(zets(nz-2)-zets(nz-4))+1.d0/(zets(nz-2)-zets(nz-3))
     &     +1.d0/(zets(nz-2)-zets(nzm))+1.d0/(zets(nz-2)-zets(nz)))
      a3 = (zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz-4))
     &     *(zets(nzm)-zets(nz-3))
     &     +(zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz-4))
     &     *(zets(nzm)-zets(nz))
     &     +(zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz-3))
     &     *(zets(nzm)-zets(nz))
     &     +(zets(nzm)-zets(nz-4))*(zets(nzm)-zets(nz-3))
     &     *(zets(nzm)-zets(nz))
      a3 = 2.d0*a3
      CA2_drdr(ind+1,nz-2) = 0.5d0*b*a3/a2-a1/a2

      ku  = 5
      ind = ku+1
      a4  =(zets(nz-2)-zets(nz-5))*(zets(nz-2)-zets(nz-4))
     &     *(zets(nz-2)-zets(nz-3))
     &     +(zets(nz-2)-zets(nz-5))*(zets(nz-2)-zets(nz-4))
     &     *(zets(nz-2)-zets(nzm))
     &     +(zets(nz-2)-zets(nz-5))*(zets(nz-2)-zets(nz-4))
     &     *(zets(nz-2)-zets(nz))
     &     +(zets(nz-2)-zets(nz-5))*(zets(nz-2)-zets(nz-3))
     &     *(zets(nz-2)-zets(nzm))
     &     +(zets(nz-2)-zets(nz-5))*(zets(nz-2)-zets(nz-3))
     &     *(zets(nz-2)-zets(nz))
     &     +(zets(nz-2)-zets(nz-5))*(zets(nz-2)-zets(nzm))
     &     *(zets(nz-2)-zets(nz))
     &     +(zets(nz-2)-zets(nz-4))*(zets(nz-2)-zets(nz-3))
     &     *(zets(nz-2)-zets(nzm))
     &     +(zets(nz-2)-zets(nz-4))*(zets(nz-2)-zets(nz-3))
     &     *(zets(nz-2)-zets(nz))
     &     +(zets(nz-2)-zets(nz-4))*(zets(nz-2)-zets(nzm))
     &     *(zets(nz-2)-zets(nz))
     &     +(zets(nz-2)-zets(nz-3))*(zets(nz-2)-zets(nzm))
     &     *(zets(nz-2)-zets(nz))
      a4 = 2.d0*a4
      a5 = (zets(nz-2)-zets(nz-5))*(zets(nz-2)-zets(nz-4))
     &     *(zets(nz-2)-zets(nz-3))*(zets(nz-2)-zets(nzm))
     &     *(zets(nz-2)-zets(nz))
      CB2_drdr(ind+1,nz-2) = a3/a5+0.5d0*b*a3*a4/(a2*a5)-a1*a4/(a2*a5)
CC$OMP SECTION
      ku  = 5
      ind = ku+1
      a1 = 1.d0/(zets(nzm)-zets(nz-5))+1.d0/(zets(nzm)-zets(nz-4))
     &     +1.d0/(zets(nzm)-zets(nz-3))+1.d0/(zets(nzm)-zets(nz))
      a2 = (zets(nz-2)-zets(nz-5))*(zets(nz-2)-zets(nz-4))
     &     *(zets(nz-2)-zets(nz-3))*(zets(nz-2)-zets(nz))
     &     /((zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz-4))
     &     *(zets(nzm)-zets(nz-3))*(zets(nzm)-zets(nz)))
      a3 = a2*(1.d0/(zets(nz-2)-zets(nz-5))+1.d0/(zets(nz-2)-zets(nz-4))
     &     +1.d0/(zets(nz-2)-zets(nz-3))+1.d0/(zets(nz-2)-zets(nz)))
      a4 = (zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz-4))
     &     +(zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz-3))
     &     +(zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz))
     &     +(zets(nzm)-zets(nz-4))*(zets(nzm)-zets(nz-3))
     &     +(zets(nzm)-zets(nz-4))*(zets(nzm)-zets(nz))
     &     +(zets(nzm)-zets(nz-3))*(zets(nzm)-zets(nz))
      a4 = 2.d0*a4/((zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz-4))
     &     *(zets(nzm)-zets(nz-3))*(zets(nzm)-zets(nz)))
      CB2_drdr(ind,nzm) = a4-2.d0*a1*(a2+2.d0*b*a3)/(b*(a2+b*a3))
     &     +2.d0*a3/(b*a2+b**2*a3)
CC$OMP SECTION
      ku  = 5
      ind = ku+1
      b2 = zets(nz-2)-zets(nz)
      a1 = (zets(nz)-zets(nz-5))*(zets(nz)-zets(nz-4))
     &     *(zets(nz)-zets(nz-3))*(zets(nz)-zets(nzm))
      a2 = (zets(nz-2)-zets(nz-5))*(zets(nz-2)-zets(nz-4))
     &     *(zets(nz-2)-zets(nz-3))*(zets(nz-2)-zets(nzm))/a1
      a3 = a2*(1.d0/(zets(nz-2)-zets(nz-5))+1.d0/(zets(nz-2)-zets(nz-4))
     &     +1.d0/(zets(nz-2)-zets(nz-3))+1.d0/(zets(nz-2)-zets(nzm)))
      a4 = (zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz-4))
     &     *(zets(nzm)-zets(nz-3))/a1
      a5 = 2.d0*((zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz-4))
     &     +(zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz-3))
     &     +(zets(nzm)-zets(nz-4))*(zets(nzm)-zets(nz-3)))/a1
      CB2_drdr(ind-1,nz)=a5*(-b*a2-b**2*a3)/(-b2*a2-b2**2*a3)
     &     -2.d0*a4*(a2+2.d0*b*a3)/(a2+b2*a3)/b2
CC$OMP SECTION      
      ku  = 5
      ind = ku+1
      b2 = zets(nz-2)-zets(nz-3)
      a1 = (zets(nz-3)-zets(nz-5))*(zets(nz-3)-zets(nz-4))
     &     *(zets(nz-3)-zets(nzm))*(zets(nz-3)-zets(nz))
      a2 = (zets(nz-2)-zets(nz-5))*(zets(nz-2)-zets(nz-4))
     &     *(zets(nz-2)-zets(nzm))*(zets(nz-2)-zets(nz))/a1
      a3 = a2*(1.d0/(zets(nz-2)-zets(nz-5))+1.d0/(zets(nz-2)-zets(nz-4))
     &     +1.d0/(zets(nz-2)-zets(nzm))+1.d0/(zets(nz-2)-zets(nz)))
      a4 = (zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz-4))
     &     *(zets(nzm)-zets(nz))/a1
      a5 = 2.d0*((zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz-4))
     &     +(zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz))
     &     +(zets(nzm)-zets(nz-4))*(zets(nzm)-zets(nz)))/a1
      CB2_drdr(ind+2,nz-3) = a5*(-b*a2-b**2*a3)/(-b2*a2-b2**2*a3)
     &     -2.d0*a4*(a2+2.d0*b*a3)/(a2+b2*a3)/b2
CC$OMP SECTION      
      ku  = 5
      ind = ku+1
      b2 = zets(nz-2)-zets(nz-4)
      a1 = (zets(nz-4)-zets(nz-5))*(zets(nz-4)-zets(nz-3))
     &     *(zets(nz-4)-zets(nzm))*(zets(nz-4)-zets(nz))
      a2 = (zets(nz-2)-zets(nz-5))*(zets(nz-2)-zets(nz-3))
     &     *(zets(nz-2)-zets(nzm))*(zets(nz-2)-zets(nz))/a1
      a3 = a2*(1.d0/(zets(nz-2)-zets(nz-5))+1.d0/(zets(nz-2)-zets(nz-3))
     &     +1.d0/(zets(nz-2)-zets(nzm))+1.d0/(zets(nz-2)-zets(nz)))
      a4 = (zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz-3))
     &     *(zets(nzm)-zets(nz))/a1
      a5 = 2.d0*((zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz-3))
     &     +(zets(nzm)-zets(nz-5))*(zets(nzm)-zets(nz))
     &     +(zets(nzm)-zets(nz-3))*(zets(nzm)-zets(nz)))/a1
      CB2_drdr(ind+3,nz-4) = a5*(-b*a2-b**2*a3)/(-b2*a2-b2**2*a3)
     &     -2.d0*a4*(a2+2.d0*b*a3)/(a2+b2*a3)/b2
CC$OMP SECTION  
      ku  = 5
      ind = ku+1
      b2 = zets(nz-2)-zets(nz-5)
      a1 = (zets(nz-5)-zets(nz-4))*(zets(nz-5)-zets(nz-3))
     &     *(zets(nz-5)-zets(nzm))*(zets(nz-5)-zets(nz))
      a2 = (zets(nz-2)-zets(nz-4))*(zets(nz-2)-zets(nz-3))
     &     *(zets(nz-2)-zets(nzm))*(zets(nz-2)-zets(nz))/a1
      a3 = a2*(1.d0/(zets(nz-2)-zets(nz-4))+1.d0/(zets(nz-2)-zets(nz-3))
     &     +1.d0/(zets(nz-2)-zets(nzm))+1.d0/(zets(nz-2)-zets(nz)))
      a4 = (zets(nzm)-zets(nz-4))*(zets(nzm)-zets(nz-3))
     &     *(zets(nzm)-zets(nz))/a1
      a5 = 2.d0*((zets(nzm)-zets(nz-4))*(zets(nzm)-zets(nz-3))
     &     +(zets(nzm)-zets(nz-4))*(zets(nzm)-zets(nz))
     &     +(zets(nzm)-zets(nz-3))*(zets(nzm)-zets(nz)))/a1
      CB2_drdr(ind+4,nz-5) = a5*(-b*a2-b**2*a3)/(-b2*a2-b2**2*a3)
     &     -2.d0*a4*(a2+2.d0*b*a3)/(a2+b2*a3)/b2
CC$OMP END PARALLEL SECTIONS
      
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,k)
      do k=1,nz
         do i=1,3
            CA2_nono(i,k) = CA2_drdr(i,k)
            CA2_nodr(i,k) = CA2_drdr(i,k)
            CA2_mimi(i,k) = CA2_drdr(i,k)
         enddo
      enddo
C$OMP END PARALLEL DO
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,k)
      do k=1,nz
         do i=1,11
            CB2_nono(i,k) = CB2_drdr(i,k)
            CB2_nodr(i,k) = CB2_drdr(i,k)
            CB2_mimi(i,k) = CB2_drdr(i,k)
         enddo
      enddo
C$OMP END PARALLEL DO
      
C     For Dirichlet-Dirichlet boundary condition
      CA2_drdr(2,1) = 1.d0
      CB2_drdr(6,1) = 1.d0
      CA2_drdr(2,nz)= 1.d0
      CB2_drdr(6,nz)= 1.d0

C     For Neuman-Neuman boundary condition
      CA2_nono(2,1) = 1.d0
      CA2_nono(2,nz)= 1.d0
      ku  = 5
      ind = ku+1
C     For i=1
      CB2_nono(ind,1)   =1.d0/(zets(1)-zets(2))+1.d0/(zets(1)-zets(3))
     &     +1.d0/(zets(1)-zets(4))+1.d0/(zets(1)-zets(5))
     &     +1.d0/(zets(1)-zets(6))
      CB2_nono(ind-1,2) = (zets(1)-zets(3))*(zets(1)-zets(4))
     &     *(zets(1)-zets(5))*(zets(1)-zets(6))/((zets(2)-zets(1))
     &     *(zets(2)-zets(3))*(zets(2)-zets(4))*(zets(2)-zets(5))
     &     *(zets(2)-zets(6)))
      CB2_nono(ind-2,3) = (zets(1)-zets(2))*(zets(1)-zets(4))
     &     *(zets(1)-zets(5))*(zets(1)-zets(6))/((zets(3)-zets(1))
     &     *(zets(3)-zets(2))*(zets(3)-zets(4))*(zets(3)-zets(5))
     &     *(zets(3)-zets(6)))
      CB2_nono(ind-3,4) = (zets(1)-zets(2))*(zets(1)-zets(3))
     &     *(zets(1)-zets(5))*(zets(1)-zets(6))/((zets(4)-zets(1))
     &     *(zets(4)-zets(2))*(zets(4)-zets(3))*(zets(4)-zets(5))
     &     *(zets(4)-zets(6)))
      CB2_nono(ind-4,5) = (zets(1)-zets(2))*(zets(1)-zets(3))
     &     *(zets(1)-zets(4))*(zets(1)-zets(6))/((zets(5)-zets(1))
     &     *(zets(5)-zets(2))*(zets(5)-zets(3))*(zets(5)-zets(4))
     &     *(zets(5)-zets(6)))
      CB2_nono(ind-5,6) = (zets(1)-zets(2))*(zets(1)-zets(3))
     &     *(zets(1)-zets(4))*(zets(1)-zets(5))/((zets(6)-zets(1))
     &     *(zets(6)-zets(2))*(zets(6)-zets(3))*(zets(6)-zets(4))
     &     *(zets(6)-zets(5)))
C     For i = N
      CB2_nono(ind,nz)     = 1.d0/(zets(nz)-zets(nzm))
     &     +1.d0/(zets(nz)-zets(nz-2))+1.d0/(zets(nz)-zets(nz-3))
     &     +1.d0/(zets(nz)-zets(nz-4))+1.d0/(zets(nz)-zets(nz-5))
      CB2_nono(ind+1,nzm)  = (zets(nz)-zets(nz-2))*(zets(nz)-zets(nz-3))
     &     *(zets(nz)-zets(nz-4))*(zets(nz)-zets(nz-5))
     &     /((zets(nzm)-zets(nz))*(zets(nzm)-zets(nz-2))
     &     *(zets(nzm)-zets(nz-3))*(zets(nzm)-zets(nz-4))
     &     *(zets(nzm)-zets(nz-5)))
      CB2_nono(ind+2,nz-2) = (zets(nz)-zets(nzm))*(zets(nz)-zets(nz-3))
     &     *(zets(nz)-zets(nz-4))*(zets(nz)-zets(nz-5))
     &     /((zets(nz-2)-zets(nz))*(zets(nz-2)-zets(nzm))
     &     *(zets(nz-2)-zets(nz-3))*(zets(nz-2)-zets(nz-4))
     &     *(zets(nz-2)-zets(nz-5)))
      CB2_nono(ind+3,nz-3) = (zets(nz)-zets(nzm))*(zets(nz)-zets(nz-2))
     &     *(zets(nz)-zets(nz-4))*(zets(nz)-zets(nz-5))
     &     /((zets(nz-3)-zets(nz))*(zets(nz-3)-zets(nzm))
     &     *(zets(nz-3)-zets(nz-2))*(zets(nz-3)-zets(nz-4))
     &     *(zets(nz-3)-zets(nz-5)))
      CB2_nono(ind+4,nz-4) = (zets(nz)-zets(nzm))*(zets(nz)-zets(nz-2))
     &     *(zets(nz)-zets(nz-3))*(zets(nz)-zets(nz-5))
     &     /((zets(nz-4)-zets(nz))*(zets(nz-4)-zets(nzm))
     &     *(zets(nz-4)-zets(nz-2))*(zets(nz-4)-zets(nz-3))
     &     *(zets(nz-4)-zets(nz-5)))
      CB2_nono(ind+5,nz-5) = (zets(nz)-zets(nzm))*(zets(nz)-zets(nz-2))
     &     *(zets(nz)-zets(nz-3))*(zets(nz)-zets(nz-4))
     &     /((zets(nz-5)-zets(nz))*(zets(nz-5)-zets(nzm))
     &     *(zets(nz-5)-zets(nz-2))*(zets(nz-5)-zets(nz-3))
     &     *(zets(nz-5)-zets(nz-4)))

c      print *,CB2_nono(ind+4,nz-4),CB2_nono(ind-4,5)
C     For mixed-mixed boundary condition
C     mixed-mixed bc for concentration
c$$$      ttat(1) = re*sc*abs(w_settling*gz(1,1)/sqrt(gy(1,1)**2
c$$$     &     +gz(1,1)**2))
      ku  = 5
      ind = ku+1
      ttat(1) = re*sc*abs(w_settling*gz(1,1)/sqrt(gy(1,1)**2
     &          +gz(1,1)**2))
      ttbt(1) = 1.d0
      ttgt(1) = 0.d0
      ttab(1) = ttat(1)
      ttbb(1) = ttbt(1)
      ttgb(1) = ttgt(1)
      CA2_mimi(2,nz) = 1.d0
      CA2_mimi(2,1)  = 1.d0
c     For i=1
      CB2_mimi(ind,1)   = CB2_nono(ind,1)*ttbt(1)+ttat(1)
      CB2_mimi(ind-1,2) = CB2_nono(ind-1,2)*ttbt(1)
      CB2_mimi(ind-2,3) = CB2_nono(ind-2,3)*ttbt(1)
      CB2_mimi(ind-3,4) = CB2_nono(ind-3,4)*ttbt(1)
      CB2_mimi(ind-4,5) = CB2_nono(ind-4,5)*ttbt(1)
      CB2_mimi(ind-5,6) = CB2_nono(ind-5,6)*ttbt(1)
c     For i=N
      CB2_mimi(ind,nz)     = CB2_nono(ind,nz)*ttbb(1)+ttab(1)
      CB2_mimi(ind+1,nzm)  = CB2_nono(ind+1,nzm)*ttbb(1)
      CB2_mimi(ind+2,nz-2) = CB2_nono(ind+2,nz-2)*ttbb(1)
      CB2_mimi(ind+3,nz-3) = CB2_nono(ind+3,nz-3)*ttbb(1)
      CB2_mimi(ind+4,nz-4) = CB2_nono(ind+4,nz-4)*ttbb(1)
      CB2_mimi(ind+5,nz-5) = CB2_nono(ind+5,nz-5)*ttbb(1)

C     For Neuman-Dirichlet boundary condition
      CA2_nodr(2,nz) = 1.d0
      CB2_nodr(6,nz) = 1.d0
      CA2_nodr(2,1)  = 1.d0
      ku = 5
      ind = ku+1  
C     For i=1
      CB2_nodr(ind,1)   = CB2_nono(ind,1)
      CB2_nodr(ind-1,2) = CB2_nono(ind-1,2)
      CB2_nodr(ind-2,3) = CB2_nono(ind-2,3)
      CB2_nodr(ind-3,4) = CB2_nono(ind-3,4)
      CB2_nodr(ind-4,5) = CB2_nono(ind-4,5)
      CB2_nodr(ind-5,6) = CB2_nono(ind-5,6)
      return
      
      end

c$$$      subroutine ddz(val,dval)
c$$$#include "common.inc"
c$$$      real*8 val(nx0,ny0,nz),dval(nx0,ny0,nz)
c$$$      real*8 tmp(nz),tmp2(nz)
c$$$C     This subroutine computes ddz in both real and Fourire spaces
c$$$      integer i,j,k
c$$$C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
c$$$      do j=1,ny
c$$$         do i=1,nx
c$$$c     Copy data into tmp
c$$$            call dcopy(nz,val(i,j,1),nx0y0,tmp,1)
c$$$c     B*tmp
c$$$            call dgbmv('n',nz,nz,5,5,1.d0,CB1,11,tmp,1,0.d0,
c$$$     &           tmp2,1)
c$$$c     solve A*tmp2=tmp2
c$$$            call dgtsv(nz,1,d1l,d1,d1u,tmp2,nz,k)
c$$$c     copy back to dval
c$$$            call dcopy(nz,tmp2,1,dval(i,j,1),nx0y0)
c$$$         enddo 
c$$$      enddo
c$$$C$OMP END PARALLEL DO    
c$$$      return
c$$$      end



      subroutine GaussEliminate(n,A,kl,ku,lda,b,x,info)
      implicit none
      integer n,ku,kl,lda,info
      real*8 A(lda,*),b(*),x(*)
c     Guassian Elimination method to solve nonsymmetric system Ax=b
c     without pivoting
C     A is a general banded matrix
c     parameter list
c     n    : A is an n*n matrix
c     A    : A general banded matrix
c     ku   : number of superdiagonals
c     kl   : number of subdiagonals
c     x    : return unknowns
c     b    : Right hand side
c     info : On output, indicate ternimation status
c            -1 : zero on diagonal
c             0 : Ternimate normally

      integer i,j,k
      integer ind,ind1,ind2
      integer upper,lower,ju
      real*8  alpha,beta
      
c     Initialization
      ind  = ku+1
      info = 0
c     Start loop
      do i=1,n
         lower = max(i,i-kl) 
         upper = min(i+ku,n)
         ju    = min(i+kl,n)
         do j= i+1,ju
            ind1  = ind
            ind2  = ind + (j-i)
            if(A(ind1,i).eq.0.d0) then
               info = -1
               print *, 'Zero diagonal, failed!'
               return
            endif
            alpha = A(ind2,i)/A(ind1,i)
            b(j)  = b(j) - alpha*b(i)
            do k= lower,upper
               ind1 = ind + (i-k)
               ind2 = ind + (j-k)
               A(ind2,k) = A(ind2,k) - alpha*A(ind1,k)
            enddo     
         enddo
      enddo
c     Debug use only
c      do i=1,n
c         do j=1,11
c            write(79,*) A(j,i)
c         enddo
c      enddo
C     Use backsubstitution to solve x
      x(n) = b(n)/A(ind,n)
      do i=n-1,1,-1
         k = min(n-i,ku)
         do j=i+1,i+k
            ind1 = ind + i-j
            b(i) = b(i) - A(ind1,j)*x(j)
         enddo
         x(i) = b(i)/A(ku+1,i)
      enddo
         
      return
      end
