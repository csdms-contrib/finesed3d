      SUBROUTINE initial(dt0,ch_fin,m_fin)
#include "common.inc"
      REAL*8 dt0
      INTEGER ch_fin,m_fin
c
c     set up the matrices and ffts, initialize the field variables
c
      INTEGER i,j,k,pt
      REAL*8 pin,hpin,sum,c100,c200,xk,yk,d00,d0n,dn0,dnn,iden
     &      ,inzm,tinzm,hinzm
      REAL*8 condo,condi,emax,emin
      REAL*8 aux_rr,h0,sumtemp,theta
      REAL*8 alpha,epsilon,tmp
c     for compact scheme
      INTEGER ind,ind1,ind2,ku,kl

      if (FLTIMER) call begin_timer(10)

c  set write and read format
#ifdef COPPER
C      open(21,file=name,form='unformatted')
      format_mode='unformatted'
#endif
#ifdef COBALT
C      open(21,file=name,form='binary')
      format_mode='binary'
#endif

c  set some variables
      pi = 4.d0*atan(1.d0)
      ipin = dble(nzm)/pi
c
c  horizontal box size for channel flow 
      rlenx = 1.d0    !pi      !2.d0*pi
      rleny = 2.d0    !2.d0*pi !4.d0*pi
c
      deltax = rlenx/dble(nx)
      deltay = rleny/dble(ny)
      idx = 1.d0/deltax
      idy = 1.d0/deltay
      dt0 = 0.d0
c
c for density currents specify ra=-1
      if(ri.eq.-1)ri=1.0
c
c set erosion coefficients
cz not used in this version
      cferos1=1.d0
      cferos2=0.6d0
      if(rey_p.le.2.34)then
         cferos1=0.586d0
         cferos2=1.23d0
      endif
c
c     zero multidimensional arrays
c---------------------------------------------------------------------
c
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP DO PRIVATE(i,j,k)
      do k=1,nz
        do j=1,ny0
          do i=1,nx0
            u(  i,j,k) = 0.d0
            v(  i,j,k) = 0.d0
            w(  i,j,k) = 0.d0
            p(  i,j,k) = 0.d0
            ru( i,j,k) = 0.d0
            rv( i,j,k) = 0.d0
            rw( i,j,k) = 0.d0
            rp( i,j,k) = 0.d0
            h1( i,j,k) = 0.d0
            h2( i,j,k) = 0.d0
            h3( i,j,k) = 0.d0
            rtt(i,j,k) = 0.d0
            extmx(i,j,k) = 0.d0
	    extmy(i,j,k) = 0.d0
            extmz(i,j,k) = 0.d0
          enddo
        enddo
      enddo
C$OMP END DO NOWAIT
      if (FLTHRM) then
C$OMP DO PRIVATE(i,j,k)
         do k=1,nz
            do j=1,ny0
               do i=1,nx0
                  tt( i,j,k) = 0.d0
                  htt(i,j,k) = 0.d0
               enddo
            enddo
         enddo
C$OMP END DO NOWAIT
      endif
      if (FLNBSQ.or.FLINERTIA) then
C$OMP DO PRIVATE(i,j,k)
         do k=1,nz
            do j=1,ny0
               do i=1,nx0
                  ou1(i,j,k) = 0.d0
                  ov1(i,j,k) = 0.d0
                  ow1(i,j,k) = 0.d0
                  ou2(i,j,k) = 0.d0
                  ov2(i,j,k) = 0.d0
                  ow2(i,j,k) = 0.d0
                  ou3(i,j,k) = 0.d0
                  ov3(i,j,k) = 0.d0
                  ow3(i,j,k) = 0.d0
                  DuDt(i,j,k) = 0.d0
                  DvDt(i,j,k) = 0.d0
                  DwDt(i,j,k) = 0.d0
                  up(i,j,k) = 0.d0
                  vp(i,j,k) = 0.d0
                  wp(i,j,k) = 0.d0
               enddo
            enddo
         enddo
C$OMP END DO NOWAIT
      endif
C$OMP SECTIONS PRIVATE(i,j)
      do j=1,nyh
        do i=1,nx0
          pbxnb(i,j) = 0.d0
          pcxnb(i,j) = 0.d0
          pbynb(i,j) = 0.d0
          pcynb(i,j) = 0.d0
          pbxob(i,j) = 0.d0
          pcxob(i,j) = 0.d0
          pbyob(i,j) = 0.d0
          pcyob(i,j) = 0.d0
          pbxnt(i,j) = 0.d0
          pcxnt(i,j) = 0.d0
          pbynt(i,j) = 0.d0
          pcynt(i,j) = 0.d0
          pbxot(i,j) = 0.d0
          pcxot(i,j) = 0.d0
          pbyot(i,j) = 0.d0
          pcyot(i,j) = 0.d0
        enddo
      enddo
c
c      if(FLAVER)then
         do k=1,nz
            um(   k) = 0.d0
            vm(   k) = 0.d0
            wm(   k) = 0.d0
            ttm(  k) = 0.d0
            u2m(  k) = 0.d0
            v2m(  k) = 0.d0
            w2m(  k) = 0.d0
            tt2m( k) = 0.d0
            u3m(  k) = 0.d0
            v3m(  k) = 0.d0
            w3m(  k) = 0.d0
            tt3m( k) = 0.d0
            u4m(  k) = 0.d0
            v4m(  k) = 0.d0
            w4m(  k) = 0.d0
            tt4m( k) = 0.d0
            uvm(  k) = 0.d0
            uwm(  k) = 0.d0
            uttm( k) = 0.d0
            vwm(  k) = 0.d0
            vttm( k) = 0.d0
            wttm( k) = 0.d0
c
            uxm(  k) = 0.d0
            uym(  k) = 0.d0
            uzm(  k) = 0.d0
            vxm(  k) = 0.d0
            vym(  k) = 0.d0
            vzm(  k) = 0.d0
            ttxm( k) = 0.d0
            ttym( k) = 0.d0
            ttzm( k) = 0.d0
c
            ux2m( k) = 0.d0
            uy2m( k) = 0.d0
            uz2m( k) = 0.d0
            vx2m( k) = 0.d0
            vy2m( k) = 0.d0
            vz2m( k) = 0.d0
            ttx2m(k) = 0.d0
            tty2m(k) = 0.d0
            ttz2m(k) = 0.d0
c
            ux3m( k) = 0.d0
            uy3m( k) = 0.d0
            uz3m( k) = 0.d0
            vx3m( k) = 0.d0
            vy3m( k) = 0.d0
            vz3m( k) = 0.d0
            ttx3m(k) = 0.d0
            tty3m(k) = 0.d0
            ttz3m(k) = 0.d0
c
            ux4m( k) = 0.d0
            uy4m( k) = 0.d0
            uz4m( k) = 0.d0
            vx4m( k) = 0.d0
            vy4m( k) = 0.d0
            vz4m( k) = 0.d0
            ttx4m(k) = 0.d0
            tty4m(k) = 0.d0
            ttz4m(k) = 0.d0
         enddo
c      endif
c     
      do j=1,ny0
         do i=1,nx0
            ubct(   i,j)=0.d0   !=ugt if non-homog, see uvci_mimi
            ubcb(   i,j)=0.d0   !=ugt if non-homog, see uvci_mimi
            vbct(   i,j)=0.d0   !=vgt if non-homog, see uvci_mimi
            vbcb(   i,j)=0.d0   !=vgt if non-homog, see uvci_mimi
            wbct(   i,j)=0.d0   !=wgt if non-homog, see  wci_mimi
            wbcb(   i,j)=0.d0   !=wgt if non-homog, see  wci_mimi
            ttbct(  i,j)=ttgt(1)
            ttbcb(  i,j)=ttgb(1)
            ttot(i,j)   = 0.d0
            ttob(i,j)   = 0.d0
            ttotn(i,j)  = 0.d0
            ttobn(i,j)  = 0.d0
            ttoto(i,j)  = 0.d0
            ttobo(i,j)  = 0.d0
            dwbct(i,j)  = 0.d0
            dwbcb(i,j)  = 0.d0
            wsh(i,j) = 0.d0
            ttwb(i,j) = 0.d0
            auxdep(i,j) = 0.d0
            deposit(i,j)=0.d0
            erosion(i,j)=0.d0
         enddo
      enddo
C$OMP SECTION
      do j=1,nz
        do i=1,nz0
          chb(    i,j) = 0.d0
          chbd1(  i,j) = 0.d0
          chbd2(  i,j) = 0.d0
          cn_drdr(i,j) = 0.d0
          ci_drdr(i,j) = 0.d0
          cn_nodr(i,j) = 0.d0
          ci_nodr(i,j) = 0.d0
        enddo
      enddo
      do j=1,nzmm
        do i=1,nzmm0
            cn_nono(i,j) = 0.d0
            ci_nono(i,j) = 0.d0
          uvcn_mimi(i,j) = 0.d0
          uvci_mimi(i,j) = 0.d0
           wcn_mimi(i,j) = 0.d0
           wci_mimi(i,j) = 0.d0
          ttcn_mimi(i,j) = 0.d0
          ttci_mimi(i,j) = 0.d0
        enddo
      enddo
C$OMP END SECTIONS
C$OMP END PARALLEL
c
c     define the coefficients for the r-k time stepping
c
c     coefficient of old advection term
c
      cfnl1v(1) = 0.d0
      cfnl1v(2) = -5.d0/9.d0
      cfnl1v(3) = -153.d0/128.d0
c
c     coefficient of advection term
c
      cfnl2v(1) = 1.d0/3.d0
      cfnl2v(2) = 15.d0/16.d0
      cfnl2v(3) = 8.d0/15.d0
c
c     coefficient of diffusion term
c
      cfdifv(1) = dt/6.d0
      cfdifv(2) = dt*5.d0/24.d0
      cfdifv(3) = dt/8.d0
c
c     sub-step sizes for RK3 stepping
c
      cfdelv(1) = 1.d0/3.d0
      cfdelv(2) = 5.d0/12.d0
      cfdelv(3) = 1.d0/4.d0
c
c     accumulation of sub-step sizes
c
      cfcumv(1) = 1.d0/3.d0
      cfcumv(2) = 3.d0/4.d0
      cfcumv(3) = 1.d0
c
c     compute coeficients for particle velocity
      if(FLINERTIA)then
       cfup1 = tau_p*(1.d0-beta_p)
c       cfup2 = tau_p*beta_p/2.d0
       cfup2 = 0.d0
       cfup3 = tau_p*w_settling
       if(FLOTAU)cfup3 = 0.d0      
      endif
c
c     set up fourier apparatus
c     initialize fft's
c      call initfft
c
      c100 = 2.d0*pi/rlenx
      c200 = 2.d0*pi/rleny
c
c     wavex is the wavenumber, for a complex DFT
c     wavexs is the square of the wavenumber
c     the subarray 1:nxhp is used for a real DFT
c
c     wave numbers above the cutoff are set to zero,
c     and are hence able to enforce cutoff (though poor style)
c
      do i=1,nx
        wavex(i) = 0.d0
        wavexs(i) = 0.d0
      enddo
      do i=2,nxh
        xk = c100*dble(i-1)
        wavex(i) = xk
        wavexs(i) = xk**2
        wavex(nx+2-i) = -xk
        wavexs(nx+2-i) = xk**2
      enddo
c
      do j=1,ny
        wavey(j) = 0.d0
        waveys(j) = 0.d0
      enddo
      do j=2,nyh
        yk = c200*dble(j-1)
        wavey(j) = yk
        waveys(j) = yk**2
        wavey(ny+2-j) = -yk
        waveys(ny+2-j) = yk**2
      enddo
c
      write(26,107)(i,wavex(i),i=1,nx)
      write(26,108)(i,wavey(i),i=1,ny)
 107  FORMAT('   wavex  ',/,129(2x,i4,2x,1pe15.7,/))
 108  FORMAT('   wavey  ',/,129(2x,i4,2x,1pe15.7,/))
c
c     compute the stretched grid in the vertical direction
c
c     wzet are the weights used to compute the average over [-1,1]
c
      pin  = pi/dble(nzm)
      hpin = 0.5d0*pin
      epsilon = 0.1d0
      alpha   = asinh(1.d0/epsilon)
      do k = 1,nz/2
c        zets(k) = epsilon*sinh(cos(dble(k-1)*pin)*alpha)
        zets(k) = cos(dble(k-1)*pin)
        zets(nz+1-k) = -zets(k)
        sum = 1.d0
        if (mod(nz,2).eq.1) sum = 1.d0 + (-1)**k / dble(nzm**2-1)
        do i = 1,nz/2-1
          sum = sum - 2.d0*cos(2.d0*pin*dble(i*(k-1)))/dble(4*i**2-1)
        enddo
        sum = sum/dble(nzm)
        if (k.eq.1) sum = 0.5d0*sum
        wzet(k) = sum
        wzet(nz+1-k) = sum
      enddo
      if (mod(nz,2).eq.1) then
        k = nz/2 + 1
        zets(k) = 0.d0
        sum = 1.d0 + (-1)**k / dble(nzm**2-1)
        do i = 1,nz/2-1
          sum = sum - 2.d0*(-1)**i / dble(4*i**2-1)
        enddo
        wzet(k) = sum/dble(nzm)
      endif
      do k = 1,nzm/2
        deltaz(k) = zets(k)-zets(k+1)
        deltaz(nz-k) = deltaz(k)
      enddo
c      if (mod(nz,2).eq.0) deltaz(nz/2) = 2.d0*sin(hpin)*
c     &  sin(dble(nzm)*hpin)
      do k=1,nzm
        ndeltaz(k) = -deltaz(k)
        idz(k) = 1.d0/deltaz(k)
        nidz(k) = -idz(k)
      enddo
c
      write(26,119)
      do k = 1,nzm
        write(26,121) k,zets(k),wzet(k),deltaz(k)
      enddo
      write(26,121) nz,zets(nz),wzet(k)
 119  FORMAT(/,'   i   gauss lobatto pts     weights         spacings')
 121  FORMAT(2x,i4,3(2x,1pe15.7,2x,1pe15.7,2x,1pe15.7))
c
c     compute anisotropic spectral viscosity
c     x-direction
      do i=1,nx
         svx(i)=0.d0
      enddo
      if(epssvx.gt.0)then
         write(26,*)'Using spectral viscosity in x-direction'
         if(msvx.lt.nxh)then
            do i=msvx+1,nxh
               svx(i)=epssvx*exp((dble(i)**2-dble(nxh )**2)
     *                                   /(dble(i)**2-dble(msvx)**2))
               svx(nx+2-i)=svx(i)
            enddo
         else
          write(26,*)'Error in Spectral viscosity:  msvx >= nxh'
          write(*,*)'Error in Spectral viscosity:  msvx >= nxh'
          stop
         endif
      endif
c     y-direction
      do j=1,nyh
         svy(j)=0.d0
      enddo
      if(epssvy.gt.0)then
         write(26,*)'Using spectral viscosity in y-direction'
         if(msvy.lt.nyh)then
            do j=msvy+1,nyh
               svy(j)=epssvy*exp((dble(j)**2-dble(nyh )**2)
     *                                   /(dble(j)**2-dble(msvy)**2))
            enddo
         else
          write(26,*)'Error in Spectral viscosity:  msvy >= nyh'
          write(*,*)'Error in Spectral viscosity:  msvy >= nyh'
          stop
         endif
      endif
c
c compute gravity vector
      FLGVAR=.false.
      theta=0.d0 !degrees
      theta=0.0 !2.d0*pi/3.6d2*theta !rad
      if(FLGVAR)then
         do k=1,nz
            do j=1,nyl0
               gy(j,k)=dsin(theta)
c     *                *derf(deltay/1.5d0*dble(j-nylhp)*sqrt(re*pr)*3.d0)
     *                *derf(dble(j-nylhp)*sqrt(re*sc)*1.d-1)
               gz(j,k)=-dcos(theta)
            enddo
         enddo
      else
c this change is done by EMRE everything is set to -1 & 0 in z and y direction      
         gy(1,1)=0.0    !1.d0                     !dsin(theta) !dsin(theta)
         gz(1,1)=-1.0       !-dcos(theta)/dsin(theta) !-dcos(theta)!-dcos(theta)
      endif
c
c---------------------------------------------------------------------
c       compute chb matrix for Chebyshev spectrum
c---------------------------------------------------------------------
c
c      pin=pi/dble(nzm) !pin is computed before for zets
      inzm=1.d0/dble(nzm)
      tinzm=2.d0*inzm
      hinzm=inzm/2.d0
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP DO PRIVATE(i,j)
      do i=2,nzm
       do j=2,nzm
        chb(i,j)=tinzm*cos(pin*dble(j-1)*dble(i-1))
       enddo 
       chb(i,1) =inzm
       chb(i,nz)=inzm*cos(pi*dble(i-1))
      enddo
C$OMP END DO
C$OMP DO PRIVATE(j)
      do j=2,nzm
       chb(1,j) =inzm
       chb(nz,j)=inzm*cos(pi*dble(j-1))
      enddo
C$OMP END DO
C$OMP END PARALLEL
      chb(1,1)  =hinzm
      chb(1,nz) =hinzm
      chb(nz,1) =hinzm
      chb(nz,nz)=hinzm*cos(pi*dble(nzm))
c
c----------------------------------------------------------------------
c     compute coefficient matrix for Compact finite difference schemes
c----------------------------------------------------------------------
      call compact
c
c---------------------------------------------------------------------
c       compute sm_d1 matrix
c---------------------------------------------------------------------
c
      ku = 5
      kl = 5
      do j=1,nz
         ind1 = j-ku
         if(ind1.le.1) ind1 = 1
         ind2 = j+kl
         if(ind2.ge.nz) ind2 = nz
         do i=ind1,ind2
            ind = ku+1+i-j
            chbd1(i,j) = CB1(ind,j)
         enddo
      enddo
      write(90,*) chbd1
C     Using dgtsv
      call dgtsv(nz,nz,d1l,d1,d1u,chbd1,nz0,i)
      if(i.ne.0) then
         print *, 'Fatal error in computing chbd1'
         stop
      endif
c
c---------------------------------------------------------------------
c     compute sm_d2 matrix by forming the product sm_d1*sm_d1
c---------------------------------------------------------------------
c       
      call dgemm('N','N',nz,nz,nz,1.d0,chbd1,nz0,chbd1,nz0,0.d0,
     &  chbd2,nz0)
c
c---------------------------------------------------------------------
c     compute the sm_d2 matrix with neumann-dirichlet bc's
c---------------------------------------------------------------------
c
      do j = 1,nz
        ci_nodr(1,j)  = chbd1(1,j)
        do i = 2,nzm
          ci_nodr(i,j) = chbd2(i,j)
        enddo
      enddo
      ci_nodr(nz,nz) = 1.d0
c
c---------------------------------------------------------------------
c     diagonalize the ci_nodr matrix
c---------------------------------------------------------------------
c
      call e_decomp(nz,ci_nodr,nz0,cn_nodr,e_nodr,condo,condi,emax,emin)
      write(26,*)'cn_nodr: ', condo, condi, emax, emin
c
c
c---------------------------------------------------------------------
c     compute the sm_d2 matrix with neumann-neumann bc's 
c---------------------------------------------------------------------
c
      d00   = chbd1(1,1)
      d0n   = chbd1(1,nz)
      dn0   = chbd1(nz,1)
      dnn   = chbd1(nz,nz)
      iden  = 1.d0/(dn0*d0n-d00*dnn)
      do i=1,nzmm
        prbc1(i) = (dnn*chbd1(1, i+1)-d0n*chbd1(nz,i+1))*iden
        prbcn(i) = (d00*chbd1(nz,i+1)-dn0*chbd1(1, i+1))*iden
      enddo
c
      do j = 2,nzm
        do i = 2,nzm
          ci_nono(i-1,j-1) = chbd2(i,j) + chbd2(i, 1)*prbc1(j-1)
     &                                  + chbd2(i,nz)*prbcn(j-1)
        enddo
      enddo
c
c---------------------------------------------------------------------
c     diagonalize the ci_nono matrix
c---------------------------------------------------------------------
c
      call e_decomp(nzmm,ci_nono,nzmm0,cn_nono,e_nono,condo,condi,
     &  emax,emin)
      write(26,*)'cn_nono: ', condo, condi, emax, emin

C     Set value for homogeneous dirichlet-dirichlet bc's
      uvat = 1.d0
      uvab = 1.d0
      uvbt = 1.d0-uvat
      uvbb = 1.d0-uvab
      ugt = 0.d0
      ugb = 0.d0
      vgt = 0.d0
      vgb = 0.d0
c---------------------------------------------------------------------
c     read in data
c---------------------------------------------------------------------
c
      rtime = 0.d0
      if (ch_fin.eq.0) then
c     initialize field for zero start run for ch_fin = 0
        call init_uvw
c        call courant(dt0,u,v,w)
      else if (ch_fin.eq.1) then
c     dt_old is read from input file
c     read field data as written by routine "output" for ch_fin = 1
        call input
      else
        write(26,*) 'invalid value for ch_fin'
        stop
      endif
c
      if(m_fin.eq.0)then
         ust2=0.d0
         usb2=0.d0
         vst2=0.d0
         vsb2=0.d0
         NT=0
      elseif(m_fin.eq.1)then
c     read mean value and NT
c         call input_means
      else
         write(26,*) 'invalid value for m_fin'
         stop
      endif
c
      if (FLTIMER) call end_timer(10)
      return
      END
c
c     Eigendecomposition
c
c     Input:  the matrix a
c     Output: eval contains the eigenvalues of a
c             evec contains the eigenvector matrix of a
c             a    now contains the inverse of evec
c
c     It is assumed that n <= nz and n0 <= nz0
c
      SUBROUTINE e_decomp(n,a,n0,evec,eval,condo,condi,evalmax,evalmin)
#include "common.inc"
      INTEGER n,n0
      REAL*8 a(n0,n),evec(n0,n),eval(n),condo,condi,evalmax,evalmin
c
      INTEGER nwork
      REAL*8 eps
      PARAMETER (nwork=nz*(nz+4),eps=1.d-9)
c
      INTEGER i0,i1,i2,i3,i,j,ipiv(nz)
      REAL*8 val,evi(nz),work(nwork)
#ifdef COPPER
c for ESSL
      COMPLEX*16 ceval(nz),cevec(nz0,nz)
      LOGICAL alltrue(nz)
      REAL*8 dnorm2
      EXTERNAL dnorm2
#endif
#ifdef COBALT
      REAL*8 dnrm2
      EXTERNAL dnrm2
#endif
      INTEGER idamax,idmax,idmin
      REAL*8 dasum
      EXTERNAL idamax,idmax,idmin,dasum
c
      i0=0
      write(26,*) 'E_DECOMP'
#ifdef COPPER
c for ESSL
      do i=1,n
        alltrue(i) = .true.
      enddo
      call dgeev(2,a,n0,ceval,cevec,nz0,alltrue,n,work,nwork)
      do i=1,n
        eval(i) = dreal(ceval(i))
        evi(i)  = dimag(ceval(i))
      enddo
      do j=1,n
        do i=1,n
          evec(i,j) = dreal(cevec(i,j))
        enddo
        val = 1.d0/dnorm2(n,evec(1,j),1)
        call dscal(n,val,evec(1,j),1)
      enddo
#endif
#ifdef COBALT
c for SCSL
      call dgeev('N','V',n,a,n0,eval,evi,evec,n0,evec,n0,work,nwork,i0)
c      write(26,*) 'DGEEV called'
      do j=1,n
        val = 1.d0/dnrm2(n,evec(1,j),1)
        call dscal(n,val,evec(1,j),1)
      enddo
#endif
c
      i1 = 0
      val = abs(evi(idamax(n,evi,1)))
      if (val.gt.eps) i1 = 1
c
      evalmax = eval(idmax(n,eval,1))
      evalmin = eval(idmin(n,eval,1))
c
      do i=1,n
        work(i) = dasum(n,evec(1,i), 1)
        evi( i) = dasum(n,evec(i,1),n0)
      enddo
      condo = work(idmax(n,work,1))
      condi = evi( idmax(n, evi,1))
c
      call dcopy(n*n0,evec,1,a,1)
      call dgetrf(n,n,a,n0,ipiv,i2)
      call dgetri(n,a,n0,ipiv,work,nwork,i3)
c
      do i=1,n
        work(i) = dasum(n,a(1,i), 1)
        evi( i) = dasum(n,a(i,1),n0)
      enddo
      condo = condo*work(idmax(n,work,1))
      condi = condi*evi( idmax(n, evi,1))
c
      if (i0.eq.0.and.i1.eq.0.and.i2.eq.0.and.i3.eq.0) return
      write(26,*) 'error in e_decomp'
      write(26,*) 'i0, i1, i2, i3:',i0,i1,i2,i3
      print *, 'error in e_decomp'
      stop
      END
c
c     Bad for parallel
c
      SUBROUTINE init_uvw
#include "common.inc"

      INTEGER i,j,k
      INTEGER jran
      REAL*8  rnum
      REAL*8  fac,Vz,sum,sum1
      INTEGER idum
      REAL*8  dum

C     Initialize the velocity profile using laminar profile + random fluctuation
      fac  = 1.d-12
c     Setting the seed for random number generator
      jran = 123
C     The laminar profile for channel flow is
c     v = 0.5*gradpbar*Re*(1-zets(k)**2)
cz      do k = 1,nz
cz         do j = 1,ny
cz            do i = 1,nx
cz               call randomgen(jran,rnum)
cz               u(i,j,k) = 0.4d0*(rnum-0.5d0)
cz               call randomgen(jran,rnum)
cz               v(i,j,k) = 0.5d0*re*(1.d0-zets(k)**2)
cz     &                   *(1.d0+0.3d0*(rnum-0.5d0))
cz               call randomgen(jran,rnum)
cz               w(i,j,k) = 0.3d0*(rnum-0.5d0)
cz            enddo
cz         enddo
cz      enddo
      
c     No-slip boundary condition for u,v,w
      do j = 1,ny
         do i = 1,nx
            u(i,j,1)  = 0.d0
            u(i,j,nz) = 0.d0
            v(i,j,1)  = 0.d0
            v(i,j,nz) = 0.d0
            w(i,j,1)  = 0.d0
            w(i,j,nz) = 0.d0
         enddo
      enddo
     
cz      open(21,file='conc.0',form=format_mode)
cz      read (21) dum,idum,dum,idum,idum,idum,dum,dum
cz      read(21) tt
cz      close(21)
      open(22,file='vel.0',form=format_mode)
      read(22) dum,idum,dum,idum,idum,idum,dum,dum
      read(22) u,v,w
      read(22) dum,dum,dum,dum,dum,dum,dum,dum
      read(22) dum,dum,dum,dum,dum,dum,dum,dum
      close(22)
c to start with sediment-free 
      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               tt(i,j,k) = 0.00d0
            enddo
         enddo
      enddo

c     The laminar profile + Orr-Sommerfeld equation solution
C     The most unstable mode for Re = 180 is
      return
      END

c
c This randomgen subroutine can be changed with any other subroutine we want
      subroutine randomgen(iseed,rnum)
      integer iseed
      real*8 rnum
      iseed  = mod(iseed*1366+150889,714025)
      rnum   = (dble(iseed)/714025.d0) 
      return
      end
