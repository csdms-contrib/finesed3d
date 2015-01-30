      IMPLICIT none
#include "global.inc"
#include "fft.inc"
#include "flags.inc"
#include "timers.inc"
c
      INTEGER istart,itime
      INTEGER oi_chan,out_press, limit
      INTEGER oi_timer,oi_cfl,oi_spec,oi_mean,oi_gbal,oi_1d
      INTEGER msvx, msvy
      REAL*8 re,yenir,rlenx,rleny,deltax,deltay,idx,idy
      REAL*8 dt,dt_old,rtime,pi,ipin,cflmax,gradpbar
      REAL*8 sc,ri,Tmax,RR
      REAL*8 uvab,   uvat,   uvbb,   uvbt,    ugb,    ugt,vgb,vgt
      REAL*8  wab,    wat,    wbb,    wbt,    wgb,    wgt
      REAL*8 ttab(m),ttat(m),ttbb(m),ttbt(m),ttgb(m),ttgt(m)
      REAL*8 w_settling,tau_p,beta_p,rey_p,phiave,phind,Ero,ushearc
      REAL*8 cferos1,cferos2,cfup1,cfup2,cfup3
      REAL*8 cfnl1v(3),cfnl2v(3),cfdifv(3),cfdelv(3),cfcumv(3)
      REAL*8 zets(nz),deltaz(nzm),idz(nzm),ndeltaz(nzm),nidz(nzm),
     &  wzet(nz)
      REAL*8 wavex(nx),wavexs(nx),wavey(ny),waveys(ny)
      REAL*8 chb(nz0,nz)
      REAL*8 chbd1(nz0,nz),chbd2(nz0,nz)
      REAL*8 cn_drdr(nz0,nz),ci_drdr(nz0,nz),e_drdr(nz)
      REAL*8 uvcn_mimi(nzmm0,nzmm),uvci_mimi(nzmm0,nzmm),uve_mimi(nzmm)
      REAL*8 wcn_mimi(nzmm0,nzmm),wci_mimi(nzmm0,nzmm),we_mimi(nzmm)
      REAL*8 ttcn_mimi(nzmm0,nzmm),ttci_mimi(nzmm0,nzmm),tte_mimi(nzmm)
      REAL*8 cn_nodr(nz0,nz),ci_nodr(nz0,nz),e_nodr(nz)
      REAL*8 cn_nono(nzmm0,nzmm),ci_nono(nzmm0,nzmm),e_nono(nzmm)
      REAL*8 mlt(nz0,nx0,nyh)
      REAL*8 pbxnb(nx0,nyh),pbynb(nx0,nyh)
      REAL*8 pbxnt(nx0,nyh),pbynt(nx0,nyh)
      REAL*8 pbxob(nx0,nyh),pbyob(nx0,nyh)
      REAL*8 pbxot(nx0,nyh),pbyot(nx0,nyh)
      REAL*8 pcxnb(nx0,nyh),pcynb(nx0,nyh)
      REAL*8 pcxnt(nx0,nyh),pcynt(nx0,nyh)
      REAL*8 pcxob(nx0,nyh),pcyob(nx0,nyh)
      REAL*8 pcxot(nx0,nyh),pcyot(nx0,nyh)
      REAL*8   prbc1(nzmm),  prbcn(nzmm),
     &       uvprbc1(nzmm), uvprbcn(nzmm),
     &        wprbc1(nzmm), wprbcn(nzmm),
     &       ttprbc1(nzmm),ttprbcn(nzmm)
      REAL*8 u(nx0,ny0,nz),v(nx0,ny0,nz),w(nx0,ny0,nz),p(nx0,ny0,nz)
      REAL*8 up(nx0,ny0,nz),vp(nx0,ny0,nz),wp(nx0,ny0,nz)
      REAL*8 ru(nx0,ny0,nz),rv(nx0,ny0,nz),rw(nx0,ny0,nz),rp(nx0,ny0,nz)
      REAL*8 h1(nx0,ny0,nz),h2(nx0,ny0,nz),h3(nx0,ny0,nz)
c      REAL*8 g(nx0,ny0,nz)
      REAL*8 ou1(nx0,ny0,nz),ov1(nx0,ny0,nz),ow1(nx0,ny0,nz)
      REAL*8 ou2(nx0,ny0,nz),ov2(nx0,ny0,nz),ow2(nx0,ny0,nz)
      REAL*8 ou3(nx0,ny0,nz),ov3(nx0,ny0,nz),ow3(nx0,ny0,nz)
      REAL*8 DuDt(nx0,ny0,nz),DvDt(nx0,ny0,nz),DwDt(nx0,ny0,nz)
      REAL*8 extmx(nx0,ny0,nz),extmy(nx0,ny0,nz),extmz(nx0,ny0,nz)
      REAL*8 auxt1(nz0,nx0,nyh),auxt2(nz0,nx0,nyh)
      REAL*8 tt(nx0,ny0,nz),rtt(nx0,ny0,nz),htt(nx0,ny0,nz)
      REAL*8 ul(nxl0,nyl0,nz),vl(nxl0,nyl0,nz),wl(nxl0,nyl0,nz)
      REAL*8 ddxul(nxl0,nyl0,nz),ddyul(nxl0,nyl0,nz),ddzul(nxl0,nyl0,nz)
      REAL*8 svx(nx),svy(nyh)
      REAL*8 epssvx,epssvy
      REAL*8 deposit(nx0,ny0),erosion(nx0,ny0)
      REAL*8 ttotn(nx0,ny0),ttobn(nx0,ny0)
      REAL*8 ttoto(nx0,ny0),ttobo(nx0,ny0)
      REAL*8 ttot(nx0,ny0),ttob(nx0,ny0)
      REAL*8 ttotl(nxl0,nyl0),ttobl(nxl0,nyl0)
      REAL*8 dwbct(nx0,ny0),dwbcb(nx0,ny0)	
      REAL*8 dwbctl(nxl0,nyl0),dwbcbl(nxl0,nyl0)
      REAL*8 wsh(nx0,ny0),ttwb(nx0,ny0),auxdep(nx0,ny0)
      REAL*8 wshl(nxl0,nyl0),ttwbl(nxl0,nyl0),auxdepl(nxl0,nyl0)
      REAL*8 ubct(nx0,ny0),ubcb(nx0,ny0),
     &       ubctax(nx0,ny0),ubcbax(nx0,ny0)
      REAL*8 vbct(nx0,ny0),vbcb(nx0,ny0),
     &       vbctax(nx0,ny0),vbcbax(nx0,ny0),
     &       uvmbct1,uvmbcb1,uvmbctn,uvmbcbn
      REAL*8 wbct(nx0,ny0),wbcb(nx0,ny0),
     &       wmbct1,wmbcb1,wmbctn,wmbcbn
      REAL*8 ttbct(nx0,ny0),ttbcb(nx0,ny0),
     &       ttbcber(nx0,ny0),
     &       ttmbct1,ttmbcb1,ttmbctn,ttmbcbn
      REAL*8 gy(nyl0,nz),gz(nyl0,nz)
      CHARACTER*32 ch_file,tt_file,format_mode
c
      INTEGER NT
      REAL*8  ust2,usb2,vst2,vsb2
      REAL*8  um(nz),vm(nz),wm(nz),ttm(nz)
     &       ,u2m(nz),v2m(nz),w2m(nz),tt2m(nz)
     &       ,u3m(nz),v3m(nz),w3m(nz),tt3m(nz)
     &       ,u4m(nz),v4m(nz),w4m(nz),tt4m(nz)
     &       ,uvm(nz),uwm(nz),uttm(nz)
     &       ,vwm(nz),vttm(nz)
     &       ,wttm(nz)
      REAL*8  uxm( nz),uym( nz),uzm( nz)
     &       ,ux2m(nz),uy2m(nz),uz2m(nz)
     &       ,ux3m(nz),uy3m(nz),uz3m(nz)
     &       ,ux4m(nz),uy4m(nz),uz4m(nz)
      REAL*8  vxm( nz),vym( nz),vzm( nz)
     &       ,vx2m(nz),vy2m(nz),vz2m(nz)
     &       ,vx3m(nz),vy3m(nz),vz3m(nz)
     &       ,vx4m(nz),vy4m(nz),vz4m(nz)
      REAL*8  wxm( nz),wym( nz),wzm( nz)
     &       ,wx2m(nz),wy2m(nz),wz2m(nz)
     &       ,wx3m(nz),wy3m(nz),wz3m(nz)
     &       ,wx4m(nz),wy4m(nz),wz4m(nz)
      REAL*8  ttxm( nz),ttym( nz),ttzm( nz)
     &       ,ttx2m(nz),tty2m(nz),ttz2m(nz)
     &       ,ttx3m(nz),tty3m(nz),ttz3m(nz)
     &       ,ttx4m(nz),tty4m(nz),ttz4m(nz)
      REAL*8 CA1(3,nz),CB1(11,nz),CA2(3,nz),CB2(11,nz)
      REAL*8 CA2_drdr(3,nz),CB2_drdr(11,nz)
      REAL*8 CA2_nono(3,nz),CB2_nono(11,nz)
      REAL*8 CA2_nodr(3,nz),CB2_nodr(11,nz)
      REAL*8 CA2_mimi(3,nz),CB2_mimi(11,nz)
      REAL*8 d1l(nzm),d1(nz),d1u(nzm)
      REAL*8 d2l(nzm),d2(nz),d2u(nzm)
c
      COMMON /iparm/   istart,itime
      COMMON /oiparm/  oi_chan,oi_timer,out_press
     &                ,oi_cfl,oi_spec,oi_mean,oi_gbal,oi_1d
      COMMON /rparm1/  re,rlenx,rleny,deltax,deltay,idx,idy
      COMMON /rparm2/  dt,dt_old,rtime,pi,ipin,cflmax,gradpbar
      COMMON /rparmt/  sc,ri,Tmax,RR
      COMMON /abg/     uvab,uvat,uvbb,uvbt, ugb, ugt,
     &                                      vgb, vgt,
     &                  wab, wat, wbb, wbt, wgb, wgt,
     &                 ttab,ttat,ttbb,ttbt,ttgb,ttgt
      COMMON /rparmw/  w_settling,tau_p,beta_p,rey_p,phiave,phind
     &                 ,cferos1,cferos2,cfup1,cfup2,cfup3,Ero,ushearc
      COMMON /rkcv/    cfnl1v,cfnl2v,cfdifv,cfdelv,cfcumv
      COMMON /zet/     zets,deltaz,idz,ndeltaz,nidz,wzet
      COMMON /wv/      wavex,wavexs,wavey,waveys
      COMMON /chb/     chb
      COMMON /chbdx/   chbd1,chbd2
      COMMON /chndd/     cn_drdr,  ci_drdr,  e_drdr,
     &                   cn_nodr,  ci_nodr,  e_nodr,
     &                   cn_nono,  ci_nono,  e_nono,
     &                 uvcn_mimi,uvci_mimi,uve_mimi,
     &                  wcn_mimi, wci_mimi, we_mimi,
     &                 ttcn_mimi,ttci_mimi,tte_mimi
      COMMON /mltm/    mlt
      COMMON /pbot/    pbxnb,pcxnb,pbynb,pcynb,pbxob,pcxob,pbyob,pcyob
      COMMON /ptop/    pbxnt,pcxnt,pbynt,pcynt,pbxot,pcxot,pbyot,pcyot
      COMMON /prbc/      prbc1,  prbcn,
     &                 uvprbc1,uvprbcn,
     &                  wprbc1, wprbcn,
     &                 ttprbc1,ttprbcn
      COMMON /veluvw/  u,v,w
      COMMON /press/   p,rp
      COMMON /veluvwp/ up,vp,wp
      COMMON /velruvw/ ru,rv,rw
      COMMON /adam/    h1,h2,h3
c      COMMON /large1g/ g
      COMMON /oldvelu/ ou1,ou2,ou3
      COMMON /oldvelv/ ov1,ov2,ov3
      COMMON /oldvelw/ ow1,ow2,ow3
      COMMON /extratms/ extmx,extmy,extmz
      COMMON /matdev/  DuDt,DvDt,DwDt
      COMMON /transp/  auxt1,auxt2
      COMMON /thrm0/   tt,rtt,htt
      COMMON /files/   ch_file,tt_file,format_mode
      COMMON /velul/   ul
      COMMON /velvl/   vl
      COMMON /velwl/   wl
      COMMON /veldxul/ ddxul
      COMMON /veldyul/ ddyul
      COMMON /veldzul/ ddzul
      COMMON /spv1/    svx,svy
      COMMON /spv2/    msvx,msvy
      COMMON /spv3/    epssvx,epssvy
      COMMON /dep/     deposit,erosion
      COMMON /bc_v/    ubct,ubcb,ubctax,ubcbax,
     &                 vbct,vbcb,vbctax,vbcbax,
     &                 uvmbct1,uvmbcb1,uvmbctn,uvmbcbn
      COMMON /bc_w/    wbct,wbcb,wmbct1,wmbcb1,wmbctn,wmbcbn
      COMMON /bc_tt/   ttbct,ttbcb,ttbcber,ttmbct1,ttmbcb1,ttmbctn,
     &                 ttmbcbn
      COMMON /grav/    gy,gz
      COMMON /avgtime/ NT
      COMMON /ustar/   ust2,usb2,vst2,vsb2
      COMMON /avgvar/  um,vm,wm,ttm
     &                ,u2m,v2m,w2m,tt2m
     &                ,u3m,v3m,w3m,tt3m
     &                ,u4m,v4m,w4m,tt4m
     &                ,uvm,uwm,uttm
     &                ,vwm,vttm
     &                ,wttm
      COMMON /avgvar1/ uxm,uym,uzm
     &                ,ux2m,uy2m,uz2m
     &                ,ux3m,uy3m,uz3m
     &                ,ux4m,uy4m,uz4m
      COMMON /avgvar2/ vxm,vym,vzm
     &                ,vx2m,vy2m,vz2m
     &                ,vx3m,vy3m,vz3m
     &                ,vx4m,vy4m,vz4m
      COMMON /avgvar3/ wxm,wym,wzm
     &                ,wx2m,wy2m,wz2m
     &                ,wx3m,wy3m,wz3m
     &                ,wx4m,wy4m,wz4m
      COMMON /avgvar4/ ttxm,ttym,ttzm
     &                ,ttx2m,tty2m,ttz2m
     &                ,ttx3m,tty3m,ttz3m
     &                ,ttx4m,tty4m,ttz4m
      COMMON /cpt/ CA1,CA2,CB1,CB2,CA2_drdr,CB2_drdr,CA2_nono,
     &	           CB2_nono,CA2_mimi,CB2_mimi,CA2_nodr,CB2_nodr,
     &		   d1l,d1,d1u,d2l,d2,d2u
      COMMON /bct/ ttotn,ttoto,ttobn,ttobo,dwbct,dwbcb,
     &	     ttotl,ttobl,dwbctl,dwbcbl,wshl,ttwbl,auxdepl,
     &       ttot,ttob, wsh,ttwb,auxdep
