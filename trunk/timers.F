      SUBROUTINE init_timers
      IMPLICIT none
#include "timers.inc"
c
      INTEGER k
c
      do k=1,max_timers
        it_cnt(k) = 0
        tm_min(k) = 1.d12
        tm_max(k) = -1.d0
        tm_sum(k) = 0.d0
      enddo
c
      return
      END
c
      SUBROUTINE begin_timer(k)
      IMPLICIT none
#include "timers.inc"
      INTEGER k
c
c  rtc must be an 8-byte real
      REAL*8 rtc
      EXTERNAL rtc
c
      time1(k) = rtc()
c
      return
      END
c
      SUBROUTINE end_timer(k)
      IMPLICIT none
#include "timers.inc"
      INTEGER k
c
      REAL*8 tdiff
c  rtc must be an 8-byte real
      REAL*8 rtc
      EXTERNAL rtc
c
      it_cnt(k) = it_cnt(k) + 1
      tdiff = rtc() - time1(k)
      if (tdiff.lt.tm_min(k)) tm_min(k) = tdiff
      if (tdiff.gt.tm_max(k)) tm_max(k) = tdiff
      tm_sum(k) = tm_sum(k) + tdiff
c
      return
      END
c
      SUBROUTINE print_timers
      IMPLICIT none
#include "timers.inc"
c
      INTEGER k,n
c
      do k=1,max_timers
        n = it_cnt(k)
        if (n.gt.0) print 1000,
     &    k, n, tm_sum(k)/dble(n),tm_min(k),tm_max(k)
      enddo
c
      return
 1000 FORMAT('Timer',i3,':',i8,' calls, mean =',f8.3,', min =',f8.3,
     &  ', max =',f8.3)
      END
