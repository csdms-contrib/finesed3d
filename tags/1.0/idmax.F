*     BLAS-Extension function IDMAX, IDMIN, both double precision


      INTEGER FUNCTION IDMAX(N,DX,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)

*     ..
*
*     Purpose
*     =======
*     Find the position i of the first or last occurence of a vector
*     element having the maximum value.
*     You get the position of the first or last occurence of an element
*     by specifying positive or negative stride, respectively, for vector
*     x. Regardless of the stride, the position i is always relative to
*     the location specified in the calling sequence for vector x.
*
*     Syntax
*     ======
*     IDMAX(n,x,inc)
*
*     On Entry
*     n : the number of elements in vector x. specified as: an integer;
*         n greater than or equal to 0.
*     x : the vector x of length n. Specified as: a one-dimensional array
*         of (at least) length 1+(n-1)|incx|
*     incx : the stride for vector x.

*     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX

*     ..
      IDMAX = 0
      IF(N.LT.1.OR.INCX.EQ.0) RETURN
      IDMAX = 1
      IF(N.EQ.1) RETURN
      IF(INCX.EQ.1) GO TO 20
      IF(INCX.LT.0) GO TO 30

*     
*     Code for increment not equal to 1 and positive INCX
*
      IX = 1
      DMAX = DX(1)
      IX = IX + INCX
      DO 10 I=2,N
         IF (DX(IX).LE.DMAX) GO TO 5
         IDMAX = I
         DMAX  = DX(IX)
 5       IX    = IX + INCX
 10   CONTINUE
      RETURN

*
*     Code for increment equal to 1
*     
 20   DMAX = DX(1)
      DO 25 I=2,N
         IF(DX(I).LE.DMAX) GO TO 25
         IDMAX = I
         DMAX  = DX(I)
 25   CONTINUE
      RETURN

*
*     Code for negative stride
*

 30   IF(INCX.EQ.-1) GO TO 50

*
*     Code for increment not euqal to -1
*
      IX = 1
      DMAX = DX(1)
      IX = IX - INCX
      DO 40 I = 2,N
         IF(DX(IX).LT.DMAX) GO TO 45
         IDMAX = I
         DMAX  = DX(IX)
 45      IX    = IX - INCX
 40   CONTINUE
      RETURN

*
*     Code for increment equal to -1
*

 50   DMAX = DX(1)
      DO 60 I = 2,N
         IF(DX(I).LT.DMAX) GO TO 60
         IDMAX = I
         DMAX  = DX(I)
 60   CONTINUE
      RETURN
      END



      INTEGER FUNCTION IDMIN(N,DX,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)

*     ..
*
*     Purpose
*     =======
*     Find the position i of the first or last occurence of a vector
*     element having the minimum value.
*     You get the position of the first or last occurence of an element
*     by specifying positive or negative stride, respectively, for vector
*     x. Regardless of the stride, the position i is always relative to
*     the location specified in the calling sequence for vector x.
*
*     Syntax
*     ======
*     IDMIN(n,x,inc)
*
*     On Entry
*     n : the number of elements in vector x. specified as: an integer;
*         n greater than or equal to 0.
*     x : the vector x of length n. Specified as: a one-dimensional array
*         of (at least) length 1+(n-1)|incx|
*     incx : the stride for vector x.

*     .. Local Scalars ..
      DOUBLE PRECISION DMIN
      INTEGER I,IX

*     ..
      IDMIN = 0
      IF(N.LT.1.OR.INCX.EQ.0) RETURN
      IDMIN = 1
      IF(N.EQ.1) RETURN
      IF(INCX.EQ.1) GO TO 20
      IF(INCX.LT.0) GO TO 30

*     
*     Code for increment not equal to 1 and positive INCX
*
      IX = 1
      DMIN = DX(1)
      IX = IX + INCX
      DO 10 I=2,N
         IF (DX(IX).GE.DMIN) GO TO 5
         IDMIN = I
         DMIN  = DX(IX)
 5       IX    = IX + INCX
 10   CONTINUE
      RETURN

*
*     Code for increment equal to 1
*     
 20   DMIN = DX(1)
      DO 25 I=2,N
         IF(DX(I).GE.DMIN) GO TO 25
         IDMIN = I
         DMIN  = DX(I)
 25   CONTINUE
      RETURN

*
*     Code for negative stride
*

 30   IF(INCX.EQ.-1) GO TO 50

*
*     Code for increment not euqal to -1
*
      IX = 1
      DMIN = DX(1)
      IX = IX - INCX
      DO 40 I = 2,N
         IF(DX(IX).GT.DMIN) GO TO 45
         IDMIN = I
         DMIN  = DX(IX)
 45      IX    = IX - INCX
 40   CONTINUE
      RETURN

*
*     Code for increment equal to -1
*

 50   DMIN = DX(1)
      DO 60 I = 2,N
         IF(DX(I).GT.DMIN) GO TO 60
         IDMIN = I
         DMIN  = DX(I)
 60   CONTINUE
      RETURN
      END
