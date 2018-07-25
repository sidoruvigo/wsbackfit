     MODULE lsq

     !  Module for unconstrained linear least-squares calculations.
     !  The algorithm is suitable for updating LS calculations as more
     !  data are added.   This is sometimes called recursive estimation.
     !  Only one dependent variable is allowed.
     !  Based upon Applied Statistics algorithm AS 274.
     !  Translation from Fortran 77 to Fortran 90 by Alan Miller.
     !  A function, VARPRD, has been added for calculating the variances
     !  of predicted values, and this uses a subroutine BKSUB2.

     !  Version 1.14, 19 August 2002 - ELF90 compatible version
     !  Author: Alan Miller
     !  e-mail : amiller @ bigpond.net.au
     !  WWW-pages: http://www.ozemail.com.au/~milleraj
     !             http://users.bigpond.net.au/amiller/

     !  Bug fixes:
     !  1. In REGCF a call to TOLSET has been added in case the user had
     !     not set tolerances.
     !  2. In SING, each time a singularity is detected, unless it is in the
     !     variables in the last position, INCLUD is called.   INCLUD assumes
     !     that a new observation is being added and increments the number of
     !     cases, NOBS.   The line:  nobs = nobs - 1 has been added.
     !  3. row_ptr was left out of the DEALLOCATE statement in routine startup
     !     in version 1.07.
     !  4. In COV, now calls SS if rss_set = .FALSE.  29 August 1997
     !  5. In TOLSET, correction to accomodate negative values of D.  19 August 2002

     !  Other changes:
     !  1. Array row_ptr added 18 July 1997.   This points to the first element
     !     stored in each row thus saving a small amount of time needed to
     !     calculate its position.
     !  2. Optional parameter, EPS, added to routine TOLSET, so that the user
     !     can specify the accuracy of the input data.
     !  3. Cosmetic change of lsq_kind to dp (`Double precision')
     !  4. Change to routine SING to use row_ptr rather than calculate the position
     !     of first elements in each row.

     !  The PUBLIC variables are:
     !  dp       = a KIND parameter for the floating-point quantities calculated
     !             in this module.   See the more detailed explanation below.
     !             This KIND parameter should be used for all floating-point
     !             arguments passed to routines in this module.

     !  nobs    = the number of observations processed to date.
     !  ncol    = the total number of variables, including one for the constant,
     !            if a constant is being fitted.
     !  r_dim   = the dimension of array r = ncol*(ncol-1)/2
     !  vorder  = an integer vector storing the current order of the variables
     !            in the QR-factorization.   The initial order is 0, 1, 2, ...
     !            if a constant is being fitted, or 1, 2, ... otherwise.
     !  initialized = a logical variable which indicates whether space has
     !                been allocated for various arrays.
     !  tol_set = a logical variable which is set when subroutine TOLSET has
     !            been called to calculate tolerances for use in testing for
     !            singularities.
     !  rss_set = a logical variable indicating whether residual sums of squares
     !            are available and usable.
     !  d()     = array of row multipliers for the Cholesky factorization.
     !            The factorization is X = Q.sqrt(D).R where Q is an ortho-
     !            normal matrix which is NOT stored, D is a diagonal matrix
     !            whose diagonal elements are stored in array d, and R is an
     !            upper-triangular matrix with 1's as its diagonal elements.
     !  rhs()   = vector of RHS projections (after scaling by sqrt(D)).
     !            Thus Q'y = sqrt(D).rhs
     !  r()     = the upper-triangular matrix R.   The upper triangle only,
     !            excluding the implicit 1's on the diagonal, are stored by
     !            rows.
     !  tol()   = array of tolerances used in testing for singularities.
     !  rss()   = array of residual sums of squares.   rss(i) is the residual
     !            sum of squares with the first i variables in the model.
     !            By changing the order of variables, the residual sums of
     !            squares can be found for all possible subsets of the variables.
     !            The residual sum of squares with NO variables in the model,
     !            that is the total sum of squares of the y-values, can be
     !            calculated as rss(1) + d(1)*rhs(1)^2.   If the first variable
     !            is a constant, then rss(1) is the sum of squares of
     !            (y - ybar) where ybar is the average value of y.
     !  sserr   = residual sum of squares with all of the variables included.
     !  row_ptr() = array of indices of first elements in each row of R.
     !
     !--------------------------------------------------------------------------

     !     General declarations

     IMPLICIT NONE

     INTEGER, SAVE                :: nobs, ncol, r_dim
     INTEGER, ALLOCATABLE, SAVE   :: vorder(:), row_ptr(:)
     LOGICAL, SAVE                :: initialized = .false.,                  &
                                     tol_set = .false., rss_set = .false.

     ! Note. dp is being set to give at least 12 decimal digit
     !       representation of floating point numbers.   This should be adequate
     !       for most problems except the fitting of polynomials.   dp is
     !       being set so that the same code can be run on PCs and Unix systems,
     !       which will usually represent floating-point numbers in `double
     !       precision', and other systems with larger word lengths which will
     !       give similar accuracy in `single precision'.

     INTEGER, PARAMETER           :: dp = SELECTED_REAL_KIND(12,60)
     double precision, ALLOCATABLE, SAVE :: d(:), rhs(:), r(:), tol(:), rss(:)
     double precision, SAVE              :: zero = 0.0_dp, one = 1.0_dp, vsmall
     double precision, SAVE              :: sserr, toly

     PUBLIC                       :: dp, nobs, ncol, r_dim, vorder, row_ptr, &
                                     initialized, tol_set, rss_set,          &
                                     d, rhs, r, tol, rss, sserr
     PRIVATE                      :: zero, one, vsmall


     CONTAINS

     SUBROUTINE startup(nvar, fit_const)

     !     Allocates dimensions for arrays and initializes to zero
     !     The calling program must set nvar = the number of variables, and
     !     fit_const = .true. if a constant is to be included in the model,
     !     otherwise fit_const = .false.
     !
     !--------------------------------------------------------------------------

     IMPLICIT NONE
     INTEGER, INTENT(IN)  :: nvar
     LOGICAL, INTENT(IN)  :: fit_const

     !     Local variable
     INTEGER   :: i

     vsmall = 10. * TINY(zero)

     nobs = 0
     IF (fit_const) THEN
       ncol = nvar + 1
     ELSE
       ncol = nvar
     END IF

     IF (initialized) DEALLOCATE(d, rhs, r, tol, rss, vorder, row_ptr)
     r_dim = ncol * (ncol - 1)/2
     ALLOCATE( d(ncol), rhs(ncol), r(r_dim), tol(ncol), rss(ncol), vorder(ncol),  &
               row_ptr(ncol) )

     d = zero
     rhs = zero
     r = zero
     sserr = zero

     IF (fit_const) THEN
       DO i = 1, ncol
         vorder(i) = i-1
       END DO
     ELSE
       DO i = 1, ncol
         vorder(i) = i
       END DO
     END IF ! (fit_const)

     ! row_ptr(i) is the position of element R(i,i+1) in array r().

     row_ptr(1) = 1
     DO i = 2, ncol-1
       row_ptr(i) = row_ptr(i-1) + ncol - i + 1
     END DO
     row_ptr(ncol) = 0

     initialized = .true.
     tol_set = .false.
     rss_set = .false.

     RETURN
     END SUBROUTINE startup

     !     ***********************************************************************************
     !     ***********************************************************************************
     SUBROUTINE endup()
          IF (initialized) DEALLOCATE(d, rhs, r, tol, rss, &
                                       vorder, row_ptr)
          initialized = .FALSE.
     RETURN
     END SUBROUTINE endup


     SUBROUTINE includ(weight, xrow, yelem)

     !     ALGORITHM AS75.1  APPL. STATIST. (1974) VOL.23, NO. 3

     !     Calling this routine updates D, R, RHS and SSERR by the
     !     inclusion of xrow, yelem with the specified weight.

     !     *** WARNING  Array XROW is overwritten.

     !     N.B. As this routine will be called many times in most applications,
     !          checks have been eliminated.
     !
     !--------------------------------------------------------------------------


     IMPLICIT NONE
     double precision,INTENT(IN)                    :: weight, yelem
     double precision, DIMENSION(:), INTENT(IN OUT) :: xrow

     !     Local variables

     INTEGER     :: i, k, nextr
     double precision   :: w, y, xi, di, wxi, dpi, cbar, sbar, xk

     nobs = nobs + 1
     w = weight
     y = yelem
     rss_set = .false.
     nextr = 1
     DO i = 1, ncol

     !     Skip unnecessary transformations.   Test on exact zeroes must be
     !     used or stability can be destroyed.

       IF (ABS(w) < vsmall) RETURN
       xi = xrow(i)
       IF (ABS(xi) < vsmall) THEN
         nextr = nextr + ncol - i
       ELSE
         di = d(i)
         wxi = w * xi
         dpi = di + wxi*xi
         cbar = di / dpi
         sbar = wxi / dpi
         w = cbar * w
         d(i) = dpi
         DO k = i+1, ncol
           xk = xrow(k)
           xrow(k) = xk - xi * r(nextr)
           r(nextr) = cbar * r(nextr) + sbar * xk
           nextr = nextr + 1
         END DO
         xk = y
         y = xk - xi * rhs(i)
         rhs(i) = cbar * rhs(i) + sbar * xk
       END IF
     END DO ! i = 1, ncol

     !     Y * SQRT(W) is now equal to the Brown, Durbin & Evans recursive
     !     residual.

     sserr = sserr + w * y * y

     RETURN
     END SUBROUTINE includ



     SUBROUTINE regcf(beta, nreq, ifault)

     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

     !     Modified version of AS75.4 to calculate regression coefficients
     !     for the first NREQ variables, given an orthogonal reduction from
     !     AS75.1.
     !
     !--------------------------------------------------------------------------

     IMPLICIT NONE
     INTEGER, INTENT(IN)                  :: nreq
     INTEGER, INTENT(OUT)                 :: ifault
     double precision, DIMENSION(:), INTENT(OUT) :: beta

     !     Local variables

     INTEGER   :: i, j, nextr

     !     Some checks.

     ifault = 0
     IF (nreq < 1 .OR. nreq > ncol) ifault = ifault + 4
     IF (ifault /= 0) RETURN

     IF (.NOT. tol_set) CALL tolset()

     DO i = nreq, 1, -1
       IF (SQRT(d(i)) < tol(i)) THEN
         beta(i) = zero
         d(i) = zero
         ifault = -i
       ELSE
         beta(i) = rhs(i)
         nextr = row_ptr(i)
         DO j = i+1, nreq
           beta(i) = beta(i) - r(nextr) * beta(j)
           nextr = nextr + 1
         END DO ! j = i+1, nreq
       END IF
     END DO ! i = nreq, 1, -1

     RETURN
     END SUBROUTINE regcf



     SUBROUTINE tolset(eps)

     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

     !     Sets up array TOL for testing for zeroes in an orthogonal
     !     reduction formed using AS75.1.

     double precision, INTENT(IN), OPTIONAL :: eps

     !     Unless the argument eps is set, it is assumed that the input data are
     !     recorded to full machine accuracy.   This is often not the case.
     !     If, for instance, the data are recorded to `single precision' of about
     !     6-7 significant decimal digits, then singularities will not be detected.
     !     It is suggested that in this case eps should be set equal to
     !     10.0 * EPSILON(1.0)
     !     If the data are recorded to say 4 significant decimals, then eps should
     !     be set to 1.0E-03
     !     The above comments apply to the predictor variables, not to the
     !     dependent variable.

     !     Correction - 19 August 2002
     !     When negative weights are used, it is possible for an alement of D
     !     to be negative.

     !     Local variables.
     !
     !--------------------------------------------------------------------------

     !     Local variables

     INTEGER    :: col, row, pos
     double precision  :: eps1, ten = 10.0, total, work(ncol)

     !     EPS is a machine-dependent constant.

     IF (PRESENT(eps)) THEN
       eps1 = MAX(ABS(eps), ten * EPSILON(ten))
     ELSE
       eps1 = ten * EPSILON(ten)
     END IF

     !     Set tol(i) = sum of absolute values in column I of R after
     !     scaling each element by the square root of its row multiplier,
     !     multiplied by EPS1.

     work = SQRT(ABS(d))
     DO col = 1, ncol
       pos = col - 1
       total = work(col)
       DO row = 1, col-1
         total = total + ABS(r(pos)) * work(row)
         pos = pos + ncol - row - 1
       END DO
       tol(col) = eps1 * total
     END DO

     tol_set = .TRUE.
     RETURN
     END SUBROUTINE tolset




     SUBROUTINE sing(lindep, ifault)

     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

     !     Checks for singularities, reports, and adjusts orthogonal
     !     reductions produced by AS75.1.

     !     Correction - 19 August 2002
     !     When negative weights are used, it is possible for an alement of D
     !     to be negative.

     !     Auxiliary routines called: INCLUD, TOLSET
     !
     !--------------------------------------------------------------------------

     INTEGER, INTENT(OUT)                :: ifault
     LOGICAL, DIMENSION(:), INTENT(OUT)  :: lindep

     !     Local variables

     double precision  :: temp, x(ncol), work(ncol), y, weight
     INTEGER    :: pos, row, pos2

     ifault = 0

     work = SQRT(ABS(d))
     IF (.NOT. tol_set) CALL tolset()

     DO row = 1, ncol
       temp = tol(row)
       pos = row_ptr(row)         ! pos = location of first element in row

     !     If diagonal element is near zero, set it to zero, set appropriate
     !     element of LINDEP, and use INCLUD to augment the projections in
     !     the lower rows of the orthogonalization.

       lindep(row) = .FALSE.
       IF (work(row) <= temp) THEN
         lindep(row) = .TRUE.
         ifault = ifault - 1
         IF (row < ncol) THEN
           pos2 = pos + ncol - row - 1
           x = zero
           x(row+1:ncol) = r(pos:pos2)
           y = rhs(row)
           weight = d(row)
           r(pos:pos2) = zero
           d(row) = zero
           rhs(row) = zero
           CALL includ(weight, x, y)
                                  ! INCLUD automatically increases the number
                                  ! of cases each time it is called.
           nobs = nobs - 1
         ELSE
           sserr = sserr + d(row) * rhs(row)**2
         END IF ! (row < ncol)
       END IF ! (work(row) <= temp)
     END DO ! row = 1, ncol

     RETURN
     END SUBROUTINE sing



     SUBROUTINE ss()

     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

     !     Calculates partial residual sums of squares from an orthogonal
     !     reduction from AS75.1.
     !
     !--------------------------------------------------------------------------

     !     Local variables

     INTEGER    :: i
     double precision  :: total

     total = sserr
     rss(ncol) = sserr
     DO i = ncol, 2, -1
       total = total + d(i) * rhs(i)**2
       rss(i-1) = total
     END DO

     rss_set = .TRUE.
     RETURN
     END SUBROUTINE ss



     SUBROUTINE cov(nreq, var, covmat, dimcov, sterr, ifault)

     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

     !     Calculate covariance matrix for regression coefficients for the
     !     first nreq variables, from an orthogonal reduction produced from
     !     AS75.1.

     !     Auxiliary routine called: INV
     !
     !--------------------------------------------------------------------------

     INTEGER, INTENT(IN)                   :: nreq, dimcov
     INTEGER, INTENT(OUT)                  :: ifault
     double precision, INTENT(OUT)                :: var
     double precision, DIMENSION(:), INTENT(OUT)  :: covmat, sterr

     !     Local variables.

     INTEGER                :: dim_rinv, pos, row, start, pos2, col, pos1, k
     double precision              :: total
     double precision, ALLOCATABLE :: rinv(:)

     !     Check that dimension of array covmat is adequate.

     IF (dimcov < nreq*(nreq+1)/2) THEN
       ifault = 1
       RETURN
     END IF

     !     Check for small or zero multipliers on the diagonal.

     ifault = 0
     DO row = 1, nreq
       IF (ABS(d(row)) < vsmall) ifault = -row
     END DO
     IF (ifault /= 0) RETURN

     !     Calculate estimate of the residual variance.

     IF (nobs > nreq) THEN
       IF (.NOT. rss_set) CALL ss()
       var = rss(nreq) / (nobs - nreq)
     ELSE
       ifault = 2
       RETURN
     END IF

     dim_rinv = nreq*(nreq-1)/2
     ALLOCATE ( rinv(dim_rinv) )

     CALL INV(nreq, rinv)
     pos = 1
     start = 1
     DO row = 1, nreq
       pos2 = start
       DO col = row, nreq
         pos1 = start + col - row
         IF (row == col) THEN
           total = one / d(col)
         ELSE
           total = rinv(pos1-1) / d(col)
         END IF
         DO K = col+1, nreq
           total = total + rinv(pos1) * rinv(pos2) / d(k)
           pos1 = pos1 + 1
           pos2 = pos2 + 1
         END DO ! K = col+1, nreq
         covmat(pos) = total * var
         IF (row == col) sterr(row) = SQRT(covmat(pos))
         pos = pos + 1
       END DO ! col = row, nreq
       start = start + nreq - row
     END DO ! row = 1, nreq

     DEALLOCATE(rinv)
     RETURN
     END SUBROUTINE cov



     SUBROUTINE inv(nreq, rinv)

     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

     !     Invert first nreq rows and columns of Cholesky factorization
     !     produced by AS 75.1.
     !
     !--------------------------------------------------------------------------

     INTEGER, INTENT(IN)                  :: nreq
     double precision, DIMENSION(:), INTENT(OUT) :: rinv

     !     Local variables.

     INTEGER    :: pos, row, col, start, k, pos1, pos2
     double precision  :: total

     !     Invert R ignoring row multipliers, from the bottom up.

     pos = nreq * (nreq-1)/2
     DO row = nreq-1, 1, -1
       start = row_ptr(row)
       DO col = nreq, row+1, -1
         pos1 = start
         pos2 = pos
         total = zero
         DO k = row+1, col-1
           pos2 = pos2 + nreq - k
           total = total - r(pos1) * rinv(pos2)
           pos1 = pos1 + 1
         END DO ! k = row+1, col-1
         rinv(pos) = total - r(pos1)
         pos = pos - 1
       END DO ! col = nreq, row+1, -1
     END DO ! row = nreq-1, 1, -1

     RETURN
     END SUBROUTINE inv



     SUBROUTINE partial_corr(in, cormat, dimc, ycorr, ifault)

     !     Replaces subroutines PCORR and COR of:
     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

     !     Calculate partial correlations after the variables in rows
     !     1, 2, ..., IN have been forced into the regression.
     !     If IN = 1, and the first row of R represents a constant in the
     !     model, then the usual simple correlations are returned.

     !     If IN = 0, the value returned in array CORMAT for the correlation
     !     of variables Xi & Xj is:
     !       sum ( Xi.Xj ) / Sqrt ( sum (Xi^2) . sum (Xj^2) )

     !     On return, array CORMAT contains the upper triangle of the matrix of
     !     partial correlations stored by rows, excluding the 1's on the diagonal.
     !     e.g. if IN = 2, the consecutive elements returned are:
     !     (3,4) (3,5) ... (3,ncol), (4,5) (4,6) ... (4,ncol), etc.
     !     Array YCORR stores the partial correlations with the Y-variable
     !     starting with YCORR(IN+1) = partial correlation with the variable in
     !     position (IN+1).
     !
     !--------------------------------------------------------------------------

     INTEGER, INTENT(IN)                  :: in, dimc
     INTEGER, INTENT(OUT)                 :: ifault
     double precision, DIMENSION(:), INTENT(OUT) :: cormat, ycorr

     !     Local variables.

     INTEGER    :: base_pos, pos, row, col, col1, col2, pos1, pos2
     double precision  :: rms(in+1:ncol), sumxx, sumxy, sumyy, work(in+1:ncol)

     !     Some checks.

     ifault = 0
     IF (in < 0 .OR. in > ncol-1) ifault = ifault + 4
     IF (dimc < (ncol-in)*(ncol-in-1)/2) ifault = ifault + 8
     IF (ifault /= 0) RETURN

     !     Base position for calculating positions of elements in row (IN+1) of R.

     base_pos = in*ncol - (in+1)*(in+2)/2

     !     Calculate 1/RMS of elements in columns from IN to (ncol-1).

     IF (d(in+1) > zero) rms(in+1) = one / SQRT(d(in+1))
     DO col = in+2, ncol
       pos = base_pos + col
       sumxx = d(col)
       DO row = in+1, col-1
         sumxx = sumxx + d(row) * r(pos)**2
         pos = pos + ncol - row - 1
       END DO ! row = in+1, col-1
       IF (sumxx > zero) THEN
         rms(col) = one / SQRT(sumxx)
       ELSE
         rms(col) = zero
         ifault = -col
       END IF ! (sumxx > zero)
     END DO ! col = in+1, ncol-1

     !     Calculate 1/RMS for the Y-variable

     sumyy = sserr
     DO row = in+1, ncol
       sumyy = sumyy + d(row) * rhs(row)**2
     END DO ! row = in+1, ncol
     IF (sumyy > zero) sumyy = one / SQRT(sumyy)

     !     Calculate sums of cross-products.
     !     These are obtained by taking dot products of pairs of columns of R,
     !     but with the product for each row multiplied by the row multiplier
     !     in array D.

     pos = 1
     DO col1 = in+1, ncol
       sumxy = zero
       work(col1+1:ncol) = zero
       pos1 = base_pos + col1
       DO row = in+1, col1-1
         pos2 = pos1 + 1
         DO col2 = col1+1, ncol
           work(col2) = work(col2) + d(row) * r(pos1) * r(pos2)
           pos2 = pos2 + 1
         END DO ! col2 = col1+1, ncol
         sumxy = sumxy + d(row) * r(pos1) * rhs(row)
         pos1 = pos1 + ncol - row - 1
       END DO ! row = in+1, col1-1

     !     Row COL1 has an implicit 1 as its first element (in column COL1)

       pos2 = pos1 + 1
       DO col2 = col1+1, ncol
         work(col2) = work(col2) + d(col1) * r(pos2)
         pos2 = pos2 + 1
         cormat(pos) = work(col2) * rms(col1) * rms(col2)
         pos = pos + 1
       END DO ! col2 = col1+1, ncol
       sumxy = sumxy + d(col1) * rhs(col1)
       ycorr(col1) = sumxy * rms(col1) * sumyy
     END DO ! col1 = in+1, ncol-1

     ycorr(1:in) = zero

     RETURN
     END SUBROUTINE partial_corr




     SUBROUTINE vmove(from, to, ifault)

     !     ALGORITHM AS274 APPL. STATIST. (1992) VOL.41, NO. 2

     !     Move variable from position FROM to position TO in an
     !     orthogonal reduction produced by AS75.1.
     !
     !--------------------------------------------------------------------------

     INTEGER, INTENT(IN)    :: from, to
     INTEGER, INTENT(OUT)   :: ifault

     !     Local variables

     double precision  :: d1, d2, x, d1new, d2new, cbar, sbar, y
     INTEGER    :: m, first, last, inc, m1, m2, mp1, col, pos, row

     !     Check input parameters

     ifault = 0
     IF (from < 1 .OR. from > ncol) ifault = ifault + 4
     IF (to < 1 .OR. to > ncol) ifault = ifault + 8
     IF (ifault /= 0) RETURN

     IF (from == to) RETURN

     IF (.NOT. rss_set) CALL ss()

     IF (from < to) THEN
       first = from
       last = to - 1
       inc = 1
     ELSE
       first = from - 1
       last = to
       inc = -1
     END IF

     DO m = first, last, inc

     !     Find addresses of first elements of R in rows M and (M+1).

       m1 = row_ptr(m)
       m2 = row_ptr(m+1)
       mp1 = m + 1
       d1 = d(m)
       d2 = d(mp1)

     !     Special cases.

       IF (d1 < vsmall .AND. d2 < vsmall) GO TO 40
       x = r(m1)
       IF (ABS(x) * SQRT(d1) < tol(mp1)) THEN
         x = zero
       END IF
       IF (d1 < vsmall .OR. ABS(x) < vsmall) THEN
         d(m) = d2
         d(mp1) = d1
         r(m1) = zero
         DO col = m+2, ncol
           m1 = m1 + 1
           x = r(m1)
           r(m1) = r(m2)
           r(m2) = x
           m2 = m2 + 1
         END DO ! col = m+2, ncol
         x = rhs(m)
         rhs(m) = rhs(mp1)
         rhs(mp1) = x
         GO TO 40
       ELSE IF (d2 < vsmall) THEN
         d(m) = d1 * x**2
         r(m1) = one / x
         r(m1+1:m1+ncol-m-1) = r(m1+1:m1+ncol-m-1) / x
         rhs(m) = rhs(m) / x
         GO TO 40
       END IF

     !     Planar rotation in regular case.

       d1new = d2 + d1*x**2
       cbar = d2 / d1new
       sbar = x * d1 / d1new
       d2new = d1 * cbar
       d(m) = d1new
       d(mp1) = d2new
       r(m1) = sbar
       DO col = m+2, ncol
         m1 = m1 + 1
         y = r(m1)
         r(m1) = cbar*r(m2) + sbar*y
         r(m2) = y - x*r(m2)
         m2 = m2 + 1
       END DO ! col = m+2, ncol
       y = rhs(m)
       rhs(m) = cbar*rhs(mp1) + sbar*y
       rhs(mp1) = y - x*rhs(mp1)

     !     Swap columns M and (M+1) down to row (M-1).

       40 pos = m
       DO row = 1, m-1
         x = r(pos)
         r(pos) = r(pos-1)
         r(pos-1) = x
         pos = pos + ncol - row - 1
       END DO ! row = 1, m-1

     !     Adjust variable order (VORDER), the tolerances (TOL) and
     !     the vector of residual sums of squares (RSS).

       m1 = vorder(m)
       vorder(m) = vorder(mp1)
       vorder(mp1) = m1
       x = tol(m)
       tol(m) = tol(mp1)
       tol(mp1) = x
       rss(m) = rss(mp1) + d(mp1) * rhs(mp1)**2
     END DO

     RETURN
     END SUBROUTINE vmove



     SUBROUTINE reordr(list, n, pos1, ifault)

     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

     !     Re-order the variables in an orthogonal reduction produced by
     !     AS75.1 so that the N variables in LIST start at position POS1,
     !     though will not necessarily be in the same order as in LIST.
     !     Any variables in VORDER before position POS1 are not moved.

     !     Auxiliary routine called: VMOVE
     !
     !--------------------------------------------------------------------------

     INTEGER, INTENT(IN)               :: n, pos1
     INTEGER, DIMENSION(:), INTENT(IN) :: list
     INTEGER, INTENT(OUT)              :: ifault

     !     Local variables.

     INTEGER    :: next, i, l, j

     !     Check N.

     ifault = 0
     IF (n < 1 .OR. n > ncol+1-pos1) ifault = ifault + 4
     IF (ifault /= 0) RETURN

     !     Work through VORDER finding variables which are in LIST.

     next = pos1
     i = pos1
     10 l = vorder(i)
     DO j = 1, n
       IF (l == list(j)) GO TO 40
     END DO
     30 i = i + 1
     IF (i <= ncol) GO TO 10

     !     If this point is reached, one or more variables in LIST has not
     !     been found.

     ifault = 8
     RETURN

     !     Variable L is in LIST; move it up to position NEXT if it is not
     !     already there.

     40 IF (i > next) CALL vmove(i, next, ifault)
     next = next + 1
     IF (next < n+pos1) GO TO 30

     RETURN
     END SUBROUTINE reordr



     SUBROUTINE hdiag(xrow, nreq, hii, ifault)

     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
     !
     !                         -1           -1
     ! The hat matrix H = x(X'X) x' = x(R'DR) x' = z'Dz

     !              -1
     ! where z = x'R

     ! Here we only calculate the diagonal element hii corresponding to one
     ! row (xrow).   The variance of the i-th least-squares residual is (1 - hii).
     !--------------------------------------------------------------------------

     INTEGER, INTENT(IN)                  :: nreq
     INTEGER, INTENT(OUT)                 :: ifault
     double precision, DIMENSION(:), INTENT(IN)  :: xrow
     double precision, INTENT(OUT)               :: hii

     !     Local variables

     INTEGER    :: col, row, pos
     double precision  :: total, wk(ncol)

     !     Some checks

     ifault = 0
     IF (nreq > ncol) ifault = ifault + 4
     IF (ifault /= 0) RETURN

     !     The elements of xrow.inv(R).sqrt(D) are calculated and stored in WK.

     hii = zero
     DO col = 1, nreq
       IF (SQRT(d(col)) <= tol(col)) THEN
         wk(col) = zero
       ELSE
         pos = col - 1
         total = xrow(col)
         DO row = 1, col-1
           total = total - wk(row)*r(pos)
           pos = pos + ncol - row - 1
         END DO ! row = 1, col-1
         wk(col) = total
         hii = hii + total**2 / d(col)
       END IF
     END DO ! col = 1, nreq

     RETURN
     END SUBROUTINE hdiag



     FUNCTION varprd(x, nreq) RESULT(fn_val)

     !     Calculate the variance of x'b where b consists of the first nreq
     !     least-squares regression coefficients.
     !
     !--------------------------------------------------------------------------

     INTEGER, INTENT(IN)                  :: nreq
     double precision, DIMENSION(:), INTENT(IN)  :: x
     double precision                            :: fn_val

     !     Local variables

     INTEGER    :: ifault, row
     double precision  :: var, wk(nreq)

     !     Check input parameter values

     fn_val = zero
     ifault = 0
     IF (nreq < 1 .OR. nreq > ncol) ifault = ifault + 4
     IF (nobs <= nreq) ifault = ifault + 8
     IF (ifault /= 0) THEN
       RETURN
     END IF

     !     Calculate the residual variance estimate.

     var = sserr / (nobs - nreq)

     !     Variance of x'b = var.x'(inv R)(inv D)(inv R')x
     !     First call BKSUB2 to calculate (inv R')x by back-substitution.

     CALL BKSUB2(x, wk, nreq)
     DO row = 1, nreq
       IF(d(row) > tol(row)) fn_val = fn_val + wk(row)**2 / d(row)
     END DO

     fn_val = fn_val * var

     RETURN
     END FUNCTION varprd



     SUBROUTINE bksub2(x, b, nreq)

     !     Solve x = R'b for b given x, using only the first nreq rows and
     !     columns of R, and only the first nreq elements of R.
     !
     !--------------------------------------------------------------------------

     INTEGER, INTENT(IN)                  :: nreq
     double precision, DIMENSION(:), INTENT(IN)  :: x
     double precision, DIMENSION(:), INTENT(OUT) :: b

     !     Local variables

     INTEGER    :: pos, row, col
     double precision  :: temp

     !     Solve by back-substitution, starting from the top.

     DO row = 1, nreq
       pos = row - 1
       temp = x(row)
       DO col = 1, row-1
         temp = temp - r(pos)*b(col)
         pos = pos + ncol - col - 1
       END DO
       b(row) = temp
     END DO

     RETURN
     END SUBROUTINE bksub2


     END MODULE lsq
!     *********************************************************
!     *********************************************************
     subroutine Ortogonaliza (X,W,n,p)
          implicit none
          integer i,n,j,p
          double precision X(n,p),W(n),Pred(n),B(p+1)
          do j=2,p
               call Predl (X,X(1,j),W,n,j-1,Pred)
               do i=1,n
                    X(i,j)=X(i,j)-Pred(i)
               end do
          end do
     end
!     *********************************************************
!     *********************************************************
     subroutine sbackMain(X,Y,W,n,npar,kbin,h,M,muhat,family,X0,M0,muhat0,n0,B)
     !!DEC$ ATTRIBUTES DLLEXPORT::sbackMain
     !!DEC$ ATTRIBUTES C, REFERENCE, ALIAS:'sbackmain_' :: sbackMain
     implicit none
     integer j,k,i,it,itmax,n,npar,kbin,n0
     double precision eps,media,s0,X(n,npar),Y(n),W(n),M(n,npar),h(npar),family,&
     linc,cte,muhat(n),devian,devold,num,den,p0,sumw,&
     etaold,alfa, X0(n0,npar),M0(n0,npar),etam,muhat0(n0),B(npar+1)
     double precision,external::slinc,diriv,weight,dev
     double precision,allocatable::Ml(:,:), Ml0(:,:), h0(:), Z2(:), Wz(:), Eta(:), &
     eta0(:), Predl(:), Y2(:)
     allocate (Ml(n,npar), Ml0(n0,npar),h0(npar),Z2(n),Wz(n),Eta(n),eta0(n0),Predl(n),Y2(n))

     linc = family
     eps = 0.01
     itmax = 10
     M = 0
     Ml = 0
     eta = 0
     B = 0
     Y2 = 0

     call Mean(Y,W,n,p0)
     Muhat = p0
     Eta = Slinc(p0,linc)
     Devian = dev(n,muhat,y,w,family)

     do it=1,itmax     !@@@@@@@@@@ Bucle de Iteraciones @@@@@@
          do i=1,n
               Z2(i)= Eta(i)+(Y(i)-muhat(i))*diriv(muhat(i),linc)
               Wz(i)= weight(w(i), muhat(i), family, linc)
          end do

          sumw = sum(wz(1:n))
          Wz = Wz/sumW

          h0 = h

          ! Parametric part
          !Predl = 0
          do i=1,n
               do j = 1,npar
                   Y2(i) = Z2(i) - M(i,j)  !Eliminamos componente no param?trica
               end do
          end do

          call RegL(X,Y2,Wz,n,npar,B,Predl)
          do i=1,n
               Y2(i) = Z2(i) - Predl(i)     !Eliminamos componente param?trica
                do j = 1, npar
                    Ml(i,j) = B(j+1)*X(i,j)
               end do
          end do
          ! Prediccion
          do i=1,n0
               do j = 1, npar
                    Ml0(i,j) = B(j+1)*X0(i,j)
               end do
          end do

          !call Mean(Y2,Wz,n,p0)
          !do i=1,n
          !     Y2(i) = Z2(i) - p0
          !end do
          !Ml = 0
          !Ml0 = 0
          !B = 0
          !B(1) = p0

          call SBack3(X,Y2,W,Wz,n,kbin,npar,h0,M,eta,X0,M0,eta0,n0)

          ! Obtain the additive predictor
          do j=1,npar
               do i=1,n
                    eta(i) = eta(i) + Ml(i,j)
               end do
               do i=1,n0
                    eta0(i) = eta0(i) + Ml0(i,j)
               end do
          end do

          eta = eta + B(1)
          eta0 = eta0 + B(1)

          call Linv(n,Eta,Muhat,linc)
          call Linv(n0,Eta0,Muhat0,linc)

          !if (family.eq.2) goto 1

          devold=devian
          devian = dev(n,muhat,y,w,family)

          if(abs((devold-devian)/(devold+0.001)).lt.eps) goto 1
     end do
     1   continue
     h=h0
     deallocate (Ml, Ml0,h0,Z2,Wz,Eta,eta0,Predl,Y2)
     end
!     *********************************************************
!     *********************************************************
     subroutine SBack3(X,Y,W,Wy,n,kbin,npar,h,M,Pred,X0,M0,eta0,n0)
     implicit none
     logical converge,calib
     integer j,k,i,n,npar,icont,irep,ih,nh,kbin,N0
     double precision X(n,npar),Y(n),W(n),Wy(n),h(npar),M(n,npar),NW(n,npar),&
     Pred(n),h0(npar),h2(npar),xmin,xmax,pb(kbin),rango,&
     erropt,haux,hmin(npar),hmax(npar),h3(npar),hopt(npar),&
     Xgrid(kbin,npar),ErrCV,PjGrid(kbin,npar), PredCV(n), &
     PjkGrid(kbin,kbin,npar,npar),X0(n0,npar),M0(n0,npar),eta0(n0)
     !logical isnan
     h2=h
     calib=.false.
     do i=1,npar
          if (h(i).lt.0) calib = .true.
     end do

     do j=1,npar
          call min_y_max(X(1,j),n,xmin,xmax,W)
          rango=xmax-xmin
          hmin(j)=rango/30
          hmax(j)=rango/4
     end do

     nh=10

     if (calib) then
        ! VENTANAS INICALES
          do j=1,npar
               ErrOpt=9e9
               if (h(j).lt.0) then
                    do ih = 1, nh
                         haux = hmin(j)+((hmax(j)-hmin(j))/nh)*(ih-1)
                         h3 = 0
                         h3(j) = haux
                         call SBack2(X,Y,W,Wy,n,npar,h3,M,Pred,PredCV,Xgrid,&
                                     kbin,pjgrid,pjkgrid,X0,M0,eta0,n0)

                         ErrCV=0
                         do i=1,n
                              ErrCV = ErrCV + Wy(i)*(predCv(i)-Y(i))**2
                         end do

                         if (ErrCV.le.ErrOPt) then
                              ErrOPt = ErrCV
                              hopt(j) = haux
                         end if
                    end do
               else
                    hopt(j) = h(j)
               end if
          end do

          do j=1,npar
               if (h(j).lt.0) then
                    !call SBack2(X,Y,W,Wy,n,npar,hopt,M,Pred,PredCV,Xgrid,&
                    !            kbin,pjgrid,pjkgrid,X0,M0,eta0,n0)
                    h3 = hopt
                    ErrOpt = 9e9
                    do ih=1,nh
                         haux = hmin(j)+((hmax(j)-hmin(j))/nh)*(ih-1)
                         h3(j) = haux
                         call SBack2(X,Y,W,Wy,n,npar,h3,M,Pred,PredCV,Xgrid,&
                              kbin,pjgrid,pjkgrid,X0,M0,eta0,n0)

                         ErrCV = 0
                         do i = 1,n
                              ErrCV = ErrCV+Wy(i)*(predCv(i)-Y(i))**2
                         end do

                         if (ErrCV.le.ErrOPt) then !.and.(.not.isnan(ErrCV))) then
                              ErrOPt = ErrCV
                              hopt(j) = haux
                         end if
                    end do
               end if
          end do
          h2 = hopt
     end if

     call SBack2(X,Y,W,Wy,n,npar,h2,M,Pred,PredCV,Xgrid,&
     kbin,pjgrid,pjkgrid,X0,M0,eta0,n0)
     h= h2

     end
!     *********************************************************
!     *********************************************************
     subroutine SBack2(X,Y,W,Wy,n,npar,h,M,Pred,PredCV,Xgrid,kbin,&
     pjgrid,pjkgrid,X0,M0,pred0,n0)
     implicit none
     logical converge
     integer maxit,j,k,i,it,i1,i2,n,npar,kbin,N0
     double precision Err(npar),eps,integral,num,den,med,erraux,Maux,&
     X(n,npar),Y(n),W(n),Wy(n),h(npar),M(n,npar),&
     Pred(n),xmin,xmax,delta,medy,PredCV(n),&
     mauxcv, Xgrid(kbin,npar), &
     X0(n0,npar),M0(n0,npar),pred0(n0), Var,&
     PjGrid(kbin,npar),PjkGrid(kbin,kbin,npar,npar), a
     double precision,allocatable::Mold(:,:),NW(:,:),Mcv(:,:),Mgrid(:,:),NWgrid(:,:)
     double precision,allocatable::Vaux(:),Wb(:),NWgridCV(:,:),MgridCV(:,:)
     double precision,external::pjk,pj,NadWat,integrate
     maxit=10
     eps=0.01

     allocate (Mold(n,npar),NW(n,npar),Mcv(n,npar),Mgrid(kbin,npar),NWgrid(kbin,npar),&
     Vaux(kbin),Wb(kbin),NWgridCV(kbin,npar),MgridCV(kbin,npar))

     !   ###     GRID DONDE HACER LAS ESTIMACIONES
          do j=1,npar
               call Min_y_max(X(1,j),n,xmin,xmax,Wy)
               Delta=(xmax-xmin)/(kbin-1)
               do i=1,kbin
                    XGrid(i,j)=xmin+(i-1)*Delta
               end do
          end do
     !     ########     DENSIDADES     ##########
          do j=1,npar
               if (h(j).gt.0) then
                    do i=1,kbin
                         PjGrid(i,j)= pj(X(1,j),Xgrid(i,j),Wy,h(j),n)
                    end     do
               end if
          end do

          do i1=1,kbin
               do i2=1,kbin
                    do j=1,npar
                         do k=j+1,npar
                              if (h(j).gt.0.and.h(k).gt.0) then
                                   PjkGrid(i1,i2,j,k)=pjk(X(1,j),X(1,k),&
                             Xgrid(i1,j),Xgrid(i2,k),Wy,h(j),h(k),n)
                              end if
                         end do
                    end do
               end do
          end do

          do i1=1,kbin
               do i2=1,kbin
                    do j=1,npar
                         do k=1,j-1
                              PjkGrid(i1,i2,j,k)=PjkGrid(i2,i1,k,j)
                         end do
                    end do
               end do
          end do
          !########  NADARAYA-WATSON ###############
          do j=1,npar
               if (h(j).gt.0) then
                    call R1GRID_(X(1,j),Y,n,WY,h(j),kbin,Xgrid(1,j),&
                NWGrid(1,j),NWGRIDCV(1,j),Wb)
               else
                    do i=1,kbin
                         NWgrid(i,j) = 0
                         NWGRIDCV(i,j) = 0
                    end do
               end if
          end do
          !####  Initial estimates ############
          call Mean(Y,Wy,n,a)
          Mgrid = NWGrid
          do i = 1, npar
               call Interpola (Xgrid(1,i), Mgrid(1,i),kbin,&
                         X(1,i), M(1,i), n)
               do j = 1, n
                    Mold(j,i) = M(j,i)
               end do
          end do
          !########  Iterations ###############
          do it=1,maxit
               do j=1,npar
                    if (h(j).gt.0) then
                         do i=1,kbin
                              Maux = NWGrid(i,j)
                              MauxCV = NWGridCV(i,j)
                              do k=1,npar
                                   if (h(k).gt.0.and.k.ne.j) then
                                        do i1=1,kbin
                                             if (pjgrid(i,j).gt.0) then
                                                  Vaux(i1) = Mgrid(i1,k)*&
                                            PjkGrid(i,i1,j,k)/Pjgrid(i,j)
                                             else
                                                  Vaux(i1) = Mgrid(i1,k)*&
                                            PjkGrid(i,i1,j,k)/0.00001
                                             end if
                                        end do

                                        integral = Integrate(XGrid(1,k),Vaux,kbin)
                                        Maux = Maux - Integral - a
                                        MauxCV = MauxCV - Integral - a
                                   end if
                              end do
                              Mgrid(i,j)= Maux
                              MgridCV(i,j)= MauxCV
                         end do

                         call Interpola (Xgrid(1,j),Mgrid(1,j),kbin,&
                         X(1,j),M(1,j),n)

                         call Interpola (Xgrid(1,j),MgridCV(1,j),kbin,&
                         X(1,j),MCV(1,j),n)

                         !Center the smooth effects
                         call Mean(M(1,j),W,n,Med)
                         do i = 1, n
                              M(i,j) = M(i,j) - Med
                         end do

                         if (n0.gt.0) then
                              call Interpola (Xgrid(1,j),Mgrid(1,j),kbin,X0(1,j),M0(1,j),n0)
                              ! Center the smooth effects
                              do i = 1, n0
                                   M0(i,j) = M0(i,j) - Med
                              end do
                         else
                              M0 = -1
                         end if
                    else
                         do i=1,n
                              M(i,j) = 0
                              MCV(i,j) = 0
                         end do
                    end if
               end do

               do j=1,npar  !Evaluaci?n del Error en cada Coordenada
                    Num = 0
                    Den = 0
                    do i = 1,n
                         Num = Num + (M(i,j) - Mold(i,j))**2
                         Den = Den + (Mold(i,j))**2
                    end do
                    Err(j) = Num/(Den+0.0001)
               end do

               converge=.true. !Criterio de parada
               do j = 1,npar
                    if (Err(j).gt.eps) converge=.false.
               end do

               if (converge) goto 1
               Mold = M
          end do
     !     ########  Fin Iterations   ###############
     !     ########  PREDICCI?N FINAL ###############
     1   continue

          Pred = 0
          Pred0 = 0
          PredCV = 0
          do j=1,npar
               do i=1,n
                    Pred(i) = Pred(i) + M(i,j)
                    PredCV(i)= PredCV(i) + MCV(i,j)
               end do
               do i=1,n0
                    Pred0(i) = Pred0(i) + M0(i,j)
               end do
          end do

     deallocate (Mold,NW,Mcv,Mgrid,NWgrid,&
     Vaux,Wb,NWgridCV,MgridCV)
     return
     end
!     *********************************************************
!     *********************************************************
     subroutine VcoefMain(X,Z,Y,W,n,npar,Zl,nparl,kbin,h,M,Mx,muhat,family,X0,Z0,Z0l,M0,Mx0,muhat0,n0,B)
     !!DEC$ ATTRIBUTES DLLEXPORT::vcoefMain
     !!DEC$ ATTRIBUTES C, REFERENCE, ALIAS:'vcoefmain_' :: vcoefMain

     implicit none
     integer j,k,i,it,itmax,n,npar,nparl,kbin,n0,II(npar),icont
     double precision eps,media,s0,X(n,npar),Z(n,npar),Zl(n, nparl),Y(n),W(n),M(n,npar),Mx(n,npar), &
     h(npar), family, linc, cte, muhat(n), devian, devold, num, den, p0, sumw,&
     etaold,alfa, X0(n0,npar),Z0(n0,npar),Z0l(n0, nparl),M0(n0,npar),MX0(n0,npar),muhat0(n0), &
     etam, B(nparl+1)

     double precision,allocatable:: h0(:),Z2(:),Wz(:),Eta(:),etacv(:),muhatcv(:), &
     eta0(:),Predl(:),Y2(:),Ml0(:,:),Ml(:,:)

     double precision,external::slinc,diriv,weight,dev

     allocate (h0(npar),Z2(n),Wz(n),Eta(n),etacv(n),muhatcv(n),eta0(n0),Predl(n),Y2(n),Ml0(n0,nparl),Ml(n,nparl))


     linc=family
     eps=0.01
     itmax=10
     M=0
     Mx=0
     ETA=0

     call Mean(Y,W,n,p0)
     s0=Slinc(p0,linc)
     Muhat=p0
     Eta=s0
     Devian=0

     do it=1,itmax            !@@@@@@@@@@ Bucle de Iteraciones @@@@@@
          do i=1,n
               Z2(i) = Eta(i)+(Y(i)-muhat(i))*diriv(muhat(i),linc)
               Wz(i) = weight(w(i),muhat(i),family,linc)
          end do

          sumw = sum(wz(1:n))

          Wz=Wz/sumW

          h0=h

          ! Parametric part
          Predl = 0
          do i=1,n
               do j = 1,npar
                    Y2(i) = Z2(i) - M(i,j)  !Eliminamos componente no param?trica
               end do
          end do
          B = 0
          if(nparl.eq.0) then
               call Mean(Y2,Wz,n,p0)
               do i=1,n
                    Y2(i) = Z2(i) - p0
               end do
               B(1) = p0
          else
               call RegL(Zl,Y2,Wz,n,nparl,B,Predl)
               do i=1,n
                    Y2(i) = Z2(i) - Predl(i)
                    do j = 1, nparl
                         Ml(i,j) = B(j+1)*Zl(i,j)
                    end do
               end do
               ! Prediccion
               do i=1,n0
                    do j = 1, nparl
                         Ml0(i,j) = B(j+1)*Z0l(i,j)
                    end do
               end do
          end if

          !call Mean(Z2,Wz,n,p0)
          !do i=1,n
          !     Y2(i) = Z2(i) - p0
          !end do
          !Ml = 0
          !Ml0 = 0
          !B = 0
          !B(1) = p0

          call Vcoef2(X,Z,Y2,W,Wz,n,kbin,npar,h0,M,Mx,eta,X0,Z0,M0,MX0,eta0,n0)

          ! Obtain the additive predictor
          if(nparl.gt.0) then
               do j=1,nparl
                    do i=1,n
                         eta(i) = eta(i) + Ml(i,j)
                    end do
                    do i=1,n0
                         eta0(i) = eta0(i) + Ml0(i,j)
                    end do
               end do
          end if

          eta = eta + B(1)
          eta0 = eta0 + B(1)

          call Linv(n,Eta,Muhat,linc)
          call Linv(n0,Eta0,Muhat0,linc)

          !if (family.eq.2) goto 1

          devold=devian
          devian=dev(n,muhat,y,w,family)

          if(abs((devold-devian)/(devold+0.001)).lt.eps) goto 1
     end do
     1     continue
          h=h0
          deallocate (h0,Z2,Wz,Eta,etacv,muhatcv,eta0,Predl,Y2,Ml0,Ml)
     return
     end
!     *********************************************************
!     *********************************************************
     subroutine VCoef2(X,Z,Y,W,Wy,n,kbin,npar,h,M,Mx,Pred,X0,Z0,M0,MX0,pred0,n0)
     implicit none
     logical converge
     integer maxit,j,k,i,it,i1,i2,n,npar,s,nh,ih,iopt,kbin,n0
     double precision Err(npar),eps,integral(npar),num,den,media,alpha,kaux,&
     X(n,npar),Z(n,npar),Y(n),W(n),Wy(n),h(npar),Maux,erraux,M(n,npar),&
     Pred(n),xmin,xmax,delta,aux,medy,h2(npar),Xb(kbin,npar),rango,&
     yb(kbin,npar),Wb(kbin,npar),ErrCV,errmin,haux,&
     Wb2(kbin,npar),hopt(npar),erropt,pred2,Xgrid(kbin,npar),&
     aux1(kbin,npar),aux2(kbin,npar),PZjGrid(kbin,npar),&
     PZjkGrid(kbin,kbin,npar,npar),PredCV(n),h3(npar),&
     aux1cv(kbin,npar),aux2cv(kbin,npar),Mx(n,npar),hmin(npar),hmax(npar),&
     X0(n0,npar),z0(n0,npar),m0(n0,npar),mx0(n0,npar),pred0(n0)
     logical calib !,isnan

     h2=h
     calib=.false.
     do i=1,npar
          if (h(i).lt.0) calib=.true.
     end do

     do j=1,npar
          call min_y_max(X(1,j),n,xmin,xmax,W)
          rango=xmax-xmin
          hmin(j)=rango/30
          hmax(j)=rango/4
     end do

     nh=10

     do j=1,npar
          call Min_y_max(X(1,j),n,xmin,xmax,Wy)
          Delta=(xmax-xmin)/(kbin-1)
          do i=1,kbin
               Xb(i,j)=xmin+(i-1)*Delta
          end do
          call Bin1d_(X(1,j),Y,Wy,n,Xb,Yb(1,j),Wb(1,j),kbin)
     end do

     if (calib) then
          call Mean(Y,Wy,n,medy)
          alpha=medy
          !   VENTANAS INICIALES
          do j=1,npar
               ErrOpt = 9e9
               if (h(j).lt.0) then
                    do ih=1,nh
                         haux = hmin(j)+((hmax(j)-hmin(j))/nh)*(ih-1)
                         iopt = 1
                         h3 = 0
                         h3(j) = haux
                         call VCoef3(X,Z,Y,W,Wy,n,npar,kbin,h3,M,Mx,Pred,PredCv,&
                                     Xgrid,PzjGrid,PzjkGrid,aux1,aux2,alpha,iopt,X0,Z0,M0,MX0,pred0,n0)

                         ErrCV = 0
                         do i=1,n
                              ErrCV = ErrCV + Wy(i)*(predCv(i)-Y(i))**2
                         end do

                         ! MIRAR SI LA PREDICCION ES CORRECTA

                         if (ErrCV.le.ErrOPt) then !.and.(.not.isnan(ErrCV))) then
                              ErrOPt = ErrCV
                              hopt(j) = haux
                         end if
                    end do
               else
                    hopt(j)=h(j)
               end if
          end do


          do j = 1,npar
               if (h(j).lt.0) then
                    iopt = 1
                    !call VCoef3(X,Z,Y,W,Wy,n,npar,kbin,hopt,M,Mx,Pred,PredCv,&
                    !            Xgrid,PzjGrid,PzjkGrid,aux1,aux2,alpha,iopt,X0,Z0,M0,MX0,pred0,n0)
                    h3 = hopt
                    ErrOpt = 9e9

                    do ih = 1,nh
                         haux = hmin(j)+((hmax(j)-hmin(j))/nh)*(ih-1)
                         h3(j) = haux
                         iopt = 1
                         call VCoef3(X,Z,Y,W,Wy,n,npar,kbin,h3,M,Mx,Pred,PredCv,&
                                     Xgrid,PzjGrid,PzjkGrid,aux1,aux2,alpha,iopt,X0,Z0,M0,MX0,pred0,n0)
                         ErrCV = 0
                         do i = 1,n
                              ErrCV = ErrCV+Wy(i)*(predCv(i)-Y(i))**2
                         end do

                         if (ErrCV.le.ErrOPt) then !.and.(.not.isnan(ErrCV))) then
                              ErrOPt=ErrCV
                              hopt(j)=haux
                         end if
                    end do
               end if
          end do
         h2=hopt
     end if
     iopt = 1
     call VCoef3(X,Z,Y,W,Wy,n,npar,kbin,h2,M,Mx,Pred,PredCv,&
     Xgrid,PzjGrid,PzjkGrid,aux1,aux2,alpha,iopt,X0,Z0,M0,MX0,pred0,n0)
     h = h2
     end
!     *********************************************************
!     *********************************************************
     subroutine VCoef3(X,Z,Y,W,Wy,n,npar,kbin,h,M,Mx,Pred,PredCV,&
     Xgrid,PzjGrid,PzjkGrid,aux1,aux2,alpha,iopt,X0,Z0,M0,MX0,pred0,n0)
     implicit none
     logical converge
     integer maxit,j,k,i,it,i1,i2,n,npar,kbin,s,iopt,n0
     double precision Err(npar),eps,integral(npar),num,den,media,alpha,kaux,&
     X(n,npar),Z(n,npar),Y(n),W(n),Wy(n),h(npar),Maux,erraux,M(n,npar),&
     Pred(n),xmin,xmax,delta,B(npar),b0,Predl(n),aux,medy,&
     Wb(kbin,npar),aux1(kbin,npar),aux2(kbin,npar),VAux(kbin),Mold(n,npar),&
     Mgrtid(kbin,npar),aux1cv(kbin,npar),aux2cv(kbin,npar),Xgrid(kbin,npar),&
     mauxcv,PZjGrid(kbin,npar),PZjkGrid(kbin,kbin,npar,npar),Mgrid(kbin,npar),&
     MgridCV(kbin,npar),PredCV(n),MCV(n,npar),MX(n,npar),X0(n0,npar),z0(n0,npar),m0(n0,npar),&
     mx0(n0,npar),pred0(n0),med
     double precision,external::pZjk,pZj,cuantil,kernh,integrate

     maxit = 10
     eps = 0.01

     if (iopt.eq.1) then
     !     ########       NODOS    ###############
          do j=1,npar
               call Min_y_max(X(1,j),n,xmin,xmax,Wy)
               Delta=(xmax-xmin)/(kbin-1)
               do i=1,kbin
                    XGrid(i,j)=xmin+(i-1)*Delta
               end do
          end do
     !     ########       DENSIDADES    ###############
          do j=1,npar
               do i=1,kbin
                    PZjGrid(i,j)= pZj(X(1,j),Z(1,j),Xgrid(i,j),Wy,h(j),n)
               end do
          end do

          do i1=1,kbin
               do i2=1,kbin
                    do j=1,npar
                         do k=j+1,npar
                              PZjkGrid(i1,i2,j,k)=pZjk(X(1,j),X(1,k),Z(1,j),&
                        Z(1,k),Xgrid(i1,j),Xgrid(i2,k),Wy,h(j),h(k),n)
                         end do
                    end do
               end do
          end do

          do i1 = 1,kbin
               do i2 = 1,kbin
                    do j = 1,npar
                         do k = 1,j-1
                              PZjkGrid(i1,i2,j,k)=PZjkGrid(i2,i1,k,j)
                         end do
                    end do
               end do
          end do
     end if

     do j = 1,npar
          call R1GRIDZ_(X(1,j),Y,Z(1,j),n,Wy,h(j),kbin,Xgrid(1,j),&
                           Aux1(1,j),Aux2(1,j),Aux1Cv(1,j),Aux2Cv(1,j),Wb(1,j))
     end do

     !####  Initial estimates ############
     call Mean(Y, Wy, n, alpha)
     do i = 1, npar
          do j = 1, kbin
               Mgrid(j,i) = Aux1(j,i) - alpha*Aux2(j,i)
          end do
     end do

     do i = 1, npar
          call Interpola (Xgrid(1,i), Mgrid(1,i),kbin,&
                    X(1,i), M(1,i), n)
          do j = 1, n
               Mold(j,i) = M(j,i)
          end do
     end do
     !     ########  Iteraciones ###############
     do it=1,maxit
          do k=1,npar
               if (h(k).gt.0) then
                    do i1=1,kbin
                         Vaux(i1) = Mgrid(i1,k)*aux2(i1,k)*PzjGrid(i1,k)
                    end do
                    Integral(k) = Integrate(Xgrid(1,k),Vaux,kbin)
                    alpha = alpha - Integral(k)/n
               end if
          end do
     !     ########  CICLOS DEL ALGORITMO ###############
          do j=1,npar     ! Calculamos efectos
               if (h(j).gt.0) then
                    do i=1,kbin
                         Maux = Aux1(i,j) - alpha*Aux2(i,j)
                         MauxCV = Aux1CV(i,j) - alpha*Aux2CV(i,j)
                         do k=1,npar
                              if (h(k).gt.0.and.k.ne.j) then
                                   do i1=1,kbin
                                        Vaux(i1) = Mgrid(i1,k)* &
                                  PZjkGrid(i,i1,j,k)/PZjgrid(i,j)
                                   end do
                                   Integral(k) = Integrate(XGrid(1,k),Vaux,kbin)
                                   Maux = Maux - Integral(k)
                                   MauxCV = MauxCV - Integral(k)
                              end if
                         end do
                         Mgrid(i,j) = Maux
                         MgridCV(i,j) = MauxCV
                    end do
                    call Interpola (Xgrid(1,j),Mgrid(1,j),kbin,X(1,j),Mx(1,j),n)

                    ! Center the smooth effects
                    call Mean(Mx(1,j),W,n,Med)
                    do i = 1, n
                         Mx(i,j) = Mx(i,j) - Med
                    end do

                    do i=1,n
                         M(i,j) = MX(i,j)*Z(i,j)
                    end do
                    call Interpola (Xgrid(1,j),MgridCV(1,j),kbin,X(1,j),MCV(1,j),n)

                    do i=1,n
                         MCV(i,j) = MCV(i,j)*Z(i,j)
                    end do
                    if (n0.gt.0) then
                         call Interpola (Xgrid(1,j),Mgrid(1,j),kbin,X0(1,j),Mx0(1,j),n0)
                         ! Center the smooth effects
                         do i = 1, n0
                              Mx0(i,j) = Mx0(i,j) - Med
                         end do
                         do i=1,n0
                              M0(i,j) = MX0(i,j)*Z0(i,j)
                         end do
                    else
                         Mx0 = -1
                         M0 = -1
                    end if
               else
                    do i=1,n
                         M(i,j) = 0
                         MCV(i,j) = 0
                    end do
                    do i=1,n0
                         M0(i,j) = 0
                    end do
               end if
          end do
          do j=1,npar  !Evaluaci?n del Error en cada Coordenada
               if (h(j).gt.0) then
                    Num=0
                    Den=0
                    do i=1,n
                         Num = Num + (M(i,j) - Mold(i,j))**2
                         Den = Den + (Mold(i,j))**2
                    end do
                    Err(j)=Num/(Den+0.0001)
               else
                    Err(j)=0
               end if
          end do

          converge=.true. !Criterio de parada
          do j=1,npar
               if (Err(j).gt.eps) converge=.false.
          end do

          if (converge) goto 1
          Mold=M
     end do
     1     continue
          !Pred = alpha
          !PredCV = alpha
          !Pred0 = alpha
          Pred = 0
          PredCV = 0
          Pred0 = 0
          do j=1,npar
               do i=1,n
                    Pred(i) = Pred(i) + M(i,j)
                    PredCV(i) = PredCV(i) + MCV(i,j)
               end do
               do i=1,n0
                    Pred0(i) = Pred0(i) + M0(i,j)
               end do
          end do
     end
!     *********************************************************
!     *********************************************************
     subroutine R1GRIDZ_(X,Y,Z,n,W,h,kbin,Xb,M1b,M2b,M1bCV,M2bCV,Wb)
          implicit none
          integer n,kbin,i,i0,iveces,icont
          double precision X(n),Y(n),W(n),Xb(kbin),Wb(kbin),Yb(kbin),Delta,d1,k,&
          kernel1(kbin),S0(kbin),t0(kbin),t1(kbin),Sumw,ErrCV,&
          SumCv,Z(n),s0cv,t0cv,t1cv,M1bcv(kbin),M2bcv(kbin),&
          kernel2(kbin),Ymed(kbin),M1b(kbin),aux,h,haux,num,den,&
          Zb(kbin),Yzb(kbin),Z2b(kbin),M2b(kbin)
          double precision,external::L1,ker
          !logical isnan
          if (h.gt.0) then
               sumw=0
               Sumcv=0
               do i=1,n
                    sumw=sumw+W(i)
               end do
               W=W/sumw

               S0=0
               t0=0
               t1=0

               call Bin1dZ_(X,Y,Z,W,n,Xb,YZb,Zb,Z2b,Wb,kbin)

               Delta=Xb(2)-Xb(1)
               call ker1D_(h,Delta,kbin,kernel1)
               do i=1,kbin
                    kernel2=kernel1
                    haux=h
                    iveces=0
          1          continue
                    icont=0
                    do i0=1,kbin
                         d1=(i0-i)*Delta
                         k=kernel2(abs(i0-i)+1)
                         if (k.gt.0.and.Wb(i0).gt.0) icont=icont+1
                         if (k.gt.0) then
                              if (Wb(i0).gt.0) then
                                   S0(i)=S0(i)+L1(d1,0,k)*Z2b(i0)
                                   t0(i)=t0(i)+L1(d1,0,k)*YZb(i0)
                                   t1(i)=t1(i)+L1(d1,0,k)*Zb(i0)
                              end if
                         end if
                    end do
                    M1b(i)=t0(i)/s0(i)
                    M2b(i)=t1(i)/s0(i)
                s0cv=s0(i)-kernel1(1)*Z2b(i)
                    t0cv=t0(i)-kernel1(1)*YZb(i)
                    t1cv=t1(i)-kernel1(1)*Zb(i)
                    M1bCV(i)=t0cv/s0cv
                    M2bCV(i)=t1cv/s0cv
               end do
               return
          end if

     13     continue

          M1b=-1.0
          M2b=-1.0
          M1bCV=-1.0
          M2bCV=-1.0
     end
!     *********************************************************
!     *********************************************************
     subroutine Bin1dZ_(X,Y,Z,W,n,Xb,YZb,Zb,Z2b,Wb,kbin)
          implicit none
          integer KBIN,n,i,ii,igrid(n),j
          double precision x(n),y(n),Z(n),W(n),Xb(kbin),YZb(kbin),Zb(kbin),Z2b(kbin),&
          Wb(kbin),d1,delta,Area(2)

          Wb=0
          Yzb=0
          Zb=0
          Z2b=0
          delta=Xb(2)-Xb(1)


          do i=1,n
               if (X(i).le.Xb(1)) then
                    ii=1
                    Wb(ii)=Wb(ii)+W(i)
                    YZb(ii)=YZb(ii)+Y(i)*Z(i)*W(i)
                    Zb(ii)=Zb(ii)+Z(i)*W(i)
                    Z2b(ii)=Z2b(ii)+(Z(i)**2)*W(i)

               elseif (X(i).ge.Xb(kbin)) then
                    ii=kbin
                    Wb(ii)=Wb(ii)+W(i)
                    YZb(ii)=YZb(ii)+Y(i)*Z(i)*W(i)
                    Zb(ii)=Zb(ii)+Z(i)*W(i)
                    Z2b(ii)=Z2b(ii)+(Z(i)**2)*W(i)
               else
                    do j=1,kbin-1
                         if (Xb(j).le.X(i).and.X(i).le.Xb(j+1)) then
                              ii=j
                              goto 1
                         end if
                    end do
     1               continue


                    d1=Xb(ii+1)-X(i)
                    Area(1)=d1/delta
                    Area(2)=(delta-d1)/delta

                    Wb(ii)=Wb(ii)+W(i)*Area(1)
                    Wb(ii+1)=Wb(ii+1)+W(i)*Area(2)

                    YZb(ii)=YZb(ii)+Y(i)*Z(i)*W(i)*Area(1)
                    YZb(ii+1)=YZb(ii+1)+Y(i)*Z(i)*W(i)*Area(2)

                    Zb(ii)=Zb(ii)+Z(i)*W(i)*Area(1)
                    Zb(ii+1)=Zb(ii+1)+Z(i)*W(i)*Area(2)

                    Z2b(ii)=Z2b(ii)+(Z(i)**2)*W(i)*Area(1)
                    Z2b(ii+1)=Z2b(ii+1)+(Z(i)**2)*W(i)*Area(2)
               end if
          end do
     end
!     *********************************************************
!     *********************************************************
     subroutine R1GRID_(X,Y,n,WY,h,kbin,Xb,M0Grid,M0CV,Wb)
          implicit none
          integer n,kbin,i,i0
          double precision X(n),Y(n),WY(n),Xb(kbin),Wb(kbin),Yb(kbin),Delta,d1,k,&
          kernel(kbin),S0(kbin),t0(kbin),m0CV(kbin),Ymed(kbin),M0grid(kbin),aux,&
          t0cv,h,num,den
          double precision,external::L1,ker
          !logical isnan
          S0=0
          t0=0
          call Bin1d_(X,Y,Wy,n,Xb,Yb,Wb,kbin)

          Delta=Xb(2)-Xb(1)
          call ker1D_(h,Delta,kbin,kernel)
          do i=1,kbin
               do i0=1,kbin
                    d1=(i0-i)*Delta
                    k=kernel(abs(i0-i)+1)
                    if (k.gt.0) then
                         if (Wb(i0).gt.0) then
                              S0(i)=S0(i)+L1(d1,0,k)*Wb(i0)
                              t0(i)=t0(i)+L1(d1,0,k)*Yb(i0)
                         end if
                    end if
               end do

               M0grid(i) = t0(i)/s0(i)
               if (Wb(i).gt.0) M0Cv(i) = (t0(i)-kernel(1)*Yb(i))/(s0(i)-kernel(1)*Wb(i))
          end do
          return

     13   M0grid = -1
          M0Cv=-1

     end
!     *********************************************************
!     *********************************************************
     subroutine Bin1d_(X,Y,Wy,n,Xb,Yb,Wb,kbin)
          implicit none
          integer KBIN,n,i,ii,j
          double precision x(n),y(n),Wy(n),Xb(kbin),Yb(kbin),Wb(kbin),d1,delta,Area(2)

          Wb=0
          Yb=0
          delta=Xb(2)-Xb(1)
          do i=1,n
               if (X(i).le.Xb(1)) then
                    ii=1
                    Wb(ii)=Wb(ii)+Wy(i)
                    Yb(ii)=Yb(ii)+Wy(i)*Y(i)
               elseif (X(i).ge.Xb(kbin)) then
                    ii=kbin
                    Wb(ii)=Wb(ii)+Wy(i)
                    Yb(ii)=Yb(ii)+Wy(i)*Y(i)
               else

                    do j=1,kbin-1
                         if (Xb(j).le.X(i).and.X(i).le.Xb(j+1)) then
                              ii=j
                              goto 1
                         end if
                    end do
     1               continue
                    d1=Xb(ii+1)-X(i)
                    Area(1)=d1/delta
                    Area(2)=(delta-d1)/delta
                    Wb(ii)=Wb(ii)+Wy(i)*Area(1)
                    Wb(ii+1)=Wb(ii+1)+Wy(i)*Area(2)
                    Yb(ii)=Yb(ii)+Y(i)*Wy(i)*Area(1)
                    Yb(ii+1)=Yb(ii+1)+Y(i)*Wy(i)*Area(2)
               end if
          end do
     end
!     *********************************************************
!     *********************************************************
     subroutine WRegresion(X,Y,W,n,nvar,beta,sterr,se,r2,iopt)
          USE lsq
          IMPLICIT NONE
          INTEGER             :: i, ier, iostatus, j, m, n,nvar,iopt
          double precision          :: x(n,nvar), y(n),W(n), xrow(0:nvar+1),&
           wt = 1.0_dp, beta(0:nvar+1),var, covmat(231), sterr(0:nvar+1), &
           totalSS,se,r2
          LOGICAL             :: fit_const = .TRUE., lindep(0:20), xfirst



          ! Least-squares calculations
          m=nvar
          CALL startup(m, fit_const)
          DO i = 1, n
            xrow(0) = 1.0_dp
            DO j = 1, m
               xrow(j) = x(i,j)
            END DO
            CALL includ(W(i), xrow, y(i))
          END DO


          if (iopt.gt.0) then
               CALL sing(lindep, ier)
          end if

          ! Calculate progressive residual sums of squares
          CALL ss()
          var = rss(m+1) / (n - m - 1)

          ! Calculate least-squares regn. coeffs.
          CALL regcf(beta, m+1, ier)

          if (iopt.gt.0) then
          ! Calculate covariance matrix, and hence std. errors of coeffs.
          CALL cov(m+1, var, covmat, 231, sterr, ier)

          se=SQRT(var)
          totalSS = rss(1)
          !WRITE(*, '(a, g20.12)') ' R^2 = ', (totalSS - rss(m+1))/totalSS

          r2=(totalSS - rss(m+1))/totalSS

          end if
     END
!     *********************************************************
!     *********************************************************
     subroutine Predl (X,Y,W,n,p,Pred)
          implicit none
          integer i,n,p,iopt,j
          double precision X(n,p),Y(n),W(n),Pred(n),beta(p+1),&
          sterr(p+1),se,r2
          iopt=0
          call WRegresion(X,Y,W,n,p,beta,sterr,se,r2,iopt)

          Pred=beta(1)
          do i=1,n
               do j=1,p
                    Pred(i)=Pred(i)+Beta(j+1)*X(i,j)
               end do
          end do
     end
!     *********************************************************
!     *********************************************************
     subroutine ker1D_(h,Delta,kbin,ker)
          implicit none
          integer kbin,i
          double precision Delta,ker(kbin),Dis,h
          ker=0
          do i=1,kbin
               Dis=(i-1)*Delta/h
               Dis=-.5*Dis**2
               if (Dis.gt.-5.and.h.gt.0) ker(i)=exp(dis)/(h*sqrt(2*3.1415))
          end do
     end
!     *********************************************************
!     *********************************************************
     double precision function L1(d1,r,Ker)
          implicit none
          integer r
          double precision d1,ker
          L1=(d1**r)*ker
     end
!     *********************************************************
!     *********************************************************
     double precision function kernh(X,x0,h)
          implicit none
          double precision X,x0,dis,h,aux
          dis=(x-x0)/h
          kernh=0
          aux=-0.5*(Dis**2)
          if (aux.gt.-5.and.h.gt.0) kernh=exp(aux)/(h*sqrt(2*3.1415))
     end
!     *********************************************************
!     *********************************************************
     double precision function pZj(X,Z,x0,W,h,n)
          implicit none
          integer n,i
          double precision X(n),Z(n),W(n),x0,h
          double precision,external::kernh
          pZj=0
          if (h.gt.0) then
               do i=1,n
                    pZj=pZj+(Z(i)**2)*kernh(X(i),x0,h)*W(i)
               end do
          end if
     end
!     *********************************************************
!     *********************************************************
     double precision function pZjk(X1,X2,Z1,Z2,x01,x02,W,h1,h2,n)
          implicit none
          integer n,i
          double precision X1(n),X2(n),W(n),Z1(n),Z2(n),x01,x02,h1,h2
          double precision,external::kernh
          pZjk=0

          if (h1.le.0.or.h2.le.0) then
               pZjk=0
          else
               do i=1,n
                    pZjk=pZjk+Z1(i)*Z2(i)*kernh(X1(i),x01,h1)*&
                kernh(X2(i),x02,h2)*W(i)
               end do
          end if
     end
!     *********************************************************
!     *********************************************************
     double precision function pj(X,x0,W,h,n)
          implicit none
          integer n,i
          double precision X(n),W(n),x0,h
          double precision,external::kernh
          pj=0
          if (h.gt.0) then
               do i=1,n
                    pj=pj+kernh(X(i),x0,h)*W(i)
               end do
          end if
     end
!     *********************************************************
!     *********************************************************
     double precision function pjk(X1,X2,x01,x02,W,h1,h2,n)
          implicit none
          integer n,i
          double precision X1(n),X2(n),W(n),x01,x02,h1,h2
          double precision,external::kernh
          pjk=0

          if (h1.le.0.or.h2.le.0) then
               pjk=0
          else
               do i=1,n
                    pjk=pjk+kernh(X1(i),x01,h1)*kernh(X2(i),x02,h2)*W(i)
               end do
          end if
     end
!     *********************************************************
!     *********************************************************
     subroutine min_y_max(x,n,min,max,W)
          double precision x(n),W(n),min,max
          do i=1,n
               if (W(i).gt.0) then
                    min=x(i)
                    max=x(i)
                    goto 1
               end if
          end do
          1 continue
          do i=1,n
               if (W(i).gt.0) then
                    if (x(i).lt.min) min=x(i)
                    if (x(i).gt.max) max=x(i)
                end if
          end do
     end
!     *********************************************************
!     *********************************************************
     subroutine Mean(vector,w,n,media)
          integer n
          double precision vector(n),media,w(n),wtot
          media=0
          wtot=0
          do i=1,n
               media=media+w(i)*vector(i)
               wtot=wtot+w(i)
          end do
          media=media/wtot
     end
!     *********************************************************
!     *********************************************************
     double precision function Integrate(X,FX,n)
          implicit none
          integer n, it
          double precision FX(n), X(n)
          Integrate=0
          do it=1,(n-3)/2
               Integrate = Integrate+ 2*FX(2*it+1)
          end do
          do it =1,(n-1)/2
               Integrate = Integrate + 4*FX(2*it)
          end do
          Integrate = (((X(n)-X(1))/(n-1))/3)*(Integrate + FX(1) + FX(n))
     end
!     *********************************************************
!     *********************************************************
     subroutine Interpola(Xgrid,Pgrid,kbin,X0,P0,n)

     ! Fit a quintic spline with user control of knot positions.
     ! If the knots are at tk1, tk2, ..., then the fitted spline is
     ! b0 + b1.t + b2.t^2 + b3.t^3 + b4.t^4 + b5.t^5    for t <= tk1
     ! b0 + ... + b5.t^5 + b6.(t-tk1)^5                 for tk1 < t <= tk2
     ! b0 + ... + b5.t^5 + b6.(t-tk1)^5 + b7.(t-tk2)^5  for tk2 < t <= tk3
     ! b0 + ... + b5.t^5 + b6.(t-tk1)^5 + b7.(t-tk2)^5 + b8.(t-tk3)^5
     !                                                  for tk3 < t <= tk4, etc.

     ! In this version, the knots are evenly spaced.
     ! Also calculates first & 2nd derivatives of the spline.

     ! Uses the author's least-squares package in file lsq.f90
     ! Latest revision - 2 November 2003
     ! Alan Miller (amiller @ bigpond.net.au)

          USE lsq
          IMPLICIT NONE

          INTEGER                 :: i, ier, iostatus, j, n, nk,next_knot, pos,kbin,icont
          double precision               :: t, t1, y, dist,&
          Xgrid(kbin),Pgrid(kbin),X0(n),P0(n),P1(n),P2(n)
          double precision, PARAMETER    :: one = 1.0_dp,cero=0.0_dp
          double precision, ALLOCATABLE  :: knot(:), xrow(:), b(:)


          icont=0
          do i=1,kbin
               if (pgrid(i).ne.-1.0) icont=icont+1
          end do

          if (icont.gt.5) then
          nk=icont/5


          ALLOCATE (knot(nk),xrow(0:5+nk), b(0:5+nk))

          ! Calculate knot positions, evenly spaced.
          dist = (Xgrid(kbin) - Xgrid(1)) / (nk + 1)
          t1=Xgrid(1)
          DO i = 1, nk
            knot(i) = t1 + dist * i
          END DO

          next_knot = 1

          ! Initialize the least-squares calculations
          CALL startup(6+nk, .FALSE.)

          DO i=1,kbin
             t=Xgrid(i)
             y=Pgrid(i)
             xrow(0) = one
               xrow(1) = (t - t1)
               xrow(2) = (t - t1) * xrow(1)
               xrow(3) = (t - t1) * xrow(2)
               xrow(4) = (t - t1) * xrow(3)
               xrow(5) = (t - t1) * xrow(4)
               IF (t > knot(next_knot)) next_knot = MIN(nk, next_knot + 1)
               DO j = 1, next_knot-1
                    xrow(5+j) = (t - knot(j))**5
               END DO
               xrow(5+next_knot:5+nk) = 0.0_dp
               if (y.ne.-1.0_dp) CALL includ(one, xrow, y)

          END DO

          CALL regcf(b, 6+nk, ier)

          next_knot = 1
          DO i = 1, n
             next_knot = 1
             t=X0(i)
             xrow(0) = one
             xrow(1) = (t - t1)
             xrow(2) = (t - t1) * xrow(1)
             xrow(3) = (t - t1) * xrow(2)
             xrow(4) = (t - t1) * xrow(3)
             xrow(5) = (t - t1) * xrow(4)
             if (i.eq.45) then
               continue
             end if
          55 continue
             IF (t > knot(next_knot)) THEN
               next_knot = next_knot + 1
               IF (next_knot <= nk) THEN
                  goto 55
               ELSE
                 next_knot = nk + 1
                 goto 56
               END IF
            END IF

          56 continue
            DO j = 1, next_knot-1
               xrow(5+j) = (t - knot(j))**5
            END DO
            p0(i) = DOT_PRODUCT( b(0:5+next_knot-1), xrow(0:5+next_knot-1) )
            p2(i) = ((20*b(5)*(t-t1) + 12*b(4))*(t-t1) + 6*b(3))*(t-t1) + 2*b(2)
            p1(i) = (((5*b(5)*(t-t1) + 4*b(4))*(t-t1) + 3*b(3))*(t-t1) + 2*b(2))*(t-t1) + b(1)
            DO j = 1, next_knot-1
               p1(i) = p1(i) + 5*b(j+5)*(t - knot(j))**4
               p2(i) = p2(i) + 20*b(j+5)*(t - knot(j))**3
            END DO
          END DO
          deallocate ( knot,xrow, b )
          else
          p0=-1
          p1=-1
          p2=-1
          end if
     end
!     *********************************************************
!     *********************************************************
     double precision function DEV(n,fits,y,w,family)
          implicit none
          integer n
          double precision fits(n),y(n),w(n),family
          double precision,external::DEVGAM,DEVPOI,DEVB,DEVG
          if (family.eq.1) then          !BINOMIAL
               dev=devb(n,fits,y,w)
          elseif(family.eq.2) then     !NORMAL
               dev=devg(n,fits,y,w)
          elseif(family.eq.3) then     !POISSON
               dev=devpoi(n,fits,y,w)
          end if
     end
!     *********************************************************
!     *********************************************************
     double precision function DEVG(n,fits,y,w)
          implicit none
          integer i,n
          double precision fits(n),y(n),w(n),rss
          rss=0.0
          do i=1,n
               rss=rss+w(i)*(y(i)-fits(i))*(y(i)-fits(i))
          end do
          devg=rss
     end
!     *********************************************************
!     *********************************************************
     double precision function DEVB(n,fits,y,w)
          implicit none
          integer i,n
          double precision fits(n),y(n),w(n)
          double precision pr,entrop,entadd,dev
          dev=0.0
          do i=1,n
               pr=fits(i)
               if(pr.lt.0.0001) pr=0.0001
               if(pr.gt.0.9999) pr=0.9999
               if((1.0-y(i))*y(i).le.0.0) then
                    entrop=0.0
               else
                    entrop=2.0*w(i)*(y(i)*log(y(i))+&
                (1.0-y(i))*log(1.0-y(i)))
               end if
               entadd=2.0*w(i)*(y(i)*log(pr)+(1-y(i))*log(1.0-pr))
               dev=dev+entrop-entadd
          end do
          devb=dev
     end
!     *********************************************************
!     *********************************************************
     double precision function DEVPOI(n,fits,y,w)
          implicit none
          integer i,n
          double precision fits(n),y(n),w(n),tempf,alog
          devpoi=0.0
          do i=1,n
               tempf=fits(i)
               if(tempf.lt.0.0001) tempf=.0001
               devpoi=devpoi+2*w(i)*(-y(i)*log(tempf)- (y(i)-fits(i)))
               if(y(i).gt.0.0) devpoi=devpoi+2*w(i)*y(i)*log(y(i))
          end do
     end
!     *********************************************************
!     *********************************************************
     double precision function WEIGHT(w,muhat,family,linc)
          implicit none
          double precision w,muhat,family,linc,temp1,temp,maux,aux
          double precision,external::diriv

          if(family.eq.1) then          !BINOMIAL
               temp=DIRIV(muhat,linc)
               temp1=muhat*(1.0-muhat)
               aux=temp1*temp**2
               if (aux.lt.0.0001) aux=0.0001
               weight=w/aux
          elseif(family.eq.2) then     !NORMAL
                    WEIGHT=w
          elseif(family.eq.3) then     !POISSON
               if (muhat.le.0.05)then
                    weight=0
               else
                    weight=w/(muhat*DIRIV(muhat,linc)**2)
               end if
          end if
     end
!     *********************************************************
!     *********************************************************
     double precision function DIRIV(muhat,linc)
          implicit none
          double precision muhat,linc
          double precision,external::dirvlo,dirvlt,dirvin,dirvid
          if(linc.eq.1) then          !LOGIT
               DIRIV=dirvlt(muhat)
          elseif(linc.eq.2) then  !IDENTIDAD
               DIRIV=1
          elseif(linc.eq.3) then    !LOGARITMO
               DIRIV=dirvlo(muhat)
          end if
     end
!     *********************************************************
!     *********************************************************
     double precision function DIRVLT(muhat)     !Derivada Logit
          implicit none
          double precision muhat,pr
          pr=muhat
          if (pr.ge..9999) pr=.9999
          if (pr.le..0001) pr=.0001
          pr=pr*(1.0-pr)
          dirvlt=1.0/pr
     end
!     *********************************************************
!     *********************************************************
     double precision function LINCLT(muhat)     !Link Logit
          implicit none
          double precision muhat,logit,d
          logit=muhat
          d=1.0-logit
          if(d.lt.0.0001) d=0.0001
          if(d.ge..9999)  d=0.9999
          logit=logit/d
          LINCLT=log(logit)
     end
!     *********************************************************
!     *********************************************************
     double precision function DIRVLO(muhat)     !Derivada Logaritmo
          implicit none
          double precision muhat,pr
          pr=muhat
          if (pr.le.0.0001) pr=0.0001
          DIRVLO=1.0/pr
     end
!     *********************************************************
!     *********************************************************
     double precision function SLINC(muhat,linc)
          implicit none
          double precision muhat,linc
          double precision,external::lincid,linclt,linclo,lincin
          if(linc.eq.1) then !Logit
               SLINC=LINCLT(muhat)
          elseif(linc.eq.2) then  !Identidad
               SLINC=LINCID(muhat)
          elseif(linc.eq.3) then !Logar
               SLINC=LINCLO(muhat)
          end if
     end
!     *********************************************************
!     *********************************************************
     subroutine linv(n,etahat,muhat,linc)
         implicit none
          integer n
         double precision linc,etahat(n),muhat(n)
          if(linc.eq.1) then !LOGIT
               call linvlt(n,etahat,muhat)
          elseif(linc.eq.2) then !IDENTIDAD
               call linvid(n,etahat,muhat)
          elseif(linc.eq.3) then !LOGAR
               call linvlo(n,etahat,muhat)
          end if
      end
!     *********************************************************
!     *********************************************************
     subroutine linvid(n,etahat,muhat)
          implicit none
          integer n,i
          double precision muhat(n),etahat(n)
          muhat=etahat
     end
!     *********************************************************
!     *********************************************************
     subroutine linvlt(n,etahat,muhat)
          implicit none
          integer n,i
          double precision muhat(n),etahat(n), pr
          do i=1,n
               pr=etahat(i)
               if (pr.gt.300)  pr=300
               if (pr.lt.-300) pr=-300
               pr=exp(pr)
               muhat(i)=pr/(1+pr)
          end do
     end
!     *********************************************************
!     *********************************************************
     subroutine linvlo(n,etahat,muhat)
          implicit none
          integer i,n
          double precision muhat(n),etahat(n),pr

          do i=1,n
               pr=etahat(i)
               if (pr.gt.300)  pr=300
               muhat(i)=exp(pr)
          end do
     end
!     *********************************************************
!     *********************************************************
     double precision function LINCID(muhat)  !Link identidad
          double precision muhat
          LINCID=muhat
     end
!     *********************************************************
!     *********************************************************
     double precision function LINCLO(muhat)     !Linc Logaritmo
          double precision muhat,pr
          pr=muhat
          if (pr.le.0.0001) pr=0.0001
          linclo=log(pr)
     end
!     *********************************************************
!     *********************************************************
     double precision function LINCIN(muhat) !Link inverso
          double precision muhat
          LINCIN=1.0/muhat
     end
!     *********************************************************
!     Multidimensional linear regression
!     *********************************************************
     subroutine Regl(X,Y,W,n,p,Beta,Pred)
          implicit none
          integer i,n,j,p,iopt
          double precision X(n,p),Y(n),W(n),Pred(n),beta(p+1), &
          sterr(p+1),se,r2
          iopt=0
          call WRegresion(X,Y,W,n,p,beta,sterr,se,r2,iopt)
          call PredLineal(X,n,p,Beta,Pred)
     end
!     *********************************************************
!     Prediction
!     *********************************************************
     subroutine PredLineal (X,n,p,B,Pred)
          implicit none
          integer i,n,j,p
          double precision X(n,p),B(p+1),Pred(n)

          Pred=0
          do i=1,n
               Pred(i) = B(1)
               do j = 1,p
                    Pred(i) = Pred(i)+B(j+1)*X(i,j)
               end do
     end do
     end
!     *********************************************************
!     *********************************************************
     subroutine Reglineal(X,Y,W,n,p,B,pred)
          implicit none
          integer i,n,j,p
          double precision X(n),Y(n),W(n),Pred(n),B(p+1)
          double precision, allocatable:: X2(:,:)
          allocate (X2(n,p))
          do i=1,n
               do j=1,p
                    X2(i,j)=X(i)**j
               end do
          end do
          call Regl(X2,Y,W,n,p,B,Pred)
          deallocate (X2)
     end
!     *********************************************************
!     *********************************************************
     subroutine RegLinealPred(X,Y,W,n,p,F,Xp,Yp,np)
          implicit none
          integer n,np,i,p,n2,p2,j
          double precision X(n), Y(n), W(n), Xp(np), Yp(np), F(n)
          double precision, allocatable:: Xp2(:,:), B(:), pred(:)
          allocate (B(p+1),Xp2(np,p),pred(n))

          call Reglineal (X,Y,W,n,p,B,F)

          ! Prediccion
          do i=1,np
               Yp(i)=B(1)
               do j=1,p
                    Xp2(i,j)=Xp(i)**j
                    Yp(i)=Yp(i)+B(j+1)*Xp2(i,j)
               end do
          end do
          deallocate(B,Xp2,pred)
     end subroutine
     !######################################################################################################################
     !######################################################################################################################
     !######################################################################################################################
     !######################################################################################################################
     subroutine Interpola2(Xgrid, Pgrid, kbin, X0, P0, n)
     ! Fit a quintic spline with user control of knot positions.
     ! If the knots are at tk1, tk2, ..., then the fitted spline is
     ! b0 + b1.t + b2.t^2 + b3.t^3 + b4.t^4 + b5.t^5    for t <= tk1
     ! b0 + ... + b5.t^5 + b6.(t-tk1)^5                 for tk1 < t <= tk2
     ! b0 + ... + b5.t^5 + b6.(t-tk1)^5 + b7.(t-tk2)^5  for tk2 < t <= tk3
     ! b0 + ... + b5.t^5 + b6.(t-tk1)^5 + b7.(t-tk2)^5 + b8.(t-tk3)^5
     !                                                  for tk3 < t <= tk4, etc.

     ! In this version, the knots are evenly spaced.
     ! Also calculates first & 2nd derivatives of the spline.

     ! Uses the author's least-squares package in file lsq.f90
     ! Latest revision - 2 November 2003
     ! Alan Miller (amiller @ bigpond.net.au)
     use lsq
     IMPLICIT NONE
     INTEGER :: i, ier, iostatus, j, n, nk,next_knot, pos,kbin,icont
     double precision :: t, t1, y, dist, &
     Xgrid(kbin),Pgrid(kbin),X0(n),P0(n)
     double precision, PARAMETER    :: one = 1.0_dp, cero=0.0_dp
     double precision, ALLOCATABLE  :: knot(:), xrow(:), b(:)

     icont=0
     do i=1,kbin
          if (pgrid(i).ne.-1.0) icont=icont+1
     end do

     if (icont.gt.5) then
          nk=icont/5
          ALLOCATE (knot(nk),xrow(0:5+nk), b(0:5+nk))

          !numero de nodos
          ! Calculate knot positions, evenly spaced.
          dist = (Xgrid(kbin) - Xgrid(1)) / (nk + 1)
          t1=Xgrid(1)
          DO i = 1, nk
               knot(i) = t1 + dist * i
          END DO
          next_knot = 1
          ! Initialize the least-squares calculations
          CALL startup(6+nk, .FALSE.)
          DO i=1,kbin
             t=Xgrid(i)
             y=Pgrid(i)
             xrow(0) = one
               xrow(1) = (t - t1)
               xrow(2) = (t - t1) * xrow(1)
               xrow(3) = (t - t1) * xrow(2)
               xrow(4) = (t - t1) * xrow(3)
               xrow(5) = (t - t1) * xrow(4)
               IF (t > knot(next_knot)) next_knot = MIN(nk, next_knot + 1)
               DO j = 1, next_knot-1
                    xrow(5+j) = (t - knot(j))**5
               END DO
               xrow(5+next_knot:5+nk) = 0.0_dp
               if (y.ne.-1.0_dp) CALL includ(one, xrow, y)
          END DO
          CALL regcf(b, 6+nk, ier)
          next_knot = 1
          DO i = 1, n
               next_knot = 1
               t=X0(i)
               xrow(0) = one
               xrow(1) = (t - t1)
               xrow(2) = (t - t1) * xrow(1)
               xrow(3) = (t - t1) * xrow(2)
               xrow(4) = (t - t1) * xrow(3)
               xrow(5) = (t - t1) * xrow(4)
               if (i.eq.45) then
                    continue
               end if
     55               continue
               IF (t > knot(next_knot)) THEN
                    next_knot = next_knot + 1
                    IF (next_knot <= nk) THEN
                         goto 55
                    ELSE
                         next_knot = nk + 1
                         goto 56
                    END IF
               END IF
     56               continue
               DO j = 1, next_knot-1
                    xrow(5+j) = (t - knot(j))**5
               END DO
               p0(i) = DOT_PRODUCT( b(0:5+next_knot-1), xrow(0:5+next_knot-1) )
          END DO
          deallocate (knot,xrow,b)
     else
          p0=-1
     end if
     call endup()
     end
