#include "sofa.h"
#include "sofam.h"

double iauS06(double date1, double date2, double x, double y)
/*
**  - - - - - - -
**   i a u S 0 6
**  - - - - - - -
**
**  The CIO locator s, positioning the Celestial Intermediate Origin on
**  the equator of the Celestial Intermediate Pole, given the CIP's X,Y
**  coordinates.  Compatible with IAU 2006/2000A precession-nutation.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date1,date2   double    TT as a 2-part Julian Date (Note 1)
**     x,y           double    CIP coordinates (Note 3)
**
**  Returned (function value):
**                   double    the CIO locator s in radians (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The CIO locator s is the difference between the right ascensions
**     of the same point in two systems:  the two systems are the GCRS
**     and the CIP,CIO, and the point is the ascending node of the
**     CIP equator.  The quantity s remains below 0.1 arcsecond
**     throughout 1900-2100.
**
**  3) The series used to compute s is in fact for s+XY/2, where X and Y
**     are the x and y components of the CIP unit vector;  this series
**     is more compact than a direct series for s would be.  This
**     function requires X,Y to be supplied by the caller, who is
**     responsible for providing values that are consistent with the
**     supplied date.
**
**  4) The model is consistent with the "P03" precession (Capitaine et
**     al. 2003), adopted by IAU 2006 Resolution 1, 2006, and the
**     IAU 2000A nutation (with P03 adjustments).
**
**  Called:
**     iauFal03     mean anomaly of the Moon
**     iauFalp03    mean anomaly of the Sun
**     iauFaf03     mean argument of the latitude of the Moon
**     iauFad03     mean elongation of the Moon from the Sun
**     iauFaom03    mean longitude of the Moon's ascending node
**     iauFave03    mean longitude of Venus
**     iauFae03     mean longitude of Earth
**     iauFapa03    general accumulated precession in longitude
**
**  References:
**
**     Capitaine, N., Wallace, P.T. & Chapront, J., 2003, Astron.
**     Astrophys. 432, 355
**
**     McCarthy, D.D., Petit, G. (eds.) 2004, IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG
**
**  This revision:  2021 May 11
**
**  SOFA release 2023-10-11
**
**  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
*/
{
/* Time since J2000.0, in Julian centuries */
   double t;

/* Miscellaneous */
   int i, j;
   double a, w0, w1, w2, w3, w4, w5;

/* Fundamental arguments */
   double fa[8];

/* Returned value */
   double s;

/* --------------------- */
/* The series for s+XY/2 */
/* --------------------- */

   typedef struct {
      int nfa[8];      /* coefficients of l,l',F,D,Om,LVe,LE,pA */
      double s, c;     /* sine and cosine coefficients */
   } TERM;

/* Polynomial coefficients */
   static const double sp[] = {

   /* 1-6 */
          94.00e-6,
        3808.65e-6,
        -122.68e-6,
      -72574.11e-6,
          27.98e-6,
          15.62e-6
   };

/* Terms of order t^0 */
   static const TERM s0[] = {

   /* 1-10 */
      {{ 0,  0,  0,  0,  1,  0,  0,  0}, -2640.73e-6,   0.39e-6 },
      {{ 0,  0,  0,  0,  2,  0,  0,  0},   -63.53e-6,   0.02e-6 },
      {{ 0,  0,  2, -2,  3,  0,  0,  0},   -11.75e-6,  -0.01e-6 },
      {{ 0,  0,  2, -2,  1,  0,  0,  0},   -11.21e-6,  -0.01e-6 },
      {{ 0,  0,  2, -2,  2,  0,  0,  0},     4.57e-6,   0.00e-6 },
      {{ 0,  0,  2,  0,  3,  0,  0,  0},    -2.02e-6,   0.00e-6 },
      {{ 0,  0,  2,  0,  1,  0,  0,  0},    -1.98e-6,   0.00e-6 },
      {{ 0,  0,  0,  0,  3,  0,  0,  0},     1.72e-6,   0.00e-6 },
      {{ 0,  1,  0,  0,  1,  0,  0,  0},     1.41e-6,   0.01e-6 },
      {{ 0,  1,  0,  0, -1,  0,  0,  0},     1.26e-6,   0.01e-6 },

   /* 11-20 */
      {{ 1,  0,  0,  0, -1,  0,  0,  0},     0.63e-6,   0.00e-6 },
      {{ 1,  0,  0,  0,  1,  0,  0,  0},     0.63e-6,   0.00e-6 },
      {{ 0,  1,  2, -2,  3,  0,  0,  0},    -0.46e-6,   0.00e-6 },
      {{ 0,  1,  2, -2,  1,  0,  0,  0},    -0.45e-6,   0.00e-6 },
      {{ 0,  0,  4, -4,  4,  0,  0,  0},    -0.36e-6,   0.00e-6 },
      {{ 0,  0,  1, -1,  1, -8, 12,  0},     0.24e-6,   0.12e-6 },
      {{ 0,  0,  2,  0,  0,  0,  0,  0},    -0.32e-6,   0.00e-6 },
      {{ 0,  0,  2,  0,  2,  0,  0,  0},    -0.28e-6,   0.00e-6 },
      {{ 1,  0,  2,  0,  3,  0,  0,  0},    -0.27e-6,   0.00e-6 },
      {{ 1,  0,  2,  0,  1,  0,  0,  0},    -0.26e-6,   0.00e-6 },

   /* 21-30 */
      {{ 0,  0,  2, -2,  0,  0,  0,  0},     0.21e-6,   0.00e-6 },
      {{ 0,  1, -2,  2, -3,  0,  0,  0},    -0.19e-6,   0.00e-6 },
      {{ 0,  1, -2,  2, -1,  0,  0,  0},    -0.18e-6,   0.00e-6 },
      {{ 0,  0,  0,  0,  0,  8,-13, -1},     0.10e-6,  -0.05e-6 },
      {{ 0,  0,  0,  2,  0,  0,  0,  0},    -0.15e-6,   0.00e-6 },
      {{ 2,  0, -2,  0, -1,  0,  0,  0},     0.14e-6,   0.00e-6 },
      {{ 0,  1,  2, -2,  2,  0,  0,  0},     0.14e-6,   0.00e-6 },
      {{ 1,  0,  0, -2,  1,  0,  0,  0},    -0.14e-6,   0.00e-6 },
      {{ 1,  0,  0, -2, -1,  0,  0,  0},    -0.14e-6,   0.00e-6 },
      {{ 0,  0,  4, -2,  4,  0,  0,  0},    -0.13e-6,   0.00e-6 },

   /* 31-33 */
      {{ 0,  0,  2, -2,  4,  0,  0,  0},     0.11e-6,   0.00e-6 },
      {{ 1,  0, -2,  0, -3,  0,  0,  0},    -0.11e-6,   0.00e-6 },
      {{ 1,  0, -2,  0, -1,  0,  0,  0},    -0.11e-6,   0.00e-6 }
   };

/* Terms of order t^1 */
   static const TERM s1[] = {

   /* 1 - 3 */
      {{ 0,  0,  0,  0,  2,  0,  0,  0},    -0.07e-6,   3.57e-6 },
      {{ 0,  0,  0,  0,  1,  0,  0,  0},     1.73e-6,  -0.03e-6 },
      {{ 0,  0,  2, -2,  3,  0,  0,  0},     0.00e-6,   0.48e-6 }
   };

/* Terms of order t^2 */
   static const TERM s2[] = {

   /* 1-10 */
      {{ 0,  0,  0,  0,  1,  0,  0,  0},   743.52e-6,  -0.17e-6 },
      {{ 0,  0,  2, -2,  2,  0,  0,  0},    56.91e-6,   0.06e-6 },
      {{ 0,  0,  2,  0,  2,  0,  0,  0},     9.84e-6,  -0.01e-6 },
      {{ 0,  0,  0,  0,  2,  0,  0,  0},    -8.85e-6,   0.01e-6 },
      {{ 0,  1,  0,  0,  0,  0,  0,  0},    -6.38e-6,  -0.05e-6 },
      {{ 1,  0,  0,  0,  0,  0,  0,  0},    -3.07e-6,   0.00e-6 },
      {{ 0,  1,  2, -2,  2,  0,  0,  0},     2.23e-6,   0.00e-6 },
      {{ 0,  0,  2,  0,  1,  0,  0,  0},     1.67e-6,   0.00e-6 },
      {{ 1,  0,  2,  0,  2,  0,  0,  0},     1.30e-6,   0.00e-6 },
      {{ 0,  1, -2,  2, -2,  0,  0,  0},     0.93e-6,   0.00e-6 },

   /* 11-20 */
      {{ 1,  0,  0, -2,  0,  0,  0,  0},     0.68e-6,   0.00e-6 },
      {{ 0,  0,  2, -2,  1,  0,  0,  0},    -0.55e-6,   0.00e-6 },
      {{ 1,  0, -2,  0, -2,  0,  0,  0},     0.53e-6,   0.00e-6 },
      {{ 0,  0,  0,  2,  0,  0,  0,  0},    -0.27e-6,   0.00e-6 },
      {{ 1,  0,  0,  0,  1,  0,  0,  0},    -0.27e-6,   0.00e-6 },
      {{ 1,  0, -2, -2, -2,  0,  0,  0},    -0.26e-6,   0.00e-6 },
      {{ 1,  0,  0,  0, -1,  0,  0,  0},    -0.25e-6,   0.00e-6 },
      {{ 1,  0,  2,  0,  1,  0,  0,  0},     0.22e-6,   0.00e-6 },
      {{ 2,  0,  0, -2,  0,  0,  0,  0},    -0.21e-6,   0.00e-6 },
      {{ 2,  0, -2,  0, -1,  0,  0,  0},     0.20e-6,   0.00e-6 },

   /* 21-25 */
      {{ 0,  0,  2,  2,  2,  0,  0,  0},     0.17e-6,   0.00e-6 },
      {{ 2,  0,  2,  0,  2,  0,  0,  0},     0.13e-6,   0.00e-6 },
      {{ 2,  0,  0,  0,  0,  0,  0,  0},    -0.13e-6,   0.00e-6 },
      {{ 1,  0,  2, -2,  2,  0,  0,  0},    -0.12e-6,   0.00e-6 },
      {{ 0,  0,  2,  0,  0,  0,  0,  0},    -0.11e-6,   0.00e-6 }
   };

/* Terms of order t^3 */
   static const TERM s3[] = {

   /* 1-4 */
      {{ 0,  0,  0,  0,  1,  0,  0,  0},     0.30e-6, -23.42e-6 },
      {{ 0,  0,  2, -2,  2,  0,  0,  0},    -0.03e-6,  -1.46e-6 },
      {{ 0,  0,  2,  0,  2,  0,  0,  0},    -0.01e-6,  -0.25e-6 },
      {{ 0,  0,  0,  0,  2,  0,  0,  0},     0.00e-6,   0.23e-6 }
   };

/* Terms of order t^4 */
   static const TERM s4[] = {

   /* 1-1 */
      {{ 0,  0,  0,  0,  1,  0,  0,  0},    -0.26e-6,  -0.01e-6 }
   };

/* Number of terms in the series */
   static const int NS0 = (int) (sizeof s0 / sizeof (TERM));
   static const int NS1 = (int) (sizeof s1 / sizeof (TERM));
   static const int NS2 = (int) (sizeof s2 / sizeof (TERM));
   static const int NS3 = (int) (sizeof s3 / sizeof (TERM));
   static const int NS4 = (int) (sizeof s4 / sizeof (TERM));

/* ------------------------------------------------------------------ */

/* Interval between fundamental epoch J2000.0 and current date (JC). */
   t = ((date1 - DJ00) + date2) / DJC;

/* Fundamental Arguments (from IERS Conventions 2003) */

/* Mean anomaly of the Moon. */
   fa[0] = iauFal03(t);

/* Mean anomaly of the Sun. */
   fa[1] = iauFalp03(t);

/* Mean longitude of the Moon minus that of the ascending node. */
   fa[2] = iauFaf03(t);

/* Mean elongation of the Moon from the Sun. */
   fa[3] = iauFad03(t);

/* Mean longitude of the ascending node of the Moon. */
   fa[4] = iauFaom03(t);

/* Mean longitude of Venus. */
   fa[5] = iauFave03(t);

/* Mean longitude of Earth. */
   fa[6] = iauFae03(t);

/* General precession in longitude. */
   fa[7] = iauFapa03(t);

/* Evaluate s. */
   w0 = sp[0];
   w1 = sp[1];
   w2 = sp[2];
   w3 = sp[3];
   w4 = sp[4];
   w5 = sp[5];

   for (i = NS0-1; i >= 0; i--) {
   a = 0.0;
   for (j = 0; j < 8; j++) {
      a += (double)s0[i].nfa[j] * fa[j];
   }
   w0 += s0[i].s * sin(a) + s0[i].c * cos(a);
   }

   for (i = NS1-1; i >= 0; i--) {
      a = 0.0;
      for (j = 0; j < 8; j++) {
         a += (double)s1[i].nfa[j] * fa[j];
      }
      w1 += s1[i].s * sin(a) + s1[i].c * cos(a);
   }

   for (i = NS2-1; i >= 0; i--) {
      a = 0.0;
      for (j = 0; j < 8; j++) {
         a += (double)s2[i].nfa[j] * fa[j];
      }
      w2 += s2[i].s * sin(a) + s2[i].c * cos(a);
   }

   for (i = NS3-1; i >= 0; i--) {
      a = 0.0;
      for (j = 0; j < 8; j++) {
         a += (double)s3[i].nfa[j] * fa[j];
      }
      w3 += s3[i].s * sin(a) + s3[i].c * cos(a);
   }

   for (i = NS4-1; i >= 0; i--) {
      a = 0.0;
      for (j = 0; j < 8; j++) {
         a += (double)s4[i].nfa[j] * fa[j];
      }
      w4 += s4[i].s * sin(a) + s4[i].c * cos(a);
   }

   s = (w0 +
       (w1 +
       (w2 +
       (w3 +
       (w4 +
        w5 * t) * t) * t) * t) * t) * DAS2R - x*y/2.0;

   return s;

/* Finished. */

/*----------------------------------------------------------------------
**
**  Copyright (C) 2023
**  Standards of Fundamental Astronomy Board
**  of the International Astronomical Union.
**
**  =====================
**  SOFA Software License
**  =====================
**
**  NOTICE TO USER:
**
**  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING SIX TERMS AND
**  CONDITIONS WHICH APPLY TO ITS USE.
**
**  1. The Software is owned by the IAU SOFA Board ("SOFA").
**
**  2. Permission is granted to anyone to use the SOFA software for any
**     purpose, including commercial applications, free of charge and
**     without payment of royalties, subject to the conditions and
**     restrictions listed below.
**
**  3. You (the user) may copy and distribute SOFA source code to others,
**     and use and adapt its code and algorithms in your own software,
**     on a world-wide, royalty-free basis.  That portion of your
**     distribution that does not consist of intact and unchanged copies
**     of SOFA source code files is a "derived work" that must comply
**     with the following requirements:
**
**     a) Your work shall be marked or carry a statement that it
**        (i) uses routines and computations derived by you from
**        software provided by SOFA under license to you; and
**        (ii) does not itself constitute software provided by and/or
**        endorsed by SOFA.
**
**     b) The source code of your derived work must contain descriptions
**        of how the derived work is based upon, contains and/or differs
**        from the original SOFA software.
**
**     c) The names of all routines in your derived work shall not
**        include the prefix "iau" or "sofa" or trivial modifications
**        thereof such as changes of case.
**
**     d) The origin of the SOFA components of your derived work must
**        not be misrepresented;  you must not claim that you wrote the
**        original software, nor file a patent application for SOFA
**        software or algorithms embedded in the SOFA software.
**
**     e) These requirements must be reproduced intact in any source
**        distribution and shall apply to anyone to whom you have
**        granted a further right to modify the source code of your
**        derived work.
**
**     Note that, as originally distributed, the SOFA software is
**     intended to be a definitive implementation of the IAU standards,
**     and consequently third-party modifications are discouraged.  All
**     variations, no matter how minor, must be explicitly marked as
**     such, as explained above.
**
**  4. You shall not cause the SOFA software to be brought into
**     disrepute, either by misuse, or use for inappropriate tasks, or
**     by inappropriate modification.
**
**  5. The SOFA software is provided "as is" and SOFA makes no warranty
**     as to its use or performance.   SOFA does not and cannot warrant
**     the performance or results which the user may obtain by using the
**     SOFA software.  SOFA makes no warranties, express or implied, as
**     to non-infringement of third party rights, merchantability, or
**     fitness for any particular purpose.  In no event will SOFA be
**     liable to the user for any consequential, incidental, or special
**     damages, including any lost profits or lost savings, even if a
**     SOFA representative has been advised of such damages, or for any
**     claim by any third party.
**
**  6. The provision of any version of the SOFA software under the terms
**     and conditions specified herein does not imply that future
**     versions will also be made available under the same terms and
**     conditions.
*
**  In any published work or commercial product which uses the SOFA
**  software directly, acknowledgement (see www.iausofa.org) is
**  appreciated.
**
**  Correspondence concerning SOFA software should be addressed as
**  follows:
**
**      By email:  sofa@ukho.gov.uk
**      By post:   IAU SOFA Center
**                 HM Nautical Almanac Office
**                 UK Hydrographic Office
**                 Admiralty Way, Taunton
**                 Somerset, TA1 2DN
**                 United Kingdom
**
**--------------------------------------------------------------------*/
}
