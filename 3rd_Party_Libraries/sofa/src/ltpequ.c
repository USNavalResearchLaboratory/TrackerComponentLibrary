#include "sofa.h"
#include "sofam.h"

void iauLtpequ(double epj, double veq[3])
/*
**  - - - - - - - - - -
**   i a u L t p e q u
**  - - - - - - - - - -
**
**  Long-term precession of the equator.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     epj     double         Julian epoch (TT)
**
**  Returned:
**     veq     double[3]      equator pole unit vector
**
**  Notes:
**
**  1) The returned vector is with respect to the J2000.0 mean equator
**     and equinox.
**
**  2) The Vondrak et al. (2011, 2012) 400 millennia precession model
**     agrees with the IAU 2006 precession at J2000.0 and stays within
**     100 microarcseconds during the 20th and 21st centuries.  It is
**     accurate to a few arcseconds throughout the historical period,
**     worsening to a few tenths of a degree at the end of the
**     +/- 200,000 year time span.
**
**  References:
**
**    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
**    expressions, valid for long time intervals, Astron.Astrophys. 534,
**    A22
**
**    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
**    expressions, valid for long time intervals (Corrigendum),
**    Astron.Astrophys. 541, C1
**
**  This revision:  2021 May 11
**
**  SOFA release 2023-10-11
**
**  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
*/
{
/* Polynomial coefficients */
   enum { NPOL = 4 };
   static const double xypol[2][NPOL] = {
      {  5453.282155,
            0.4252841,
           -0.00037173,
           -0.000000152},
      {-73750.930350,
           -0.7675452,
           -0.00018725,
            0.000000231}
   };

/* Periodic coefficients */
   static const double xyper[][5] = {
      { 256.75, -819.940624,75004.344875,81491.287984, 1558.515853},
      { 708.15,-8444.676815,  624.033993,  787.163481, 7774.939698},
      { 274.20, 2600.009459, 1251.136893, 1251.296102,-2219.534038},
      { 241.45, 2755.175630,-1102.212834,-1257.950837,-2523.969396},
      {2309.00, -167.659835,-2660.664980,-2966.799730,  247.850422},
      { 492.20,  871.855056,  699.291817,  639.744522, -846.485643},
      { 396.10,   44.769698,  153.167220,  131.600209,-1393.124055},
      { 288.90, -512.313065, -950.865637, -445.040117,  368.526116},
      { 231.10, -819.415595,  499.754645,  584.522874,  749.045012},
      {1610.00, -538.071099, -145.188210,  -89.756563,  444.704518},
      { 620.00, -189.793622,  558.116553,  524.429630,  235.934465},
      { 157.87, -402.922932,  -23.923029,  -13.549067,  374.049623},
      { 220.30,  179.516345, -165.405086, -210.157124, -171.330180},
      {1200.00,   -9.814756,    9.344131,  -44.919798,  -22.899655}
   };
   static const int NPER = (int) ( sizeof xyper / 5 / sizeof (double) );

/* Miscellaneous */
   int i;
   double t, x, y, w, a, s, c;


/* Centuries since J2000. */
   t  = ( epj - 2000.0 ) / 100.0;

/* Initialize X and Y accumulators. */
   x = 0.0;
   y = 0.0;

/* Periodic terms. */
   w = D2PI * t;
   for ( i = 0; i < NPER; i++ ) {
      a = w / xyper[i][0];
      s = sin(a);
      c = cos(a);
      x += c*xyper[i][1] + s*xyper[i][3];
      y += c*xyper[i][2] + s*xyper[i][4];
   }

/* Polynomial terms. */
   w = 1.0;
   for ( i = 0; i < NPOL; i++ ) {
      x += xypol[0][i]*w;
      y += xypol[1][i]*w;
      w *= t;
   }

/* X and Y (direction cosines). */
   x *= DAS2R;
   y *= DAS2R;

/* Form the equator pole vector. */
   veq[0] = x;
   veq[1] = y;
   w = 1.0 - x*x - y*y;
   veq[2] = w < 0.0 ? 0.0 : sqrt(w);

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
