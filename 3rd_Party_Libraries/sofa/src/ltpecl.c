#include "sofa.h"

void iauLtpecl(double epj, double vec[3])
/*
**  - - - - - - - - - -
**   i a u L t p e c l
**  - - - - - - - - - -
**
**  Long-term precession of the ecliptic.
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
**     vec     double[3]      ecliptic pole unit vector
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
**  This revision:  2016 February 9
**
**  SOFA release 2016-05-03
**
**  Copyright (C) 2016 IAU SOFA Board.  See notes at end.
*/
{
/* Obliquity at J2000.0 (radians). */
   static const double eps0 = 84381.406 * DAS2R;

/* Polynomial coefficients */
   enum { NPOL = 4 };
   static const double pqpol[2][NPOL] = {
      { 5851.607687,
          -0.1189000,
          -0.00028913,
           0.000000101},
      {-1600.886300,
           1.1689818,
          -0.00000020,
          -0.000000437}
   };

/* Periodic coefficients */
   static const double pqper[][5] = {
      { 708.15,-5486.751211,-684.661560,  667.666730,-5523.863691},
      {2309.00,  -17.127623,2446.283880,-2354.886252, -549.747450},
      {1620.00, -617.517403, 399.671049, -428.152441, -310.998056},
      { 492.20,  413.442940,-356.652376,  376.202861,  421.535876},
      {1183.00,   78.614193,-186.387003,  184.778874,  -36.776172},
      { 622.00, -180.732815,-316.800070,  335.321713, -145.278396},
      { 882.00,  -87.676083, 198.296701, -185.138669,  -34.744450},
      { 547.00,   46.140315, 101.135679, -120.972830,   22.885731}
   };
   static const int NPER = (int) ( sizeof pqper / 5 / sizeof (double) );

/* Miscellaneous */
   int i;
   double t, p, q, w, a, s, c;


/* Centuries since J2000. */
   t  = ( epj - 2000.0 ) / 100.0;

/* Initialize P_A and Q_A accumulators. */
   p = 0.0;
   q = 0.0;

/* Periodic terms. */
   w = D2PI*t;
   for ( i = 0; i < NPER; i++ ) {
      a = w/pqper[i][0];
      s = sin(a);
      c = cos(a);
      p += c*pqper[i][1] + s*pqper[i][3];
      q += c*pqper[i][2] + s*pqper[i][4];
   }

/* Polynomial terms. */
   w = 1.0;
   for ( i = 0; i < NPOL; i++ ) {
      p += pqpol[0][i]*w;
      q += pqpol[1][i]*w;
      w *= t;
   }

/* P_A and Q_A (radians). */
   p *= DAS2R;
   q *= DAS2R;

/* Form the ecliptic pole vector. */
   w = 1.0 - p*p - q*q;
   w = w < 0.0 ? 0.0 : sqrt(w);
   s = sin(eps0);
   c = cos(eps0);
   vec[0] = p;
   vec[1] = - q*c - w*s;
   vec[2] = - q*s + w*c;

/*----------------------------------------------------------------------
**
**  Copyright (C) 2016
**  Standards Of Fundamental Astronomy Board
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
