#include "sofa.h"
#include "sofam.h"

void iauFk45z(double r1950, double d1950, double bepoch,
              double *r2000, double *d2000)
/*
**  - - - - - - - - -
**   i a u F k 4 5 z
**  - - - - - - - - -
**
**  Convert a B1950.0 FK4 star position to J2000.0 FK5, assuming zero
**  proper motion in the FK5 system.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  This function converts a star's catalog data from the old FK4
**  (Bessel-Newcomb) system to the later IAU 1976 FK5 (Fricke) system,
**  in such a way that the FK5 proper motion is zero.  Because such a
**  star has, in general, a non-zero proper motion in the FK4 system,
**  the function requires the epoch at which the position in the FK4
**  system was determined.
**
**  Given:
**     r1950,d1950    double   B1950.0 FK4 RA,Dec at epoch (rad)
**     bepoch         double   Besselian epoch (e.g. 1979.3)
**
**  Returned:
**     r2000,d2000    double   J2000.0 FK5 RA,Dec (rad)
**
**  Notes:
**
**  1) The epoch bepoch is strictly speaking Besselian, but if a
**     Julian epoch is supplied the result will be affected only to a
**     negligible extent.
**
**  2) The method is from Appendix 2 of Aoki et al. (1983), but using
**     the constants of Seidelmann (1992).  See the function iauFk425
**     for a general introduction to the FK4 to FK5 conversion.
**
**  3) Conversion from equinox B1950.0 FK4 to equinox J2000.0 FK5 only
**     is provided for.  Conversions for different starting and/or
**     ending epochs would require additional treatment for precession,
**     proper motion and E-terms.
**
**  4) In the FK4 catalog the proper motions of stars within 10 degrees
**     of the poles do not embody differential E-terms effects and
**     should, strictly speaking, be handled in a different manner from
**     stars outside these regions.  However, given the general lack of
**     homogeneity of the star data available for routine astrometry,
**     the difficulties of handling positions that may have been
**     determined from astrometric fields spanning the polar and non-
**     polar regions, the likelihood that the differential E-terms
**     effect was not taken into account when allowing for proper motion
**     in past astrometry, and the undesirability of a discontinuity in
**     the algorithm, the decision has been made in this SOFA algorithm
**     to include the effects of differential E-terms on the proper
**     motions for all stars, whether polar or not.  At epoch J2000.0,
**     and measuring "on the sky" rather than in terms of RA change, the
**     errors resulting from this simplification are less than
**     1 milliarcsecond in position and 1 milliarcsecond per century in
**     proper motion.
**
**  References:
**
**     Aoki, S. et al., 1983, "Conversion matrix of epoch B1950.0
**     FK4-based positions of stars to epoch J2000.0 positions in
**     accordance with the new IAU resolutions".  Astron.Astrophys.
**     128, 263-267.
**
**     Seidelmann, P.K. (ed), 1992, "Explanatory Supplement to the
**     Astronomical Almanac", ISBN 0-935702-68-7.
**
**  Called:
**     iauAnp       normalize angle into range 0 to 2pi
**     iauC2s       p-vector to spherical
**     iauEpb2jd    Besselian epoch to Julian date
**     iauEpj       Julian date to Julian epoch
**     iauPdp       scalar product of two p-vectors
**     iauPmp       p-vector minus p-vector
**     iauPpsp      p-vector plus scaled p-vector
**     iauPvu       update a pv-vector
**     iauS2c       spherical to p-vector
**
**  This revision:   2023 March 4
**
**  SOFA release 2023-10-11
**
**  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
*/
{
/* Radians per year to arcsec per century */
   const double PMF = 100.0*DR2AS;

/* Position and position+velocity vectors */
   double r0[3], p[3], pv[2][3];

/* Miscellaneous */
   double w, djm0, djm;
   int i, j, k;

/*
** CANONICAL CONSTANTS (Seidelmann 1992)
*/

/* Vectors A and Adot (Seidelmann 3.591-2) */
   static double a[3]  = { -1.62557e-6, -0.31919e-6, -0.13843e-6 };
   static double ad[3] = { +1.245e-3,   -1.580e-3,   -0.659e-3   };

/* 3x2 matrix of p-vectors (cf. Seidelmann 3.591-4, matrix M) */
   static double em[2][3][3] = {
         { { +0.9999256782, -0.0111820611, -0.0048579477 },
           { +0.0111820610, +0.9999374784, -0.0000271765 },
           { +0.0048579479, -0.0000271474, +0.9999881997 } },
         { { -0.000551,     -0.238565,     +0.435739     },
           { +0.238514,     -0.002667,     -0.008541     },
           { -0.435623,     +0.012254,     +0.002117     } }
                               };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Spherical coordinates to p-vector. */
   iauS2c(r1950, d1950, r0);

/* Adjust p-vector A to give zero proper motion in FK5. */
   w  = (bepoch - 1950) / PMF;
   iauPpsp(a, w, ad, p);

/* Remove E-terms. */
   iauPpsp(p, -iauPdp(r0,p), r0, p);
   iauPmp(r0, p, p);

/* Convert to Fricke system pv-vector (cf. Seidelmann 3.591-3). */
   for ( i = 0; i < 2; i++ ) {
      for ( j = 0; j < 3; j++ ) {
         w = 0.0;
         for ( k = 0; k < 3; k++ ) {
            w += em[i][j][k] * p[k];
         }
         pv[i][j] = w;
      }
   }

/* Allow for fictitious proper motion. */
   iauEpb2jd(bepoch, &djm0, &djm);
   w = (iauEpj(djm0,djm)-2000.0) / PMF;
   iauPvu(w, pv, pv);

/* Revert to spherical coordinates. */
   iauC2s(pv[0], &w, d2000);
   *r2000 = iauAnp(w);

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
