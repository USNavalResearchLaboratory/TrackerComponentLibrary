#include "sofa.h"
#include "sofam.h"

void iauApcs(double date1, double date2, double pv[2][3],
             double ebpv[2][3], double ehp[3],
             iauASTROM *astrom)
/*
**  - - - - - - - -
**   i a u A p c s
**  - - - - - - - -
**
**  For an observer whose geocentric position and velocity are known,
**  prepare star-independent astrometry parameters for transformations
**  between ICRS and GCRS.  The Earth ephemeris is supplied by the
**  caller.
**
**  The parameters produced by this function are required in the space
**  motion, parallax, light deflection and aberration parts of the
**  astrometric transformation chain.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1  double       TDB as a 2-part...
**     date2  double       ...Julian Date (Note 1)
**     pv     double[2][3] observer's geocentric pos/vel (m, m/s)
**     ebpv   double[2][3] Earth barycentric PV (au, au/day)
**     ehp    double[3]    Earth heliocentric P (au)
**
**  Returned:
**     astrom iauASTROM*   star-independent astrometry parameters:
**      pmt    double       PM time interval (SSB, Julian years)
**      eb     double[3]    SSB to observer (vector, au)
**      eh     double[3]    Sun to observer (unit vector)
**      em     double       distance from Sun to observer (au)
**      v      double[3]    barycentric observer velocity (vector, c)
**      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
**      bpn    double[3][3] bias-precession-nutation matrix
**      along  double       unchanged
**      xpl    double       unchanged
**      ypl    double       unchanged
**      sphi   double       unchanged
**      cphi   double       unchanged
**      diurab double       unchanged
**      eral   double       unchanged
**      refa   double       unchanged
**      refb   double       unchanged
**
**  Notes:
**
**  1) The TDB date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TDB)=2450123.7 could be expressed in any of these ways, among
**     others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in cases
**     where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 method is best matched to the way the
**     argument is handled internally and will deliver the optimum
**     resolution.  The MJD method and the date & time methods are both
**     good compromises between resolution and convenience.  For most
**     applications of this function the choice will not be at all
**     critical.
**
**     TT can be used instead of TDB without any significant impact on
**     accuracy.
**
**  2) All the vectors are with respect to BCRS axes.
**
**  3) Providing separate arguments for (i) the observer's geocentric
**     position and velocity and (ii) the Earth ephemeris is done for
**     convenience in the geocentric, terrestrial and Earth orbit cases.
**     For deep space applications it maybe more convenient to specify
**     zero geocentric position and velocity and to supply the
**     observer's position and velocity information directly instead of
**     with respect to the Earth.  However, note the different units:
**     m and m/s for the geocentric vectors, au and au/day for the
**     heliocentric and barycentric vectors.
**
**  4) In cases where the caller does not wish to provide the Earth
**     ephemeris, the function iauApcs13 can be used instead of the
**     present function.  This computes the Earth ephemeris using the
**     SOFA function iauEpv00.
**
**  5) This is one of several functions that inserts into the astrom
**     structure star-independent parameters needed for the chain of
**     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
**
**     The various functions support different classes of observer and
**     portions of the transformation chain:
**
**          functions         observer        transformation
**
**       iauApcg iauApcg13    geocentric      ICRS <-> GCRS
**       iauApci iauApci13    terrestrial     ICRS <-> CIRS
**       iauApco iauApco13    terrestrial     ICRS <-> observed
**       iauApcs iauApcs13    space           ICRS <-> GCRS
**       iauAper iauAper13    terrestrial     update Earth rotation
**       iauApio iauApio13    terrestrial     CIRS <-> observed
**
**     Those with names ending in "13" use contemporary SOFA models to
**     compute the various ephemerides.  The others accept ephemerides
**     supplied by the caller.
**
**     The transformation from ICRS to GCRS covers space motion,
**     parallax, light deflection, and aberration.  From GCRS to CIRS
**     comprises frame bias and precession-nutation.  From CIRS to
**     observed takes account of Earth rotation, polar motion, diurnal
**     aberration and parallax (unless subsumed into the ICRS <-> GCRS
**     transformation), and atmospheric refraction.
**
**  6) The context structure astrom produced by this function is used by
**     iauAtciq* and iauAticq*.
**
**  Called:
**     iauCp        copy p-vector
**     iauPm        modulus of p-vector
**     iauPn        decompose p-vector into modulus and direction
**     iauIr        initialize r-matrix to identity
**
**  This revision:   2021 February 24
**
**  SOFA release 2023-10-11
**
**  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
*/
{
/* au/d to m/s */
   const double AUDMS = DAU/DAYSEC;

/* Light time for 1 au (day) */
   const double CR = AULT/DAYSEC;

   int i;
   double dp, dv, pb[3], vb[3], ph[3], v2, w;


/* Time since reference epoch, years (for proper motion calculation). */
   astrom->pmt = ( (date1 - DJ00) + date2 ) / DJY;

/* Adjust Earth ephemeris to observer. */
   for (i = 0; i < 3; i++) {
      dp = pv[0][i] / DAU;
      dv = pv[1][i] / AUDMS;
      pb[i] = ebpv[0][i] + dp;
      vb[i] = ebpv[1][i] + dv;
      ph[i] = ehp[i] + dp;
   }

/* Barycentric position of observer (au). */
   iauCp(pb, astrom->eb);

/* Heliocentric direction and distance (unit vector and au). */
   iauPn(ph, &astrom->em, astrom->eh);

/* Barycentric vel. in units of c, and reciprocal of Lorenz factor. */
   v2 = 0.0;
   for (i = 0; i < 3; i++) {
      w = vb[i] * CR;
      astrom->v[i] = w;
      v2 += w*w;
   }
   astrom->bm1 = sqrt(1.0 - v2);

/* Reset the NPB matrix. */
   iauIr(astrom->bpn);

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
