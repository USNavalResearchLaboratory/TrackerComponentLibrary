#include "sofa.h"
#include "sofam.h"

void iauApco(double date1, double date2,
             double ebpv[2][3], double ehp[3],
             double x, double y, double s, double theta,
             double elong, double phi, double hm,
             double xp, double yp, double sp,
             double refa, double refb,
             iauASTROM *astrom)
/*
**  - - - - - - - -
**   i a u A p c o
**  - - - - - - - -
**
**  For a terrestrial observer, prepare star-independent astrometry
**  parameters for transformations between ICRS and observed
**  coordinates.  The caller supplies the Earth ephemeris, the Earth
**  rotation information and the refraction constants as well as the
**  site coordinates.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1  double       TDB as a 2-part...
**     date2  double       ...Julian Date (Note 1)
**     ebpv   double[2][3] Earth barycentric PV (au, au/day, Note 2)
**     ehp    double[3]    Earth heliocentric P (au, Note 2)
**     x,y    double       CIP X,Y (components of unit vector)
**     s      double       the CIO locator s (radians)
**     theta  double       Earth rotation angle (radians)
**     elong  double       longitude (radians, east +ve, Note 3)
**     phi    double       latitude (geodetic, radians, Note 3)
**     hm     double       height above ellipsoid (m, geodetic, Note 3)
**     xp,yp  double       polar motion coordinates (radians, Note 4)
**     sp     double       the TIO locator s' (radians, Note 4)
**     refa   double       refraction constant A (radians, Note 5)
**     refb   double       refraction constant B (radians, Note 5)
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
**      along  double       adjusted longitude (radians)
**      xpl    double       polar motion xp wrt local meridian (radians)
**      ypl    double       polar motion yp wrt local meridian (radians)
**      sphi   double       sine of geodetic latitude
**      cphi   double       cosine of geodetic latitude
**      diurab double       magnitude of diurnal aberration vector
**      eral   double       "local" Earth rotation angle (radians)
**      refa   double       refraction constant A (radians)
**      refb   double       refraction constant B (radians)
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
**  2) The vectors eb, eh, and all the astrom vectors, are with respect
**     to BCRS axes.
**
**  3) The geographical coordinates are with respect to the WGS84
**     reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN
**     CONVENTION:  the longitude required by the present function is
**     right-handed, i.e. east-positive, in accordance with geographical
**     convention.
**
**     The adjusted longitude stored in the astrom array takes into
**     account the TIO locator and polar motion.
**
**  4) xp and yp are the coordinates (in radians) of the Celestial
**     Intermediate Pole with respect to the International Terrestrial
**     Reference System (see IERS Conventions), measured along the
**     meridians 0 and 90 deg west respectively.  sp is the TIO locator
**     s', in radians, which positions the Terrestrial Intermediate
**     Origin on the equator.  For many applications, xp, yp and
**     (especially) sp can be set to zero.
**
**     Internally, the polar motion is stored in a form rotated onto the
**     local meridian.
**
**  5) The refraction constants refa and refb are for use in a
**     dZ = A*tan(Z)+B*tan^3(Z) model, where Z is the observed
**     (i.e. refracted) zenith distance and dZ is the amount of
**     refraction.
**
**  6) It is advisable to take great care with units, as even unlikely
**     values of the input parameters are accepted and processed in
**     accordance with the models used.
**
**  7) In cases where the caller does not wish to provide the Earth
**     Ephemeris, the Earth rotation information and refraction
**     constants, the function iauApco13 can be used instead of the
**     present function.  This starts from UTC and weather readings etc.
**     and computes suitable values using other SOFA functions.
**
**  8) This is one of several functions that inserts into the astrom
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
**  9) The context structure astrom produced by this function is used by
**     iauAtioq, iauAtoiq, iauAtciq* and iauAticq*.
**
**  Called:
**     iauIr        initialize r-matrix to identity
**     iauRz        rotate around Z-axis
**     iauRy        rotate around Y-axis
**     iauRx        rotate around X-axis
**     iauAnpm      normalize angle into range +/- pi
**     iauC2ixys    celestial-to-intermediate matrix, given X,Y and s
**     iauPvtob     position/velocity of terrestrial station
**     iauTrxpv     product of transpose of r-matrix and pv-vector
**     iauApcs      astrometry parameters, ICRS-GCRS, space observer
**     iauCr        copy r-matrix
**
**  This revision:   2021 February 24
**
**  SOFA release 2023-10-11
**
**  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
*/
{
   double r[3][3], a, b, eral, c, pvc[2][3], pv[2][3];


/* Form the rotation matrix, CIRS to apparent [HA,Dec]. */
   iauIr(r);
   iauRz(theta+sp, r);
   iauRy(-xp, r);
   iauRx(-yp, r);
   iauRz(elong, r);

/* Solve for local Earth rotation angle. */
   a = r[0][0];
   b = r[0][1];
   eral = ( a != 0.0 || b != 0.0 ) ?  atan2(b, a) : 0.0;
   astrom->eral = eral;

/* Solve for polar motion [X,Y] with respect to local meridian. */
   a = r[0][0];
   c = r[0][2];
   astrom->xpl = atan2(c, sqrt(a*a+b*b));
   a = r[1][2];
   b = r[2][2];
   astrom->ypl = ( a != 0.0 || b != 0.0 ) ? -atan2(a, b) : 0.0;

/* Adjusted longitude. */
   astrom->along = iauAnpm(eral - theta);

/* Functions of latitude. */
   astrom->sphi = sin(phi);
   astrom->cphi = cos(phi);

/* Refraction constants. */
   astrom->refa = refa;
   astrom->refb = refb;

/* Disable the (redundant) diurnal aberration step. */
   astrom->diurab = 0.0;

/* CIO based BPN matrix. */
   iauC2ixys(x, y, s, r);

/* Observer's geocentric position and velocity (m, m/s, CIRS). */
   iauPvtob(elong, phi, hm, xp, yp, sp, theta, pvc);

/* Rotate into GCRS. */
   iauTrxpv(r, pvc, pv);

/* ICRS <-> GCRS parameters. */
   iauApcs(date1, date2, pv, ebpv, ehp, astrom);

/* Store the CIO based BPN matrix. */
   iauCr(r, astrom->bpn );

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
