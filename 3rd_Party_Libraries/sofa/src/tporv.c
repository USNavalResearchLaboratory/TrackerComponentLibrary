#include "sofa.h"

int iauTporv(double xi, double eta, double v[3],
             double v01[3], double v02[3])
/*
**  - - - - - - - - -
**   i a u T p o r v
**  - - - - - - - - -
**
**  In the tangent plane projection, given the rectangular coordinates
**  of a star and its direction cosines, determine the direction
**  cosines of the tangent point.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     xi,eta   double    rectangular coordinates of star image (Note 2)
**     v        double[3] star's direction cosines (Note 3)
**
**  Returned:
**     v01      double[3] tangent point's direction cosines, Solution 1
**     v02      double[3] tangent point's direction cosines, Solution 2
**
**  Returned (function value):
**                int     number of solutions:
**                        0 = no solutions returned (Note 4)
**                        1 = only the first solution is useful (Note 5)
**                        2 = both solutions are useful (Note 5)
**
**  Notes:
**
**  1) The tangent plane projection is also called the "gnomonic
**     projection" and the "central projection".
**
**  2) The eta axis points due north in the adopted coordinate system.
**     If the direction cosines represent observed (RA,Dec), the tangent
**     plane coordinates (xi,eta) are conventionally called the
**     "standard coordinates".  If the direction cosines are with
**     respect to a right-handed triad, (xi,eta) are also right-handed.
**     The units of (xi,eta) are, effectively, radians at the tangent
**     point.
**
**  3) The vector v must be of unit length or the result will be wrong.
**
**  4) Cases where there is no solution can arise only near the poles.
**     For example, it is clearly impossible for a star at the pole
**     itself to have a non-zero xi value, and hence it is meaningless
**     to ask where the tangent point would have to be.
**
**  5) Also near the poles, cases can arise where there are two useful
**     solutions.  The return value indicates whether the second of the
**     two solutions returned is useful;  1 indicates only one useful
**     solution, the usual case.
**
**  6) The basis of the algorithm is to solve the spherical triangle
**     PSC, where P is the north celestial pole, S is the star and C is
**     the tangent point.  Calling the celestial spherical coordinates
**     of the star and tangent point (a,b) and (a0,b0) respectively, and
**     writing rho^2 = (xi^2+eta^2) and r^2 = (1+rho^2), and
**     transforming the vector v into (a,b) in the normal way, side c is
**     then (pi/2-b), side p is sqrt(xi^2+eta^2) and side s (to be
**     found) is (pi/2-b0), while angle C is given by sin(C) = xi/rho
**     and cos(C) = eta/rho;  angle P (to be found) is (a-a0).  After
**     solving the spherical triangle, the result (a0,b0) can be
**     expressed in vector form as v0.
**
**  7) This function is a member of the following set:
**
**         spherical      vector         solve for
**
**         iauTpxes      iauTpxev         xi,eta
**         iauTpsts      iauTpstv          star
**         iauTpors    > iauTporv <       origin
**
**  References:
**
**     Calabretta M.R. & Greisen, E.W., 2002, "Representations of
**     celestial coordinates in FITS", Astron.Astrophys. 395, 1077
**
**     Green, R.M., "Spherical Astronomy", Cambridge University Press,
**     1987, Chapter 13.
**
**  This revision:   2018 January 2
**
**  SOFA release 2023-10-11
**
**  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
*/
{
   double x, y, z, rxy2, xi2, eta2p1, r, rsb, rcb, w2, w, c;


   x = v[0];
   y = v[1];
   z = v[2];
   rxy2 = x*x + y*y;
   xi2 = xi*xi;
   eta2p1 = eta*eta + 1.0;
   r = sqrt(xi2 + eta2p1);
   rsb = r*z;
   rcb = r*sqrt(x*x + y*y);
   w2 = rcb*rcb - xi2;
   if ( w2 > 0.0 ) {
      w = sqrt(w2);
      c = (rsb*eta + w) / (eta2p1*sqrt(rxy2*(w2+xi2)));
      v01[0] = c * (x*w + y*xi);
      v01[1] = c * (y*w - x*xi);
      v01[2] = (rsb - eta*w) / eta2p1;
      w = - w;
      c = (rsb*eta + w) / (eta2p1*sqrt(rxy2*(w2+xi2)));
      v02[0] = c * (x*w + y*xi);
      v02[1] = c * (y*w - x*xi);
      v02[2] = (rsb - eta*w) / eta2p1;
      return (fabs(rsb) < 1.0) ? 1 : 2;
   } else {
      return 0;
   }

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
