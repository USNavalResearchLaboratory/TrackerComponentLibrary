#include "sofa.h"
#include "sofam.h"
#include <stdlib.h>

void iauMoon98 ( double date1, double date2, double pv[2][3] )
/*
**  - - - - - - - - - -
**   i a u M o o n 9 8
**  - - - - - - - - - -
**
**  Approximate geocentric position and velocity of the Moon.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  n.b. Not IAU-endorsed and without canonical status.
**
**  Given:
**     date1  double         TT date part A (Notes 1,4)
**     date2  double         TT date part B (Notes 1,4)
**
**  Returned:
**     pv     double[2][3]   Moon p,v, GCRS (au, au/d, Note 5)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways, among
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
**     good compromises between resolution and convenience.  The limited
**     accuracy of the present algorithm is such that any of the methods
**     is satisfactory.
**
**  2) This function is a full implementation of the algorithm
**     published by Meeus (see reference) except that the light-time
**     correction to the Moon's mean longitude has been omitted.
**
**  3) Comparisons with ELP/MPP02 over the interval 1950-2100 gave RMS
**     errors of 2.9 arcsec in geocentric direction, 6.1 km in position
**     and 36 mm/s in velocity.  The worst case errors were 18.3 arcsec
**     in geocentric direction, 31.7 km in position and 172 mm/s in
**     velocity.
**
**  4) The original algorithm is expressed in terms of "dynamical time",
**     which can either be TDB or TT without any significant change in
**     accuracy.  UT cannot be used without incurring significant errors
**     (30 arcsec in the present era) due to the Moon's 0.5 arcsec/sec
**     movement.
**
**  5) The result is with respect to the GCRS (the same as J2000.0 mean
**     equator and equinox to within 23 mas).
**
**  6) Velocity is obtained by a complete analytical differentiation
**     of the Meeus model.
**
**  7) The Meeus algorithm generates position and velocity in mean
**     ecliptic coordinates of date, which the present function then
**     rotates into GCRS.  Because the ecliptic system is precessing,
**     there is a coupling between this spin (about 1.4 degrees per
**     century) and the Moon position that produces a small velocity
**     contribution.  In the present function this effect is neglected
**     as it corresponds to a maximum difference of less than 3 mm/s and
**     increases the RMS error by only 0.4%.
**
**  References:
**
**     Meeus, J., Astronomical Algorithms, 2nd edition, Willmann-Bell,
**     1998, p337.
**
**     Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G. & Laskar, J., Astron.Astrophys., 1994, 282, 663
**
**  Defined in sofam.h:
**     DAU           astronomical unit (m)
**     DJC           days per Julian century
**     DJ00          reference epoch (J2000.0), Julian Date
**     DD2R          degrees to radians
**
**  Called:
**     iauS2pv      spherical coordinates to pv-vector
**     iauPfw06     bias-precession F-W angles, IAU 2006
**     iauIr        initialize r-matrix to identity
**     iauRz        rotate around Z-axis
**     iauRx        rotate around X-axis
**     iauRxpv      product of r-matrix and pv-vector
**
**  This revision:  2023 March 20
**
**  SOFA release 2023-10-11
**
**  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
*/
{
/*
**  Coefficients for fundamental arguments:
**
**  . Powers of time in Julian centuries
**  . Units are degrees.
*/

/* Moon's mean longitude (wrt mean equinox and ecliptic of date) */
   static double elp0 = 218.31665436,        /* Simon et al. (1994). */
                 elp1 = 481267.88123421,
                 elp2 = -0.0015786,
                 elp3 = 1.0 / 538841.0,
                 elp4 = -1.0 / 65194000.0;
   double elp, delp;

/* Moon's mean elongation */
   static double d0 = 297.8501921,
                 d1 = 445267.1114034,
                 d2 = -0.0018819,
                 d3 = 1.0 / 545868.0,
                 d4 = 1.0 / 113065000.0;
   double d, dd;

/* Sun's mean anomaly */
   static double em0 = 357.5291092,
                 em1 = 35999.0502909,
                 em2 = -0.0001536,
                 em3 = 1.0 / 24490000.0,
                 em4 = 0.0;
   double em, dem;

/* Moon's mean anomaly */
   static double emp0 = 134.9633964,
                 emp1 = 477198.8675055,
                 emp2 = 0.0087414,
                 emp3 = 1.0 / 69699.0,
                 emp4 = -1.0 / 14712000.0;
   double emp, demp;

/* Mean distance of the Moon from its ascending node */
   static double f0 = 93.2720950,
                 f1 = 483202.0175233,
                 f2 = -0.0036539,
                 f3 = 1.0 / 3526000.0,
                 f4 = 1.0 / 863310000.0;
   double f, df;

/*
** Other arguments
*/

/* Meeus A_1, due to Venus (deg) */
   static double a10 = 119.75,
                 a11 = 131.849;
   double a1, da1;

/* Meeus A_2, due to Jupiter (deg) */
   static double a20 = 53.09,
                 a21 = 479264.290;
   double a2, da2;

/* Meeus A_3, due to sidereal motion of the Moon in longitude (deg) */
   static double a30 = 313.45,
                 a31 = 481266.484;
   double a3, da3;

/* Coefficients for Meeus "additive terms" (deg) */
   static double al1 =  0.003958,
                 al2 =  0.001962,
                 al3 =  0.000318;
   static double ab1 = -0.002235,
                 ab2 =  0.000382,
                 ab3 =  0.000175,
                 ab4 =  0.000175,
                 ab5 =  0.000127,
                 ab6 = -0.000115;

/* Fixed term in distance (m) */
   static double r0 = 385000560.0;

/* Coefficients for (dimensionless) E factor */
   static double e1 = -0.002516,
                 e2 = -0.0000074;
   double e, de, esq, desq;

/*
** Coefficients for Moon longitude and distance series
*/
   struct termlr {
      int nd;           /* multiple of D  in argument           */
      int nem;          /*     "    "  M   "    "               */
      int nemp;         /*     "    "  M'  "    "               */
      int nf;           /*     "    "  F   "    "               */
      double coefl;     /* coefficient of L sine argument (deg) */
      double coefr;     /* coefficient of R cosine argument (m) */
   };

static struct termlr tlr[] = {{0,  0,  1,  0,  6.288774, -20905355.0},
                              {2,  0, -1,  0,  1.274027,  -3699111.0},
                              {2,  0,  0,  0,  0.658314,  -2955968.0},
                              {0,  0,  2,  0,  0.213618,   -569925.0},
                              {0,  1,  0,  0, -0.185116,     48888.0},
                              {0,  0,  0,  2, -0.114332,     -3149.0},
                              {2,  0, -2,  0,  0.058793,    246158.0},
                              {2, -1, -1,  0,  0.057066,   -152138.0},
                              {2,  0,  1,  0,  0.053322,   -170733.0},
                              {2, -1,  0,  0,  0.045758,   -204586.0},
                              {0,  1, -1,  0, -0.040923,   -129620.0},
                              {1,  0,  0,  0, -0.034720,    108743.0},
                              {0,  1,  1,  0, -0.030383,    104755.0},
                              {2,  0,  0, -2,  0.015327,     10321.0},
                              {0,  0,  1,  2, -0.012528,         0.0},
                              {0,  0,  1, -2,  0.010980,     79661.0},
                              {4,  0, -1,  0,  0.010675,    -34782.0},
                              {0,  0,  3,  0,  0.010034,    -23210.0},
                              {4,  0, -2,  0,  0.008548,    -21636.0},
                              {2,  1, -1,  0, -0.007888,     24208.0},
                              {2,  1,  0,  0, -0.006766,     30824.0},
                              {1,  0, -1,  0, -0.005163,     -8379.0},
                              {1,  1,  0,  0,  0.004987,    -16675.0},
                              {2, -1,  1,  0,  0.004036,    -12831.0},
                              {2,  0,  2,  0,  0.003994,    -10445.0},
                              {4,  0,  0,  0,  0.003861,    -11650.0},
                              {2,  0, -3,  0,  0.003665,     14403.0},
                              {0,  1, -2,  0, -0.002689,     -7003.0},
                              {2,  0, -1,  2, -0.002602,         0.0},
                              {2, -1, -2,  0,  0.002390,     10056.0},
                              {1,  0,  1,  0, -0.002348,      6322.0},
                              {2, -2,  0,  0,  0.002236,     -9884.0},
                              {0,  1,  2,  0, -0.002120,      5751.0},
                              {0,  2,  0,  0, -0.002069,         0.0},
                              {2, -2, -1,  0,  0.002048,     -4950.0},
                              {2,  0,  1, -2, -0.001773,      4130.0},
                              {2,  0,  0,  2, -0.001595,         0.0},
                              {4, -1, -1,  0,  0.001215,     -3958.0},
                              {0,  0,  2,  2, -0.001110,         0.0},
                              {3,  0, -1,  0, -0.000892,      3258.0},
                              {2,  1,  1,  0, -0.000810,      2616.0},
                              {4, -1, -2,  0,  0.000759,     -1897.0},
                              {0,  2, -1,  0, -0.000713,     -2117.0},
                              {2,  2, -1,  0, -0.000700,      2354.0},
                              {2,  1, -2,  0,  0.000691,         0.0},
                              {2, -1,  0, -2,  0.000596,         0.0},
                              {4,  0,  1,  0,  0.000549,     -1423.0},
                              {0,  0,  4,  0,  0.000537,     -1117.0},
                              {4, -1,  0,  0,  0.000520,     -1571.0},
                              {1,  0, -2,  0, -0.000487,     -1739.0},
                              {2,  1,  0, -2, -0.000399,         0.0},
                              {0,  0,  2, -2, -0.000381,     -4421.0},
                              {1,  1,  1,  0,  0.000351,         0.0},
                              {3,  0, -2,  0, -0.000340,         0.0},
                              {4,  0, -3,  0,  0.000330,         0.0},
                              {2, -1,  2,  0,  0.000327,         0.0},
                              {0,  2,  1,  0, -0.000323,      1165.0},
                              {1,  1, -1,  0,  0.000299,         0.0},
                              {2,  0,  3,  0,  0.000294,         0.0},
                              {2,  0, -1, -2,  0.000000,      8752.0}};

   static int NLR = ( sizeof tlr / sizeof ( struct termlr ) );

/*
** Coefficients for Moon latitude series
*/
   struct termb {
      int nd;           /* multiple of D  in argument           */
      int nem;          /*     "    "  M   "    "               */
      int nemp;         /*     "    "  M'  "    "               */
      int nf;           /*     "    "  F   "    "               */
      double coefb;     /* coefficient of B sine argument (deg) */
   };

static struct termb tb[] = {{0,  0,  0,  1,  5.128122},
                            {0,  0,  1,  1,  0.280602},
                            {0,  0,  1, -1,  0.277693},
                            {2,  0,  0, -1,  0.173237},
                            {2,  0, -1,  1,  0.055413},
                            {2,  0, -1, -1,  0.046271},
                            {2,  0,  0,  1,  0.032573},
                            {0,  0,  2,  1,  0.017198},
                            {2,  0,  1, -1,  0.009266},
                            {0,  0,  2, -1,  0.008822},
                            {2, -1,  0, -1,  0.008216},
                            {2,  0, -2, -1,  0.004324},
                            {2,  0,  1,  1,  0.004200},
                            {2,  1,  0, -1, -0.003359},
                            {2, -1, -1,  1,  0.002463},
                            {2, -1,  0,  1,  0.002211},
                            {2, -1, -1, -1,  0.002065},
                            {0,  1, -1, -1, -0.001870},
                            {4,  0, -1, -1,  0.001828},
                            {0,  1,  0,  1, -0.001794},
                            {0,  0,  0,  3, -0.001749},
                            {0,  1, -1,  1, -0.001565},
                            {1,  0,  0,  1, -0.001491},
                            {0,  1,  1,  1, -0.001475},
                            {0,  1,  1, -1, -0.001410},
                            {0,  1,  0, -1, -0.001344},
                            {1,  0,  0, -1, -0.001335},
                            {0,  0,  3,  1,  0.001107},
                            {4,  0,  0, -1,  0.001021},
                            {4,  0, -1,  1,  0.000833},
                            {0,  0,  1, -3,  0.000777},
                            {4,  0, -2,  1,  0.000671},
                            {2,  0,  0, -3,  0.000607},
                            {2,  0,  2, -1,  0.000596},
                            {2, -1,  1, -1,  0.000491},
                            {2,  0, -2,  1, -0.000451},
                            {0,  0,  3, -1,  0.000439},
                            {2,  0,  2,  1,  0.000422},
                            {2,  0, -3, -1,  0.000421},
                            {2,  1, -1,  1, -0.000366},
                            {2,  1,  0,  1, -0.000351},
                            {4,  0,  0,  1,  0.000331},
                            {2, -1,  1,  1,  0.000315},
                            {2, -2,  0, -1,  0.000302},
                            {0,  0,  1,  3, -0.000283},
                            {2,  1,  1, -1, -0.000229},
                            {1,  1,  0, -1,  0.000223},
                            {1,  1,  0,  1,  0.000223},
                            {0,  1, -2, -1, -0.000220},
                            {2,  1, -1, -1, -0.000220},
                            {1,  0,  1,  1, -0.000185},
                            {2, -1, -2, -1,  0.000181},
                            {0,  1,  2,  1, -0.000177},
                            {4,  0, -2, -1,  0.000176},
                            {4, -1, -1, -1,  0.000166},
                            {1,  0,  1, -1, -0.000164},
                            {4,  0,  1, -1,  0.000132},
                            {1,  0, -1, -1, -0.000119},
                            {4, -1,  0, -1,  0.000115},
                            {2, -2,  0,  1,  0.000107}};

   static int NB = ( sizeof tb / sizeof ( struct termb ) );

/* Miscellaneous */
   int n, i;
   double t, elpmf, delpmf, vel, vdel, vr, vdr, a1mf, da1mf, a1pf,
          da1pf, dlpmp, slpmp, vb, vdb, v, dv, emn, empn, dn, fn, en,
          den, arg, darg, farg, coeff, el, del, r, dr, b, db, gamb,
          phib, psib, epsa, rm[3][3];

/* ------------------------------------------------------------------ */

/* Centuries since J2000.0 */
   t = ((date1 - DJ00) + date2) / DJC;

/* --------------------- */
/* Fundamental arguments */
/* --------------------- */

/* Arguments (radians) and derivatives (radians per Julian century)
   for the current date. */

/* Moon's mean longitude. */
   elp = DD2R * fmod ( elp0
                   + ( elp1
                   + ( elp2
                   + ( elp3
                   +   elp4 * t ) * t ) * t ) * t, 360.0 );
   delp = DD2R * (     elp1
                   + ( elp2 * 2.0
                   + ( elp3 * 3.0
                   +   elp4 * 4.0 * t ) * t ) * t );

/* Moon's mean elongation. */
   d = DD2R * fmod ( d0
                 + ( d1
                 + ( d2
                 + ( d3
                 +   d4 * t ) * t ) * t ) * t, 360.0 );
   dd = DD2R * (     d1
                 + ( d2 * 2.0
                 + ( d3 * 3.0
                 +   d4 * 4.0 * t ) * t ) * t );

/* Sun's mean anomaly. */
   em = DD2R * fmod ( em0
                  + ( em1
                  + ( em2
                  + ( em3
                  +   em4 * t ) * t ) * t ) * t, 360.0 );
   dem = DD2R * (     em1
                  + ( em2 * 2.0
                  + ( em3 * 3.0
                  +   em4 * 4.0 * t ) * t ) * t );

/* Moon's mean anomaly. */
   emp = DD2R * fmod ( emp0
                   + ( emp1
                   + ( emp2
                   + ( emp3
                   +   emp4 * t ) * t ) * t ) * t, 360.0 );
   demp = DD2R * (     emp1
                   + ( emp2 * 2.0
                   + ( emp3 * 3.0
                   +   emp4 * 4.0 * t ) * t ) * t );

/* Mean distance of the Moon from its ascending node. */
   f = DD2R * fmod ( f0
                 + ( f1
                 + ( f2
                 + ( f3
                 +   f4 * t ) * t ) * t ) * t, 360.0 );
   df = DD2R * (     f1
                 + ( f2 * 2.0
                 + ( f3 * 3.0
                 +   f4 * 4.0 * t ) * t ) * t );

/* Meeus further arguments. */
   a1 = DD2R * ( a10 + a11*t );
   da1 = DD2R * al1;
   a2 = DD2R * ( a20 + a21*t );
   da2 = DD2R * a21;
   a3 = DD2R * ( a30 + a31*t );
   da3 = DD2R * a31;

/* E-factor, and square. */
   e = 1.0 + ( e1 + e2*t ) * t;
   de = e1 + 2.0*e2*t;
   esq = e*e;
   desq = 2.0*e*de;

/* Use the Meeus additive terms (deg) to start off the summations. */
   elpmf = elp - f;
   delpmf = delp - df;
   vel = al1 * sin(a1)
       + al2 * sin(elpmf)
       + al3 * sin(a2);
   vdel = al1 * cos(a1) * da1
        + al2 * cos(elpmf) * delpmf
        + al3 * cos(a2) * da2;

   vr = 0.0;
   vdr = 0.0;

   a1mf = a1 - f;
   da1mf = da1 - df;
   a1pf = a1 + f;
   da1pf = da1 + df;
   dlpmp = elp - emp;
   slpmp = elp + emp;
   vb = ab1 * sin(elp)
      + ab2 * sin(a3)
      + ab3 * sin(a1mf)
      + ab4 * sin(a1pf)
      + ab5 * sin(dlpmp)
      + ab6 * sin(slpmp);
   vdb = ab1 * cos(elp) * delp
       + ab2 * cos(a3) * da3
       + ab3 * cos(a1mf) * da1mf
       + ab4 * cos(a1pf) * da1pf
       + ab5 * cos(dlpmp) * (delp-demp)
       + ab6 * cos(slpmp) * (delp+demp);

/* ----------------- */
/* Series expansions */
/* ----------------- */

/* Longitude and distance plus derivatives. */
   for ( n = NLR-1; n >= 0; n-- ) {
      dn = (double) tlr[n].nd;
      emn = (double) ( i = tlr[n].nem );
      empn = (double) tlr[n].nemp;
      fn = (double) tlr[n].nf;
      switch ( abs(i) ) {
      case 1:
         en = e;
         den = de;
         break;
      case 2:
         en = esq;
         den = desq;
         break;
      default:
         en = 1.0;
         den = 0.0;
      }
      arg = dn*d + emn*em + empn*emp + fn*f;
      darg = dn*dd + emn*dem + empn*demp + fn*df;
      farg = sin(arg);
      v = farg * en;
      dv = cos(arg)*darg*en + farg*den;
      coeff = tlr[n].coefl;
      vel += coeff * v;
      vdel += coeff * dv;
      farg = cos(arg);
      v = farg * en;
      dv = -sin(arg)*darg*en + farg*den;
      coeff = tlr[n].coefr;
      vr += coeff * v;
      vdr += coeff * dv;
   }
   el = elp + DD2R*vel;
   del = ( delp + DD2R*vdel ) / DJC;
   r = ( vr + r0 ) / DAU;
   dr = vdr / DAU / DJC;

/* Latitude plus derivative. */
   for ( n = NB-1; n >= 0; n-- ) {
      dn = (double) tb[n].nd;
      emn = (double) ( i = tb[n].nem );
      empn = (double) tb[n].nemp;
      fn = (double) tb[n].nf;
      switch ( abs(i) ) {
      case 1:
         en = e;
         den = de;
         break;
      case 2:
         en = esq;
         den = desq;
         break;
      default:
         en = 1.0;
         den = 0.0;
      }
      arg = dn*d + emn*em + empn*emp + fn*f;
      darg = dn*dd + emn*dem + empn*demp + fn*df;
      farg = sin(arg);
      v = farg * en;
      dv = cos(arg)*darg*en + farg*den;
      coeff = tb[n].coefb;
      vb += coeff * v;
      vdb += coeff * dv;
   }
   b = vb * DD2R;
   db = vdb * DD2R / DJC;

/* ------------------------------ */
/* Transformation into final form */
/* ------------------------------ */

/* Longitude, latitude to x, y, z (au). */
   iauS2pv ( el, b, r, del, db, dr, pv );

/* IAU 2006 Fukushima-Williams bias+precession angles. */
   iauPfw06 ( date1, date2, &gamb, &phib, &psib, &epsa );

/* Mean ecliptic coordinates to GCRS rotation matrix. */
   iauIr ( rm );
   iauRz ( psib, rm );
   iauRx ( -phib, rm );
   iauRz ( -gamb, rm );

/* Rotate the Moon position and velocity into GCRS (Note 6). */
   iauRxpv ( rm, pv, pv );

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
