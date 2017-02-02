%-Abstract
%
%   CSPICE_DRDPGR computes the Jacobian matrix of the transformation
%   from planetographic to rectangular coordinates.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED
%   "AS-IS" TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING
%   ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR
%   A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY,
%   OR NASA BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING,
%   BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
%   ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY
%   AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR
%   NASA BE ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL
%   KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE
%   OF THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO
%   INDEMNIFY CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING
%   FROM THE ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      body   name of the body with which the planetographic coordinate system
%             is associated.
%
%             [1,m] = size(body); char = class(body)
%
%             `body' is used by this routine to look up from the
%             kernel pool the prime meridian rate coefficient giving
%             the body's spin sense.
%
%      lon    planetographic longitude of the input point. This is the angle
%             between the prime meridian and the meridian containing the input
%             point. For bodies having prograde (aka direct) rotation, the
%             direction of increasing longitude is positive west:  from the +X
%             axis of the rectangular coordinate system toward the -Y axis. For
%             bodies having retrograde rotation, the direction of increasing
%             longitude is positive east: from the +X axis toward the +Y axis.
%
%             [1,n] = size(lon); double = class(lon)
%
%             The earth, moon, and sun are exceptions:
%             planetographic longitude is measured positive east for
%             these bodies.
%
%             The default interpretation of longitude by this
%             and the other planetographic coordinate conversion
%             routines can be overridden; see the discussion in
%             Particulars below for details.
%
%             Longitude is measured in radians. On input, the range
%             of longitude is unrestricted.
%
%      lat    planetographic latitude of the input point.  For a point P on the
%             reference spheroid, this is the angle  between the XY plane and
%             the outward normal vector at P. For a point P not on the
%             reference spheroid, the planetographic latitude is that of the
%             closest point to P on the spheroid.
%
%             [1,n] = size(lat); double = class(lat)
%
%             Latitude is measured in radians.  On input, the range of
%             latitude is unrestricted.
%
%      alt    Altitude of point above the reference spheroid. Units of `alt'
%             must match those of `re'.
%
%             [1,n] = size(alt); double = class(alt)
%
%      re     equatorial radius of a reference spheroid. This spheroid is a
%             volume of revolution: its horizontal cross sections are circular.
%              The shape of the spheroid is defined by an equatorial radius
%             `re' and a polar radius `rp'. Units of `re' must match those of
%             `alt'.
%
%             [1,1] = size(re); double = class(re)
%
%      f      the flattening coefficient
%
%             [1,1] = size(f); double = class(f)
%
%                f = (re-rp) / re
%
%             where rp is the polar radius of the spheroid. (More importantly
%             rp = re*(1-f).) The units of `rp' match those of `re'.
%
%   the call:
%
%      jacobi = cspice_drdpgr( body, lon, lat, alt, re, f)
%
%   returns:
%
%      jacobi   the matrix of partial derivatives of the conversion from
%               planetographic to rectangular coordinates evaluated at the
%               input coordinates. This matrix has the form
%
%               If [1,1] = size(lon) then [3,3]   = size(jacobi)
%               If [1,n] = size(lon) then [3,3,n] = size(jacobi)
%                                          double = class(jacobi)
%
%                   -                              -
%                  |  dx/dlon   dx/dlat   dx/dalt   |
%                  |                                |
%                  |  dy/dlon   dy/dlat   dy/dalt   |
%                  |                                |
%                  |  dz/dlon   dz/dlat   dz/dalt   |
%                   -                              -
%
%               evaluated at the input values of 'lon', 'lat' and 'alt'.
%
%-Examples
%
%   None.
%
%-Particulars
%
%   It is often convenient to describe the motion of an object in the
%   planetographic coordinate system.  However, when performing
%   vector computations it's hard to beat rectangular coordinates.
%
%   To transform states given with respect to planetographic
%   coordinates to states with respect to rectangular coordinates,
%   one makes use of the Jacobian of the transformation between the
%   two systems.
%
%   Given a state in planetographic coordinates
%
%      ( lon, lat, alt, dlon, dlat, dalt )
%
%   the velocity in rectangular coordinates is given by the matrix
%   equation:
%
%                  t          |                                  t
%      (dx, dy, dz)   = jacobi|              * (dlon, dlat, dalt)
%                             |(lon,lat,alt)
%
%
%   This routine computes the matrix
%
%            |
%      jacobi|
%            |(lon,lat,alt)
%
%
%   In the planetographic coordinate system, longitude is defined
%   using the spin sense of the body.  Longitude is positive to the
%   west if the spin is prograde and positive to the east if the spin
%   is retrograde.  The spin sense is given by the sign of the first
%   degree term of the time-dependent polynomial for the body's prime
%   meridian Euler angle "W":  the spin is retrograde if this term is
%   negative and prograde otherwise.  For the sun, planets, most
%   natural satellites, and selected asteroids, the polynomial
%   expression for W may be found in a SPICE PCK kernel.
%
%   The earth, moon, and sun are exceptions: planetographic longitude
%   is measured positive east for these bodies.
%
%   If you wish to override the default sense of positive longitude
%   for a particular body, you can do so by defining the kernel
%   variable
%
%      BODY<body ID>_PGR_POSITIVE_LON
%
%   where <body ID> represents the NAIF ID code of the body. This
%   variable may be assigned either of the values
%
%      'WEST'
%      'EAST'
%
%   For example, you can have this routine treat the longitude
%   of the earth as increasing to the west using the kernel
%   variable assignment
%
%      BODY399_PGR_POSITIVE_LON = 'WEST'
%
%   Normally such assignments are made by placing them in a text
%   kernel and loading that kernel via furnsh_c.
%
%   The definition of this kernel variable controls the behavior of
%   the CSPICE planetographic routines
%
%      cspice_pgrrec
%      cspice_recpgr
%      cspice_dpgrdr
%      cspice_drdpgr
%
%   It does not affect the other SPICE coordinate conversion
%   routines.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine drdpgr_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 09-NOV-2012, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   Jacobian of rectangular w.r.t. planetographic coordinates
%
%-&

function [jacobi] = cspice_drdpgr( body, lon, lat, alt, re, f)

   switch nargin
      case 6

         body = zzmice_str(body);
         lon  = zzmice_dp(lon);
         lat  = zzmice_dp(lat);
         alt  = zzmice_dp(alt);
         re   = zzmice_dp(re);
         f    = zzmice_dp(f);

      otherwise

         error( ['Usage: [_jacobi(3,3)_] = '...
                 'cspice_drdpgr( `body`, _lon_, _lat_, _alt_, re, f)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [jacobi] = mice('drdpgr_c', body, lon, lat, alt, re, f);
   catch
      rethrow(lasterror)
   end




