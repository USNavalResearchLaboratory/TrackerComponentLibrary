%-Abstract
%
%   CSPICE_DPGRDR computes the Jacobian matrix of the transformation
%   from rectangular to planetographic coordinates.
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
%      x      [1,n] = size(x); double = class(x)
%
%      y      [1,n] = size(y); double = class(y)
%
%      z      [1,n] = size(z); double = class(z)
%
%             the rectangular coordinates of the point at which the Jacobian of
%             the map from rectangular to planetographic coordinates is
%             desired.
%
%      re     equatorial radius of a reference spheroid. This spheroid is a
%             volume of revolution: its horizontal cross sections are circular.
%             The shape of the spheroid is defined by an equatorial radius `re'
%             and a polar radius `rp'.
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
%             rp = re*(1-f).)
%
%   the call:
%
%      jacobi = cspice_dpgrdr( body, x, y, z, re, f)
%
%   returns:
%
%      jacobi   the matrix of partial derivatives of the conversion from
%               rectangular to planetographic coordinates. It has the form
%
%               If [1,1] = size(x) then [3,3]   = size(jacobi).
%               If [1,n] = size(x) then [3,3,n] = size(jacobi).
%               double = class(jacobi)
%
%                   -                               -
%                  |  dlon/dx    dlon/dy   dlon/dz   |
%                  |                                 |
%                  |  dlat/dx    dlat/dy   dlat/dz   |
%                  |                                 |
%                  |  dalt/dx    dalt/dy   dalt/dz   |
%                   -                               -
%
%               evaluated at the input values of 'x', 'y', and 'z'.
%
%-Examples
%
%   None.
%
%-Particulars
%
%   When performing vector calculations with velocities it is usually
%   most convenient to work in rectangular coordinates. However, once
%   the vector manipulations have been performed, it is often
%   desirable to convert the rectangular representations into
%   planetographic coordinates to gain insights about phenomena in
%   this coordinate frame.
%
%   To transform rectangular velocities to derivatives of coordinates
%   in a planetographic system, one uses the Jacobian of the
%   transformation between the two systems.
%
%   Given a state in rectangular coordinates
%
%      ( x, y, z, dx, dy, dz )
%
%   the velocity in planetographic coordinates is given by the matrix
%   equation:
%                        t          |                     t
%      (dlon, dlat, dalt)   = jacobi|       * (dx, dy, dz)
%                                   |(x,y,z)
%
%   This routine computes the matrix
%
%            |
%      jacobi|
%            |(x, y, z)
%
%
%   The planetographic definition of latitude is identical to the
%   planetodetic (also called "geodetic" in SPICE documentation)
%   definition. In the planetographic coordinate system, latitude is
%   defined using a reference spheroid.  The spheroid is
%   characterized by an equatorial radius and a polar radius. For a
%   point P on the spheroid, latitude is defined as the angle between
%   the X-Y plane and the outward surface normal at P.  For a point P
%   off the spheroid, latitude is defined as the latitude of the
%   nearest point to P on the spheroid.  Note if P is an interior
%   point, for example, if P is at the center of the spheroid, there
%   may not be a unique nearest point to P.
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
%   It does not affect the other CSPICE coordinate conversion
%   routines.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dpgrdr_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 11-NOV-2013, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   Jacobian of planetographic  w.r.t. rectangular coordinates
%
%-&

function [jacobi] = cspice_dpgrdr( body, x, y, z, re, f)

   switch nargin
      case 6

         body = zzmice_str(body);
         x    = zzmice_dp(x);
         y    = zzmice_dp(y);
         z    = zzmice_dp(z);
         re   = zzmice_dp(re);
         f    = zzmice_dp(f);

      otherwise

         error( ['Usage: [_jacobi(3,3)_] = ' ...
                 'cspice_dpgrdr( `body`, _x_, _y_, _z_, re, f)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [jacobi] = mice('dpgrdr_c', body, x, y, z, re, f);
   catch
      rethrow(lasterror)
   end




