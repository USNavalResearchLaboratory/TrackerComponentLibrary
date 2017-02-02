%-Abstract
%
%   CSPICE_PGRREC converts planetographic coordinates to
%   rectangular coordinates.
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
%             is associated, optionally, you may supply the integer ID code for
%             the object as an integer string, i.e. both 'MOON' and '301' are
%             legitimate strings that indicate the Moon is the target body.
%
%             [1,m] = size(body); char = class(body)
%
%       lon    the planetographic longitude of the input point. This is the
%             angle between the prime meridian and the meridian containing the
%             input point. For bodies having prograde (aka direct) rotation,
%             the direction of increasing longitude is positive west:  from the
%             +X axis of the rectangular coordinate system toward the -Y axis.
%             For bodies having retrograde rotation, the direction of
%             increasing longitude is positive east: from the +X axis toward
%             the +Y axis.
%
%             [1,n] = size(lon); double = class(lon)
%
%             The earth, moon, and sun are exceptions: planetographic
%             longitude is measured positive east for these bodies.
%
%             The default interpretation of longitude by this
%             and the other planetographic coordinate conversion
%             routines can be overridden; see the discussion in
%             Particulars below for details.
%
%             'lon' is measured in radians. On input, the range
%             of longitude is unrestricted.
%
%       lat    the planetographic latitude of the input point.  For a point P
%             on the reference spheroid, this is the angle between the XY plane
%             and the outward normal vector at P. For a point P not on the
%             reference spheroid, the planetographic latitude is that of the
%             closest point to P on the spheroid.
%
%             [1,n] = size(lat); double = class(lat)
%
%             'lat' is measured in radians.  On input, the
%             range of latitude is unrestricted.
%
%       alt    the altitude above the reference spheroid.
%
%             [1,n] = size(alt); double = class(alt)
%
%             Units of 'alt' must match those of  're'.
%
%       re    equatorial radius of the body of interest.
%
%             [1,1] = size(re); double = class(re)
%
%       f     flattening coefficient of the body, a dimensionless value defined
%             as:
%
%                    equatorial_radius - polar_radius
%                    --------------------------------
%                           equatorial_radius
%
%             [1,1] = size(f); double = class(f)
%
%   the call:
%
%      rectan = cspice_pgrrec( body, lon, lat, alt, re, f)
%
%   returns:
%
%      rectan   the rectangular body-fixed coordinates of the position or set
%               of positions.
%
%               [3,n] = size(rectan); double = class(rectan)
%
%               'rectan' returns with the same units associated with
%               'alt' and 're'.
%
%               'rectan' returns with the same vectorization measure
%                (N) as 'lon', 'lat', and 'alt'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Load a PCK file containing a triaxial
%      % ellipsoidal shape model and orientation
%      % data for Mars.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Example 1: convert a single set of planetographic
%      %            coordinates to rectangular bodyfixed
%      %            coordinates.
%      %
%      % Look up the radii for Mars.  Although we
%      % omit it here, we could check the kernel pool
%      % to make sure the variable BODY499_RADII
%      % has three elements and numeric data type.
%      % If the variable is not present in the kernel
%      % pool, cspice_bodvrd will signal an error.
%      %
%      body = 'MARS';
%      radii = cspice_bodvrd( body, 'RADII', 3 );
%
%      %
%      %
%      % Calculate the flatness coefficient. Set a bodyfixed
%      % position vector, 'x'.
%      %
%      re   = radii(1);
%      rp   = radii(3);
%      flat = ( re - rp ) / re;
%
%      % Set a longitude, latitude, altitude position.
%      % Note that we must provide longitude and
%      % latitude in radians.
%      %
%      lon  = 90. * cspice_rpd;
%      lat  = 45.  * cspice_rpd;
%      alt  = 3.d2;
%
%      %
%      % Do the conversion.
%      %
%      x = cspice_pgrrec( body, lon, lat, alt, re, flat );
%
%      %
%      % Output.
%      %
%      disp( 'Scalar:' )
%      disp(' ')
%
%      disp( 'Rectangular coordinates in km (x, y, z)' )
%      fprintf( '%9.3f   %9.3f   %9.3f\n', x' )
%
%      disp( 'Planetographic coordinates in degs and km (lon, lat, alt)' )
%      fprintf( '%9.3f   %9.3f   %9.3f\n', lon *cspice_dpr() ...
%                                        , lat *cspice_dpr() ...
%                                        , alt               )
%      disp(' ')
%
%
%      %
%      % Example 2: convert a vectorized set of planetographic coordinates
%      %            to rectangular bodyfixed coordinates.
%      %
%      % Define 1xN arrays of longitudes, latitudes, and altitudes.
%      %
%      lon = [ 0.,   ...
%              180., ...
%              180., ...
%              180., ...
%              90.,  ...
%              270., ...
%              0.,   ...
%              0.,   ...
%              0. ];
%
%      lat = [ 0.,  ...
%              0.,  ...
%              0.,  ...
%              0.,  ...
%              0.,  ...
%              0.,  ...
%              90., ...
%             -90., ...
%              90. ];
%
%      alt = [ 0., ...
%              0., ...
%              10., ...
%              10., ...
%              0., ...
%              0., ...
%              0., ...
%              0., ...
%             -3376.200 ];
%
%      %
%      % Convert angular measures to radians.
%      %
%      lon = lon*cspice_rpd;
%      lat = lat*cspice_rpd;
%
%      %
%      % Using the same Mars parameters, convert the 'lon', 'lat', 'alt'
%      % vectors to bodyfixed rectangular coordinates.
%      %
%      x = cspice_pgrrec( body, lon, lat, alt, re, flat);
%
%      disp('Vector:')
%      disp(' ')
%
%      disp( ['rectan(1)   rectan(2)   rectan(3)' ...
%             '         lon         lat         alt'] )
%      disp( ['---------------------------------' ...
%             '------------------------------------'] )
%
%      %
%      % Create an array of values for output.
%      %
%      output = [  x(1,:);         x(2,:);         x(3,:); ...
%                  lon*cspice_dpr; lat*cspice_dpr; alt ];
%
%      txt = sprintf( '%9.3f   %9.3f   %9.3f   %9.3f   %9.3f   %9.3f\n', ...
%                                                                   output);
%      disp( txt )
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Scalar:
%
%      Rectangular coordinates in km (x, y, z)
%          0.000   -2620.679    2592.409
%      Planetographic coordinates in degs and km (lon, lat, alt)
%         90.000      45.000     300.000
%
%      Vector:
%
%      rectan(1)   rectan(2)   rectan(3)         lon         lat         alt
%      ---------------------------------------------------------------------
%       3396.190      -0.000       0.000       0.000       0.000       0.000
%      -3396.190      -0.000       0.000     180.000       0.000       0.000
%      -3406.190      -0.000       0.000     180.000       0.000      10.000
%      -3406.190      -0.000       0.000     180.000       0.000      10.000
%          0.000   -3396.190       0.000      90.000       0.000       0.000
%         -0.000    3396.190       0.000     270.000       0.000       0.000
%          0.000      -0.000    3376.200       0.000      90.000       0.000
%          0.000      -0.000   -3376.200       0.000     -90.000       0.000
%          0.000       0.000       0.000       0.000      90.000   -3376.200
%
%-Particulars
%
%   Given the planetographic coordinates of a point, this routine
%   returns the body-fixed rectangular coordinates of the point.  The
%   body-fixed rectangular frame is that having the X-axis pass
%   through the 0 degree latitude 0 degree longitude direction, the
%   Z-axis pass through the 90 degree latitude direction, and the
%   Y-axis equal to the cross product of the unit Z-axis and X-axis
%   vectors.
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
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine pgrrec_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%      Corrected misspellings.
%
%   -Mice Version 1.0.0, 22-JAN-2008, EDW (JPL)
%
%-Index_Entries
%
%   convert planetographic to rectangular coordinates
%
%-&

function [rectan] = cspice_pgrrec( body, lon, lat, alt, re, f)

   switch nargin
      case 6

         body = zzmice_str(body);
         lon  = zzmice_dp(lon);
         lat  = zzmice_dp(lat);
         alt  = zzmice_dp(alt);
         re   = zzmice_dp(re);
         f    = zzmice_dp(f);

      otherwise

         error ( ['Usage: [_rectan(3)_] = ' ...
                  'cspice_pgrrec( body, _lon_, _lat_, _alt_, re, f)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [rectan] = mice( 'pgrrec_c', body, lon, lat, alt, re, f);
   catch
      rethrow(lasterror)
   end
