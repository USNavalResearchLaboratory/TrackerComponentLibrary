%-Abstract
%
%   CSPICE_GEOREC converts geodetic coordinates to rectangular
%   coordinates.
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
%       lon   a double precision scalar or 1xN array describing
%             the geodetic longitude measured in radians
%
%       lat   a double precision scalar or 1xN array describing
%             the geodetic latitude measured in radians
%
%       alt   a double precision scalar or 1xN array describing
%             the altitude above the reference spheroid
%
%       re    the scalar, double precision equatorial radius of
%             the body of interest
%
%       f     the scalar, double precision flattening coefficient
%             of the body, a dimensionless value defined as:
%
%                    equatorial_radius - polar_radius
%                    --------------------------------
%                           equatorial_radius
%
%   the call:
%
%      rectan = cspice_georec( lon, lat, alt, re, f)
%
%   returns:
%
%      rectan   a double precision 3x1 array or double precision
%               3xN array containing the rectangular coordinates of the
%               position or set of positions
%
%               'rectan' returns with the same units associated with
%               'alt' and 're'
%
%               'rectan' returns with the same vectorization measure
%                (N) as 'lon', 'lat', and 'alt'
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Load the standard kernel set.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Retrieve the triaxial radii of the earth
%      %
%      radii = cspice_bodvrd( 'EARTH', 'RADII', 3 );
%
%      %
%      % Calculate the flatness coefficient.
%      %
%      flat = (radii(1) - radii(3))/radii(1);
%
%      %
%      % Set a latitude,
%      % longitude, altitude coordinate at 118 west,
%      % 32 North, 0 altitude (convert the angular measures
%      % to radians).
%      %
%      lon  = 118. * cspice_rpd;
%      lat  = 32.  * cspice_rpd;
%      alt  = 0.;
%
%      x = cspice_georec( lon, lat, alt, radii(1), flat );
%
%      disp( 'Scalar:' )
%      txt = sprintf( '%14.6f   %14.6f   %14.6f', x );
%      disp( txt )
%
%      disp( ' ' )
%
%      %
%      % Using the equatorial radius of the Clark66 spheroid
%      % (CLARKR = 6378.2064 km) and the Clark 66 flattening
%      % factor (CLARKF = 1.0 / 294.9787 ) convert to
%      % body fixed rectangular coordinates.
%      %
%      CLARKR = 6378.2064;
%      CLARKF = 1./294.9787;
%
%      %
%      % Define a vector of scalar longitudes, latitudes, and altitudes.
%      % This is a vector of scalars (1xN), NOT the same as an N-vector (Nx1).
%      %
%      lon = [  0., ...
%               0., ...
%              90., ...
%               0., ...
%             180., ...
%             -90., ...
%               0., ...
%              45., ...
%               0., ...
%              90., ...
%              45. ];
%
%      lat = [ 90.      , ...
%              88.677225, ...
%              88.677225, ...
%              90.      , ...
%              88.677225, ...
%              88.677225, ...
%              -90.     , ...
%              88.129144, ...
%              88.707084, ...
%              88.707084, ...
%              88.171393 ];
%
%      alt = [ -6356.5838  , ...
%              -6356.572258, ...
%              -6356.572258, ...
%              -6355.5838,   ...
%              -6356.572258, ...
%              -6356.572258, ...
%              -6355.5838,   ...
%              -6356.560715, ...
%              -6355.572518, ...
%              -6355.572518, ...
%              -6355.561236  ];
%
%      %
%      % Convert angular measures to radians.
%      %
%      lon = lon*cspice_rpd;
%      lat = lat*cspice_rpd;
%
%      %
%      % Calculate then output the rectangular coordinates.
%      %
%      x = cspice_georec( lon, lat, alt, CLARKR, CLARKF);
%
%      disp( 'Vector:' )
%
%      %
%      % Create an array of values for output.
%      %
%      output = [  x(1,:);  x(2,:);  x(3,:) ];
%      txt    = sprintf( '%14.6f   %14.6f   %14.6f\n', output);
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
%        -2541.748162      4780.333036      3360.428190
%
%      Vector:
%            0.000000         0.000000         0.000000
%            1.000000         0.000000        -0.000000
%            0.000000         1.000000        -0.000000
%            0.000000         0.000000         1.000000
%           -1.000000         0.000000        -0.000000
%            0.000000        -1.000000        -0.000000
%            0.000000         0.000000        -1.000000
%            1.000000         1.000000         0.000000
%            1.000000         0.000000         1.000000
%            0.000000         1.000000         1.000000
%            1.000000         1.000000         1.000000
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine georec_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   geodetic to rectangular coordinates
%
%-&

function [rectan] = cspice_georec( lon, lat, alt, re, f)

   switch nargin
      case 5

         lon  = zzmice_dp(lon);
         lat  = zzmice_dp(lat);
         alt  = zzmice_dp(alt);
         re   = zzmice_dp(re);
         f    = zzmice_dp(f);

      otherwise

         error ( ['Usage: [_rectan(3)_] = ' ...
                  'cspice_georec( _lon_, _lat_, _alt_, re, f)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [rectan] = mice( 'georec_c', lon, lat, alt, re, f);
   catch
      rethrow(lasterror)
   end



