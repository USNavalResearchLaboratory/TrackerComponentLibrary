%-Abstract
%
%   CSPICE_RECGEO converts rectangular coordinates to geodetic
%   coordinates.
%%
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
%      rectan   a double precision 3x1 array or double precision
%               3xN array containing the rectangular coordinates of the
%               position or set of positions
%
%       re      the scalar, double precision equatorial radius of
%               the body of interest
%
%       f       the scalar, double precision flattening coefficient
%               of the body, a dimensionless value defined as:
%
%                    equatorial_radius - polar_radius
%                    --------------------------------
%                           equatorial_radius
%
%   the call:
%
%      [ lon, lat, alt ] = cspice_recgeo( rectan, re, f)
%
%   returns:
%
%       lon   a double precision scalar or 1XN-vector describing
%             the geodetic longitude measured in radians.
%
%       lat   a double precision scalar or 1XN-vector describing
%             the geodetic latitude measured in radians.
%
%       alt   a double precision scalar or 1XN-vector describing
%             the altitude above the reference spheroid.
%
%             'lon', 'lat', and 'alt' return with the same vectorization
%             measure (N) as 'rectan'.
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
%      % Calculate the flatness coefficient. Set a bodyfixed
%      % position.
%      %
%      flat = (radii(1) - radii(3))/radii(1);
%      x    = [ -2541.748162; 4780.333036; 3360.428190];
%
%      [ lon, lat, alt] = cspice_recgeo( x, radii(1), flat );
%
%      %
%      % Output, convert the angular values to degrees.
%      %
%      lon = lon * cspice_dpr;
%      lat = lat * cspice_dpr;
%
%      disp('Scalar:')
%      txt = sprintf( '%12.8f   %12.8f   %12.8f', lon , lat , alt );
%      disp( txt )
%
%      disp(' ')
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
%      x = [ [ 0, 1, 0, 0, -1,  0,  0, 1, 1, 0, 1];
%            [ 0, 0, 1, 0,  0, -1,  0, 1, 0, 1, 1];
%            [ 0, 0, 0, 1,  0,  0, -1, 0, 1, 1, 1] ];
%
%      [ lon, lat, alt] = cspice_recgeo(  x, CLARKR, CLARKF );
%
%      %
%      % Output, convert the angular values to degrees.
%      %
%      lon = lon * cspice_dpr;
%      lat = lat * cspice_dpr;
%
%      disp('Vector:')
%
%      %
%      % Output banner.
%      %
%    disp('    lon        lat          alt         x          y          z    ')
%    disp('  --------   --------   ----------   --------   --------   --------')
%
%      output = [ lon; lat; alt; x(1,:); x(2,:); x(3,:) ];
%      txt    = sprintf( '%10.4f %10.4f %12.6f %10.4f %10.4f %10.4f\n',output);
%      disp(txt)
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
%         118.00000000    32.00000000     0.00000024
%
%      118 degrees west, 32 north, 0.24 mm altitude.
%
%   Vector:
%       lon        lat          alt         x          y          z
%     --------   --------   ----------   --------   --------   --------
%       0.0000    90.0000 -6356.583800     0.0000     0.0000     0.0000
%       0.0000    88.6772 -6356.572258     1.0000     0.0000     0.0000
%      90.0000    88.6772 -6356.572258     0.0000     1.0000     0.0000
%       0.0000    90.0000 -6355.583800     0.0000     0.0000     1.0000
%     180.0000    88.6772 -6356.572258    -1.0000     0.0000     0.0000
%     -90.0000    88.6772 -6356.572258     0.0000    -1.0000     0.0000
%       0.0000   -90.0000 -6355.583800     0.0000     0.0000    -1.0000
%      45.0000    88.1291 -6356.560715     1.0000     1.0000     0.0000
%       0.0000    88.7071 -6355.572518     1.0000     0.0000     1.0000
%      90.0000    88.7071 -6355.572518     0.0000     1.0000     1.0000
%      45.0000    88.1714 -6355.561236     1.0000     1.0000     1.0000
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine recgeo_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   rectangular to geodetic
%
%-&

function [lon, lat, alt] = cspice_recgeo(rectan, re, f)

   switch nargin
      case 3

         rectan = zzmice_dp(rectan);
         re     = zzmice_dp(re);
         f      = zzmice_dp(f);

      otherwise
         error ( ['Usage: [_lon_, _lat_, _alt_] = '...
                  'cspice_recgeo(_rectan(3)_, re, f)' ] )
   end

   %
   % Call the MEX library.
   %
   try
      [lon, lat, alt] = mice( 'recgeo_c', rectan, re, f);
   catch
      rethrow(lasterror)
   end






