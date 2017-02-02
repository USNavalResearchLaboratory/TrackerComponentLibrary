%-Abstract
%
%   CSPICE_LATSPH converts latitudinal coordinates to spherical
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
%      radius   a double precision scalar or double precision 1xN array
%               describing the distance of the position from origin
%
%      lon      a double precision scalar or double precision 1xN array
%               the angle of the position from the XZ plane
%               measured in radians
%
%      lat      a double precision scalar or double precision 1xN array
%               the angle of the position from the XY plane
%               measured in radians
%
%   the call:
%
%      [rho, colat, lons] = cspice_latsph( radius, lon, lat)
%
%   returns:
%
%      rho     a double precision scalar or double precision 1XN
%              array describing the distance of the position from origin
%
%      colat   a double precision scalar or double precision 1XN
%              array describing the angle between the point and the
%              positive z-axis, measured in radians (also referred to
%              as the polar angle)
%
%      lons    a double precision scalar or double precision 1XN array
%              describing the angle of the projection of the point to the XY
%              plane from the positive X-axis, measured in radians,
%              with range:
%
%                  -pi < lons <= pi
%
%              The positive Y-axis is at longitude PI/2 radians.
%
%              The argument 'rho' returns in the same units associated
%              with 'radius'.
%
%              'rho', 'colat', and 'lons' return with the same vectorization
%              measure (N) as 'radius', 'lon', and 'lat'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%    Example (1):
%
%      %
%      % Load an SPK, leapseconds, and PCK kernel set.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Create a vector of scalar times.
%      %
%      et = [0:2]*2.*cspice_spd;
%
%      %
%      % Retrieve the position of the moon seen from earth at 'et'
%      % in the J2000 frame without aberration correction.
%      %
%      [pos, et] = cspice_spkpos( 'MOON', et, 'J2000', 'NONE', 'EARTH' );
%
%      %
%      % Convert the array of position vectors 'pos' to latitudinal
%      % coordinates.
%      %
%      [radius, longitude, latitude] = cspice_reclat(pos);
%
%      %
%      % Convert the latitudinal coords to spherical.
%      %
%      [rho, colat, lon] = cspice_latsph( radius, longitude, latitude);
%
%      %
%      % Convert the spherical to rectangular.
%      %
%      [rectan] = cspice_sphrec(radius, colat, lon);
%
%      %
%      % Calculate the relative error against the original position
%      % vectors.
%      %
%      (rectan-pos) ./ pos
%
%   MATLAB outputs:
%
%      1.0e-14 *
%
%                     0  -0.03701547067225   0.63783453323816
%      0.02182376758148   0.01641520435413  -0.01531271963894
%     -0.01912147275010  -0.04855217028457   0.02039513446643
%
%      The relative error between the original array of position vectors
%      and those that resulted from the various coordinate conversion
%      has magnitude on the order of 10^(-14).
%
%    Example (2):
%
%      %
%      % Define six sets of latitudinal coordinates, 'lon' and 'lat'
%      % expressed in degrees - converted to radians by use
%      % of cspice_rpd.
%      %
%      rad = [ 1.,  1., sqrt(2.), sqrt(2.),   1.,  0. ];
%      lon = [ 0., 90.,     180.,     180., 180., 33. ] * cspice_rpd;
%      lat = [ 0.,  0.,      45.,      -45., 90.,  0. ] * cspice_rpd;
%
%      %
%      % ...convert the latitudinal coordinates to spherical coordinates
%      %
%      [rho, colat, slon] = cspice_latsph(rad, lon, lat);
%
%      %
%      % ...convert angular measure to degrees.
%      %
%      colat = colat * cspice_dpr;
%      lon   = lon   * cspice_dpr;
%      slon  = slon  * cspice_dpr;
%      lat   = lat   * cspice_dpr;
%
%      %
%      % Output banner.
%      %
%      disp('    rho       colat       slon       rad        lon        lat   ')
%      disp('  --------   --------   --------   --------   --------   --------')
%
%      %
%      % Create an array of values for output.
%      %
%      output = [ rho; colat; slon; rad; lon; lat];
%      txt    = sprintf( '%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', ...
%                        output );
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
%          rho       colat       slon       rad        lon        lat
%        --------   --------   --------   --------   --------   --------
%          1.0000    90.0000     0.0000     1.0000     0.0000     0.0000
%          1.0000    90.0000    90.0000     1.0000    90.0000     0.0000
%          1.4142    45.0000   180.0000     1.4142   180.0000    45.0000
%          1.4142   135.0000   180.0000     1.4142   180.0000   -45.0000
%          1.0000     0.0000   180.0000     1.0000   180.0000    90.0000
%          0.0000    90.0000    33.0000     0.0000    33.0000     0.0000
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine latsph_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   latitudinal to spherical coordinates
%
%-&

function [radius, colat, lon] = cspice_latsph( radius, longitude, latitude)

   switch nargin
      case 3

         radius    = zzmice_dp(radius);
         longitude = zzmice_dp(longitude);
         latitude  = zzmice_dp(latitude);

      otherwise

         error ( ['Usage: [_radius_, _colat_, _lon_] = '...
                  'cspice_latsph( _radius_, _longitude_, _latitude_)' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [radius, colat, lon] = mice('latsph_c', radius, longitude, latitude);
   catch
      rethrow(lasterror)
   end





