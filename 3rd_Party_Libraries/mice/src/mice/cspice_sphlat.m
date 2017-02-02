%-Abstract
%
%   CSPICE_SPHLAT converts spherical coordinates to latitudinal
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
%      r       a double precision scalar or double precision 1XN
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
%   the call:
%
%      [radius, lon, lat] = cspice_sphlat(r, colat, lons)
%
%   returns:
%
%      radius   a double precision scalar or double precision 1xN array
%               describing the distance of the position from origin
%
%      lon      a double precision scalar or double precision 1xN array
%               describing the angle of the position from the XZ plane
%               measured in radians
%
%      lat      a double precision scalar or double precision 1xN array
%               describing the angle of the position from the XY plane
%               measured in radians
%
%               The argument 'radius' returns in the same units associated
%               with 'r'.
%
%               'radius', 'lon', and 'lat' return with the same
%                vectorization measure (N) as the 'r', 'colat',
%                and 'lons'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
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
%      [r, colat, lons] = cspice_recsph(pos);
%
%      %
%      % Convert the latitudinal coords to spherical.
%      %
%      [ radius, lon, lat] = cspice_sphlat(r, colat, lons);
%
%      %
%      % Convert the spherical coords to rectangular.
%      %
%      [rectan] = cspice_latrec( radius, lon, lat);
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
%                     0  -0.05552320600838   0.63783453323816
%      0.02182376758148                  0  -0.01531271963894
%      0.01912147275010  -0.02427608514229   0.02039513446643
%
%   Example(2):
%
%      %
%      % Define six sets of spherical coordinates, 'lon' and 'colat'
%      % expressed in degrees - converted to radians by use of cspice_rpd.
%      %
%      r     = [  1.,  1., 1.4142, 1.4142, 1.  , 0. ];
%      colat = [ 90., 90., 45.   , 135.  , 0.  , 0. ] * cspice_rpd;
%      lons  = [  0., 90., 180.  , 180.  , 180., 33.] * cspice_rpd;
%
%      %
%      % ...convert the latitudinal coordinates to spherical coordinates
%      %
%      [rad, lon, lat] = cspice_sphlat(r, colat, lons);
%
%      %
%      % ...convert angular measure to degrees.
%      %
%      colat = colat * cspice_dpr;
%      lon   = lon   * cspice_dpr;
%      lons  = lons  * cspice_dpr;
%      lat   = lat   * cspice_dpr;
%
%      %
%      % Output banner.
%      %
%      disp('    r         colat       lons       rad        lon        lat   ')
%      disp('  --------   --------   --------   --------   --------   --------')
%
%      %
%      % Create an array of values for output.
%      %
%      output = [ r; colat; lons; rad; lon; lat];
%      txt    = sprintf( '%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', output);
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
%       r         colat       lons       rad        lon        lat
%     --------   --------   --------   --------   --------   --------
%       1.0000    90.0000     0.0000     1.0000     0.0000     0.0000
%       1.0000    90.0000    90.0000     1.0000    90.0000     0.0000
%       1.4142    45.0000   180.0000     1.4142   180.0000    45.0000
%       1.4142   135.0000   180.0000     1.4142   180.0000   -45.0000
%       1.0000     0.0000   180.0000     1.0000   180.0000    90.0000
%       0.0000     0.0000    33.0000     0.0000    33.0000    90.0000
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine sphlat_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 12-DEC-2005, EDW (JPL)
%
%-Index_Entries
%
%   spherical to latitudinal coordinates
%
%-&

function [radius, lon, lat] = cspice_sphlat(r, colat, lons)

   switch nargin
      case 3

         r     = zzmice_dp(r);
         colat = zzmice_dp(colat);
         lons  = zzmice_dp(lons);

      otherwise

         error ( [ 'Usage: [_radius_, _lon_, _lat_] = ' ...
                   'cspice_sphlat(_r_, _colat_, _lons_)' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [radius, lon, lat] = mice('sphlat_c', r, colat, lons);
   catch
      rethrow(lasterror)
   end

