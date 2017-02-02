%-Abstract
%
%   CSPICE_RECSPH converts rectangular (Cartesian) coordinates to
%   spherical coordinates. All coordinates are expressed as
%   double precision values.
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
%      rectan   a double precision 3x1 array or double precision
%               3xN array containing the rectangular coordinates of the
%               position or set of positions
%
%   the call:
%
%      [r, colat, lon] = cspice_recsph(rectan)
%
%   returns:
%
%      r       a double precision scalar or double precision 1XN
%              array describing the distance of the position from origin
%
%      colat   a double precision scalar or double precision 1XN
%              array describing the angle between the point and the
%              positive z-axis, measured in radians (also referred to
%              as the polar angle)
%
%      lon     a double precision scalar or double precision 1XN array
%              describing the angle of the projection of the point to the XY
%              plane from the positive X-axis, measured in radians,
%              with range:
%
%                  -pi < lon <= pi
%
%              The positive Y-axis is at longitude PI/2 radians.
%
%              The argument 'r' returns with the same units associated
%              with 'rectan'.
%
%              'r', 'colat', and 'lon' return with the same vectorization
%              measure (N) as 'rectan'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example (1):
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
%      % Convert the array of position vectors 'pos' to spherical
%      % coordinates.
%      %
%      [r, colat, lon] = cspice_recsph(pos);
%
%      %
%      % Convert the spherical to rectangular.
%      %
%      [rectan] = cspice_sphrec(r, colat, lon);
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
%   Example (2):
%
%      %
%      % Define eleven sets of rectangular coordinates.
%      %
%      rec = [ [ 0., 1., 0., 0., -1., 0., 0., 1., 1., 0., 1. ]; ...
%              [ 0., 0., 1., 0., 0., -1., 0., 1., 0., 1., 1. ]; ...
%              [ 0., 0., 0., 1., 0., 0., -1., 0., 1., 1., 1. ]    ];
%
%      %
%      % ...convert the rectangular coordinates to spherical coordinates
%      %
%      [r, colat, lon] = cspice_recsph(rec);
%
%      %
%      % Convert 'colat' and 'lon to degrees...
%      %
%      colat = colat * cspice_dpr;
%      lon   = lon   * cspice_dpr;
%
%      %
%      % Output banner.
%      %
%      disp('     r        colat       lon          x         y           z   ')
%      disp('  --------   --------   --------   --------   --------   --------')
%
%      %
%      % Create an array of values for output.
%      %
%      output = [ r; colat; lon; rec(1,:); rec(2,:); rec(3,:) ];
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
%        r        colat       lon          x         y           z
%     --------   --------   --------   --------   --------   --------
%       0.0000     0.0000     0.0000     0.0000     0.0000     0.0000
%       1.0000    90.0000     0.0000     1.0000     0.0000     0.0000
%       1.0000    90.0000    90.0000     0.0000     1.0000     0.0000
%       1.0000     0.0000     0.0000     0.0000     0.0000     1.0000
%       1.0000    90.0000   180.0000    -1.0000     0.0000     0.0000
%       1.0000    90.0000   -90.0000     0.0000    -1.0000     0.0000
%       1.0000   180.0000     0.0000     0.0000     0.0000    -1.0000
%       1.4142    90.0000    45.0000     1.0000     1.0000     0.0000
%       1.4142    45.0000     0.0000     1.0000     0.0000     1.0000
%       1.4142    45.0000    90.0000     0.0000     1.0000     1.0000
%       1.7321    54.7356    45.0000     1.0000     1.0000     1.0000
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine recsph_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   rectangular to spherical coordinates
%
%-&

function [r, colat, lon] = cspice_recsph(rectan)

   switch nargin
      case 1
         rectan = zzmice_dp(rectan);
      otherwise
         error ( 'Usage: [_r_, _colat_, _lon_] = cspice_recsph(_rectan(3)_)' )
   end

   %
   % Call the MEX library.
   %
   try
      [r, colat, lon] = mice('recsph_c',rectan);
   catch
      rethrow(lasterror)
   end

