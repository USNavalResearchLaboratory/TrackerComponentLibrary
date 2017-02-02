%-Abstract
%
%   CSPICE_RECLAT converts rectangular (Cartesian) coordinates to
%   latitudinal coordinates. All coordinates are expressed as
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
%      [radius, longitude, latitude] = cspice_reclat(rectan)
%
%   returns:
%
%      radius      a double precision scalar or double precision 1xN array
%                  describing the distance of the position from origin
%
%      longitude   a double precision scalar or double precision 1xN array
%                  describing the angle of the position from the XZ plane
%                  measured in radians
%
%      latitude    a double precision scalar or double precision 1xN array
%                  describing the angle of the position from the XY plane
%                  measured in radians
%
%                  The argument 'radius' returns in the same units associated
%                  with 'rectan'.
%
%                  'radius', 'longitude', and 'latitude' return with
%                  the same vectorization measure (N) as 'rectan'.
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
%      % Convert the array of position vectors 'pos' to latitudinal
%      % coordinates.
%      %
%      [radius, longitude, latitude] = cspice_reclat(pos);
%
%      %
%      % Convert the cylindrical to rectangular.
%      %
%      [rectan] = cspice_latrec(radius, longitude, latitude);
%
%      %
%      % Calculate the relative error against the original position
%      % vectors.
%      %
%      (rectan-pos) ./ pos
%
%    MATLAB outputs:
%
%      1.0e-14 *
%
%     -0.01996090072080  -0.05552320600838   0.63783453323816
%      0.02182376758148                  0  -0.01531271963894
%      0.01912147275010   0.01213804257114   0.02039513446643
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
%      % ...convert the rectangular coordinates to latitudinal coordinates
%      %
%      [radius, longitude, latitude] = cspice_reclat(rec);
%
%      %
%      % Convert 'longitude' and 'latitude' to degrees.
%      %
%      longitude = longitude * cspice_dpr;
%      latitude  = latitude  * cspice_dpr;
%
%      %
%      % Output banner.
%      %
%      disp('   radius    longitude  latitude       x         y           z   ')
%      disp('  --------   --------   --------   --------   --------   --------')
%
%      %
%      % Create an array of values for output.
%      %
%      output = [radius; longitude; latitude; rec(1,:); rec(2,:); rec(3,:) ];
%      txt    = sprintf( '%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', ...
%                        output );
%
%      disp( txt );
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%    MATLAB outputs:
%
%      radius    longitude  latitude       x         y           z
%     --------   --------   --------   --------   --------   --------
%       0.0000     0.0000     0.0000     0.0000     0.0000     0.0000
%       1.0000     0.0000     0.0000     1.0000     0.0000     0.0000
%       1.0000    90.0000     0.0000     0.0000     1.0000     0.0000
%       1.0000     0.0000    90.0000     0.0000     0.0000     1.0000
%       1.0000   180.0000     0.0000    -1.0000     0.0000     0.0000
%       1.0000   -90.0000     0.0000     0.0000    -1.0000     0.0000
%       1.0000     0.0000   -90.0000     0.0000     0.0000    -1.0000
%       1.4142    45.0000     0.0000     1.0000     1.0000     0.0000
%       1.4142     0.0000    45.0000     1.0000     0.0000     1.0000
%       1.4142    90.0000    45.0000     0.0000     1.0000     1.0000
%       1.7321    45.0000    35.2644     1.0000     1.0000     1.0000
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine reclat_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   rectangular to latitudinal coordinates
%
%-&

function [radius, longitude, latitude] = cspice_reclat(rectan)

   switch nargin
      case 1

         rectan = zzmice_dp(rectan);

      otherwise
         error ( ['Usage: [_radius_, _longitude_, _latitude_] = ' ...
                  'cspice_reclat(_rectan(3)_)'] )
   end

   %
   % Call the MEX library.
   %
   try
      [radius, longitude, latitude] = mice('reclat_c',rectan);
   catch
      rethrow(lasterror)
   end

