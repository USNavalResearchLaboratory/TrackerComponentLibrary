%-Abstract
%
%   CSPICE_RECCYL converts rectangular (Cartesian) coordinates to
%   cylindrical coordinates.
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
%       [r, lonc, z] = cspice_reccyl( rectan)
%
%   returns:
%
%      r      a double precision scalar or double precision 1xN array
%             describing the distance of the point of interest from z axis
%
%      lonc   a double precision scalar or double precision 1xN array
%             describing the cylindrical angle of the point of interest
%             from the XZ plane measured in radians
%
%      z      a double precision scalar or double precision 1xN array
%             describing the height of the point above the XY plane
%
%             The arguments 'r' and 'z' return in the same units associated
%             with 'rectan'.
%
%             'r', 'lonc', and 'z' return with the same vectorization
%             measure as 'rectan'.
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
%      % Convert the array of position vectors 'pos' to cylindrical
%      % coordinates.
%      %
%      [r, lonc, z] = cspice_reccyl(pos);
%
%      %
%      % Convert the cylindrical to rectangular.
%      %
%      [rectan] = cspice_cylrec(r, lonc, z);
%
%      %
%      % Calculate the relative error against the original position
%      % vectors.
%      %
%      (rectan-pos) ./ pos
%
%   MATLAB outputs:
%
%      1.0e-13 *
%
%      0.00199609007208                  0  -0.25513381329527
%     -0.00218237675815                  0  -0.00153127196389
%                     0                  0                  0
%
%       The relative error between the original array of position vectors
%       and those that resulted from the various coordinate conversion
%       has magnitude on the order of 10^(-13).  A numerical
%       demonstration of equality.
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
%      % ...convert the rectangular coordinates to cylindrical coordinates
%      %
%      [r, lonc, z] = cspice_reccyl(rec);
%
%      %
%      % Convert 'lonc' to degrees...
%      %
%      lonc = lonc * cspice_dpr;
%
%      %
%      % Output banner.
%      %
%      disp('     r         lonc        z           x         y           z   ')
%      disp('  --------   --------   --------   --------   --------   --------')
%
%      %
%      % Create an array of values for output.
%      %
%      output = [ r; lonc; z; rec(1,:); rec(2,:); rec(3,:) ];
%      txt    = sprintf( '%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', output);
%      disp( txt );
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%
%   MATLAB outputs:
%
%        r         lonc        z           x         y           z
%     --------   --------   --------   --------   --------   --------
%       0.0000     0.0000     0.0000     0.0000     0.0000     0.0000
%       1.0000     0.0000     0.0000     1.0000     0.0000     0.0000
%       1.0000    90.0000     0.0000     0.0000     1.0000     0.0000
%       0.0000     0.0000     1.0000     0.0000     0.0000     1.0000
%       1.0000   180.0000     0.0000    -1.0000     0.0000     0.0000
%       1.0000   270.0000     0.0000     0.0000    -1.0000     0.0000
%       0.0000     0.0000    -1.0000     0.0000     0.0000    -1.0000
%       1.4142    45.0000     0.0000     1.0000     1.0000     0.0000
%       1.0000     0.0000     1.0000     1.0000     0.0000     1.0000
%       1.0000    90.0000     1.0000     0.0000     1.0000     1.0000
%       1.4142    45.0000     1.0000     1.0000     1.0000     1.0000
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine reccyl_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   rectangular to cylindrical coordinates
%
%-&

function [r, lonc, z] = cspice_reccyl(rectan)

   switch nargin
      case 1

         rectan = zzmice_dp(rectan);

      otherwise
         error ( 'Usage: [_r_, _lonc_, _z_] = cspice_reccyl(_rectan(3)_)' )
   end

   %
   % Call the MEX library.
   %
   try
      [r, lonc, z] = mice('reccyl_c',rectan);
   catch
      rethrow(lasterror)
   end



