%-Abstract
%
%   CSPICE_SPHCYL converts spherical coordinates to cylindrical
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
%      radius   a double precision scalar or double precision 1XN
%               array describing the distance of the position from origin
%
%      colat    a double precision scalar or double precision 1XN array
%               describing the angle of the projection of the point to the XY
%               positive z-axis, measured in radians (also referred to
%               as the polar angle)
%
%      slon     a double precision scalar or 1XN array describing
%               the angle of the projection of the point to the XY
%               plane from the positive X-axis, measured in radians,
%               with range:
%
%                   -pi < slon <= pi
%
%               The positive Y-axis is at longitude PI/2 radians.
%
%   the call:
%
%      [ r, lonc, z] = cspice_sphcyl( radius, colat, slon)
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
%             with 'radius'.
%
%             'r', 'lonc', and 'z' return with the same vectorization
%             measure(N) as 'radius', 'colat', and 'slon'.
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
%      [radius, colat, lon] = cspice_recsph(pos);
%
%      %
%      % Convert the spherical coords to cylindrical.
%      %
%      [r, lon, z] = cspice_sphcyl(radius, colat, lon);
%
%      %
%      % Convert the cylindrical to rectangular.
%      %
%      [rectan] = cspice_cylrec(r, lon, z);
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
%   Example (2):
%
%      %
%      % Define six sets of spherical coordinates, 'slon' and 'colat'
%      % expressed in degrees - converted to radians by use of cspice_rpd.
%      %
%      radius = [  1.,  1., 1.4142, 1.4142, 1.  , 0. ];
%      colat  = [ 90., 90., 45.   , 135.  , 0.  , 0. ] * cspice_rpd;
%      slon   = [  0., 90., 180.  , 180.  , 180., 33.] * cspice_rpd;
%
%      %
%      % ...convert the spherical coordinates to cylindrical coordinates
%      %
%      [r, lonc, z] = cspice_sphcyl(radius, colat, slon);
%
%     %
%      % ...convert angular measure to degrees.
%      %
%      colat = colat * cspice_dpr;
%      lonc = lonc   * cspice_dpr;
%      slon = slon   * cspice_dpr;
%
%      %
%      % Output banner.
%      %
%      disp('     r         lonc        z        radius     slon       colat  ')
%      disp('  --------   --------   --------   --------   --------   --------')
%
%      %
%      % Create an array of values for output.
%      %
%      output = [ r; lonc; z; radius; slon; colat ];
%      txt   = sprintf( '%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', output);
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
%        r         lonc        z        radius      lon       colat
%     --------   --------   --------   --------   --------   --------
%       1.0000     0.0000     0.0000     1.0000     0.0000    90.0000
%       1.0000    90.0000     0.0000     1.0000    90.0000    90.0000
%       1.0000   180.0000     1.0000     1.4142   180.0000    45.0000
%       1.0000   180.0000    -1.0000     1.4142   180.0000   135.0000
%       0.0000   180.0000     1.0000     1.0000   180.0000     0.0000
%       0.0000    33.0000     0.0000     0.0000    33.0000     0.0000
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine sphcyl_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 12-DEC-2005, EDW (JPL)
%
%-Index_Entries
%
%   spherical to cylindrical coordinates
%
%-&

function [ r, lonc, z] = cspice_sphcyl(radius, colat, slon)

   switch nargin
      case 3

         radius = zzmice_dp(radius);
         colat  = zzmice_dp(colat);
         slon   = zzmice_dp(slon);

      otherwise

         error ( ['Usage: [ _r_, _lonc_, _z_] = '...
                  'cspice_sphcyl(_radius_, _colat_, _slon_)' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [ r, lonc, z] = mice('sphcyl_c', radius, colat, slon );
   catch
      rethrow(lasterror)
   end


