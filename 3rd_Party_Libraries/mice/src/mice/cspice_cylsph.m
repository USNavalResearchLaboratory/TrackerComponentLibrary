%-Abstract
%
%   CSPICE_CYLSPH converts cylindrical coordinates to spherical
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
%      r      a double precision scalar or double precision 1xN array
%             describing the distance of the point of interest from z axis.
%
%      lonc   a double precision scalar or double precision 1xN array
%             describing the cylindrical angle of the point of interest
%             from the XZ plane measured in radians
%
%      z      a double precision scalar or double precision 1xN array
%             describing the height of the point above the XY plane
%
%   the call:
%
%      [radius, colat, lon] = cspice_cylsph( r, lonc, z)
%
%   returns:
%
%      radius   a double precision scalar or double precision 1XN
%               array describing the distance of the position from origin
%
%      colat    a double precision scalar or double precision 1XN
%               array describing the angle between the point and the
%               positive z-axis, measured in radians (also referred to
%               as the polar angle)
%
%      lon      a double precision scalar or double precision 1XN array
%               describing the angle of the projection of the point to the XY
%               plane from the positive X-axis, measured in radians,
%               with range:
%
%                   -pi < lon <= pi
%
%               The positive Y-axis is at longitude PI/2 radians.
%
%               The argument 'radius' returns in the same units associated
%               with 'r' and 'z'.
%
%               'radius', 'colat', and 'lon' return with the vectorization
%               measure (N) as 'r', 'lonc', and 'z'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example (1):
%
%         %
%         % Load an SPK, leapseconds, and PCK kernel set.
%         %
%         cspice_furnsh( 'standard.tm' )
%
%         %
%         % Create a vector of scalar times.
%         %
%         et = [0:2]*2.*cspice_spd;
%
%         %
%         % Retrieve the position of the moon seen from earth at 'et'
%         % in the J2000 frame without aberration correction.
%         %
%         [pos, et] = cspice_spkpos( 'MOON', et, 'J2000', 'NONE', 'EARTH' );
%
%         %
%         % Convert the array of position vectors 'pos' to cylindrical
%         % coordinates.
%         %
%         [r, lon, z] = cspice_reccyl(pos);
%
%         %
%         % Convert the cylindrical coords to spherical.
%         %
%         [radius, colat, lon] = cspice_cylsph(r, lon, z);
%
%         %
%         % Convert the spherical to rectangular.
%         %
%         [rectan] = cspice_sphrec(radius, colat, lon);
%
%         %
%         % Calculate the relative error against the original position
%         % vectors.
%         %
%         (rectan-pos) ./ pos
%
%   MATLAB outputs:
%
%      1.0e-13 *
%
%      0.00199609007208   0.00185077353361  -0.25349834013312
%     -0.00436475351630                  0  -0.00153127196389
%     -0.00382429455002  -0.00485521702846   0.00203951344664
%
%      The relative error between the original array of position vectors
%      and those resulting from the coordinate conversions
%      has magnitude on the order of 10^(-13).
%
%   Example (2):
%
%      %
%      % Define six sets of cylindrical coordinates, 'lonc' expressed
%      % in degrees - converted to radians by use of cspice_rpd.
%      %
%      r     = [ 1.,  1.,   1.,   1.,   0.,  0. ];
%      lonc  = [ 0., 90., 180., 180., 180., 33. ] * cspice_rpd;
%      z     = [ 0.,  0.,   1.,  -1.,   1.,  0. ];
%
%      %
%      % ...convert the cylindrical coordinates to spherical coordinates
%      %
%      [radius, colat, lon] = cspice_cylsph(r, lonc, z);
%
%      %
%      % ...convert angular measure to degrees.
%      %
%      colat = colat * cspice_dpr;
%      lon  = lon    * cspice_dpr;
%      lonc = lonc   * cspice_dpr;
%
%      %
%      % Output banner.
%      %
%      disp('     r         lonc        z        radius      lon       colat  ')
%      disp('  --------   --------   --------   --------   --------   --------')
%
%      %
%      % Create an array of values for output.
%      %
%      output = [ r; lonc; z; radius; lon; colat ];
%
%      txt = sprintf( '%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', output );
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
%   the CSPICE routine cylsph_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 01-DEC-2005, EDW (JPL)
%
%-Index_Entries
%
%   cylindrical to spherical
%
%-&

function [radius, colat, lon] = cspice_cylsph(r, lonc, z)

   switch nargin
      case 3

         r    = zzmice_dp(r);
         lonc = zzmice_dp(lonc);
         z    = zzmice_dp(z);

      otherwise

         error ( ['Usage: [_radius_, _colat_, _lon_] = '...
                  'cspice_cylsph(_r_, _lonc_, _z_)' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [radius, colat, lon] = mice('cylsph_c', r, lonc, z);
   catch
      rethrow(lasterror)
   end


