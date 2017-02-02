%-Abstract
%
%   CSPICE_NEARPT calculates the point on the surface of an
%   ellipsoid nearest to a specified off-ellipsoid position.
%   The routine also returns the altitude of the position
%   above the ellipsoid
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
%      positn   a double precision position 3-vector or 3XN array defining the
%               position of a point with respect to the center of an
%               ellipsoid. The vector is expressed in a body-fixed reference
%               frame. The semi-axes of the ellipsoid are aligned with the
%               x, y, and z-axes of the body-fixed frame.
%
%      a,       the scalar double precision values of the ellipsoid's
%      b,       triaxial radii ellipsoid, where:
%      c
%                  'a' is length in kilometers of the semi-axis of the ellipsoid
%                   parallel to the x-axis of the body-fixed reference frame
%
%                  'b' is length in kilometers of the semi-axis of the ellipsoid
%                   parallel to the y-axis of the body-fixed reference frame
%
%                  'c' is length in kilometers of the semi-axis of the ellipsoid
%                   parallel to the z-axis of the body-fixed reference frame
%
%   the call:
%
%      [ npoint, alt ] = cspice_nearpt( positn, a, b, c )
%
%   returns:
%
%      npoint   the double precision 3-vector or 3xN array defining the
%               location on the ellipsoid closest to 'positn'
%
%      alt      the double precision scalar or 1xN array of altitudes of
%              'positn' above the ellipsoid.  If 'positn' is inside the
%               ellipsoid, 'alt' will be negative and have magnitude
%               equal to the distance between 'pos' and 'positn'.
%
%               'npoint' and 'alt' return with the same vectorization
%               measure (N) as 'positn'
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      Example(1):
%
%      %
%      % Define the radii of an ellipsoid.
%      %
%      a  =  1.;
%      b  =  2.;
%      c  =  3.;
%
%      %
%      % Use point on the X axis, outside the ellipsoid.
%      %
%      point        = [ 3.5; 0.; 0. ];
%      [pnear, alt] = cspice_nearpt( point, a, b, c)
%
%   MATLAB outputs:
%
%      pnear =
%
%           1
%           0
%           0
%
%      alt =
%
%         2.50000000000000
%
%      Example(2):
%
%      %
%      % Load a meta kernel containing SPK and leapseconds kernels.
%      %
%      cspice_furnsh( 'standard.tm')
%
%      %
%      % Retrieve the position of the Moon wrt the Earth at
%      % ephemeris time 0.d (Jan 1 2000 12:00 TDB) in the Earth-fixed
%      % reference frame.
%      %
%      epoch       = 0.;
%      abcorr      = 'LT+S';
%      [pos,ltime] = cspice_spkpos( 'moon', epoch, 'IAU_EARTH', ...
%                                    abcorr, 'earth');
%
%      %
%      % Retrieve the triaxial radii for Earth (body ID 399).
%      %
%      radii = cspice_bodvrd( 'EARTH', 'RADII', 3);
%
%      %
%      % Now calculate the point on the Earth nearest to the Moon
%      % given LT+S aberration correction at the epoch time.
%      %
%      [npoint,alt]= cspice_nearpt( pos, radii(1), radii(2), radii(3) )
%
%
%   MATLAB outputs:
%
%      npoint =
%
%         1.0e+03 *
%
%         3.34708386495926
%        -5.29453888129091
%        -1.19828126398311
%
%      alt =
%
%         3.960372197033597e+05
%
%-Particulars
%
%   A sister version of this routine exists named mice_nearpt that returns
%   the output arguments as fields in a single structure.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine nearpt_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 21-DEC-2005, EDW (JPL)
%
%-Index_Entries
%
%   distance from point to ellipsoid
%   nearest point on an ellipsoid
%
%-&

function [ npoint, alt ] = cspice_nearpt( positn, a, b, c )

   switch nargin
      case 4

         positn = zzmice_dp(positn);
         a      = zzmice_dp(a);
         b      = zzmice_dp(b);
         c      = zzmice_dp(c);

      otherwise

         error ( ['Usage: [_npoint(3)_, _alt_] = ' ...
                  'cspice_nearpt( _positn(3)_, a, b, c )'] )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [nearpt] = mice( 'nearpt_s', positn, a, b, c  );
      npoint   = reshape( [nearpt.pos], 3, [] );
      alt      = reshape( [nearpt.alt], 1, [] );
   catch
      rethrow(lasterror)
   end


