%-Abstract
%
%   CSPICE_SURFNM computes the double precision, outward-pointing
%   normal unit 3-vector at a point defined on the surface of an
%   ellipsoid.
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
%      point   a double precision 3-vector or a 3xN array of 3-vectors defining
%              some location(s) on the ellipsoid
%
%   the call:
%
%      normal = cspice_surfnm( a, b, c, point)
%
%   returns:
%
%      normal   a double precision unit 3-vector or 3XN array of unit vectors
%               normal to the ellipsoid at 'point' in the direction away
%               from the ellipsoid
%
%               'normal' returns with the same vectorization measure (N)
%                as 'point'.
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Define the radii of an ellipsoid.
%      %
%      a  =  1.;
%      b  =  2.;
%      c  =  3.;
%
%      %
%      % Select a set of locations, three 3-vectors.
%      %
%      point = [ [ 0.; 0.; 3.], [ 0.; 2.; 0.], [-1; 0; 0] ];
%
%      %
%      % Calculate the surface normal to the ellipsoid at 'point'.
%      %
%      out_norm = cspice_surfnm( a, b, c, point)
%
%   MATLAB outputs:
%
%      out_norm =
%
%           0     0    -1
%           0     1     0
%           1     0     0
%
%      Three 3-vectors:
%         the normal at (0,0,3) equals (0,0,1)
%         the normal at (0,2,0) equals (0,0,1)
%         the normal at (-1,0,0) equals (-1,0,0)
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine surfnm_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 15-JUN-2006, EDW (JPL)
%
%-Index_Entries
%
%   surface normal vector on an ellipsoid
%
%-&

function [normal] = cspice_surfnm(a, b, c, point)

   switch nargin
      case 4

         a     = zzmice_dp(a);
         b     = zzmice_dp(b);
         c     = zzmice_dp(c);
         point = zzmice_dp(point);

      otherwise

         error ( ['Usage: [_normal(3)_] = '             ...
                  'cspice_surfnm( a, b, c, _point(3)_ )'] )

   end

   %
   % Call the MEX library.
   %
   try
      [normal] = mice( 'surfnm_c',  a, b, c, point);
   catch
      rethrow(lasterror)
   end


