%-Abstract
%
%   CSPICE_VROTV rotates a double precision 3-vector about a specified
%   axis vector by a specified angle (measured in radians) then
%   returns the rotated vector.
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
%      v       a double precision 3x1 array to rotate
%
%      axis    a double precision 3x1 array defining the axis about which
%              to rotate 'v'
%
%      theta   a double precision scalar angle measured in radians through which
%              which rotate 'v' about 'axis'
%
%   the call:
%
%      r = cspice_vrotv( v, axis, theta )
%
%   returns:
%
%      r   a double precision 3x1 array, the result of rotating 'v' about
%          'axis' through an angle of 'theta'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Given an axis of rotation and angle of rotation.
%      %
%      axis  = [ 0.; 0.; 1.];
%      theta = cspice_halfpi;
%
%      %
%      % Perform rotations on various vectors...
%      %
%
%   Example(1):
%
%      v1 = [ 1.; 2.; 3. ];
%
%      r1 = cspice_vrotv( v1, axis, theta )
%
%   MATLAB outputs:
%
%      r1 =
%
%         -2.0000
%          1.0000
%          3.0000
%
%   Example(2):
%
%      v2 = [ 1.; 0.; 0. ];
%
%      r2 = cspice_vrotv( v2, axis, theta )
%
%   MATLAB outputs:
%
%      r2 =
%
%          0.0000
%          1.0000
%               0
%
%   Example(3):
%
%      v3 = [ 0.; 1.; 0. ];
%
%      r3 = cspice_vrotv( v3, axis, theta )
%
%   MATLAB outputs:
%
%      r3 =
%
%         -1.0000
%          0.0000
%               0
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine vrotv_c.
%
%   MICE.REQ
%   ROTATION.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 14-JUL-2010, EDW (JPL)
%
%      Corrected minor typo in header.
%
%   -Mice Version 1.0.0, 17-APR-2008, EDW (JPL)
%
%-Index_Entries
%
%   vector rotation about an axis
%
%-&

function [r] = cspice_vrotv( v, axis, theta)

   switch nargin
      case 3

         v     = zzmice_dp(v);
         axis  = zzmice_dp(axis);
         theta = zzmice_dp(theta);

      otherwise

         error ( ['Usage: [r(3)] = cspice_vrotv( v(3), axis(3), theta)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [r] = mice('vrotv_c', v, axis, theta);
   catch
      rethrow(lasterror)
   end


