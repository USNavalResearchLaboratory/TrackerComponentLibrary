%-Abstract
%
%   CSPICE_RAXISA computes the axis of the rotation given by an input matrix
%   and the angle of the rotation about that axis.
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
%      matrix   a double precision 3x3 array defining a rotation
%
%   the call:
%
%      [axis, angle] = cspice_raxisa( matrix )
%
%   returns:
%
%      axis   a double precision 3x1 unit array pointing along the axis
%             of the rotation. In other words, 'axis' is a unit eigenvector
%             of the input matrix, corresponding to the eigenvalue 1. If
%             the input matrix is the identity matrix, 'axis' will be the
%             vector (0, 0, 1). If the input rotation is a rotation by pi
%             radians, both 'axis' and -axis may be regarded as the axis
%             of the rotation.
%
%      angle  a double precision scalar defining the angle
%             between 'v' and matrix*'v' for any non-zero vector 'v'
%             orthogonal to 'axis'.  'angle' is given in radians.
%             The angle returned will be in the range from 0 to
%             pi radians.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      Example:
%      %
%      % Define an axis and an angle for rotation.
%      %
%      axis = [ 1.; 2.; 3. ];
%      angle = .1 * cspice_twopi;
%
%      %
%      % Determine the rotation matrix.
%      %
%      rot_mat = cspice_axisar( axis, angle );
%
%      %
%      % Now calculate the rotation axis and angle based on the
%      % matrix as input.
%      %
%      [ axout, angout ] = cspice_raxisa( rot_mat);
%
%      %
%      % Now input the axout and angout to cspice_axisar to
%      % compare against the original rotation matrix rot_mat.
%      %
%      rot_out = cspice_axisar( axout, angout );
%      rot_mat - rot_out
%
%   MATLAB outputs:
%
%      1.0e-15 *
%
%                     0  -0.11102230246252   0.05551115123126
%      0.11102230246252                  0                  0
%     -0.05551115123126   0.02775557561563                  0
%
%   The zero matrix accurate to round-off error. A numerical
%   demonstration of equality.
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine raxisa_c.
%
%   MICE.REQ
%   ROTATION.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 29-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   rotation axis of a matrix
%
%-&

function [axis, angle] = cspice_raxisa(matrix)

   switch nargin
      case 1

         matrix = zzmice_dp(matrix);

      otherwise

         error ( [ 'Usage: [axis(3), angle] = ' ...
                   'cspice_raxisa(matrix(3,3))' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [axis, angle] = mice('raxisa_c', matrix);
   catch
      rethrow(lasterror)
   end


