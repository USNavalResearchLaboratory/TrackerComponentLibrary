%-Abstract
%
%   CSPICE_AXISAR returns a 3x3 double rotation matrix that rotates
%   vectors by a specified angle about a specified axis.
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
%      axis    an arbitrary, non-zero, double precision 3x1 array
%              defining a rotation axis
%
%      angle   the double precision scalar angle in radians
%              defining the measure of rotation about 'axis'
%
%   the call:
%
%      r = cspice_axisar( axis, angle)
%
%   returns:
%
%      r   a double precision 3x3 array representing the coordinate
%          transformation determined by 'axis' and 'angle', i.e. the
%          application of 'r' to a 3x1 array returns the result
%          of rotating the vector about 'axis' through 'angle' radians
%
%   Please note cspice_raxisa is not guaranteed to invert the
%   operation of cspice_axisar.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      Example:
%
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
%         1.0e-15 *
%
%                        0  -0.11102230246252   0.05551115123126
%         0.11102230246252                  0                  0
%        -0.05551115123126   0.02775557561563                  0
%
%     The zero matrix accurate to round-off error. A numerical
%     demonstration of equality.
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine axisar_c.
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
%   axis and angle to rotation
%
%-&

function [r] = cspice_axisar( axis, angle)

   switch nargin
      case 2

         axis  = zzmice_dp(axis);
         angle = zzmice_dp(angle);

      otherwise

         error ( 'Usage: [r(3,3)] = cspice_axisar( axis(3), angle)' )

   end

   %
   % Call the MEX library.
   %
   try
      [r] = mice( 'axisar_c', axis, angle );
   catch
      rethrow(lasterror)
   end


