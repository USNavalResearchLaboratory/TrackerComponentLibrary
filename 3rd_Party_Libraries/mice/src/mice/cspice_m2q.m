%-Abstract
%
%   CSPICE_M2Q calculates a unit quaternion corresponding to a
%   specified 3x3 double precision, rotation matrix.
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
%      r   a double precision 3x3 or 3x3xN array of rotation
%          matrices
%
%   the call:
%
%      q = cspice_m2q(r)
%
%   returns:
%
%      q   a double precision 4-vector or 4xN array, the quaternion
%          representation of the matrix 'r'
%
%          Note that multiple styles of quaternions are in use.
%          This routine returns a quaternion that conforms to
%          the SPICE convention. See the Particulars section
%          for details.
%
%          'q' returns with the same vectorization measure (N)
%          as 'r' .
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Create a rotation matrix of 90 degrees about the Z axis.
%      %
%      r = cspice_rotate( cspice_halfpi, 3)
%
%   MATLAB outputs:
%
%      r =
%
%         0.00000000000000   1.00000000000000                  0
%        -1.00000000000000   0.00000000000000                  0
%                        0                  0   1.00000000000000
%
%      q = cspice_m2q( r )
%
%   MATLAB outputs:
%
%      q =
%
%         0.70710678118655
%                        0
%                        0
%        -0.70710678118655
%
%      %            _
%      % Confirm || q || = 1.
%      %
%      q'  * q
%
%   MATLAB outputs:
%
%      ans =
%
%           1
%
%-Particulars
%
%   About SPICE quaternions
%   =======================
%
%   There are (at least) two popular "styles" of quaternions; these
%   differ in the layout of the quaternion elements, the definition
%   of the multiplication operation, and the mapping between the set
%   of unit quaternions and corresponding rotation matrices.
%
%   SPICE-style quaternions have the scalar part in the first
%   component and the vector part in the subsequent components. The
%   SPICE convention, along with the multiplication rules for SPICE
%   quaternions, are those used by William Rowan Hamilton, the
%   inventor of quaternions.
%
%   Another common quaternion style places the scalar component
%   last.  This style is often used in engineering applications.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine m2q_c.
%
%   MICE.REQ
%   ROTATION.REQ
%
%-Version
%
%    -Mice Version 1.0.0, 10-JAN-2006, EDW (JPL)
%
%-Index_Entries
%
%   matrix to quaternion
%
%-&

function [q] = cspice_m2q(r)

   switch nargin
      case 1

         r = zzmice_dp(r);

      otherwise

         error ( 'Usage: [_q(4)_] = cspice_m2q( _r(3,3)_ )' )

   end

   %
   % Call the MEX library.
   %
   try
      [q] = mice('m2q_c', r);
   catch
      rethrow(lasterror)
   end



