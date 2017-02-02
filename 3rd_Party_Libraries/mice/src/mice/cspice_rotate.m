%-Abstract
%
%   CSPICE_ROTATE calculates the 3x3 rotation matrix generated
%   by a rotation of a specified angle about a specified axis.
%   This rotation operates as a rotation of the coordinate
%   system.
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
%      angle  the double precision scalar or 1XN-vector of rotation
%             angles measured in radians
%
%      iaxis   the integer ID of the axis of rotation where
%              X=1, Y=2, Z=3
%
%   the call:
%
%      mout = cspice_rotate( angle, iaxis)
%
%   returns:
%
%      mout   a double precision 3x3 or 3x3xN array of rotation matrices
%             that describe a rotation of 'angle' radians about 'iaxis'
%
%             'mout' return with the same vectorization measure
%             (N) as 'angle'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % A Pi/10 rotation about the Z axis.
%      %
%      rot_mat = cspice_rotate( 0.1*cspice_pi, 3 )
%
%      %
%      % Apply the coordinate rotation to a vector.
%      %
%      vec = [ 1.2; 3.4; 4.5 ];
%
%      vec1 = rot_mat * vec
%
%   MATLAB outputs:
%
%      rot_mat =
%
%          0.9511    0.3090         0
%         -0.3090    0.9511         0
%               0         0    1.0000
%
%      vec1 =
%
%          2.1919
%          2.8628
%          4.5000
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine rotate_c.
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
%   generate a rotation matrix
%
%-&

function [mout] = cspice_rotate( angle, iaxis )

   switch nargin

      case 2

         angle = zzmice_dp(angle);
         iaxis = zzmice_int( iaxis );

      otherwise
         error ( 'Usage: [_mout(3,3)_] = cspice_rotate( _angle_, iaxis )' )

   end

   %
   % Call the MEX library.
   %
   try
      [mout] = mice('rotate_c', angle, iaxis );
   catch
      rethrow(lasterror)
   end


