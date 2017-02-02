%-Abstract
%
%   CSPICE_calculates the transformation matrix to the
%   right-handed reference frame having an input vector as a
%   specified axis and a second input vector lying in a
%   define coordinate plane.
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
%      axdef    a double precision 3x1 array defining the principal
%               axes of a coordinate frame
%
%      indexa   the integer scalar signifying which of the three coordinate
%               axes contains 'axdef' (1, 2 or 3)
%
%                  If 'indexa' is 1 then axdef defines the X axis of the
%                  coordinate frame.
%
%                  If 'indexa' is 2 then axdef defines the Y axis of the
%                  coordinate frame.
%
%                  If 'indexa' is 3 then axdef defines the Z axis of the
%                  coordinate frame
%
%      plndef   a double precision 3x1 array defining a vector in the same
%               plane as 'axdef'
%
%      indexp   the integer scalar signifying the second principle axis,
%               orthogonal to 'axdef' (1, 2 or 3)
%
%                  If 'indexp' is 1, the second axis of the principal
%                  plane is the X-axis.
%
%                  If 'indexp' is 2, the second axis of the principal
%                  plane is the Y-axis.
%
%                  If 'indexp' is 3, the second axis of the principal plane
%                  is the Z-axis.
%
%   the call:
%
%      mout = cspice_twovec( axdef, indexa, plndef, indexp)
%
%   returns:
%
%      mout   a double precision 3x3 array defining a rotation matrix from
%             the frame of the original vectors to the new frame
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example:
%
%      %
%      % A trivial example. Define the reference vectors...
%      %
%      %  The i unit vector
%      %
%      axdef  = [ 1.; 0; 0.];
%      indexa = 1 ;
%
%      %
%      %  The -j unit vector. For this example, any vector
%      %  in the x-y plane linearly independent of 'axdef'
%      %  will suffice.
%      %
%      plndef = [ 0.; -1.; 0.];
%      indexp = 2;
%
%      %
%      % Calculate the transformation matrix. The new frame
%      % has 'axdef' as axis 'indexa', with 'plndef' in the same
%      % plane, the direction axis 'indexp' in that plane
%      % and orthogonal to 'axdef'. A third direction vector
%      % completes the right handed frame.
%      %
%      mout = cspice_twovec( axdef, indexa, plndef, indexp )
%
%   MATLAB outputs:
%
%      mout =
%
%           1     0     0
%           0    -1     0
%           0     0    -1
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine twovec_c.
%
%   MICE.REQ
%
%-Version
%
%    -Mice Version 1.0.0, 10-JAN-2006, EDW (JPL)
%
%-Index_Entries
%
%   define an orthonormal frame from two vectors
%
%-&

function [mout] = cspice_twovec(axdef, indexa, plndef, indexp)

   switch nargin
      case 4

         axdef  = zzmice_dp(axdef);
         indexa = zzmice_int(indexa);
         indexp = zzmice_int(indexp);
         plndef = zzmice_dp(plndef);

      otherwise

         error ( ['Usage: [mout(3,3)] = cspice_twovec( axdef(3), ' ...
                                       'indexa, plndef(3), indexp)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [mout] = mice('twovec_c', axdef, indexa, plndef, indexp);
   catch
      rethrow(lasterror)
   end


