%-Abstract
%
%   CSPICE_UNORM normalizes a double precision 3-vector and
%   returns its magnitude.
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
%      v1   any arbitrary 3-vector(s).
%
%           [3,n] = size(v1); double = class(v1)
%
%   the call:
%
%      [vout, vmag] = cspice_unorm(v1)
%
%   returns:
%
%      vout   unit vector(s) in the direction of 'v1'. If 'v1'
%             represents the zero vector, then 'vout' will also be
%             the zero vector.
%
%             vout =   v1
%                    ------
%                    ||v1||
%
%             [3,n] = size(vout); double = class(vout)
%
%      vmag   the positive definite magnitude(s) of 'v1', ||v1||, calculated
%              in a numerically stable way.
%
%             [1,n] = size(vmag); double = class(vmag)
%
%             'vout' and 'vmag' return with the same measure of vectorization
%              (N) as 'v1'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example(1):
%
%      >> v1 = [ 5; 12; 0]
%
%      v1 =
%
%           5
%          12
%           0
%
%      >> [vout, vmag] = cspice_unorm(v1 )
%
%   MATLAB outputs:
%
%      vout =
%
%          0.3846
%          0.9231
%               0
%
%      vmag =
%
%          13
%
%   Example(2):
%
%      >> v2 = [ 1D-7; 2D-7; 2D-7]
%
%      v2 =
%
%         1.0e-06 *
%
%          0.1000
%          0.2000
%          0.2000
%
%      >> [vout, vmag] = cspice_unorm(v2)
%
%   MATLAB outputs:
%
%      vout =
%
%          0.3333
%          0.6667
%          0.6667
%
%      vmag =
%
%         3.0000e-07
%
%   Example(3):
%
%      >> v = [v1, v2 ]
%
%      v =
%
%          5.0000    0.0000
%         12.0000    0.0000
%               0    0.0000
%
%      >> [vout, vmag] = cspice_unorm(v )
%
%   MATLAB outputs:
%
%      vout =
%
%          0.3846    0.3333
%          0.9231    0.6667
%               0    0.6667
%
%      vmag =
%
%         13.0000    0.0000
%
%   The second element of 'vmag' displays as 0.0 due to the eight
%   orders of magnitude difference between that element and the first.
%   Confirm the expected value:
%
%      >> vmag(2)
%
%   MATLAB outputs:
%
%      ans =
%
%         3.0000e-07
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine unorm_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.2, 09-NOV-2012, EDW (JPL)
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.1, 21-APR-2010, EDW (JPL)
%
%      Corrected description of call example to list all output arguments,
%      from:
%
%         vout = cspice_unorm(v1)
%
%      to
%
%         [vout, vmag] = cspice_unorm(v1)
%
%   -Mice Version 1.0.0, 25-APR-2006, EDW (JPL)
%
%-Index_Entries
%
%   3-dimensional unit vector and norm
%
%-&

function [vout, vmag] = cspice_unorm(v1)

   switch nargin
      case 1

         v1 = zzmice_dp(v1);

      otherwise

         error ( 'Usage: [_vout(3)_, _vmag_] = cspice_unorm(_v1(3)_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [vout, vmag] = mice('unorm_c',v1);
   catch
      rethrow(lasterror)
   end



