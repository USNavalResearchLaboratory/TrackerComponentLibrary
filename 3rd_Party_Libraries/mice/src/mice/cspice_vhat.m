%-Abstract
%
%   CSPICE_VHAT returns the unit vector along a double precision
%   3-dimensional vector.
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
%      v1   an arbitrary double precision, 3-dimensional vector
%           or 3xN array
%
%   the call:
%
%      vout = cspice_vhat(v1)
%
%   returns:
%
%      vout   contains the unit 3-vector or 3xN array of unit vectors
%             in the direction of 'v1'.
%
%                   ^       --
%                 vhat =    v1
%                        --------
%                           --
%                        || v1 ||
%
%                      _                                               _
%             where || x || indicates the Euclidean norm of the vector x.
%
%             If 'v1' represents the zero vector, then 'vout' will
%             also be the zero vector.
%
%             'vout' returns with the same vectorization measure as 'v1'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      >> v1 = [ 5; 12; 0]
%
%      v1 =
%
%           5
%          12
%           0
%
%      >> cspice_vhat(v1)
%
%   MATLAB outputs:
%
%      ans =
%
%          0.3846
%          0.9231
%               0
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
%      >> cspice_vhat(v2)
%
%   MATLAB outputs:
%
%      ans =
%
%          0.3333
%          0.6667
%          0.6667
%
%      >> v = [v1, v2 ]
%
%      v =
%
%          5.0000    0.0000
%         12.0000    0.0000
%               0    0.0000
%
%      >> cspice_vhat(v)
%
%   MATLAB outputs:
%
%      ans =
%
%          0.3846    0.3333
%          0.9231    0.6667
%               0    0.6667
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine vhat_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 30-DEC-2008, EDW (JPL)
%
%      Corrected misspellings.
%
%   -Mice Version 1.0.0, 25-APR-2006, EDW (JPL)
%
%-Index_Entries
%
%   unitize a 3-dimensional vector
%
%-&

function [vout] = cspice_vhat(v1)

   switch nargin
      case 1

         v1 = zzmice_dp(v1);

      otherwise

         error ( 'Usage: [_vout(3)_] = cspice_vhat(_v1(3)_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [vout] = mice('vhat_c',v1);
   catch
      rethrow(lasterror)
   end



