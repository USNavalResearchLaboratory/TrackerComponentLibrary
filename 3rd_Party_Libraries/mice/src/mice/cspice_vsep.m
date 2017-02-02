%-Abstract
%
%   CSPICE_VSEP returns the scalar double precision
%   separation angle in radians between two double precision,
%   3-vectors.  This angle is defined as zero if either vector
%   is zero.
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
%      v1   is an arbitrary double precision, 3-dimensional vector
%           or 3xN array
%
%      v2   is also an arbitrary double precision, 3-dimensional
%           vector or 3xN array
%
%              'v1' or 'v2' or both may be the zero vector.
%
%      An implicit assumption exists that 's1' and 's2' are specified
%      in the same reference frame. If this is not the case, the numerical
%      result has no meaning.
%
%   the call:
%
%      vsep = cspice_vsep(v1, v2)
%
%   returns:
%
%      vsep   the double precision, positive definite, scalar or
%             1xN array of the angular separation(s) between 'v1' and 'v2'
%             expressed in radians. If either 'v1' or 'v2' is the zero
%             vector, then cspice_vsep return value has value 0 radians.
%
%             'sep' returns with the same measure of vectorization (N) as
%             'v1' and 'v2'
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Define a set of vectors, calculate angular
%      % separation as measured in radians.
%      %
%      v1 = [1; 0; 0]
%      v2 = [1; 0; 0]
%
%      sep = cspice_vsep( v1, v2 )
%
%   MATLAB outputs:
%
%      v1 =
%
%           1
%           0
%           0
%
%      v2 =
%
%           1
%           0
%           0
%
%      sep =
%
%           0
%
%      %
%      % Another vector set.
%      %
%
%      v1 = [1; 0; 0]
%      v2 = [0; 1; 0]
%
%      sep = cspice_vsep( v1, v2 )
%
%   MATLAB outputs:
%
%      v1 =
%
%           1
%           0
%           0
%
%      v2 =
%
%           0
%           1
%           0
%
%      sep =
%
%         1.57079632679490
%
%      %
%      % Instead of two calls with 3-vectors,
%      % vectorize the input as two 3X2 array.
%      %
%      v1 = [ [1; 0; 0], [1; 0; 0] ]
%      v2 = [ [1; 0; 0], [0; 1; 0] ]
%
%      sep = cspice_vsep( v1, v2 )
%
%   MATLAB outputs:
%
%      v1 =
%
%           1     1
%           0     0
%           0     0
%
%      v2 =
%
%           1     0
%           0     1
%           0     0
%
%      sep =
%
%         0   1.57079632679490
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine vsep_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 10-APR-2010, EDW (JPL)
%
%      Edits to header I/O section.
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   angular separation of 3-dimensional vectors
%
%-&

function [sep] = cspice_vsep(v1, v2)

   switch nargin
      case 2

         v1 = zzmice_dp(v1);
         v2 = zzmice_dp(v2);

      otherwise

         error ( 'Usage: [_sep_] = cspice_vsep(_v1(3)_, _v2(3)_)' )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      [sep] = mice('vsep_c',v1, v2);
   catch
      rethrow(lasterror)
   end



