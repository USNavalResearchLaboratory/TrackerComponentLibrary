%-Abstract
%
%   CSPICE_INVORT returns the inverse of a 3x3 matrix with orthogonal
%   columns and non-zero norms using a numerical stable algorithm.
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
%      m     a 3x3 matrix.
%
%            [3,3] = size(m); double = class(m)
%
%   the call:
%
%      mit = cspice_invort(m)
%
%   returns:
%
%      mit   the matrix obtained by transposing 'm' and dividing the rows by
%            squares of their norms.
%
%            [3,3] = size(mit); double = class(mit)
%
%-Examples
%
%   Suppose that you have a matrix 'm' whose columns are orthogonal
%   and have non-zero norm (but not necessarily norm 1).  Then the
%   routine cspice_invort can be used to construct the inverse of 'm':
%
%        invers = cspice_invort( m )
%
%-Particulars
%
%   Suppose that m is the matrix
%
%           -                      -
%          |   A*u    B*v     C*w   |
%          |      1      1       1  |
%          |                        |
%          |   A*u    B*v     C*w   |
%          |      2      2       2  |
%          |                        |
%          |   A*u    B*v     C*w   |
%          |      3      3       3  |
%           -                      -
%
%   where the vectors (u , u , u ),  (v , v , v ),  and (w , w , w )
%                       1   2   3      1   2   3          1   2   3
%
%   are unit vectors. This routine produces the matrix:
%
%
%           -                      -
%          |   a*u    a*u     a*u   |
%          |      1      2       3  |
%          |                        |
%          |   b*v    b*v     b*v   |
%          |      1      2       3  |
%          |                        |
%          |   c*w    c*w     c*w   |
%          |      1      2       3  |
%           -                      -
%
%   where a = 1/A, b = 1/B, and c = 1/C.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine invort_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 14-NOV-2013, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   Transpose a matrix and invert the lengths of the rows
%   Invert a pseudo orthogonal matrix
%
%-&

function [mit] = cspice_invort( m )

   switch nargin
      case 1

         m = zzmice_dp(m);

      otherwise

         error( 'Usage: [mit(3,3)] = cspice_rotmat( m(3,3) )' )

   end

   %
   % Call the MEX library.
   %
   try
      [mit] = mice('invort_c', m );
   catch
      rethrow(lasterror)
   end


