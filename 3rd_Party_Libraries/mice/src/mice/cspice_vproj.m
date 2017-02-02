%-Abstract
%
%   CSPICE_VPROJ calculates the projection of a set of 3-vectors onto
%   another set of 3-vectors.
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
%      a   the vector(s) to project onto the vector(s) 'b'.
%
%          [3,n] = size(a); double = class(a)
%
%      b   the vector(s) to receive the projection(s).
%
%          [3,n] = size(b); double = class(b)
%
%      An implicit assumption exists that 'a' and 'b' are specified
%      in the same reference frame. If this is not the case, the numerical
%      result has no meaning.
%
%   the call:
%
%      vproj = cspice_vproj( a, b )
%
%   returns:
%
%      vproj   vector containing the projection(s) of 'a' onto 'b' ('vproj' is
%              necessarily parallel to 'b'.)  If 'b' equals the zero vector
%              then the zero vector will return as 'vproj'.
%
%              [3,n] = size(vproj); double = class(vproj)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Define two vector sets.
%      %
%      a = [ [ 6, 6, 6]', ...
%            [ 6, 6, 6]', ...
%            [ 6, 6, 0]', ...
%            [ 6, 0, 0]' ]
%
%      b = [ [ 2, 0, 0]', ...
%            [-3, 0, 0]', ...
%            [ 0, 7, 0]', ...
%            [ 0, 0, 9]' ]
%
%      %
%      % Calculate the projection.
%      %
%      p = cspice_vproj( a, b )
%
%   MATLAB outputs:
%
%      a =
%
%           6     6     6     6
%           6     6     6     0
%           6     6     0     0
%
%
%      b =
%
%           2    -3     0     0
%           0     0     7     0
%           0     0     0     9
%
%
%      p =
%
%           6     6     0     0
%           0     0     6     0
%           0     0     0     0
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine vproj_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   3-vector projection
%
%-&

function [vproj] = cspice_vproj( a, b)

   switch nargin
      case 2

         a = zzmice_dp(a);
         b = zzmice_dp(b);

      otherwise

         error ( 'Usage: [_vproj(3)_] = cspice_vproj(_a(3)_, _b(3)_)' )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      [vproj] = mice('vproj_c', a, b);
   catch
      rethrow(lasterror)
   end

