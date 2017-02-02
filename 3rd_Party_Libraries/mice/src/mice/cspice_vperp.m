%-Abstract
%
%   CSPICE_VPERP calculates the component of a vector perpendicular to a
%   second vector.
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
%      a   the 3-vector(s) whose component orthogonal to 'b' is sought.
%
%          [3,n] = size(a); double = class(a)
%
%          (There is a unique decomposition of a into a sum v + p, where v is
%          parallel to b and p is orthogonal to b.  We want the component p.)
%
%      b   the second 3-vector(s) used as a reference for the decomposition
%          of 'a'.
%
%          [3,n] = size(b); double = class(b)
%
%      An implicit assumption exists that 'a' and 'b' are specified
%      in the same reference frame. If this is not the case, the numerical
%      result has no meaning.
%
%   the call:
%
%      vperp = cspice_vperp( a, b )
%
%   returns:
%
%      vperp   the 3-vector(s) containing the component of 'a' orthogonal
%              to 'b'.
%
%              [3,n] = size(vperp); double = class(vperp)
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
%      % Calculate the decomposition.
%      %
%      p = cspice_vperp( a, b )
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
%           0     0     6     6
%           6     6     0     0
%           6     6     0     0
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine vperp_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 09-NOV-2012, EDW (JPL)
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.0, 22-APR-2010, EDW (JPL)
%
%-Index_Entries
%
%   perpendicular component of a 3-vector
%
%-&

function [vperp] = cspice_vperp( a, b)

   switch nargin
      case 2

         a = zzmice_dp(a);
         b = zzmice_dp(b);

      otherwise

         error ( 'Usage: [_vperp(3)_] = cspice_vperp(_a(3)_, _b(3)_)' )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      [vperp] = mice('vperp_c', a, b);
   catch
      rethrow(lasterror)
   end



