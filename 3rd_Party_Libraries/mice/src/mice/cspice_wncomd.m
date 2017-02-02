%-Abstract
%
%   CSPICE_WNCOMD returns the complement of a double precision
%   window with respect to a specified interval.
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
%      left,
%      right    values defining the left and right endpoints of the complement
%               interval.
%
%               [1,1] = size(right); double = class(right)
%
%      window   SPICE window containing zero or more intervals
%
%               [2m,1] = size(window); double = class(window)
%
%   the call:
%
%      result = cspice_wncomd( left, right, window)
%
%   returns:
%
%      result   SPICE window containing the complement of 'window' with respect
%               to the interval 'left' to 'right'
%
%               [2n,1] = size(result); double = class(result)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Let 'window' contain the intervals
%      %
%      window = [ [ 1; 3 ];  [ 7; 11 ];  [ 23; 27 ];  ];
%
%      %
%      % The floating point complement of window with respect
%      % to [2,20]
%      %
%      cspice_wncomd( 2, 20, window )
%
%   MATLAB outputs:
%
%      a =
%
%           3
%           7
%          11
%          20
%
%      Representing the intervals:
%
%      [ 3, 7 ]  [ 11, 20 ]
%
%      %
%      % The complement with respect to [ 0, 100 ]
%      %
%      cspice_wncomd( 0, 100, window )
%
%   MATLAB outputs:
%
%      b =
%
%           0
%           1
%           3
%           7
%          11
%          23
%          27
%         100
%
%      Representing the intervals:
%
%      [ 0, 1 ]  [ 3, 7 ]  [ 11, 23 ]  [ 27, 100 ]
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine wncomd_c.
%
%   MICE.REQ
%   WINDOWS.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.0, 24-JUL-2007, EDW (JPL)
%
%-Index_Entries
%
%   complement a d.p. window
%
%-&

function [result] = cspice_wncomd( left, right, window)

   switch nargin

      case 3

         left   = zzmice_dp(left);
         right  = zzmice_dp(right);
         window = zzmice_win(window);

      otherwise

         error ( 'Usage: [result] = cspice_wncomd( left, right, window)' )

   end

   %
   % Call the windows routine, add to 'a' and 'b' the space needed for
   % the control segments.
   %
   try
      [result] = mice('wncomd_c', left, right, [zeros(6,1); window] );
   catch
      rethrow(lasterror)
   end


