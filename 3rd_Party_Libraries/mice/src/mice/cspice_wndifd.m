%-Abstract
%
%   CSPICE_WNDIFD returns the difference of two double precision
%   windows.
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
%      a   SPICE window
%
%          [2l,1] = size(a); double = class(a)
%
%      b   SPICE window
%
%          [2m,1] = size(b); double = class(b)
%
%          Two SPICE windows containing zero or more intervals.
%
%   the call:
%
%      c = cspice_wndifd( a, b )
%
%   returns:
%
%      c   the window difference (in the SPICE sense) of 'a' and 'b', every
%          point contained in 'a', but not contained in 'b'.
%
%          [2n,1] = size(c); double = class(c)
%
%          'c' can overwrite 'a' or 'b'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Let 'a' contain the intervals
%      %
%      a = [ [ 1; 3 ]; [ 7; 11 ]; [ 23; 27 ]; ];
%
%      %
%      % and b contain the intervals
%      %
%      b = [ [ 2; 4 ]; [ 8; 10 ]; [ 16; 18 ]; ];
%
%      %
%      % Then the difference of 'a' and 'b', 'c':
%      %
%      c = cspice_wndifd(a, b)
%
%   MATLAB outputs:
%
%      c =
%
%           1
%           2
%           7
%           8
%          10
%          11
%          23
%          27
%
%      Representing the intervals:
%
%         [ 1, 2 ]  [ 7, 8 ]  [ 10, 11 ]  [ 23, 27 ]
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine wndifd_c.
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
%   -Mice Version 1.0.0, 23-JUL-2007, EDW (JPL)
%
%-Index_Entries
%
%   difference two d.p. windows
%
%-&

function [c] = cspice_wndifd( a, b )

   switch nargin

      case 2

         a    = zzmice_win(a);
         b    = zzmice_win(b);

      otherwise

         error ( 'Usage: [c] = cspice_wndifd( a, b )' )

   end

%
% Call the windows routine, add to 'a' and 'b' the space needed for
% the control segments.
%
   try
      [c] = mice('wndifd_c', [zeros(6,1); a], [zeros(6,1); b] );
   catch
      rethrow(lasterror)
   end



