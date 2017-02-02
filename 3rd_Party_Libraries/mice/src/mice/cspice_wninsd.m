%-Abstract
%
%   CSPICE_WNINSD inserts an interval into a double precision window.
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
%      right      left and right endpoints of the interval to insert.
%
%                 [1,1] = size(right); double = class(right)
%
%      window_i   optional input SPICE window containing zero or more
%                 intervals. Inclusion of this window argument results in an
%                 output window consisting of a union of the data in 'window_i'
%                 and the window defined as [left,right].
%
%                 [2m,1] = size(window_i); double = class(window_i)
%
%   the call:
%
%      window = cspice_wninsd( left, right )
%
%         or
%
%      window = cspice_wninsd( left, right, window_i )
%
%   returns:
%
%      window  SPICE Window consisting of either the interval [left,right] or
%              the window union of 'window_i' and [left,right]
%
%              [2n,1] = size(window); double = class(window)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Create an empty 'window' to hold the data.
%      %
%      window = zeros(0,1)
%
%      %  Let 'window' contain the intervals
%      %
%      %  [ 1, 3 ]  [ 7, 11 ]  [ 23, 27 ]
%      %
%      data = [ 1, 3, 7, 11, 23, 27];
%
%      %
%      % Add the data to 'window'.
%      %
%      for i=1:numel(data)/2
%         window = cspice_wninsd( data(2*i - 1), data(2*i), window );
%      end
%
%      %
%      % Note, the direct assignment:
%      %
%      %   window = [ [1; 3]; [7; 11]; [23; 27] ];
%      %
%      % will also perform the assignment of 'data' to 'window' but
%      % NAIF recommends using cspice_wninsd when possible.
%      %
%
%      %
%      % Perform a series of cspice_wninsd calls, adding intervals
%      % to 'window':
%      %
%
%      window = cspice_wninsd( 5,5, window )
%
%   MATLAB outputs:
%
%      window =
%
%           1
%           3
%           5
%           5
%           7
%          11
%          23
%          27
%
%      Representing the intervals:
%
%         [ 1, 3] [ 5, 5] [ 7, 11] [ 23, 27]
%
%
%      window = cspice_wninsd( 4,8, window )
%
%   MATLAB outputs:
%
%      window =
%
%           1
%           3
%           4
%          11
%          23
%          27
%
%      Representing the intervals:
%
%         [ 1, 3] [ 4, 11] [ 23, 27]
%
%
%      window = cspice_wninsd( 0,30, window )
%
%   MATLAB outputs:
%
%      window =
%
%           0
%          30
%
%      Representing the intervals:
%
%         [ 0, 30]
%
%
%      window = cspice_wninsd( 31,32, window )
%
%   MATLAB outputs:
%
%      window =
%
%           0
%          30
%          31
%          32
%
%      Representing the intervals:
%
%         [ 0, 30] [31, 32]
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine wninsd_c.
%
%   MICE.REQ
%   WINDOWS.REQ
%
%-Version
%
%   -Mice Version 1.0.2, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.1, 21-OCT-2008, EDW (JPL)
%
%      Edited example to demonstrate creation of, and loading of,
%      an empty window.
%
%   -Mice Version 1.0.0, 10-JUL-2007, EDW (JPL)
%
%-Index_Entries
%
%   insert an interval into a d.p. window
%
%-&

function [window] = cspice_wninsd( left, right, window_i )

   switch nargin
      case 2

         left   = zzmice_dp(left);
         right  = zzmice_dp(right);

      case 3

         left     = zzmice_dp(left);
         right    = zzmice_dp(right);
         window_i = zzmice_win(window_i);

      otherwise

         error ( 'Usage: [window] = cspice_wninsd( left, right, [window_i] )' )

   end

   %
   % The call passed either two or three arguments. Branch accordingly.
   %
   if ( nargin == 2 )

      try
         [window] = mice( 'wninsd_c', left, right );
      catch
         rethrow(lasterror)
      end

   else

      try
         [window] = mice( 'wninsd_c', left, right );
      catch
         rethrow(lasterror)
   end

   %
   % Perform a window union with the window defined by [left,right]
   % and the input window 'window_i'.
   %
   window = cspice_wnunid( window, window_i );

end
