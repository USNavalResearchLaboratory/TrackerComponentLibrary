%-Abstract
%
%   CSPICE_WNEXTD extracts the left or right endpoints from
%   a double precision window.
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
%      side     indicating whether to extract the left or right endpoints of
%               the intervals in 'window'
%
%               [1,m] = size(side); char = class(side)
%
%                 'L', 'l'       Left endpoints.
%                 'R', 'r'       Right endpoints.
%
%      window   SPICE window containing zero or more intervals.
%
%               [2n,1] = size(window); double = class(window)
%
%   the call:
%
%      window_f = cspice_wnextd( side, window)
%
%   returns:
%
%      window   SPICE window containing zero or more intervals, representing
%               the collection of singleton intervals containing either the
%               left or the right endpoints of the intervals in the original
%               'window'
%
%               [2n,1] = size(window); double = class(window)
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
%      window = [ [ 1; 3 ]; [ 7; 11 ]; [ 23; 27 ]; [ 29; 29 ]; ];
%
%      %
%      %
%      %
%      left  = cspice_wnextd( 'L', window )
%
%      right = cspice_wnextd( 'R', window )
%
%   MATLAB outputs:
%
%      left =
%
%           1
%           1
%           7
%           7
%          23
%          23
%          29
%          29
%
%      right =
%
%           3
%           3
%          11
%          11
%          27
%          27
%          29
%          29
%
%      %
%      % A similar operation using MATLAB native functionality,
%      % though the returned arrays are not Mice windows identical
%      % to the cspice_wnextd result.
%      %
%      index_left = 1: 2 : numel(window);
%      index_right= index_left + 1;
%
%      left  = window( index_left )
%      right = window( index_right)
%
%   MATLAB outputs:
%
%      left =
%
%           1
%           7
%          23
%          29
%
%      right =
%
%           3
%          11
%          27
%          29
%
%-Particulars
%
%   This function returns a window composed of singleton intervals
%   containing one of the endpoints of the intervals in 'window'.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine wnextd_c.
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
%   extract the endpoints from a d.p. window
%
%-&

function [window_f] = cspice_wnextd( side, window)

   switch nargin

      case 2

         side   = zzmice_str(side);
         window = zzmice_win(window);

      otherwise

         error ( 'Usage: [window_f] = cspice_wnextd( side, window)' )

   end

%
% Please note, this call does not require addition of space for the window
% control segment as needed by other windows interfaces. The interface
% copies the data in 'window' to a work variable rather than directly
% pass 'window' into a CSPICE call.
%
   try
      [window_f] = mice('wnextd_c', side, window );
   catch
      rethrow(lasterror)
   end



