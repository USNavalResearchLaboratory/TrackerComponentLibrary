%-Abstract
%
%   CSPICE_WNFETD fetches (returns) a particular interval from a
%   double precision window.
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
%      window   SPICE window containing zero or more intervals.
%
%               [2n,1] = size(window); double = class(window)
%
%      n        index of a particular interval within the window. Indices range
%               from 1 to n, where n is the number of intervals in the window:
%
%               [1,1] = size(n); int32 = class(n)
%
%                  n = size(window,1)/2
%
%   the call:
%
%      [left, right] = cspice_wnfetd( window, n )
%
%   returns:
%
%      left,
%      right   values defining the left and right endpoints of the nth interval
%              in the input 'window';
%
%              [1,1] = size(left);  int32 = class(left)
%              [1,1] = size(right); int32 = class(right)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%     %
%     % Let 'window' contain the intervals
%     %
%     window = [ [ 1; 3 ]; [ 7; 11 ]; [ 23; 27 ];  ];
%
%     %
%     % Output the intervals.  Number of intervals equals
%     % half the number of elements for the Nx1 'window'.
%     %
%     disp( 'Window contents:')
%     for i=1:numel(window)/2
%
%        [left, right] = cspice_wnfetd( window, i );
%        fprintf( '%12.5f  %12.5f\n', left, right)
%
%     end
%
%   MATLAB outputs:
%
%      Window contents:
%           1.00000       3.00000
%           7.00000      11.00000
%          23.00000      27.00000
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine wnfetd_c.
%
%   MICE.REQ
%   WINDOWS.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 11-JUN-2013, EDW (JPL)
%
%       I/O descriptions edits to conform to Mice documentation format.
%
%   -Mice Version 1.0.0, 24-JUL-2007, EDW (JPL)
%
%-Index_Entries
%
%   fetch an interval from a d.p. window
%
%-&

function [left, right] = cspice_wnfetd(window, n)

   switch nargin

      case 2

         window = zzmice_win(window);
         n      = zzmice_int(n);

      otherwise

         error( '[left, right] = cspice_wnfetd(window, n)' )

      end

   try
      [left, right] = mice( 'wnfetd_c', [zeros(6,1); window], n );
   catch
      rethrow(lasterror)
   end




