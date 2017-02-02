%-Abstract
%
%   CSPICE_WNINCD determines whether an interval is included in a
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
%      left,
%      right    values defining the endpoints of an interval, which may or may
%               not exist in one of the intervals in 'window'
%
%               [1,1] = size(right); double = class(right)
%
%      window   SPICE window containing zero or more intervals.
%
%               [2n,1] = size(window); double = class(window)
%
%   the call:
%
%      retval = cspice_wnincd( left, right, window )
%
%   returns:
%
%      A a boolean with value true if the input interval exists in
%     'window'
%
%         a(i)  <  left  <  right  <  b(i)
%               -        -         -
%
%      for some interval [ a(i), b(i) ] in 'window', false
%      otherwise.
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
%      window = [ [ 1; 3 ]; [ 7; 11 ]; [ 23; 27 ]; ];
%
%      %
%      % Define a set of test intervals.
%      %
%      t_arr = [ [  1, 3 ]; ...
%                [  9, 10]; ...
%                [  0, 20]; ...
%                [ 13, 15]; ...
%                [ 29, 30] ];
%
%      %
%      % Loop over the test intervals.
%      %
%      % The number of test intervals equals half the
%      % number of elements for the Nx2 array 't_arr'.
%      %
%      for i=1:numel(t_arr)/2
%
%        chk = cspice_wnincd( t_arr(i,1), t_arr(i,2), window);
%
%        if( chk )
%
%            fprintf(                                        ...
%               '%12.5f %12.5f - an element of the window\n', ...
%                                     t_arr(i,1), t_arr(i,2))
%
%         else
%
%            fprintf(                                             ...
%               '%12.5f %12.5f - not an element of the window\n', ...
%                                         t_arr(i,1), t_arr(i,2))
%
%         end
%
%      end
%
%   MATLAB outputs:
%
%           1.00000      3.00000 - an element of the window
%           9.00000     10.00000 - an element of the window
%           0.00000     20.00000 - not an element of the window
%          13.00000     15.00000 - not an element of the window
%          29.00000     30.00000 - not an element of the window
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine wnincd_c.
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
%      "logical" call replaced with "zzmice_logical."
%
%      Corrected version ID in 23-JUL-2009 entry, "1.0.0" to "1.0.1."
%
%   -Mice Version 1.0.1, 23-JUL-2009, EDW (JPL)
%
%      Replaced 'boolean' calls with 'logical' as 'boolean' functionally
%      aliases 'logical'.
%
%   -Mice Version 1.0.0, 17-JUL-2007, EDW (JPL)
%
%-Index_Entries
%
%   included in a d.p. window
%
%-&

function retval = cspice_wnincd( left, right, window )

   switch nargin

      case 3

         left   = zzmice_dp(left);
         right  = zzmice_dp(right);
         window = zzmice_win(window);

      otherwise

         error( 'boolean = cspice_wnincd( left, right, window )' )

      end

   try
      [retval] = mice( 'wnincd_c', left, right, [zeros(6,1); window] );
      [retval] = zzmice_logical(retval);
   catch
      rethrow(lasterror)
   end

