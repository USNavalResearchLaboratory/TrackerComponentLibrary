%-Abstract
%
%   CSPICE_WNCARD returns the cardinality (number of intervals) of a
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
%   the call:
%
%      card = cspice_wncard( window )
%
%   returns:
%
%      card   the cardinality (number of intervals) of the input 'window'.
%
%             [1,1] = size(card); int32 = class(card)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Define a window with three intervals (six values).
%      %
%      darray = [ [ 1., 3.]; [ 7., 11.]; [23., 27.] ];
%
%      %
%      % Create a window, insert data into the window.
%      %
%      win1 = eye(0,1);
%
%      for i=1:3
%
%         win1 = cspice_wninsd( darray(i,1), darray(i,2) , win1 );
%
%      end
%
%      %
%      % What's the window cardinality of 'win1'?
%      %
%      cardinality = cspice_wncard(win1);
%
%      %
%      % Print the window cardinality (this ought to show "3" ).
%      %
%      fprintf('Number of intervals in the window: %d\n', cardinality )
%
%   Matlab outputs:
%
%      Number of intervals in the window: 3
%
%-Particulars
%
%   This function returns the value numel(window)/2.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine wncard_c.
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
%   -Mice Version 1.0.0, 30-DEC-2008, EDW (JPL)
%
%-Index_Entries
%
%   cardinality of a d.p. window
%
%-&

function [card] = cspice_wncard(window)

   switch nargin

      case 1

         window = zzmice_win(window);

      otherwise

         error ( 'Usage: [card] = cspice_wncard(window)' )

   end

   %
   % Call the windows routine, add to 'window' the space needed for
   % the control segments.
   %
   try
      [card] = mice('wncard_c', [zeros(6,1); window] );
   catch
      rethrow(lasterror)
   end


