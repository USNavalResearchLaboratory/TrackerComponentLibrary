%-Abstract
%
%   CSPICE_WNVALD forms a valid double precision window from the contents
%   of an 2Mx1 array.
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
%      window_i   2M endpoints of (possibly unordered and non-disjoint)
%                 intervals.
%
%                 [2m,1] = size(window_i); double = class(window_i)
%
%   the call:
%
%      window_f = cspice_wnvald(window_i)
%
%   returns:
%
%      window_f   SPICE window containing the ordered union of the intervals in
%                 the input array 'window_i'.
%
%                 [2n,1] = size(window_f); double = class(window_f)
%
%                 The 'window_f' output may overwrite the input 'window_i'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      Define an array containing a set of unordered intervals:
%
%      >> windata = [ [0;   0]; ...
%                     [10; 12]; ...
%                     [2;   7]; ...
%                     [13; 15]; ...
%                     [1;   5]; ...
%                     [23; 29] ];
%
%      Validate the array as a SPICE window:
%
%      >> w1 = cspice_wnvald(windata)
%
%      w1 =
%
%           0
%           0
%           1
%           7
%          10
%          12
%          13
%          15
%          23
%          29
%
%   Representing the intervals:
%
%     [0, 0] [1, 7] [10, 12] [13, 15] [23, 29]
%
%-Particulars
%
%   On input, 'window' is a 2Mx1 array. During validation, the intervals
%   are ordered, and overlapping intervals are merged. The size of the output
%   window equals the number of endpoints remaining, and window is ready for
%   use with any of the window routines.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine wnvald_c.
%
%   MICE.REQ
%   WINDOWS.REQ
%
%-Version
%
%   -Mice Version 1.1.1, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.1.0, 27-JUL-2009, EDW (JPL)
%
%      'zzmice_cell' modified to include the cell type identifier 'double'.
%
%   -Mice Version 1.0.0, 17-DEC-2008, EDW (JPL)
%
%-Index_Entries
%
%   validate a d.p. window
%
%-&

function [window_f] = cspice_wnvald(window_i)

   switch nargin

      case 1

         %
         % In this case, the interface requires only a numeric
         % 'window_i' with 2Nx1 dimension.
         %

         window_i = zzmice_cell( window_i, 'double');

      otherwise

         error ( 'Usage: [window_f] = cspice_wnvald(window_i)' )

   end

   %
   % Please note, this call does not require addition of space for the window
   % control segment as needed by other windows interfaces. The interface
   % copies the data in 'window_i' to a work variable rather than directly
   % pass 'window_i' into a CSPICE call.
   %
   try
      [window_f] = mice('wnvald_c', window_i );
   catch
      rethrow(lasterror)
   end



