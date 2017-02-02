%-Abstract
%
%   CSPICE_NPLNPT calculates the location on a defined line
%   nearest to a specified point, then determines the distance
%   between the two points.
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
%      linpt    [3,1] = size(linpt); double = class(linpt)
%
%      lindir   [3,1] = size(lindir); double = class(lindir)
%
%               are, respectively, a point and a direction vector that define
%               a line.  The line is the set of vectors
%
%                     linept   +   t * linedr
%
%               where t is any real number.
%
%      point    a point in 3-dimensional space.
%
%               [3,n] = size(point); double = class(point)
%
%   the call:
%
%      [ pnear, dist ] = cspice_nplnpt( linpt, lindir, point )
%
%   returns:
%
%      pnear   the nearest point on the input line to the input 'point'.
%
%              [3,1] = size(pnear); double = class(pnear)
%
%      dist    distance between the input line and input point.
%
%              [1,n] = size(dist); double = class(dist)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Define a point on a line, a direction for the line, and
%      % an arbitrary point in space.
%      %
%      linept = [  1, 2,  3 ]';
%      linedr = [  0, 1,  1 ]';
%      point  = [ -6, 9, 10 ]';
%
%      %
%      % Calculate the location on the line nearest the point
%      % and the distance between the location and the defined
%      % point.
%      %
%      [ pnear, dist ] = cspice_nplnpt( linept, linedr, point )
%
%   MATLAB outputs:
%
%      pnear =
%
%           1
%           9
%          10
%
%
%      dist =
%
%           7
%
%-Particulars
%
%   For every line L and point P, there is a unique closest point
%   on L to P.  Call this closest point C.  It is always true that
%   P - C  is perpendicular to L, and the length of P - C is called
%   the "distance" between P and L.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine nplnpt_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 14-NOV-2013, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   distance between point and line
%   nearest point on line to point
%
%-&

function [ pnear, dist ] = cspice_nplnpt( linpt, lindir, point )

   switch nargin
      case 3

         linpt  = zzmice_dp(linpt);
         lindir = zzmice_dp(lindir);
         point  = zzmice_dp(point);

      otherwise

         error ( ['Usage: [ pnear(3), dist] = ' ...
                  'cspice_nplnpt( linpt(3), lindir(3), point(3) )'] )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [nplnpt] = mice( 'nplnpt_s', linpt, lindir, point );
      pnear    = reshape( [nplnpt.pos], 3, [] );
      dist     = reshape( [nplnpt.alt], 1, [] );
   catch
      rethrow(lasterror)
   end


