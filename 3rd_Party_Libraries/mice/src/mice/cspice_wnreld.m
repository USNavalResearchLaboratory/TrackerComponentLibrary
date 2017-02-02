%-Abstract
%
%   CSPICE_WNRELD compares two double precision windows returning
%   a scalar boolean.
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
%      a    SPICE window containing zero or more intervals.
%
%           [2l,1] = size(a); double = class(a)
%
%      b    SPICE window containing zero or more intervals.
%
%           [2m,1] = size(b); double = class(b)
%
%      op   comparison operator, indicating the way to compare the input
%           windows. 'op' may have any of the following values:
%
%           [1,m] = size(op); char = class(op)
%
%              Operator             Meaning
%              --------  -------------------------------------
%                "="     a = b is true if 'a' and 'b' are equal
%                        (contain the same intervals).
%
%                "<>"    a <> b is true if 'a' and 'b' are not
%                               equal.
%
%                "<="    a <= b is true if 'a' is a subset of 'b'.
%
%                "<"     a < b is true is 'a' is a proper subset
%                        of 'b'.
%
%                ">="    a >= b is true if 'b' is a subset of 'a'.
%
%                ">"     a > b is true if 'b' is a proper subset
%                        of 'a'.
%
%   the call:
%
%      retval = cspice_wnreld( a, op, b )
%
%   returns:
%
%      A scalar boolean with value of the comparison.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      %  Let a contain the intervals
%      %
%      a = [ [ 1; 3 ];  [ 7; 11 ];  [ 23; 27 ] ];
%
%      %
%      %  Let b and c contain the intervals
%      %
%
%      b = [ [ 1; 2 ];  [  9; 9 ];  [ 24; 27 ] ];
%      c = b;
%
%      %
%      %  Let d contain the intervals
%      %
%      d = [ [ 5; 10 ];  [ 15; 25 ] ];
%
%      %
%      %  Finally, let e and f be empty windows (containing no intervals).
%      %
%      e = zeros(0,1);
%      f = e;
%
%      %
%      % Because b and c contain the same intervals,
%      %
%      cspice_wnreld( b, '=',  c )
%      cspice_wnreld( b, '<=', c )
%      cspice_wnreld( b, '>=', c )
%
%      %
%      % are all true, while
%      %
%      cspice_wnreld( b, '<>', c )
%
%      %
%      % is false. Because neither b nor c contains any points not also
%      % contained by the other, neither is a proper subset of the other.
%      % Thus,
%      %
%      cspice_wnreld( b, '<', c )
%      cspice_wnreld( b, '>', c )
%
%      %
%      % are both false.
%      %
%      % Every point contained in b and c is also contained in a. Thus,
%      %
%      cspice_wnreld( b, '<=', a )
%      cspice_wnreld( a, '>=', c )
%
%      %
%      % are both true. In addition, a contains points not contained in
%      % b and c. (That is, the differences a-b and a-c are not empty.)
%      % Thus, b and c are proper subsets of a as well, and
%      %
%      cspice_wnreld( b, '<', a )
%      cspice_wnreld( a, '>', b )
%
%      %
%      % are both true.
%      %
%      % Although a and d have points in common, neither contains the
%      % other. Thus
%      %
%      cspice_wnreld( a, '=',  d )
%      cspice_wnreld( a, '<=', d )
%      cspice_wnreld( a, '>=', d )
%
%      %
%      % are all false.
%      %
%      % In addition, any window is equal to itself, a subset of itself,
%      % and a superset of itself. Thus,
%      %
%      cspice_wnreld( a, '=',  a )
%      cspice_wnreld( a, '<=', a )
%      cspice_wnreld( a, '>=', a )
%
%      %
%      % are always true. However, no window is a proper subset or a
%      % proper superset of itself. Thus,
%      %
%      cspice_wnreld( a, '<', a )
%      cspice_wnreld( a, '>', a )
%
%      %
%      % are always false.
%      %
%      % Finally, an empty window is a proper subset of any window
%      % except another empty window. Thus,
%      %
%      cspice_wnreld( e, '<', a )
%
%      %
%      % is true, but
%      %
%      cspice_wnreld( e, '<', f )
%
%      %
%      % is false.
%      %
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine wnreld_c.
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
%   -Mice Version 1.0.0, 22-JUL-2007, EDW (JPL)
%
%-Index_Entries
%
%   compare two d.p. windows
%
%-&

function retval = cspice_wnreld( a, op, b )

   switch nargin

      case 3

         a  = zzmice_win(a);
         op = zzmice_str(op);
         b  = zzmice_win(b);

      otherwise

         error( 'boolean = cspice_wnreld( a, `op`, b )' )

      end

   try
      [retval] = mice( 'wnreld_c', [zeros(6,1); a], op, [zeros(6,1); b] );
      [retval] = zzmice_logical(retval);
   catch
      rethrow(lasterror)
   end



