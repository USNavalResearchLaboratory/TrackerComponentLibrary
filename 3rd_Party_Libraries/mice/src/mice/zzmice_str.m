%-Abstract
%
%   ZZMICE_STR enforces an input as a string (i.e, an array of
%   characters) or string cell, converting the cell to a
%   character array.
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
%      x   an input to confirm as a character array representing
%          a scalar string, vector or strings, or a string cell.
%
%   the call:
%
%      y = zzmice_str(x)
%
%   returns:
%
%      y   the character array form of 'x'
%
%-Examples
%
%   None.
%
%-Particulars
%
%   This routine exists to support the NAIF MATLAB-CSPICE interface.
%
%-Required Reading
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.1.0, 27-JUL-2009, EDW (JPL)
%
%      Added value check on 'nargin'. Incorrect input argument type/form
%      error tag changed from "MICE(BADVAL)" to "MICE(BADARG)."
%
%   -Mice Version 1.0.1, 17-OCT-2008, EDW (JPL)
%
%      Edited error description string for clarity.
%
%      Completed header sections.
%
%   -Mice Version 1.0.0, 03-APR-2007, EDW (JPL)
%
%-Index_Entries
%
%   None.
%
%-&

function [y] = zzmice_str(x)

   if( ~isequal(nargin,1) )

      error( 'MICE(USAGE): _`y`_ = zzmice_str( _`x`_ )' )

   end

   %
   % Check the type of 'x'.
   %
   if( iscellstr(x) )

      %
      %  Covert 'x' to string(s) if a cell array.
      %
      x = char(x);

   end

   if( ischar(x) )

      %
      % 'x' is a character array. Set the return value 'y'.
      %
      y = x;

   else

      %
      % Input 'x' is not a character array or string cell.
      % Signal an error.
      %
      error( [ 'MICE(BADARG): Improper type of input ' ...
               'argument passed to function. Value ' ...
               'or values expected as string cell or ' ...
               'character array.' ] )

   end

   %
   % Do not allow an empty string as input.
   %
   if( numel(x) == 0  )

      error( 'MICE(BADARG): Attempt to use null string as input.' )

   end

