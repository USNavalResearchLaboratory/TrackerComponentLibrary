%-Abstract
%
%   ZZMICE_WIN converts an numeric input to double precision format,
%   enforcing an even array length (Nx1, with N even).
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
%      x   an input numeric to convert to double precision
%          window, i.e. an Nx1 double precision array.
%
%   the call:
%
%      y = zzmice_win(x)
%
%   returns:
%
%      y   the double precision representation of 'x' confirmed
%          having an even number of elements.
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
%   -Mice Version 1.2.0, 27-JUL-2009, EDW (JPL)
%
%      Added value check on 'nargin'. Incorrect input argument type/form
%      error tag changed from "MICE(BADVAL)" to "MICE(BADARG)."
%
%      'even' function renamed to 'zzmice_even'.
%
%   -Mice Version 1.1.0, 30-DEC-2008, EDW (JPL)
%
%      Improved type-shape checks on the input array 'x'. Added the tests:
%
%         issorted
%         isequal( size(x_dim), [1,2])
%
%      Allow users to pass '[]' as an input window. Return the zeros(0,1)
%      array as output.
%
%      Corrected misspellings.
%
%   -Mice Version 1.0.0, 26-JUN-2007, EDW (JPL)
%
%-Index_Entries
%
%   None.
%
%-&

function [y] = zzmice_win(x)

   if( ~isequal(nargin,1) )

      error( 'MICE(USAGE): y = zzmice_win( x)' )

   end

   %
   % Syntax sugar. Allow an empty array as input.
   %
   if( isequal( x, [] ) )

      %
      % Return a double precision 0x1 array of zeros.
      %
      y = zeros(0,1);
      return;

   end


   %
   % Retrieve the dimension of 'x'. Note that all Matlab variables
   % have at least two dimensions, [m,n].
   %
   x_dim     = size(x);
   [ N, col] = size(x);

   %
   % Check the validity of the window:
   %
   %    - Numerical data
   %    - Nx1 dimension
   %    - N even
   %    - Numeric values in ascending oder (sorted)
   %    - Size of size array equals [1,2]
   %
   if( isnumeric(x)        && ...
       zzmice_even(N)      && ...
       (col == 1)          && ...
       issorted(x)         && ...
       isequal(size(x_dim), [1,2])   )

      y = double(x);

   else

      error( ['MICE(BADARG): Improper type of input '          ...
              'argument passed to function. Input expected ' ...
              'as numeric Nx1 array in ascending order, with ' ...
              'either an even number of elements or an empty ' ...
              '0x1 array.' ] )

   end


%
% Simple function to test a numeric for evenness.
%
function bool = zzmice_even(n)

   bool = ( mod(n,2)==0 );


