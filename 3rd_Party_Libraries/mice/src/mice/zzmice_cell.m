%-Abstract
%
%   ZZMICE_CELL converts an numeric input to double precision format,
%   enforcing an Nx1 dimension.
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
%      x           an input numeric to convert to double precision
%                  cell, i.e. an Nx1 double precision array
%
%      cell_type   string naming the type of data for the output cell.
%                  Most likely values for 'cell_type':
%
%                     'double'
%                     'int32'
%
%   the call:
%
%      y = zzmice_cell(x cell_type)
%
%   returns:
%
%      y   the double precision representation of 'x'
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
%      Added 'cell_type' argument to identify the output type.
%
%      Added value check on 'nargin'. Incorrect input argument type/form
%      error tag changed from "MICE(BADVAL)" to "MICE(BADARG)."
%
%      Replaced "~=" with "~isequal."
%
%   -Mice Version 1.0.0, 30-DEC-2008, EDW (JPL)
%
%-Index_Entries
%
%   None.
%
%-&

function [y] = zzmice_cell(x, cell_type)

   if( ~isequal(nargin,2) )

      error( 'MICE(USAGE): _y_ = zzmice_cell( _x_, `cell_type`)' )

   end

   %
   % Syntax sugar. Allow an empty array as input.
   %
   if( isequal( x, [] ) )

      %
      % Return a 'cell_type' 0x1 array of zeros.
      %
      y = cast( zeros(0,1), cell_type);
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
   %    Numerical data
   %    Nx1 dimension
   %
   if( isnumeric(x) && (col == 1) && isequal(size(x_dim), [1,2]) )

     y = cast( x, cell_type);
     return;

   else

      error( ['MICE(BADARG): Improper type of input '        ...
              'argument passed to function. Input expected ' ...
              'as numeric Nx1 array.' ] )

   end



