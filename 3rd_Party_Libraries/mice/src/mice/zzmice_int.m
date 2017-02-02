%-Abstract
%
%   ZZMICE_INT converts an numeric input to integer 32 format,
%   truncating the value if necessary.
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
%      x       an input numeric to convert to integer 32.
%
%      range   an optional description on the range of 'x',
%              define as a 1x2 array. When used, 'range' must
%              satisfy dimension of 1x2 and range(1) < range(2).
%
%   the call:
%
%      y = zzmice_int(x)
%
%         or
%
%      y = zzmice_int(x, range )
%
%   returns:
%
%      y   the 32 bit integer representation of 'x' with the
%          property range(1) <= x <= range(2) if 'range' included
%          as an input.
%
%-Examples
%
%   None.
%
%-Particulars
%
%   This routine exists to support the NAIF MATLAB-CSPICE interface.
%
%   By default, MATLAB treats numeric values as double precision values,
%   ZZMICE_INT converts an input number to a MATLAB integer 32 value before
%   passing to the MATLAB-CSPICE interface layer.
%
%   Any Mice interface passing integers to the mex library must
%   use this routine.
%
%-Required Reading
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.3.0, 11-MAY-2011, EDW (JPL)
%
%      Edits to argument check logic to improve code structure.
%
%   -Mice Version 1.2.0, 27-JUL-2009, EDW (JPL)
%
%      Incorrect input argument type/form error tag changed from
%      "MICE(BADVAL)" to "MICE(BADARG)."
%
%      Replaced "==" with "isequal."
%
%   -Mice Version 1.1.0, 30-DEC-2008, EDW (JPL)
%
%      Function ensures all input values as finite.
%
%      Range argument as an optional input to restrict
%      value of 'x' to within a closed set.
%
%      Corrected misspellings.
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   None.
%
%-&

function [y] = zzmice_int(x, range )

   %
   % Confirm a valid type for 'x'.
   %
   is_valid = (isnumeric(x) || islogical(x)) && all( isfinite( x(:) ) );

   if( is_valid )

      switch nargin

         case 0

            error( 'MICE(USAGE): _y_ = zzmice_int( _x_, [range])' )

         case 1

            y = int32(x);

         case 2

            y = int32(x);

            if( ~isequal( size(range), [1,2] ) )

               error( ['MICE(BADVAL): The range input requires ' ...
                      'dimension 1x2.' ] )

            end

            if( range(1) >= range(2) )

               error( [ 'MICE(BADVAL): The range input requires ' ...
                       'range(1) < range(2).' ] )

            end

            if ( any( y(:) < range(1) ) || any( y(:) > range(2) ) )

               txt = sprintf( ['MICE(BADVAL): Integer input value not ' ...
                               'within required range [ %ld, %ld ]' ],  ...
                               range(1), range(2) );
               error(txt)

            end

         otherwise

            %
            % Program flow cannot reach this block (I hope), but I included it
            % for completeness.
            %
            error( 'MICE(USAGE): [_y_] = zzmice_int( _x_, [range])' )

      end

   else

      error( ['MICE(BADARG): Improper type of input ' ...
              'argument passed to function. Value ' ...
              'or values expected as a finite numeric or logical.'] )

   end


