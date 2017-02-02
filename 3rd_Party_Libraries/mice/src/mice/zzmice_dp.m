%-Abstract
%
%   ZZMICE_DP converts an numeric input to double precision format.
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
%
%   the call:
%
%      y = zzmice_dp(x)
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
%      Added value check on 'nargin'. Incorrect input argument type/form
%      error tag changed from "MICE(BADVAL)" to "MICE(BADARG)."
%
%   -Mice Version 1.0.1, 30-DEC-2008, EDW (JPL)
%
%      Function ensures all input values as finite.
%
%      Corrected misspellings.
%
%   -Mice Version 1.0.0, 30-JAN-2006, EDW (JPL)
%
%-Index_Entries
%
%   None.
%
%-&

function [y] = zzmice_dp(x)

   if( ~isequal(nargin,1) )

      error( 'MICE(USAGE): [_y_] = zzmice_dp( _x_ )' )

   end

   if( isnumeric(x) && all( isfinite( x(:) ) ) )

     y = double(x);

   else

      error( ['MICE(BADARG): Improper type of input ' ...
              'argument passed to function. Value '   ...
              'or values expected as finite double '  ...
              'precision or integer.'] )

   end


