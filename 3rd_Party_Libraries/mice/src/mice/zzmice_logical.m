%-Abstract
%
%   CSPICE_LOGICAL converts boolean or numeric input to the corresponding
%   logical (boolean) expression.
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
%      x   an input logical or numeric representation of a logical
%          to convert to logical.
%
%   the call:
%
%      y = zzmice_dp(x)
%
%   returns:
%
%      y   the logical representation of 'x'.
%
%-Examples
%
%   None.
%
%-Particulars
%
%   This routine exists to support the NAIF MATLAB-CSPICE interface.
%
%   Any input evaluating to 0 returns as a false, otherwise true.
%
%   All numeric inputs OTHER THAN 0 evaluate to TRUE.
%
%-Required Reading
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 08-MAY-2011, EDW (JPL)
%
%       Implemented a numeric-to-logical function due to the change
%       in "logical" behavior in octave 3.4.0. The 3.4.0 version
%       does not allow int32 inputs, which seems strange.
%
%-Index_Entries
%
%   None.
%
%-&

function y = zzmice_logical (x)

   if( ~isequal(nargin,1) )

      error( 'MICE(USAGE): [_y_] = zzmice_logical( _x_ )' )

   end

   %
   % Check input type.
   %
   if (islogical(x) && ~any(isnan(x)) )

      %
      % Logical in, logical out. Easy-peasy.
      %
      y = x;

   elseif( isnumeric(x) && all(isfinite(x))  )

      %
      % Accept logical or finite, non NaN, numeric variable. Lemon-squeezy.
      %

      %
      % Numeric input or'd with zero (false) array. The "or" operation
      % returns a logical.
      %
      y =  [ double(x) | zeros(size(x)) ];

   else

      %
      % Not a logical or numeric. Signal an error.
      %
      error( ['MICE(BADARG): Improper type of input ' ...
              'argument passed to function. Value ' ...
              'or values expected as finite numeric or logical.'] )

   end
