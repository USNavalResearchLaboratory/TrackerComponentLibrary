%-Abstract
%
%   ZZMICE_PLN enforces a structure as a SPICE plane.
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
%      x   a scalar structure representing a SPICE plane
%
%   the call:
%
%      y = zzmice_pln(x)
%
%   returns:
%
%      y   a copy of 'x', confirmed to have the proper fields,
%          field types, and field dimensions of a SPICE plane
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
%   -Mice Version 1.2.0, 05-AUG-2010, EDW (JPL)
%
%      Update to remove unneeded vectorization check. No Mice routines
%      allow vectorized planes as input arguments.
%
%   -Mice Version 1.1.0, 27-JUL-2009, EDW (JPL)
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

function [y] = zzmice_pln(x)


   %
   % Confirm one input argument.
   %
   if( ~isequal(nargin,1) )

      error( 'MICE(USAGE): y = zzmice_pln( x)' )

   end


   %
   % Confirm input as a struct.
   %
   if( isstruct(x) )

      %
      % Confirm proper shape for scalar expression, [1,1].
      %
      if( ~isequal( size( x ), [1,1] )  )

         error( ['MICE(BADARG): Improper type of input '    ...
                'argument passed to zzmice_pln. Incorrect ' ...
                'shape, not [1,1].'] )
      end

      fields = { 'normal'; 'constant' };
      names  = fieldnames( x );

   else

      error( ['MICE(BADARG): Improper type of input ' ...
              'argument passed to zzmice_pln. Not a structure.'] )

   end


   %
   % Confirm expected fields.
   %
   if( isequal(fields,names) )

      %
      % The 'normal' field must have dimensions 3x1.
      %
      if( ~isequal( size( x.normal ), [3,1] ) )

         error( ['MICE(BADFIELD): Improper dimensions for ' ...
                    'normal field, not [3,1].'] )
      end

      %
      % The 'constant' field must have dimensions 1x1.
      %
      if( ~isequal( size( x.constant ), [1,1] ) )

         error( ['MICE(BADFIELD): Improper dimensions for ' ...
                    'constant field, not [1,1].'] )
      end

      %
      % Cast the fields to double, check appropriate type.
      %
      y.normal   = zzmice_dp(x.normal);
      y.constant = zzmice_dp(x.constant);

   else

      error( ['MICE(BADARG): Improper type of input '    ...
              'argument passed to zzmice_pln. Improper ' ...
              'form for plane structure.'] )

   end


