%-Abstract
%
%   CSPICE_STPOOL retrieves the 'nth' string from the kernel pool variable
%   'item' , where the string may be continued across several components
%   of the kernel pool variable.
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
%      item     the scalar string name of a kernel pool variable for
%               which the caller wants to retrieve a full (potentially
%               continued) string
%
%      nth      the scalar integer index of the string to retrieve from
%               the kernel pool variable 'item' (index array base 1)
%
%      contin   a sequence of characters which (if they appear as the
%               last non-blank sequence of characters in a component of a
%               value of a kernel pool variable) act as a continuation
%               marker:  the marker indicates that the string associated
%               with the component containing it is continued into the
%               next literal component of the kernel pool variable
%
%               If 'contin' is a blank, all of the components of 'item'
%               will return as a single string.
%
%   the call:
%
%      [string, found] = cspice_stpool( item, nth, contin )
%
%   returns:
%
%      string   the 'nth' scalar string value corresponding to
%               the kernel pool variable specified by 'item'
%
%      found    a scalar boolean indicating true if the request
%               to retrieve the 'nth' string associated with 'item'
%               succeeds, false if not.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Load a kernel containing the variable assignment:
%      %
%      % LONG_VAL = ( 'This is part of the first component //'
%      %             'that needs more than one line when //'
%      %             'inserting it into the kernel pool.'
%      %             'This is the second string that is split //'
%      %             'up as several components of a kernel pool //'
%      %             'variable.' )
%      %
%      cspice_furnsh( 'pool_t.ker' )
%
%      %
%      % Retrieve the 'nth' entry for kernel pool variable
%      % 'LONG_VAL' to 'string'.
%      %
%      ITEM   = 'LONG_VAL';
%      CONTIN = '//';
%
%      for nth=1:3
%
%         [string, found] = cspice_stpool( ITEM, nth, CONTIN );
%
%         if ( found )
%
%            fprintf( ['Found index = %d component of kernel variable %s ' ...
%                     'in the kernel pool.\n\n'], nth, ITEM)
%
%            fprintf( 'String = ``%s``\n\n', string )
%
%         else
%
%            fprintf( ['No index = %d component of kernel variable %s ' ...
%                      'found in the kernel pool.\n'], nth, ITEM)
%
%         end
%
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs (approximately):
%
%      Found index = 1 component of kernel variable LONG_VAL in the
%      kernel pool.
%
%      String = ``This is part of the first component that needs more
%      than one line when inserting it into the kernel pool.``
%
%      Found index = 2 component of kernel variable LONG_VAL in the
%      kernel pool.
%
%      String = ``This is the second string that is split up as several
%      components of a kernel pool variable.``
%
%      No index = 3 component of kernel variable LONG_VAL found in the
%      kernel pool.
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine stpool_c.
%
%   MICE.REQ
%   KERNEL.REQ
%
%-Version
%
%   -Mice Version 1.1.0, 10-MAY-2011, EDW (JPL)
%
%      "logical" call replaced with "zzmice_logical."
%
%   -Mice Version 1.0.0, 26-SEP-2007, EDW (JPL)
%
%-Index_Entries
%
%   Retrieve a continued string value from the kernel pool
%
%-&

function [string, found] = cspice_stpool( item, nth, contin )

   switch nargin
      case 3

         item   = zzmice_str(item);
         nth    = zzmice_int(nth);
         contin = zzmice_str(contin);

      otherwise

         error ( ['Usage: [`string`, found] = ' ...
                   'cspice_stpool( `item`, nth, `contin` )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [string, found] = mice( 'stpool_c', item, nth, contin );

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch
      rethrow(lasterror)
   end


