%-Abstract
%
%   CSPICE_GNPOOL returns the names of kernel variables matching a
%   specified template.
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
%      name    the matchi_c template to use when searching for variable names
%              in the kernel pool. The characters '*' and '%' are used for the
%              wild card string and wild card character respectively.  For
%              details of string pattern matching see the header of the routine
%              matchi_c.
%
%              [1,m] = size(name); char = class(name)
%
%      start   value for the index indicating the first component of the data
%              vector assigned to 'name' for return (index 1 for all elements).
%
%              [1,1] = size(start); int32 = class(start)
%
%      room    value specifying the maximum number of components that can
%              return for 'name'.
%
%              [1,1] = size(room); int32 = class(room)
%
%   the call:
%
%      [kvars, found] = cspice_gnpool( name, start, room )
%
%   returns:
%
%      kvars   the values assigned to 'name' beginning at index 'start'.
%              'kvars' returns empty if variables described by 'name' does not
%              exist in the kernel pool.
%
%              [n,m] = size(kvars); char = class(kvars)
%
%      found   the flag indicating true if variables matching the 'name'
%              template exist in the kernel pool, false otherwise.
%
%              [1,1] = size(found); logical = class(found)
%
%              N <= 'room', M = length of longest string in return
%              array 'kvars'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example(1):
%
%      %
%      % Load a leapseconds kernel.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % A template for Jupiter kernel variables.
%      %
%      VAR = 'BODY599*';
%
%      %
%      % Query for the variable name, return all matches from
%      % index 1.
%      %
%      INDEX  = 1;
%      ROOM   = 10;
%
%      [kervar, found] = cspice_gnpool( VAR, INDEX, ROOM );
%
%      if( found )
%
%         n_elements = size(kervar, 1);
%
%         %
%         % Output the returned variable names.
%         %
%         for n=1: n_elements
%            txt = sprintf( 'Variable %d matching %s: %s', n, VAR, kervar(n,:));
%            disp( txt )
%         end
%
%      else
%         txt = sprintf( ['Failed to find  ' VAR ' in the kernel pool.'] );
%         disp( txt )
%      end
%
%   MATLAB outputs:
%
%      Variable 1 matching BODY599*: BODY599_POLE_RA
%      Variable 2 matching BODY599*: BODY599_POLE_DEC
%      Variable 3 matching BODY599*: BODY599_PM
%      Variable 4 matching BODY599*: BODY599_RADII
%      Variable 5 matching BODY599*: BODY599_LONG_AXIS
%
%   Example(2):
%
%      %
%      % Return to the array 'kervar' the names of the first
%      % 'ROOM' pool variables. Use the * wildcard character
%      % as a template to indicate a request for all kernel
%      % variables.
%      %
%      % Return all matches from 'INDEX' 1.
%      %
%      INDEX  = 1;
%      ROOM   = 10;
%      VAR    = '*';
%
%      [kervar, found] = cspice_gnpool( VAR, INDEX, ROOM );
%
%      if ( found )
%
%         n_elements = size(kervar, 1);
%
%         %
%         % Output the returned variable names.
%         %
%         for n=1: n_elements
%            txt = sprintf( 'Variable %d matching %s: %s', n, VAR, kervar(n,:));
%            disp( txt )
%         end
%
%      else
%         txt = sprintf( ['Failed to find  ' VAR ' in the kernel pool.'] );
%         disp( txt )
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Variable 1 matching *: BODY402_PM
%      Variable 2 matching *: BODY705_RADII
%      Variable 3 matching *: BODY618_POLE_RA
%      Variable 4 matching *: BODY610_PM
%      Variable 5 matching *: BODY710_LONG_AXIS
%      Variable 6 matching *: BODY806_PM
%      Variable 7 matching *: BODY299_PM
%      Variable 8 matching *: BODY705_NUT_PREC_DEC
%      Variable 9 matching *: BODY603_RADII
%      Variable 10 matching *: BODY805_LONG_AXIS
%
%   Note, the seemingly random order of the output list reflects the
%   order used by the SPICE kernel subsystem to store/lookup the
%   variable names.
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine gnpool_c.
%
%   MICE.REQ
%   KERNEL.REQ
%
%-Version
%
%   -Mice Version 1.2.0, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%      "logical" call replaced with "zzmice_logical."
%
%      I/O descriptions edits to conform to Mice documentation format.
%
%   -Mice Version 1.0.0, 15-DEC-2006, EDW (JPL)
%
%-Index_Entries
%
%   return names of kernel pool variables matching a template
%
%-&

function [kvars,found] = cspice_gnpool( name, start, room )

   switch nargin
      case 3

         name  = zzmice_str(name);
         start = zzmice_int(start);
         room  = zzmice_int(room);

      otherwise

         error ( ['Usage: [kvars(), found] = ' ...
                  'cspice_gnpool( `name`, start, room )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [kvars, found] = mice( 'gnpool_c', name, start, room );

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch
      rethrow(lasterror)
   end



