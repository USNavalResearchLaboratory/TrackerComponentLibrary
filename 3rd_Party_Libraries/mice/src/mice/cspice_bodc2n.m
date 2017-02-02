%-Abstract
%
%   CSPICE_BODC2N returns the body name corresponding to an input numeric
%   ID value.
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
%      code   an integer scalar or integer 1XN array of SPICE codes
%             for a set of bodies: planets, satellites, barycenters,
%             DSN stations, spacecraft, asteroids, comets, or other
%             ephemeris object.
%
%   the call:
%
%      [name, found] = cspice_bodc2n( code )
%
%   returns:
%
%      name    the scalar string or NXM character array of names associated
%              with 'code' (if 'code' has more than one translation, then
%              the most recently defined name corresponding to 'code'
%              is returned).
%
%      found   a boolean scalar or boolean 1XN array flagging if the kernel
%              subsystem translated 'code' to a corresponding name.
%
%              'found' and 'name' return with the same vectorization
%              measure (N) as 'code'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Retrieve the current body name associated to a given NAIF ID.
%      %
%      disp( 'Scalar:' )
%      naif_id = 501;
%      [name, found] = cspice_bodc2n( naif_id );
%
%      %
%      % Output the mapping if it exists.
%      %
%      if ( found )
%         fprintf( 'Body ID %i maps to name %s\n', naif_id, name );
%      end
%
%      disp( ' ' )
%
%      %
%      % Create an array of IDs. Include one unknown ID.
%      %
%      disp( 'Vector:' )
%      naif_id       = [ 502, 503, 504, 505, 5006 ];
%      [name, found] = cspice_bodc2n( naif_id );
%
%      n_elements = size(found,2);
%
%      %
%      % Loop over the output array.
%      %
%      for n=1:n_elements
%
%         %
%         % Check for a valid name/ID mapping.
%         %
%         if( found(n) )
%            fprintf( 'Body ID %i maps to name %s\n', ...
%                           naif_id(n) , name(n,:) );
%         else
%            fprintf( 'Unknown body ID %i\n', naif_id(n) );
%         end
%
%      end
%
%   MATLAB outputs:
%
%      Scalar:
%      Body ID 501 maps to name IO
%
%      Vector:
%      Body ID 502 maps to name EUROPA
%      Body ID 503 maps to name GANYMEDE
%      Body ID 504 maps to name CALLISTO
%      Body ID 505 maps to name AMALTHEA
%      Unknown body ID 5006
%
%-Particulars
%
%   A sister version of this routine exists named mice_bodc2n that returns
%   the output arguments as fields in a single structure.
%
%   cspice_bodc2n is one of five related subroutines,
%
%      cspice_bods2c      Body string to code
%      cspice_bodc2s      Body code to string
%      cspice_bodn2c      Body name to code
%      cspice_bodc2n      Body code to name
%      cspice_boddef      Body name/code definition
%
%   cspice_bods2c, cspice_bodc2s, cspice_bodn2c, and cspice_bodc2n
%   perform translations between body names and their corresponding
%   integer ID codes which are used in SPICE files and routines.
%
%   cspice_bods2c is a slightly more general version of cspice_bodn2c:
%   support for strings containing ID codes in string format enables a caller
%   to identify a body using a string, even when no name is associated with
%   that body.
%
%   cspice_bodc2s is a general version of cspice_bodc2n; the routine returns
%   either the name assigned in the body ID to name mapping or a string
%   representation of the 'code' value if no mapping exists.
%
%   cspice_boddef assigns a body name to ID mapping. The mapping has
%   priority in name-to-ID and ID-to-name translations.
%
%   Refer to NAIF_IDS.REQ for the list of name/code associations built
%   into SPICE, and for details concerning adding new name/code
%   associations at run time by loading text kernels.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine bodc2n_c.
%
%   MICE.REQ
%   NAIF_IDS.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 16-MAY-2009 (EDW)
%
%       Edit to Particulars section to document the cspice_bodc2s routine.
%       Extended argument descriptions in the I/O section.
%
%       Corrected typo in usage string.
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   body id code to name
%
%-&

function [name, found] = cspice_bodc2n(code)

   switch nargin
      case 1

         code = zzmice_int(code);

      otherwise

         error ( 'Usage: [_`name`_, _found_] = cspice_bodc2n(_code_)' )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [bodc2n] = mice('bodc2n_s', code);
      name     = char( bodc2n.name );
      found    = reshape( [bodc2n.found], 1, [] );
   catch
      rethrow(lasterror)
   end

