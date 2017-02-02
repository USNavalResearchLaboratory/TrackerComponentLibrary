%-Abstract
%
%   CSPICE_BODN2C translates the name of a body or object to the corresponding
%   SPICE integer ID code.
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
%      name   name(s) of a body or object,  such as a planet, satellite, comet,
%             asteroid, barycenter, DSN station, spacecraft, or instrument,
%             "known" to the SPICE system, whether through hard-coded
%             registration or run-time registration in the SPICE kernel pool
%
%             [n,m] = size(name); char = class(name)
%
%             Case and leading and trailing blanks in a name are not
%             significant. However when a name is made up of more than one
%             word, they must be separated by at least one blank. That is,
%             all of the following strings are equivalent names:
%
%                      'JUPITER BARYCENTER'
%                      'Jupiter Barycenter'
%                      'JUPITER BARYCENTER   '
%                      'JUPITER    BARYCENTER'
%                      '   JUPITER BARYCENTER'
%
%              However, 'JUPITERBARYCENTER' is not equivalent to the names
%              above.
%
%   the call:
%
%      [code, found] = cspice_bodn2c( name )
%
%   returns:
%
%      code    containing the SPICE code(s) assigned either by SPICE or the
%              user to 'name'.
%
%              [1,n] = size(code); int32 = class(code)
%
%      found   flag(s) indicating if the kernel subsystem translated 'name' to
%              a corresponding 'code'.
%
%              [1,n] = size(found); logical = class(found)
%
%              'found' and 'code' return with the same vectorization
%              measure (N) as 'name'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Retrieve the NAIF ID associated to a body name.
%      %
%      disp( 'Scalar:' )
%      name            = 'Hyperion';
%      [ code, found ] = cspice_bodn2c( name );
%
%      %
%      % Output the mapping if it exists.
%      %
%      if ( found )
%         txt = sprintf( 'Body code %i maps to name %s', ...
%                         code, name );
%         disp(txt)
%      end
%
%      disp(' ')
%
%      %
%      % Create an array of body names. Include one unknown name.
%      %
%      disp( 'Vector:' )
%      name           = strvcat( 'Triton', 'Mimas', ...
%                                'Oberon', 'Callisto', 'Halo' );
%      [ code, found] = cspice_bodn2c( name );
%
%      n_elements = size(code,2);
%
%      %
%      % Loop over the output array.
%      %
%      for n=1:n_elements
%
%         %
%         % Check for a valid name/code mapping.
%         %
%         if( found(n) )
%            txt = sprintf( 'Body code %i maps to name %s', ...
%                            code(n), name(n,:) );
%            disp(txt)
%         else
%            txt = sprintf( 'Unknown body name %s', name(n,:) );
%            disp(txt)
%         end
%
%      end
%
%   MATLAB outputs:
%
%      Scalar:
%      Body code 607 maps to name Hyperion
%
%      Vector:
%      Body code 801 maps to name Triton
%      Body code 601 maps to name Mimas
%      Body code 704 maps to name Oberon
%      Body code 504 maps to name Callisto
%      Unknown body name Halo
%
%-Particulars
%
%   A sister version of this routine exists named mice_bodn2c that returns
%   the output arguments as fields in a single structure.
%
%   cspice_bodn2c is one of five related subroutines,
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
%   Programmers writing user interface code should consider using the
%   Mice routine cspice_bods2c. cspice_bods2c provides more flexibility
%   in handling input strings, since it accepts both body names and
%   strings representing integer ID codes, for example '399'.
%
%   Refer to NAIF_IDS.REQ for the list of name/code associations built
%   into SPICE, and for details concerning adding new name/code
%   associations at run time by loading text kernels.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine bodn2c_c.
%
%   MICE.REQ
%   NAIF_IDS.REQ
%
%-Version
%
%   -Mice Version 1.0.2, 12-MAR-2012 (EDW), SCK (JPL)
%
%       I/O descriptions edits to conform to Mice documentation format.
%
%       Corrected minor typo in header.
%
%   -Mice Version 1.0.1, 16-MAY-2009 (EDW)
%
%       Edit to Particulars section to document the cspice_bodc2s routine.
%       Extended argument descriptions in the I/O section.
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   body name to code
%
%-&

function [code,found] = cspice_bodn2c(name)

   switch nargin
      case 1

         name = zzmice_str(name);

      otherwise

         error ( 'Usage: [_code_, _found_] = cspice_bodn2c(_`name`_)' )

   end


   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [bodn2c] = mice('bodn2c_s',name);
      code     = reshape( [bodn2c.code],  1, [] );
      found    = reshape( [bodn2c.found], 1, [] );
   catch
      rethrow(lasterror)
   end


