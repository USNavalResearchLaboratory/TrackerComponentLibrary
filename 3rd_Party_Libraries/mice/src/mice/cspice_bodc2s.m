%-Abstract
%
%   CSPICE_BODC2S translates a body ID code to either the corresponding name
%   or if no name to ID code mapping exists, the string representation of
%   the body ID value.
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
%      code   an integer scalar or integer 1XN array of SPICE ID codes
%             for a body: planet, satellite, barycenter, spacecraft,
%             asteroid, comet, or other ephemeris object.
%
%   the call:
%
%      [name] = cspice_bodc2s( code )
%
%   returns:
%
%      name   the scalar string or NXM character array of names of the
%             bodies identified by 'code' if a mapping between 'code' and
%             a body name exists within SPICE.
%
%             If 'code' has more than one translation, then the most recently
%             defined 'name' corresponding to 'code' is returned. 'name' will
%             have the exact format (case and blanks) as when the name/code pair
%             was defined.
%
%             If the input value of 'code' does not map to a body name, 'name'
%             returns with the string representation of 'code'.
%
%             'name' returns with the same vectorization measure (N) as 'code'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Assign an array of body ID codes. Not all the listed codes
%      % map to a body name.
%      %
%      code = [ 399, 0, 3, -77, 11, -1, 6000001 ];
%
%      name = cspice_bodc2s( code );
%
%      for i=1:7
%         fprintf( '%8d  %s\n', code(i), name(i,:) )
%      end
%
%   MATLAB outputs:
%
%          399  EARTH
%            0  SOLAR SYSTEM BARYCENTER
%            3  EARTH BARYCENTER
%          -77  GALILEO ORBITER
%           11  11
%           -1  GEOTAIL
%      6000001  6000001
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
%   the CSPICE routine bodc2s_c.
%
%   MICE.REQ
%   NAIF_IDS.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 01-JUN-2009, EDW (JPL)
%
%-Index_Entries
%
%   body id code to string
%
%-&

function [name] = cspice_bodc2s(code)

   switch nargin
      case 1

         code = zzmice_int(code);

      otherwise

         error ( 'Usage: [_`name`_] = cspice_bodc2s(_code_)' )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [bodc2s] = mice('bodc2s_s', code);
      name     = char( bodc2s.name );
   catch
      rethrow(lasterror)
   end

