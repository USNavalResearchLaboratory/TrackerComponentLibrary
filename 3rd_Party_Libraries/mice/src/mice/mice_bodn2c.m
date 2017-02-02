%-Abstract
%
%   MICE_BODN2C translates the name of a body or object to the
%   corresponding SPICE integer ID code.
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
%      name    the scalar string or NxM character array of names of a set
%              of bodies or objects, such as planets, satellites, comets,
%              asteroids, barycenters, DSN stations, spacecraft, or
%              instruments, "known" to the SPICE system, whether through
%              hard-coded registration or run-time registration in
%              the SPICE kernel pool
%
%   the call:
%
%      ID = mice_bodn2c( name )
%
%   returns:
%
%      ID   the scalar or 1xN array of structures associating
%           a body name with a corresponding NAIF ID. Each structure
%           contains two fields:
%
%              name   the "name" of a particular body. If a mapping
%                     does not exist, the 'name' field returns as NULL.
%
%              code   a scalar integer SPICE code assigned either
%                     by SPICE or the user to 'name'. If a mapping
%                     does not exist, the 'code' field returns as 0.
%
%              found  a scalar boolean indicating if a mapping exists (true)
%
%              'ID' returns with the same vectorization measure (N) as 'name'.
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
%      name = 'Hyperion';
%      ID   = mice_bodn2c( name );
%
%      %
%      % Output the mapping if it exists.
%      %
%      if ( ID.found )
%         txt = sprintf( 'Body ID %i maps to name %s', ...
%                         ID.code, ID.name );
%         disp(txt)
%      end
%
%      disp(' ')
%
%      %
%      % Create an array of body names. Include one unknown name.
%      %
%      disp( 'Vector:' )
%      name = strvcat( 'Triton', 'Mimas', 'Oberon', 'Callisto', 'Halo' );
%      ID   = mice_bodn2c( name );
%
%      n_elements = size(ID);
%
%      %
%      % Loop over the output array.
%      %
%      for i=1:n_elements(1)
%
%         %
%         % Check for a valid name/ID mapping.
%         %
%         if(  ID(i).name ) )
%            txt = sprintf( 'Body ID %i maps to name %s', ...
%                            ID(i).code, ID(i).name );
%            disp(txt)
%         else
%            txt = sprintf( 'Unknown body name %s', name(i,:) );
%            disp(txt)
%         end
%
%      end
%
%   MATLAB outputs:
%
%      Scalar:
%      Body ID 607 maps to name Hyperion
%
%      Vector:
%      Body ID 801 maps to name Triton
%      Body ID 601 maps to name Mimas
%      Body ID 704 maps to name Oberon
%      Body ID 504 maps to name Callisto
%      Unknown body name Halo
%
%-Particulars
%
%   A sister version of this routine exists named cspice_bodn2c that returns
%   the structure field data as separate arguments.
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
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   body name to code
%
%-&

function [code] = mice_bodn2c(name)

   switch nargin
      case 1

         name = zzmice_str(name);

      otherwise

         error ( 'Usage: [_code_] = mice_bodn2c(_`name`_)' )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [code] = mice('bodn2c_s',name);
   catch
      rethrow(lasterror)
   end


