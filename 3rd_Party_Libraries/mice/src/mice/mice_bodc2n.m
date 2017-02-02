%-Abstract
%
%   MICE_BODC2N returns the body name corresponding to an input numeric
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
%      code   a scalar integer SPICE code or 1XN array of integer SPICE
%             codes for a set of bodies: planets, satellites, barycenters,
%             DSN stations, spacecraft, asteroids, comets, or other
%             ephemeris object
%
%   the call:
%
%      ID = mice_bodc2n( code )
%
%   returns:
%
%      ID   the scalar or 1xN array of structures associating
%           a body name with a corresponding NAIF ID. Each structure
%           contains two fields:
%
%              name   the "name" of a particular body. If a mapping
%                     does not exist, the 'name' field returns as NULL
%
%              code   a scalar integer SPICE code assigned either
%                     by SPICE or the user to 'name'. If a mapping
%                     does not exist, the 'code' field returns as 0.
%
%              found  a scalar boolean indicating if a mapping exists (true)
%
%              'ID' returns with the same vectorization measure (N) as 'code'.
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
%      ID      = mice_bodc2n( naif_id );
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
%      disp( ' ' )
%
%      %
%      % Create an array of IDs. Include one unknown ID.
%      %
%      disp( 'Vector:' )
%      naif_id = [ 502, 503, 504, 505, 5006 ];
%      ID      = mice_bodc2n( naif_id );
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
%         if( ID(i).found )
%            txt = sprintf( 'Body ID %i maps to name %s', ...
%                            ID(i).code, ID(i).name );
%            disp(txt)
%         else
%            txt = sprintf( 'Unknown body ID %i', naif_id(i) );
%            disp(txt)
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
%   A sister version of this routine exists named cspice_bodc2n that returns
%   the structure field data as separate arguments.
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
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   body id code to name
%
%-&

function [name] = mice_bodc2n(code)

   switch nargin
      case 1

         code = zzmice_int(code);

      otherwise

         error ( 'Usage: [_`name`_] = mice_bodc2n(_code_)' )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [name] = mice('bodc2n_s', code);
   catch
      rethrow(lasterror)
   end



