%-Abstract
%
%   CSPICE_CIDFRM retrieves the ID code and name of the preferred 
%   frame associated with a given body ID.
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
%      cent   NAIF ID of a body for which a preferred reference frame
%             exists.
%
%             [1,n] = size(cent); int32 = class(cent)
%
%   the call:
%
%      [ frcode, frname, found] = cspice_cidfrm(cent)
%
%   returns:
%
%
%      frcode   the SPICE frame code(s) associated with 'cent'.
%
%               [1,n] = size(frcode); int32 = class(frcode)
%
%      frname   the name(s) corresponding to 'frcode'.
%
%               [n,m] = size(frname); char = class(frname)
%
%      found    the flag(s) indicating if the appropriate frame ID-code and
%               frame name can be determined from 'cent'.
%
%               [1,n] = size(found); logical = class(found)
%
%               'frcode', 'frname', and 'found' return with the same
%               vectorization measure (N) as 'cent'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Retrieve frame information for body ID 501.
%      %
%      disp('Scalar' )
%
%      bodies = 501;
%      [ frcode, frname, found ] = cspice_cidfrm( bodies );
%
%      if ( found )
%         fprintf( '%8d  %8d  %s\n', bodies, frcode, frname )
%      end
%
%      %
%      % Retrieve frame information for a vector of body IDs.
%      %
%      disp('Vector' )
%
%      bodies = [ 0, 301, 401, -501 ];
%
%      [ frcode, frname, found ] = cspice_cidfrm( bodies );
%
%      for i=1:numel( bodies)
%
%         if ( found(i) )
%            fprintf( '%8d  %8d  %s\n', bodies(i), frcode(i), frname(i,:) )
%         else
%            fprintf( 'No frame associated with ID %d\n', bodies(i) )
%         end
%
%      end
%
%   MATLAB outputs:
%
%      Scalar
%           501     10023  IAU_IO
%
%      Vector
%             0        15  DE-202
%           301     10020  IAU_MOON
%           401     10021  IAU_PHOBOS
%      No frame associated with ID -501
%
%-Particulars
%
%   This routine allows the user to determine the frame that should
%   be associated with a particular object. For example, if you
%   need the frame name and ID associated with Io, you can call 
%   cspice_cidfrm to return these values.
%
%   The preferred frame to use with an object is specified via one
%   of the kernel pool variables:
%
%       OBJECT_<cent>_FRAME
%
%   where <cent> is the NAIF ID or name of the object.
%
%   For those objects that have "built-in" frame names this
%   routine returns the corresponding "IAU" frame and frame ID code.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine cidfrm_c.
%
%   MICE.REQ
%   FRAMES.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 11-NOV-2013, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   Fetch reference frame attributes
%
%-&

function [ frcode, frname, found] = cspice_cidfrm(cent)

   switch nargin
      case 1

         cent = zzmice_int(cent);

      otherwise

         error ( ['Usage: [_frcode_, _`frname`_, _found_]' ...
                                             ' = cspice_cidfrm(_cent_)'] )

   end


   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      cidfrm = mice( 'cidfrm_s', cent ) ;
      frcode = reshape( [cidfrm(:).code],  1, [] );
      frname = char( cidfrm.name );
      found  = reshape( [cidfrm(:).found], 1, [] );
   catch
      rethrow(lasterror)
   end


