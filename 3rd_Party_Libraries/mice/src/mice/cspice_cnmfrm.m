%-Abstract
%
%   CSPICE_CNMFRM retrieves the ID code and name of the preferred
%   frame associated with a given body name.
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
%      cname   the name(s) of an object for which a preferred reference frame
%              exists.
%
%              [n,m] = size(cname); char = class(cname)
%
%                  or
%
%              [1,n] = size(cname); cell = class(cname)
%
%   the call:
%
%      [ frcode, frname, found] = cspice_cnmfrm(cname)
%
%   returns:
%
%      frcode   the SPICE frame code(s) associated with 'cname'.
%
%               [1,n] = size(frcode); int32 = class(frcode)
%
%      frname   the name(s) corresponding to 'frcode'.
%
%               [n,m] = size(frname); char = class(frname)
%
%      found    the flag(s) indicating if the appropriate frame ID-code and
%               frame name can be determined from 'cname'.
%
%               [1,n] = size(found); logical = class(found)
%
%               'frcode', 'frname', and 'found' return with the same
%               vectorization measure (N) as 'cname'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Return the body frame code and name for Io.
%      %
%
%      [ frcode, frname, found ] = cspice_cnmfrm( 'IO' );
%
%      if ( found )
%         fprintf( '%d  %s\n', frcode, frname )
%      end
%
%
%      %
%      % Return the body frame code and name for a vector of body names.
%      %
%
%      bodies = {'EARTH', 'MOON', 'HALO_DELTA'};
%
%      [ frcode, frname, found ] = cspice_cnmfrm( bodies );
%
%      for i=1:numel( bodies)
%
%         if ( found(i) )
%            fprintf( '%d  %s\n', frcode(i), frname(i,:) )
%         else
%            fprintf( 'No frame associated with body %s\n', char(bodies(i)) )
%         end
%
%      end
%
%   MATLAB outputs:
%
%      10023  IAU_IO
%      10013  IAU_EARTH
%      10020  IAU_MOON
%      No frame associated with body HALO_DELTA
%
%-Particulars
%
%   This routine allows the user to determine the frame that should
%   be associated with a particular object. For example, if you
%   need the frame name and ID associated with Io, you can call
%   cspice_cnmfrm to return these values.
%
%   The preferred frame to use with an object is specified via one
%   of the kernel pool variables:
%
%       OBJECT_<cname>_FRAME
%
%   where <cname> is the NAIF ID or name of the object.
%
%   For those objects that have "built-in" frame names this
%   routine returns the corresponding "IAU" frame and frame ID code.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine cnmfrm_c.
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

function [ frcode, frname, found] = cspice_cnmfrm(cname)

   switch nargin
      case 1

         cname = zzmice_str(cname);

      otherwise

         error ( ['Usage: [_frcode_, _`frname`_, _found_]' ...
                                             ' = cspice_cnmfrm(_`cname`_)'] )

   end


   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      cnmfrm = mice( 'cnmfrm_s', cname ) ;
      frcode = reshape( [cnmfrm(:).code],  1, [] );
      frname = char( cnmfrm.name );
      found  = reshape( [cnmfrm(:).found], 1, [] );
   catch
      rethrow(lasterror)
   end


