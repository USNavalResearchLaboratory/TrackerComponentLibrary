%-Abstract
%
%   CSPICE_FRINFO retrieves the minimal attributes associated with a
%   frame needed for converting transformations to and from it.
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
%     frcode   a SPICE ID for some reference frame.
%
%              [1,n] = size(frcode); int32 = class(frcode)
%
%   the call:
%
%      [cent, clss, clssid, found] = cspice_frinfo( frcode )
%
%   returns:
%
%      cent     the SPICE body ID for the center of the reference frame
%               (if such an ID is appropriate).
%
%               [1,n] = size(cent); int32 = class(cent)
%
%      clss     the class ID or type of the frame. This identifies which
%               subsystem will perform frame transformations.
%
%               [1,n] = size(clss); int32 = class(clss)
%
%      clssid   the ID used for the frame within its class. This may be
%               different from the frame ID.
%
%               [1,n] = size(clssid); int32 = class(clssid)
%
%      found    flag returning true if 'cent', 'frclss' and 'frcode' are
%               available, false if not.
%
%               [1,n] = size(found); logical = class(found)
%
%               'cent', 'clss', 'clssid' and 'found' return with the same
%               vectorization measure (N) as 'frcode'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Retrieve frame information for a scalar 'code'.
%      %
%      disp('Scalar' )
%      code = 13000;
%
%      [cent, clss, clssid, found] = cspice_frinfo( code );
%      fprintf(' code   center  class  class_ID  found\n' );
%      fprintf( '%d    %d      %d     %d        %d\n', ...
%               code, cent, clss, clssid, int32(found) );
%
%      %
%      % Retrieve frame information for a vector of 'codes'.
%      %
%      disp('Vector' )
%      codes = [1:5];
%
%      [cent, clss, clssid, found] = cspice_frinfo( codes );
%
%      fprintf( 'code center class class_ID found\n')
%      fprintf( '%d    %d      %d     %d        %d\n', ...
%             [codes; cent; clss; clssid; int32(found) ] );
%
%   MATLAB outputs:
%
%      Scalar
%       code   center  class  class_ID  found
%      13000    399      2     3000        1
%
%      Vector
%      code center class class_ID found
%      1    0      1     1        1
%      2    0      1     2        1
%      3    0      1     3        1
%      4    0      1     4        1
%      5    0      1     5        1
%
%-Particulars
%
%   This is a low level routine needed by state transformation
%   software to transform states and attitudes between different
%   reference frames.
%
%   The routine first examines local "hard-coded" information about
%   reference frames to see if the requested frame belongs to this
%   set.  If it does that information is returned.
%
%   If the requested information is not stored locally, the routine
%   then examines the kernel pool to see if the requested information
%   is stored there.  If it is and has the expected format, the data
%   is retrieved and returned.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine frinfo_c.
%
%   MICE.REQ
%   FRAMES.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   fetch reference frame attributes
%
%-&

function [cent, clss, clssid, found] = cspice_frinfo( frcode )

   switch nargin
      case 1

         frcode = zzmice_int(frcode);

      otherwise

         error( ['Usage: [_cent_, _clss_, _clssid_, _found_] = ' ...
                                           'cspice_frinfo(_frcode_)'] )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [frinfo] = mice('frinfo_s', frcode);

      cent   = reshape( [frinfo.center],   1, [] );
      clss   = reshape( [frinfo.class],    1, [] );
      clssid = reshape( [frinfo.class_ID], 1, [] );
      found  = reshape( [frinfo.found],    1, [] );

   catch
      rethrow(lasterror)
   end




