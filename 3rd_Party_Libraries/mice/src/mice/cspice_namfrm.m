%-Abstract
%
%   CSPICE_NAMFRM retrieves the SPICE frame ID code associated
%   with a frame name.
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
%      frname   the name of some reference frame (either inertial or
%               non-inertial).
%
%               [n,m] = size(frname); char = class(frname)
%
%                  of
%
%               [1,n] = size(frname); cell = class(frname)
%
%               Leading blanks in 'frname' are ignored as is character case.
%
%               Note that all legitimate frame names contain 32 or fewer
%               characters.
%
%   the call:
%
%      frcode = cspice_namfrm(frname)
%
%   returns:
%
%      frcode   the SPICE code(s) used for internal representation of the named
%               reference frame.
%
%               [1,n] = size(frcode); int32 = class(frcode)
%
%               If the name input through frname is not recognized, 'frcode'
%               will be returned with a value of zero.
%
%              'frcode' returns with the same vectorization measure (N)
%               as 'frname'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Retrieve frame information for a single frame.
%      %
%      disp('Scalar' )
%      name = 'ITRF93';
%
%      %
%      % Output the frame name corresponding to 'name'.
%      %
%      frcode = cspice_namfrm( name )
%
%
%      %
%      % Retrieve frame information for a vector of names.
%      % Create a vector of frame IDs, 1 to 5.
%      %
%      disp('Vector' )
%      codes = [1:5];
%
%      %
%      % Convert 'codes' to the corresponding frame name.
%      % 
%      names = cspice_frmnam( codes );
%
%      %
%      % Output the frame IDs corresponding to 'names'.
%      % The result should match the 'codes' vector'.
%      %
%      frcode = cspice_namfrm( names )
%
%  MATLAB outputs:
%
%      Scalar
%
%      frcode =
%
%             13000
%
%      Vector
%
%      frcode =
%
%                 1           2           3           4           5
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine namfrm_c.
%
%   MICE.REQ
%   FRAMES.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 14-NOV-2014, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   frame name to frame ID code translation
%
%-&

function [frcode] = cspice_namfrm(frname)

   switch nargin
      case 1

         frname = zzmice_str(frname);

      otherwise

         error ( 'Usage: [_frcode_] = cspice_namfrm(_`frname`_)' )

   end

   try
      [frcode] = mice('namfrm_c',frname);
   catch
      rethrow(lasterror)
   end


