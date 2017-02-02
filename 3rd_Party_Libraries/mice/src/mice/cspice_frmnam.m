%-Abstract
%
%   CSPICE_FRMNAM retrieves the name of a reference frame associated with a
%   SPICE frame ID code.
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
%     frcode   value defining a SPICE reference frame ID code.
%
%              [1,n] = size(frcode); int32 = class(frcode)
%
%   the call:
%
%      frmname = cspice_frmnam( frcode )
%
%   returns:
%
%      frmnam   the frame name corresponding to the 'frcode' code.
%
%               [n,m] = size(frmnam); char = class(frmnam)
%
%               If frcode is not recognized as the name of a known reference
%               frame, 'frname' will be returned as an empty string.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Retrieve frame information for a scalar code.
%      %
%      disp('Scalar' )
%      code = 13000;
%
%      %
%      % Output the frame name corresponding to 'code'.
%      %
%      frmname = cspice_frmnam( code )
%
%
%      %
%      % Retrieve frame information for a vector of codes.
%      %
%      disp('Vector' )
%      codes = [1:5];
%
%      %
%      % Output the frame names corresponding to the 'codes'.
%      %
%      frmname = cspice_frmnam( codes )
%
%   MATLAB outputs:
%
%      Scalar
%
%      frmname =
%
%      ITRF93
%
%      Vector
%
%      frmname =
%
%      J2000
%      B1950
%      FK4
%      DE-118
%      DE-96
%
%-Particulars
%
%   This routine retrieves the name of a reference frame associated
%   with a SPICE frame ID code.
%
%   The ID codes stored locally are scanned for a match with frcode.
%   If a match is found, the name stored locally will be returned
%   as the name for the frame.
%
%   If frcode is not a member of the list of internally stored
%   ID codes, the kernel pool will be examined to see if the
%   variable
%
%      FRAME_idcode_NAME
%
%   is present (where idcode is the decimal character equivalent
%   of frcode).  If the variable is located and it has both
%   character type and dimension 1, the string value of the
%   kernel pool variable is returned as the name of the reference
%   frame.
%
%   Note that because the local information is always examined
%   first and searches of the kernel pool are performed only
%   after exhausting local information, it is not possible to
%   override the local name for any reference frame that is
%   known by this routine.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine frmnam_c.
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
%   frame ID code to frame name translation
%
%-&

function [frmname] = cspice_frmnam( frcode )

   switch nargin
      case 1

         frcode = zzmice_int(frcode);

      otherwise

         error( 'Usage: [_`frmnam`_] = cspice_frmnam(_frcode_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [frmname] = mice('frmnam_c', frcode);
   catch
      rethrow(lasterror)
   end




