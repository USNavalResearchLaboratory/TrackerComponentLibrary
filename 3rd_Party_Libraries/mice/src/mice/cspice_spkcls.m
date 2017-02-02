%-Abstract
%
%   CSPICE_SPKCLS closes a SPK file opened for read or write.
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
%      handle   the file handle for an open SPK file.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%   the call:
%
%      cspice_spkcls( handle )
%
%   returns:
%
%   The routine closes the file indicated by 'handle'. The close operation
%   tests the file to ensure the presence of data segments.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%-Particulars
%
%   A cspice_spkcls call should balance every cspice_spkopn
%   call.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine spkcls_c.
%
%   MICE.REQ
%   SPK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 23-MAY-2012, EDW (JPL)
%
%-Index_Entries
%
%   close an spk file
%
%-&

function cspice_spkcls( handle)

   switch nargin
      case 1

         handle = zzmice_int( handle );

      otherwise

         error ( 'Usage: cspice_spkcls(handle)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'spkcls_c', handle);
   catch
      rethrow(lasterror)
   end


