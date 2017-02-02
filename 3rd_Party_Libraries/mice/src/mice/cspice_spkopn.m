%-Abstract
%
%   CSPICE_SPKOPN opens a new SPK file, returning the handle
%   of the opened file.
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
%      fname    the name of the SPK file to open.
%
%               [1,c1] = size(fname); char = class(fname)
%
%                  or
%
%               [1,1] = size(fname); cell = class(fname)
%
%      ifname   the descriptive internal filename for the SPK.
%
%               [1,c2] = size(ifname); char = class(ifname)
%
%                  or
%
%               [1,1] = size(ifname); cell = class(ifname)
%
%      ncomch   the scalar integer number of characters to
%               reserve for comments.
%
%               [1,1] = size(ncomch); int32 = class(ncomch)
%
%   the call:
%
%      handle = cspice_spkopn( name, ifname, ncomch )
%
%   returns:
%
%      handle   the file handle assigned to 'fname'.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%
%-Particulars
%
%   A cspice_spkcls call should balance every cspice_spkopn
%   call.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine spkopn_c.
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
%   open a new spk file
%
%-&

function [handle] = cspice_spkopn( fname, ifname, ncomch )

   switch nargin
      case 3

         fname  = zzmice_str(fname);
         ifname = zzmice_str(ifname);
         ncomch = zzmice_int(ncomch);

      otherwise

         error ( 'Usage: [handle] = cspice_spkopn(`fname`, `ifname`, ncomch)' )

   end

   %
   % Call the MEX library.
   %
   try
      [handle] = mice( 'spkopn_c', fname, ifname, ncomch );
   catch
      rethrow(lasterror)
   end




