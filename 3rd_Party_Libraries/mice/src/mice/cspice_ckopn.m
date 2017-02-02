%-Abstract
%
%   CSPICE_CKOPN opens a new CK file, returning the handle
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
%      fname    a scalar string defining the name of the CK file
%               to open
%
%      ifname   a scalar string defining the descriptive internal
%               filename for the CK
%
%      ncomch   the scalar integer number of characters to
%               reserve for comments.
%
%   the call:
%
%      handle = cspice_ckopn( name, ifname, ncomch )
%
%   returns:
%
%      handle   a scalar integer file handle assigned to 'fname'
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Define needed parameters, a name for the CK, the
%      % file internal name, and the number of characters
%      % to reserve for a comment block.
%      %
%      CK1        = 'type1.bc';
%      IFNAME     = 'CK';
%      NCOMCH     = 10;
%
%      %
%      % Open a new kernel.
%      %
%       handle = cspice_ckopn( CK1, IFNAME, NCOMCH);
%
%         ... do some writes to the open CK file ...
%
%      %
%      % SAFELY close the file
%      %
%      cspice_ckcls( handle )
%
%-Particulars
%
%   A cspice_ckcls call should balance every cspice_ckopn
%   call.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine ckopn_c.
%
%   MICE.REQ
%   CK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   open a new ck file
%
%-&

function [handle] = cspice_ckopn( fname, ifname, ncomch )

   switch nargin
      case 3

         fname  = zzmice_str(fname);
         ifname = zzmice_str(ifname);
         ncomch = zzmice_int( ncomch );

      otherwise

         error ( 'Usage: [handle] = cspice_ckopn(`fname`, `ifname`, ncomch)' )

   end

   %
   % Call the MEX library.
   %
   try
      [handle] = mice( 'ckopn_c', fname, ifname, ncomch );
   catch
      rethrow(lasterror)
   end




