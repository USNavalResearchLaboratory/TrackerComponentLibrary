%-Abstract
%
%   CSPICE_TKVRSN returns the latest version string given an item such
%   as the Toolkit or an entry point name
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
%      item   is the item for which a version string is to be
%             returned.
%
%             Currently, the only item supported is "toolkit"
%             and it will return the toolkit version number.
%
%             Any other 'item' will return "No version found."
%
%   the call:
%
%      value = cspice_tkvrsn( item )
%
%   returns:
%
%      value   a scalar string corresponding to the latest version
%              string for the specified 'item'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      >> disp( cspice_tkvrsn('toolkit') )
%
%   MATLAB outputs:
%
%      CSPICE_N0060
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine tkvrsn_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 11-JUN-2013, EDW (JPL)
%
%       I/O descriptions edits to conform to Mice documentation format.
%
%   -Mice Version 1.0.0, 26-NOV-2006, EDW (JPL)
%
%-Index_Entries
%
%   Return version strings
%-&

function [value] = cspice_tkvrsn( item )

   switch nargin
      case 1

         item = zzmice_str(item);

      otherwise

         error ( 'Usage: [`value`] = cspice_tkvrsn( `item` )' )

   end

   %
   % Call the MEX library.
   %
   try
      [value] = mice( 'tkvrsn_c', item );
   catch
      rethrow(lasterror)
   end



