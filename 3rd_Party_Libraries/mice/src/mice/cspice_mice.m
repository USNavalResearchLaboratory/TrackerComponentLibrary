%-Abstract
%
%   CSPICE_MICE returns compile time, date, and version information
%   on the Mice shared object library.
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
%      item   a scalar string identifying the compile value to return.
%             Allowed values for 'item' :
%
%                'DATE'    - the date of the Mice interface compile
%
%                'TIME'    - the time on the date of the Mice
%                            interface compile
%
%                'VERSION' - the Mice version string
%
%   the call:
%
%      value = cspice_mice( item )
%
%   returns:
%
%      value   a scalar string corresponding to the quantity identified by
%              'item'
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      >> disp( [ cspice_mice( 'version' ) ' compiled ' ...
%                 cspice_mice( 'date' )    ' '          ...
%                 cspice_mice( 'time' ) ] )
%
%   MATLAB outputs:
%
%      Mice 0.9.54 27-APR-2006 (EDW) compiled Jun 27 2006 08:48:24
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 26-JUN-2006, EDW (JPL)
%
%-Index_Entries
%
%   return compile time, date, and version
%
%-&

function [value] = cspice_mice( item )

   switch nargin
      case 1

         item = zzmice_str(item);

      otherwise

         error ( 'Usage: [`value`] = cspice_mice( `item` )' )

   end


   %
   % Call the MEX library.
   %
   try
      [value] = mice( 'cspice_mice', item );
   catch
      rethrow(lasterror)
   end



