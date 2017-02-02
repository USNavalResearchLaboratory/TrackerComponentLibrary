%-Abstract
%
%   CSPICE_DAFGN returns the name for current array in the current
%   DAF being searched
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
%      None.
%
%   the call:
%
%      name = cspice_dafgn
%
%   returns:
%
%      name     the name of the current DAF array - that array
%               found by a previous call to cspice_daffna or cspice_daffpa.
%
%               [1,c1] = size(name); char = class(name)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%
%      %
%      % Define a DAF from which to read the name for each array in
%      % the DAF.
%      %
%      DAF = 'daftest.bsp';
%      NI  = 6;
%      ND  = 2;
%
%      %
%      % Open the DAF for read
%      %
%      handle = cspice_dafopr( DAF );
%
%      %
%      % Begin a forward search on 'DAF'.
%      %
%      cspice_dafbfs( handle )
%      found = cspice_daffna;
%
%      %
%      % Loop while found
%      %
%      while ( found )
%
%         [dc, ic] = cspice_dafgs( ND, NI );
%         name = cspice_dafgn;
%
%         %
%         % Output each array name.
%         %
%         fprintf( '%s\n', name)
%
%         %
%         % Check for a next segment.
%         %
%         found = cspice_daffna;
%
%      end
%
%      %
%      % SAFELY close the file.
%      %
%      cspice_dafcls( handle )
%
%   MATLAB outputs:
%
%      PHOENIX SPACECRAFT
%      MERCURY BARYCENTER
%      VENUS BARYCENTER
%      EARTH BARYCENTER
%
%         ...
%
%      CANBERRA
%      MADRID
%      PHOBOS BASECAMP
%      TRANQUILITY BASE
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dafgn_c.
%
%   MICE.REQ
%   DAF.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 11-JUN-2013, EDW (JPL)
%
%-Index_Entries
%
%   get DAF array name
%
%-&

function [name] = cspice_dafgn

   switch nargin
      case 0
         ;
      otherwise

         error ( 'Usage: [`name`] = cspice_dafgn' )

   end


   %
   % Call the MEX library.
   %
   try
      [name] = mice( 'dafgn_c' );
   catch
      rethrow(lasterror)
   end



