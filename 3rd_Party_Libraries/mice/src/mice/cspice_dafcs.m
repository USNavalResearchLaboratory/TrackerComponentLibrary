%-Abstract
%
%   CSPICE_DAFCS sets the active DAF to search. A search must be
%   in progress for the DAF.
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
%      handle   the file handle referring to a DAF to
%               set as the "active" file for a search.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%   the call:
%
%      cspice_dafcs( handle )
%
%   causes all DAF search activity apply to the file
%   referred to by 'handle'.
%
%   returns:
%
%      None.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example(1):
%
%      %
%      % Define two SPK test files.
%      %
%      SPK1 = 'test.bsp';
%      SPK2 = 'test8.bsp';
%
%      %
%      % Open the DAFs for read
%      %
%      han1 = cspice_dafopr( SPK1 );
%      han2 = cspice_dafopr( SPK2 );
%
%      %
%      % Begin a forward search on SPK1
%      %
%      cspice_dafbfs( han1 )
%      found = cspice_daffna;
%
%      %
%      % Begin a backwards search on SPK2
%      %
%      cspice_dafbbs( han2 )
%      found2 = cspice_daffpa;
%
%      %
%      % Reinstitute the search on han1, loop
%      % so long as segment data are found.
%      %
%      cspice_dafcs( han1 )
%
%      while ( found )
%
%         segid    = cspice_dafgn;
%         found    = cspice_daffna;
%
%         %
%         % Output each segment ID.
%         %
%         fprintf( '%s\n', segid )
%
%      end
%
%      %
%      % Close the files.
%      %
%      cspice_dafcls( han1 )
%      cspice_dafcls( han2 )
%
%   Matlab outputs:
%
%      PHOENIX SPACECRAFT
%      MERCURY BARYCENTER
%      VENUS BARYCENTER
%      EARTH BARYCENTER
%      MARS BARYCENTER
%      JUPITER BARYCENTER
%      SATURN BARYCENTER
%      URANUS BARYCENTER
%      NEPTUNE BARYCENTER
%      PLUTO BARYCENTER
%      MOON
%      PHOBOS
%      DEIMOS
%      IO
%      EUROPA
%      GANYMEDE
%      CALLISTO
%      TETHYS
%      DIONE
%      RHEA
%      TITAN
%      HYPERION
%      IAPETUS
%      ARIEL
%      UMBRIEL
%      TITANIA
%      OBERON
%      MIRANDA
%      TRITON
%      NERIED
%      CHARON
%      MERCURY
%      VENUS
%      EARTH
%      MARS
%      JUPITER
%      SATURN
%      URANUS
%      NEPTUNE
%      PLUTO
%      SUN
%      GOLDSTONE
%      CANBERRA
%      MADRID
%      PHOBOS BASECAMP
%      TRANQUILITY BASE
%
%   Example(2), switch the definitions for SPK1 and SPK2:
%
%      %
%      % Define two SPK test files.
%      %
%      SPK2 = 'test.bsp';
%      SPK1 = 'test8.bsp';
%
%         ... remainder of example unchanged ..
%
%   Matlab outputs:
%
%      SPK type 8 test segment
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dafcs_c.
%
%   MICE.REQ
%   DAF.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 10-JUL-2012, EDW (JPL)
%
%-Index_Entries
%
%   select a DAF to continue searching
%
%-&

function cspice_dafcs( handle )

   switch nargin
      case 1

         handle  = zzmice_int(handle);

      otherwise

         error ( 'Usage: cspice_dafcs(handle)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'dafcs_c', handle );
   catch
      rethrow(lasterror)
   end




