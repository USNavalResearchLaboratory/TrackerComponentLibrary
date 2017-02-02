%-Abstract
%
%   CSPICE_DAFUS unpacks an array summary into its double precision and
%   integer components.
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
%      sum   an containing the DAF array summary. This identifies the
%            contents and location of a single array within a DAF.
%
%            [1,n] = size(sum); double = class(sum)
%
%      nd    the size of the return double precision array.
%
%            [1,1] = size(nd); int32 = class(nd)
%
%      ni    the size of the return integer array.
%
%            [1,1] = size(ni); int32 = class(ni)
%
%      For an SPK file, 'nd' always equals 2, 'ni' always equals 6. The precise
%      contents of the vectors depend on the type of DAF but the
%      final two elements of the 'ic' (integer) vector always contains the
%      initial and final addresses respectively of the array.
%
%   the call:
%
%      [dc, ic] = cspice_dafus( sum, nd, ni )
%
%   returns:
%
%      dc   the array of double precision components of the summary.
%
%           [1,nd] = size(dc); double = class(dc)
%
%      ic   the array of integer components of the summary.
%
%           [1,ni] = size(ic); int32 = class(ic)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Local constants
%      %
%      META   =  'standard.tm';
%      ND     =  2;
%      NI     =  6;
%
%      %
%      % Load a meta-kernel that specifies a planetary SPK file
%      % and leapseconds kernel. The contents of this meta-kernel
%      % are displayed above.
%      %
%      cspice_furnsh( META )
%
%      %
%      % Get the NAIF ID code for the Pluto system barycenter.
%      % This is a built-in ID code, so something's seriously
%      % wrong if we can't find the code.
%      %
%      [idcode, found] = cspice_bodn2c( 'PLUTO BARYCENTER' );
%
%      if ~found
%         cspice_kclear
%         error( 'SPICE(BUG)' )
%      end
%
%      %
%      % Pick a request time; convert to seconds past J2000 TDB.
%      %
%      reqtim = '2011 FEB 18 UTC';
%
%      et = cspice_str2et( reqtim );
%
%      %
%      % Find a loaded segment for the specified body and time.
%      %
%
%      [handle, descr, segid, found] = cspice_spksfs( idcode, et );
%
%      if ~found
%         cspice_kclear
%         txt = sprintf( 'No descriptor found for the body %d at time %s', ...
%                         idcode, et );
%         error( txt )
%      else
%
%         %
%         % Display the DAF file handle.
%         %
%         fprintf( 'DAF handle:    %d\n', handle )
%
%         %
%         % Display the segment ID.
%         %
%         %
%         % Unpack the descriptor. Display the contents.
%         %
%         [dc, ic] = cspice_dafus( descr, ND, NI );
%
%         fprintf( 'Segment ID:       %s\n', segid )
%         fprintf( 'Body ID code:     %d\n', ic(1) )
%         fprintf( 'Center ID code:   %d\n', ic(2) )
%         fprintf( 'Frame ID code:    %d\n', ic(3) )
%         fprintf( 'SPK data type:    %d\n', ic(4) )
%         fprintf( 'Start time (TDB): %f\n', dc(1) )
%         fprintf( 'Stop time  (TDB): %f\n', dc(2) )
%
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in IDL due to data persistence.
%      %
%      cspice_kclear
%
%   Matlab outputs:
%
%      DAF handle:    5
%      Segment ID:       DE-0421LE-0421
%      Body ID code:     9
%      Center ID code:   0
%      Frame ID code:    1
%      SPK data type:    2
%      Start time (TDB): -3169195200.000000
%      Stop time  (TDB): 1696852800.000000
%
%-Particulars
%
%   The components of array summaries are packed into double
%   precision arrays.
%
%   The total size of the summary is
%
%           (ni - 1)
%      nd + -------- + 1
%               2
%
%   double precision words (where nd, ni are nonnegative).
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dafus_c.
%
%   MICE.REQ
%   DAF.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 29-OCT-2012, EDW (JPL)
%
%-Index_Entries
%
%   get DAF summary
%
%-&

function [dc, ic] = cspice_dafus( sum, nd, ni)

   switch nargin
      case 3

         sum = zzmice_dp(sum);
         nd  = zzmice_int(nd);
         ni  = zzmice_int(ni);

      otherwise

         error ( 'Usage: [dc(nd), ic(nd)] = cspice_dafus( sum(), nd, ni)' )

   end

   %
   % Call the MEX library.
   %
   try
      [dc, ic] = mice( 'dafus_c', sum, nd, ni );
   catch
      rethrow(lasterror)
   end

