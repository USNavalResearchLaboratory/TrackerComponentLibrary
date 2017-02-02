%-Abstract
%
%   CSPICE_DAFOPW opens a DAF for subsequent write requests.
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
%      fname   the string name of a DAF to open for write access.
%
%              [1,c1] = size(fname); char = class(fname)
%
%                 or
%
%              [1,1] = size(fname); cell = class(fname)
%
%   the call:
%
%      handle = cspice_dafopw( fname )
%
%   returns:
%
%      handle   the file handle used other DAF routines
%               to refer to 'fname'.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%   Use cspice_dafcls to close files opened by this routine.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      %    Given an SPK with a comment section:
%      %
%      %      Body list
%      %       1 MERCURY BARYCENTER
%      %       2 VENUS BARYCENTER
%      %       3 EARTH BARYCENTER
%      %       4 MARS BARYCENTER
%      %       5 JUPITER BARYCENTER
%      %       6 SATURN BARYCENTER
%      %       7 URANUS BARYCENTER
%      %       8 NEPTUNE BARYCENTER
%      %
%      %             ...
%
%      %
%      % Define the SPK file from which to remove the comments section.
%      %
%      SPK = 'test.spk';
%
%      %
%      % Open for writing the 'SPK', return the corresponding
%      % file handle to 'handle'.
%      %
%      handle = cspice_dafopw( SPK );
%
%      %
%      % Remove the comments section from the DAF referred to by 'handle'.
%      %
%      cspice_dafdc( handle )
%
%      %
%      % SAFELY close the file.
%      %
%      cspice_dafcls( handle )
%
%   Examine the 'SPK' comment after the cspice_dafdc call.
%
%      $ commnt -r test.spk
%
%      There were no comments in the file 'test.spk'.
%
%   Matlab outputs:
%
%      None.
%
%-Particulars
%
%   Most DAFs require only read access. If you do not need to
%   change the contents of a file, you should open it with cspice_dafopr.
%   Use cspice_dafopw when you need to
%
%      -- change (update) one or more summaries, names, or
%         arrays within a file; or
%
%      -- add new arrays to a file.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dafopw_c.
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
%   open DAF for write
%
%-&

function [handle] = cspice_dafopw( fname )

   switch nargin
      case 1

         fname  = zzmice_str(fname);

      otherwise

         error ( 'Usage: [handle] = cspice_dafopw(`fname`)' )

   end

   %
   % Call the MEX library.
   %
   try
      [handle] = mice( 'dafopw_c', fname );
   catch
      rethrow(lasterror)
   end




