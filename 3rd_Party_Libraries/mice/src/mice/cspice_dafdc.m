%-Abstract
%
%   CSPICE_DAFDC deletes the entire comment area of a specified DAF file.
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
%      handle   the file handle referring to a DAF file opened with
%               write access. This handle refers to the DAF file from which
%               to delete the comment section.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%   the call:
%
%      cspice_dafdc( handle )
%
%   removes the comment area of the DAF file referred to by 'handle'.
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
%      % Remove the comment section from the DAF referred to by 'handle'.
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
%   A binary DAF contains an area which is reserved for storing
%   annotations or descriptive textual information about the data
%   contained in a file. This area is referred to as the ``comment
%   area'' of the file. The comment area of a DAF is a line oriented
%   medium for storing textual information. The comment area preserves
%   any leading or embedded white space in the line(s) of text which are
%   stored, so that the appearance of the of information will be
%   unchanged when it is retrieved (extracted) at some other time.
%   Trailing blanks, however, are NOT preserved, due to the way that
%   character strings are represented in standard Fortran 77.
%
%   This routine will delete the entire comment area from the binary DAF
%   attached to `handle'. The size of the binary DAF will remain
%   unchanged. The space that was used by the comment records is
%   reclaimed:  the data area of the DAF is shifted toward the beginning
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dafdc_c.
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
%    delete DAF comment area
%
%-&

function cspice_dafdc( handle )

   switch nargin
      case 1

         handle  = zzmice_int(handle);

      otherwise

         error ( 'Usage: cspice_dafdc(handle)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'dafdc_c', handle );
   catch
      rethrow(lasterror)
   end




