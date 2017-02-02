%-Abstract
%
%   CSPICE_DAFAC adds comments from a buffer of character strings to the
%   comment area of a binary DAF file, appending them to any comments which
%   are already present in the file's comment area.
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
%      handle   file handle referring to a DAF.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%      buffer   vector containing comments which to write into
%               the comment area of the binary DAF attached to 'handle'.
%
%               Each element of 'buffer' should contain one comment line.
%
%               [n,c1] = size(buffer); char = class(buffer)
%
%                  or
%
%               [1,n] = size(buffer); cell = class(buffer)
%
%   the call:
%
%      cspice_dafac( handle, buffer )
%
%   returns:
%
%      The call adds the contents of 'buffer' to the DAF referred
%      to by 'handle'.
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
%      % Create a comment text block. Use the list
%      % of bodies in an SPK for this example.
%      %
%      comments = {                                                       ...
%      '-9',                     'DEIMOS (402)',   'TITANIA (703)',       ...
%      'MERCURY BARYCENTER (1)', 'MARS (499)',     'OBERON (704)',        ...
%      'VENUS BARYCENTER (2)',   'IO (501)',       'MIRANDA (705)',       ...
%      'EARTH BARYCENTER (3)',   'EUROPA (502)',   'URANUS (799)',        ...
%      'MARS BARYCENTER (4)',    'GANYMEDE (503)', 'TRITON (801)',        ...
%      'JUPITER BARYCENTER (5)', 'CALLISTO (504)', 'NEREID (802)',        ...
%      'SATURN BARYCENTER (6)',  'JUPITER (599)',  'NEPTUNE (899)',       ...
%      'URANUS BARYCENTER (7)',  'TETHYS (603)',   'CHARON (901)',        ...
%      'NEPTUNE BARYCENTER (8)', 'DIONE (604)',    'PLUTO (999)',         ...
%      'PLUTO BARYCENTER (9)',   'RHEA (605)',     '301001*',             ...
%      'SUN (10)',               'TITAN (606)',    'GOLDSTONE (399001)*', ...
%      'MERCURY (199)',          'HYPERION (607)', 'CANBERRA (399002)*',  ...
%      'VENUS (299)',            'IAPETUS (608)',  'MADRID (399003)*',    ...
%      'MOON (301)',             'SATURN (699)',   '401001*',             ...
%      'EARTH (399)',            'ARIEL (701)',                           ...
%      'PHOBOS (401)',           'UMBRIEL (702)' };
%
%      %
%      % Define the SPK file to which to add the 'comments' text.
%      %
%      SPK = 'test.spk';
%
%      %
%      % Open the 'SPK' for writing; return the corresponding
%      % file handle to 'handle'.
%      %
%      handle = cspice_dafopw( SPK );
%
%      %
%      % Add the comments to the 'SPK', use a default line length
%      % of 80 characters.
%      %
%      cspice_dafac( handle, comments )
%
%      %
%      % SAFELY close the file.
%      %
%      cspice_dafcls( handle )
%
%   Matlab outputs:
%
%      None.
%
%   Assuming 'SPK' originally lacked comments, the file now
%   contains the comments:
%
%      -9
%      DEIMOS (402)
%      TITANIA (703)
%      MERCURY BARYCENTER (1)
%      MARS (499)
%      OBERON (704)
%      VENUS BARYCENTER (2)
%      IO (501)
%      MIRANDA (705)
%      EARTH BARYCENTER (3)
%      EUROPA (502)
%      URANUS (799)
%      MARS BARYCENTER (4)
%      GANYMEDE (503)
%      TRITON (801)
%      JUPITER BARYCENTER (5)
%      CALLISTO (504)
%      NEREID (802)
%      SATURN BARYCENTER (6)
%      JUPITER (599)
%      NEPTUNE (899)
%      URANUS BARYCENTER (7)
%      TETHYS (603)
%      CHARON (901)
%      NEPTUNE BARYCENTER (8)
%      DIONE (604)
%      PLUTO (999)
%      PLUTO BARYCENTER (9)
%      RHEA (605)
%      301001*
%      SUN (10)
%      TITAN (606)
%      GOLDSTONE (399001)*
%      MERCURY (199)
%      HYPERION (607)
%      CANBERRA (399002)*
%      VENUS (299)
%      IAPETUS (608)
%      MADRID (399003)*
%      MOON (301)
%      SATURN (699)
%      401001*
%      EARTH (399)
%      ARIEL (701)
%      PHOBOS (401)
%      UMBRIEL (702)
%
%   If 'SPK' contained comments before running the program, the comments
%   defined in 'comments' are appended to the existing comments.
%
%-Particulars
%
%   A binary DAF contains a data area which is reserved for storing
%   annotations or descriptive textual information about the data
%   contained in a file. This area is referred to as the ``comment
%   area'' of the file. The comment area of a DAF is a line oriented
%   medium for storing textual information. The comment area preserves
%   leading or embedded white space in the line(s) of text which are
%   stored so that the appearance of the information will be unchanged
%   when it is retrieved (extracted) at some other time. Trailing
%   blanks, however, are NOT preserved, due to the way that character
%   strings are represented in standard Fortran 77.
%
%   This routine will take a buffer of text lines and add (append) them
%   to the comment area of a binary DAF. If there are no comments in the
%   comment area of the file, then space will be allocated and the text
%   lines in `buffer' will be placed into the comment area. The text lines
%   may contain only printable ASCII characters (decimal values 32 -
%   126).
%
%   There is NO maximum length imposed on the significant portion of a
%   text line that may be placed into the comment area of a DAF. The
%   maximum length of a line stored in the comment area should be
%   reasonable, however, so that they may be easily extracted. A good
%   maximum value for this would be 255 characters, as this can easily
%   accommodate ``screen width'' lines as well as long lines which may
%   contain some other form of information.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dafac_c.
%
%   MICE.REQ
%   DAF.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 20-JUL-2012, EDW (JPL)
%
%-Index_Entries
%
%   add comments to a binary DAF
%   append comments to a DAF comment area
%
%-&

function cspice_dafac( handle, buffer )

   switch nargin
      case 2

         handle  = zzmice_int(handle);
         buffer  = zzmice_str(buffer);

      otherwise

         error ( 'Usage: cspice_dafac( handle, buffer )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'dafac_c', handle, buffer );
   catch
      rethrow(lasterror)
   end




