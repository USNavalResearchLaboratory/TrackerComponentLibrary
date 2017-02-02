%-Abstract
%
%   CSPICE_DAFEC reads comment text from the comment area of a DAF.
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
%      handle   file handle referring to a DAF file.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%      bufsiz   the maximum number of comment lines to copy to 'buffer'.
%
%               [1,1] = size(bufsiz); int32 = class(bufsiz)
%
%      lenout   allowed length of each string element of the output 'buffer'.
%               This length must be large enough to hold the longest output
%               string. The SPICE system imposes no limit on the length of
%               comment lines, so 'lenout' normally should be set to a
%               "generous" value that is unlikely to be exceeded.
%
%               [1,1] = size(lenout); int32 = class(lenout)
%
%   the call:
%
%      [buffer, done] = cspice_dafec( handle, bufsiz, lenout )
%
%   returns:
%
%      buffer   array containing the comment lines read from the DAF
%               associated with 'handle'.
%
%               On output, 'buffer' contains 'bufsiz' or less strings of comment
%               text, with one comment line per string ( bufsiz >= n).
%
%               [n,c1] = size(buffer); char = class(buffer)
%
%      done     logical indicating whether or not all of the comment
%               lines from the comment area of the DAF have been read. This
%               variable has value true after the last comment line has been
%               read. It will have a false value otherwise.
%
%               If no comments exist in the comment area, this variable
%               returns as true.
%
%               [1,1] = size(done); logical = class(done)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Define an SPK from which to extract the comment area.
%      %
%      SPK        = 'de421.bsp';
%
%      %
%      % Define the size of the comment area to read from the SPK.
%      % 15 lines, each with length 80 characters.
%      %
%      BUFSIZE    = 15;
%      LINLEN     = 80;
%
%      %
%      % Open the 'SPK' for reading, return the corresponding
%      % file handle to 'handle'.
%      %
%      handle = cspice_dafopr( SPK );
%
%      done = false;
%
%      [buf, done] = cspice_dafec( handle, BUFSIZE, LINLEN );
%      output = cellstr(buf);
%
%      for i=1:numel(output)
%         fprintf( '%s\n', char(output(i)) );
%      end
%
%      if done
%         fprintf( 'All comments read from file.\n' );
%      else
%         fprintf( 'Not all comments read from file.\n' );
%      end
%
%      %
%      % SAFELY close the file.
%      %
%      cspice_dafcls( handle )
%
%   Matlab outputs:
%
%      ; de421.bsp LOG FILE
%      ;
%      ; Created 2008-02-12/11:33:34.00.
%      ;
%      ; BEGIN NIOSPK COMMANDS
%
%      LEAPSECONDS_FILE    = naif0007.tls
%      SPK_FILE            = de421.bsp
%        SPK_LOG_FILE      = de421_spk_conversion.log
%        NOTE              = NIOSPK 6.1.0 Conversion
%        SOURCE_NIO_FILE   = de421.nio
%          BEGIN_TIME      = CAL-ET 1899 JUL 29 00:00:00.000
%          END_TIME        = CAL-ET 2053 OCT 09 00:00:00.000
%
%      ; END NIOSPK COMMANDS
%
%      Not all comments read from file.
%
%   The program outputs BUFSIZ (15) lines from the 'SPK' comment area.
%   Additional calls to cspice_dafec will read more comment lines
%   from the SPK in slices of BUFSIZ.
%
%   Reading all comment lines from 'SPK' requires a large value for BUFSIZ.
%   In this case, a BUFSIZ value of 50 will read all comment lines from
%   'SPK' in a single cspice_dafec.
%
%-Particulars
%
%   A binary DAF contains an area which is reserved for storing
%   annotations or descriptive textual information describing the data
%   contained in a file. This area is referred to as the ``comment
%   area'' of the file. The comment area of a DAF is a line
%   oriented medium for storing textual information. The comment
%   area preserves any leading or embedded white space in the line(s)
%   of text which are stored, so that the appearance of the of
%   information will be unchanged when it is retrieved (extracted) at
%   some other time. Trailing blanks, however, are NOT preserved,
%   due to the way that character strings are represented in
%   standard Fortran 77.
%
%   This routine will read the comments from the comment area of
%   a binary DAF, placing them into a line buffer. If the line
%   buffer is not large enough to hold the entire comment area,
%   the portion read will be returned to the caller, and the DONE
%   flag will be set to false. This allows the comment area to be
%   read in ``chunks,'' a buffer at a time. After all of the comment
%   lines have been read, the `done' flag will be set to SPICETRUE.
%
%   This routine can be used to ``simultaneously'' extract comments
%   from the comment areas of multiple binary DAFs.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dafec_c.
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
%    extract comments from a DAF
%
%-&

function [buffer, done] = cspice_dafec( handle, bufsiz, lenout )

   switch nargin
      case 3

         handle  = zzmice_int(handle);
         bufsiz  = zzmice_int(bufsiz);
         lenout  = zzmice_int(lenout);

      otherwise

         error ( ['Usage: [buffer, done] = ' ...
                          'cspice_dafec( handle, bufsiz, lenout )'] )

   end

   %
   % Call the MEX library.
   %
   try
      [buffer, done] = mice( 'dafec_c', handle, bufsiz, lenout );
   catch
      rethrow(lasterror)
      done = zzmice_logical(done);
   end




