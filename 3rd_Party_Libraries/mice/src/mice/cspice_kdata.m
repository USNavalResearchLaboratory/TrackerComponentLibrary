%-Abstract
%
%   CSPICE_KDATA returns data for the nth kernel among a list of specified
%   kernel types.
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
%      which   index of the kernel to fetch (matching the type specified by
%              kind) from the list of kernels loaded using cspice_furnsh but
%              not unloaded using cspice_unload or cleared by cspice_klear.
%
%              [1,1] = size(which); int32 = class(which)
%
%              The range of 'which' is 1 to 'count', where 'count' is
%              the number of kernels loaded via cspice_furnsh. Retrieve
%              this value from a cspice_ktotal call.  See the
%              Examples section for an illustrative code fragment.
%
%     kind     list of types of kernels to consider when fetching kernels from
%              the list of loaded kernels. 'kind' should consist of a list of
%              words of kernels to examine. Recognized types are
%
%              [1,m] = size(kind); char = class(kind)
%
%                 SPK
%                 CK
%                 PCK
%                 EK
%                 TEXT
%                 META
%                 ALL
%
%              'kind' lacks case sensitivity. The cspice_kdata algorithm
%              ignores words in 'kind' if not one of those listed above.
%
%              See the routine cspice_ktotal for example use 'kind'.
%
%   the call:
%
%      [file, filtyp, source, handle, found] = cspice_kdata(which, kind)
%
%   returns:
%
%      file     name of the file having index 'which' in the sequence of files
%               of type 'kind' currently loaded via cspice_furnsh. 'file'
%               returns empty if no loaded kernels match the specification of
%               'which' and 'kind'.
%
%               [1,m] = size(file); char = class(file)
%
%      filtyp   name of the type of kernel specified by 'file'. 'filtyp' will
%               be empty if no loaded kernels match the specification of
%               'which' and 'kind'.
%
%               [1,m] = size(filtyp); char = class(filtyp)
%
%      source   name of the source file used to specify file as one to load. If
%               file was loaded directly via a call to cspice_furnsh, 'source'
%               will be empty. If there is no file matching the specification
%               of which and kind, 'source' will be empty.
%
%               [1,m] = size(source); char = class(source)
%
%      handle   handle attached to file if it is a binary kernel.  If file is a
%               text kernel or meta-text kernel handle will be zero.  If there
%               is no file matching the specification of 'which' and 'kind',
%               'handle' will be set to zero.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%      found    flag indicating if a file matching the specification of 'which'
%               and 'kind' exists. If there no such file exists, 'found'
%               returns false (if 'found' returns as false, all return strings
%               are empty, not null).
%
%               [1,1] = size(found); logical = class(found)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%         KPL/MK
%
%         File name: standard.tm
%
%         This meta-kernel is intended to support operation of SPICE
%         example programs. The kernels shown here should not be
%         assumed to contain adequate or correct versions of data
%         required by SPICE-based user applications.
%
%         In order for an application to use this meta-kernel, the
%         kernels referenced here must be present in the user's
%         current working directory.
%
%         The names and contents of the kernels referenced
%         by this meta-kernel are as follows:
%
%            File name                     Contents
%            ---------                     --------
%            de421.bsp                     Planetary ephemeris
%            pck00009.tpc                  Planet orientation and
%                                          radii
%            naif0009.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00009.tpc',
%                                'naif0009.tls'  )
%
%         \begintext
%
%   Example:
%
%      %
%      % Load several kernel files.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Count the number of loaded kernel files.
%      %
%      count = cspice_ktotal( 'ALL' );
%
%      %
%      % Loop over the count, outputting file information as we loop.
%      % The loop tells us all files loaded via cspice_furnsh, their
%      % type, and how they were loaded.
%      %
%      for i = 1:count+1
%
%         [ file, type, source, handle, found ] = ...
%                                          cspice_kdata( i, 'ALL');
%
%         if ( found )
%            fprintf( 'Index : %d\n', i     );
%            fprintf( 'File  : %s\n', file  );
%            fprintf( 'Type  : %s\n', type  );
%            fprintf( 'Source: %s\n\n', source);
%
%         else
%
%            fprintf( 'No kernel found with index: %d\n', i );
%
%         end
%
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in Mice due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Index : 1
%      File  : standard.tm
%      Type  : META
%      Source:
%
%      Index : 2
%      File  : de421.bsp
%      Type  : SPK
%      Source: standard.tm
%
%      Index : 3
%      File  : pck00009.tpc
%      Type  : TEXT
%      Source: standard.tm
%
%      Index : 4
%      File  : naif0009.tls
%      Type  : TEXT
%      Source: standard.tm
%
%      No kernel found with index: 5
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine kdata_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.2.0, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%      "logical" call replaced with "zzmice_logical."
%
%      I/O descriptions edits to parallel to Icy version.
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%      Edits to Example section, proper description of "standard.tm"
%      meta kernel.
%
%   -Mice Version 1.0.1, 06-MAY-2009, EDW (JPL)
%
%      Added MICE.REQ reference to the Required Reading section.
%
%   -Mice Version 1.0.0, 30-MAR-2007, EDW (JPL)
%
%-Index_Entries
%
%   Retrieve information on loaded SPICE kernels
%
%-&

function [ file, filtyp, source, handle, found ] = cspice_kdata( which, kind )

   switch nargin
      case 2

         which = zzmice_int(which);
         kind  = zzmice_str(kind);

      otherwise

         error( [ 'Usage: [ `file`, `filtyp`, `source`, handle ] = ' ...
                                         'cspice_kdata( which, `kind` )']  )

   end

   %
   % Call the MEX library.
   %
   try
      [ file, filtyp, source, handle, found ]  = mice('kdata_c', which, kind);

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch
      rethrow(lasterror)
   end


