%-Abstract
%
%   CSPICE_KINFO returns information about a loaded kernel
%   specified by name.
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
%      file   the scalar string name of a kernel file for which descriptive
%             information is desired.
%
%   the call:
%
%      [ filtyp, source, handle, found] = cspice_kinfo( file)
%
%   returns:
%
%      filtyp   the scalar string type name of the kernel specified by 'file'.
%               'filtyp' will be empty if file is not on the list of kernels
%               loaded via cspice_furnsh.
%
%      source   the scalar string name of the source file used to
%               specify 'file' as one to load.  If 'file' was loaded
%               directly via a call to cspice_furnsh, 'source' will be empty.
%               If file is not on the list of kernels loaded via
%               cspice_furnsh, 'source' will be empty.
%
%      handle   the integer handle attached to 'file' if it is a binary
%               kernel.  If file is a text kernel or meta-text kernel
%               handle will be zero. If file is not on the list of
%               kernels loaded via cspice_furnsh, 'handle' has value zero.
%
%      found    returns true if the specified file exists.
%               If there is no such file, 'found' will be set to
%               false.
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
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( '/kernels/gen/lsk/naif0009.tls'
%                                '/kernels/gen/spk/de421.bsp'
%                                '/kernels/gen/pck/pck00009.tpc'
%                      )
%
%         \begintext
%
%      %
%      % Load a meta kernel listing a path to an SPK file.
%      %
%      cspice_kclear
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Use cspice_kinfo to ensure the kernel system loaded
%      % the SPK file of interest.
%      %
%      file = '/kernels/gen/spk/de421.bsp';
%
%      [ filtyp, source, handle, found ] = cspice_kinfo( file );
%
%      %
%      % Take appropriate action depending on the returned
%      % state of found. If found has value false, then
%      % 'file' is not loaded.
%      %
%      if ( found )
%         disp( [ 'File type: ' filtyp ] )
%         disp( [ 'Source   : ' source ] )
%      else
%         disp( [ 'Kernel not loaded: ' file ] )
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
%      File type: SPK
%      Source   : standard.tm
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine kinfo_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.2.0, 10-MAY-2011, EDW (JPL)
%
%      "logical" call replaced with "zzmice_logical."
%
%   -Mice Version 1.0.1, 06-MAY-2009, EDW (JPL)
%
%      Added MICE.REQ reference to the Required Reading section.
%
%   -Mice Version 1.0.0, 01-DEC-2006, EDW (JPL)
%
%-Index_Entries
%
%   Fetch information about a loaded SPICE kernel
%
%-&

function [ filtyp, source, handle, found] = cspice_kinfo(file)

   switch nargin
      case 1

         file = zzmice_str(file);

      otherwise

         error( [ 'Usage: [ `filtyp`, `source`, handle, found ] = ' ...
                                             'cspice_kinfo( `file` )']  )

   end

   %
   % Call the MEX library.
   %
   try
      [ filtyp, source, handle, found]  = mice('kinfo_c', file);

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch
      rethrow(lasterror)
   end

