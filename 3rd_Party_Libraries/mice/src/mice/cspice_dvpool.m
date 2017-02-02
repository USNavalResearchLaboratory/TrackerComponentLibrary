%-Abstract
%
%   CSPICE_DVPOOL deletes a variable from the kernel pool.
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
%      name   name(s) of a pool variable(s) to delete from the kernel pool. The
%             name and associated values are removed from the kernel pool,
%             freeing the occupied space.
%
%             [n,m] = size(name); char = class(name)
%
%             If watches are set on the variable(s) designated by 'name',
%             the corresponding agents are placed on the list of agents
%             to notify of a kernel variable update.
%
%   the call:
%
%      cspice_dvpool( name )
%
%   performs the delete operation.
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
%
%      %
%      % Load a kernel set.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % A template for the leapseconds kernel variables.
%      %
%      VAR = 'DELTET*';
%
%      %
%      % Query for the variable name, return 10 or less matches from
%      % index 1.
%      %
%      INDEX  = 1;
%      ROOM   = 10;
%
%
%      txt = sprintf( 'Kernel pool state after load.' );
%      disp( txt )
%
%      [kervar, found] = cspice_gnpool( VAR, INDEX, ROOM );
%
%      if( found )
%
%         n_elements = size(kervar, 1);
%
%         %
%         % Output the returned variable names.
%         %
%         for n=1: n_elements
%            txt = sprintf( 'Variable %d matching %s: %s', ...
%                                        n, VAR, kervar(n,:) );
%            disp( txt )
%         end
%
%      else
%         txt = sprintf( ['Failed to find  ' VAR ' in the kernel pool.'] );
%         disp( txt )
%      end
%
%
%      %
%      % Delete the kernel pool variables returned from cspice_gnpool.
%      %
%      cspice_dvpool( kervar )
%
%      txt = sprintf( '\nKernel pool state after deletion.' );
%      disp( txt )
%
%      %
%      % Confirm the variables were deleted from the pool.
%      %
%      [kervar, found] = cspice_gnpool( VAR, INDEX, ROOM );
%
%      if ( found )
%
%         n_elements = size(kervar, 1);
%
%         %
%         % Output the returned variable names.
%         %
%         for n=1: n_elements
%            txt = sprintf( 'Variable %d matching %s: %s', ...
%                                        n, VAR, kervar(n,:) );
%            disp( txt )
%         end
%
%      else
%         txt = sprintf( ['Failed to find  ' VAR ' in the kernel pool.'] );
%         disp( txt )
%      end
%
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in Matlab due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Kernel pool state after load.
%      Variable 1 matching DELTET*: DELTET/DELTA_T_A
%      Variable 2 matching DELTET*: DELTET/DELTA_AT
%      Variable 3 matching DELTET*: DELTET/K
%      Variable 4 matching DELTET*: DELTET/M
%      Variable 5 matching DELTET*: DELTET/EB
%
%      Kernel pool state after deletion.
%      Failed to find  DELTET* in the kernel pool.
%
%-Particulars
%
%   This routine enables users to programmatically remove variables
%   from the kernel pool, as opposed to having to clear the pool and
%   reload it.
%
%   Note that it is not necessary to remove kernel variables in order
%   to simply update them; this routine should be used only when
%   variables are to be removed.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dvpool_c.
%
%   MICE.REQ
%   KERNEL.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   delete a kernel pool variable
%
%-&

function cspice_dvpool(name)

   switch nargin
      case 1

         file = zzmice_str(name);

      otherwise

         error ( 'Usage: cspice_dvpool(_`name`_)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('dvpool_c', name);
   catch
      rethrow(lasterror)
   end


