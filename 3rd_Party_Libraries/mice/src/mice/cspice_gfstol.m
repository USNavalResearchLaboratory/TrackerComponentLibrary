%-Abstract
%
%   CSPICE_GFSTOL overrides the default GF convergence value used in the high
%   level GF routines.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
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
%      value   value to use as the GF subsystem convergence tolerance. This
%              value will override the default tolerance, SPICE_GF_CNVTOL,
%              defined in SpiceGF.h. Units are TDB seconds.
%
%              [1,1] = size(value); double = class(value)
%
%   the call:
%
%      cspice_gfstol( value )
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
%      In 14 A.D., the Roman princeps Tiberius sent his son Drusus to subdue
%      a revolt of a Roman Legion stationed in Pannonia. A Lunar eclipse
%      occurred during this mission.
%
%      Perform a search for occultation events of the sun by earth as
%      observed from the Moon center. Search during the interval from
%      14 A.D. SEP 1 to 14 A.D. SEP 30 (Julian).
%
%      TIMFMT  = 'YYYY ERA MON DD HR:MN:SC.#### ::JCAL';
%      MAXWIN  = 100;
%
%      %
%      % Load kernels.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Use an SPK covering year 14 AD.
%      %
%      cspice_furnsh( 'de408.bsp' )
%
%      %
%      % Store the time bounds of our search interval in
%      % the cnfine confinement window.
%      %
%      et = cspice_str2et( { '14 A.D. SEP 1  00:00:00', ...
%                            '14 A.D. SEP 30 00:00:00'} );
%
%      cnfine = cspice_wninsd( et(1), et(2) );
%
%      %
%      % Select a 3-minute step. We'll ignore any occultations
%      % lasting less than 3 minutes.
%      %
%      step    = 180.;
%
%      occtyp  = 'any';
%      front   = 'earth';
%      fshape  = 'ellipsoid';
%      fframe  = 'iau_earth';
%      back    = 'sun';
%      bshape  = 'ellipsoid';
%      bframe  = 'iau_sun';
%      obsrvr  = 'moon';
%      abcorr  = 'lt';
%
%      %
%      % Perform the search. 'et(1)' and 'et(2)' have values ~-6*10^10,
%      % SPICE_GF_CNVTOL has value 10^-6, so double precision addition or
%      % subtraction of 'et(1)' and 'et(2)' with SPICE_GF_CNVTOL returns
%      % a result indistinguishable from 'et(1)' and 'et(2)'.
%      %
%      % Reduce the GF convergence tolerance by an order of magnitude
%      % to resolve this condition.
%      %
%
%      cspice_gfstol( 1e-5 )
%
%      result = cspice_gfoclt( occtyp, front,  fshape, fframe, ...
%                              back,   bshape, bframe,         ...
%                              abcorr, obsrvr, step,   cnfine, ...
%                              MAXWIN);
%
%      %
%      % List the beginning and ending times in each interval
%      % if result contains data.
%      %
%      count = cspice_wncard(result);
%
%      for i=1:count
%
%         [left, right] = cspice_wnfetd( result, i );
%
%         output = cspice_timout( [left,right], TIMFMT );
%
%         if( isequal( left, right) )
%
%            disp( ['Event time: ' output(1,:)] )
%
%         else
%
%            disp( ['Start time :' output(1,:)] )
%            disp( ['Stop time  :' output(2,:)] )
%            disp( ' ')
%
%         end
%
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in Matlab due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Start time :  14 A.D. SEP 27 05:02:02.8250
%      Stop time  :  14 A.D. SEP 27 09:33:31.6995
%
%-Particulars
%
%   The high level GF routines (see GF.REQ for a listing) use a default
%   value for the convergence tolerance, SPICE_GF_CNVTOL, defined in
%   SpiceGF.h. It may occur that a GF search run needs a different
%   convergence tolerance. cspice_gfstol programmatically changes the
%   tolerance used by those routines.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine gfstol_c.
%
%   MICE.REQ
%   GF.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   change default convergence tolerance for GF routines
%
%-&

function cspice_gfstol( value )

  switch nargin

      case 1

         value = zzmice_dp(value);

      otherwise

         error ( 'Usage: cspice_gfstol(value)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('gfstol_c', value);
   catch
      rethrow(lasterror)
   end


