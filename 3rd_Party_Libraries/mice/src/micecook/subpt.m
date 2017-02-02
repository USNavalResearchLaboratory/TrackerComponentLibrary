%-Abstract
%
%   This "cookbook" program demonstrates the use of the CSPICE
%   Toolkit by computing the apparent sub-observer point on a target
%   body. It uses light time and stellar aberration corrections in
%   order to do this.
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
%   The user is prompted for the following:
%
%      - The name of a leapseconds kernel file.
%      - The name of a Planetary constants (PCK) kernel file.
%      - The name of a NAIF SPK Ephemeris file.
%      - The name of the observing body.
%      - The name of the target body.
%      - The name of the body-fixed reference frame
%        associated with the target body (for example, IAU_MARS).
%      - Number of evaluations to perform
%      - A UTC time interval of interest.
%
%   Output
%
%
%   The program calculates the planetocentric latitude and longitude
%   of the nearest point on the target body to the observing body
%   for a UTC epoch (see Input above). The program outputs
%
%      - The epoch of interest as supplied by the user.
%      - The planetocentric longitude of the nearest point on
%        the target body to the observing body.
%      - The planetocentric latitude of the nearest point on
%        the target body to the observing body.
%
%
%-Examples
%
%   None.
%
%-Particulars
%
%   The SPK file must contain data for both the observing body and
%   the target body during the specified time interval.
%
%   The "apparent sub-observer point" is defined in this program to
%   be the point on the target body that appears to be closest to the
%   observer. The apparent sub-observer point may also be defined as
%   the intercept on the target's surface of the ray emanating from
%   the observer and passing through the apparent target body's
%   center, but we don't demonstrate use of that definition here. See
%   the header of cspice_subpnt for details.
%
%   In order to compute the apparent location of the sub-observer
%   point, we correct the position of the sub-observer point for both
%   light time and stellar aberration, and we correct the orientation
%   of the target body for light time. We consider "light time" to be
%   the time it takes a photon to travel from the sub-observer point
%   to the observer. If the light time is given the name LT, then the
%   apparent position of the sub-observer point relative to the
%   observer is defined by the vector from the sub-observer point's
%   location (relative to the solar system barycenter) at ET-LT,
%   minus the observer's location (again, relative to the solar
%   system barycenter) at ET, where this difference vector is
%   corrected for stellar aberration.
%
%   See the header of the Mice routine cspice_spkezr for more
%   information on light time and stellar aberration corrections; see
%   the header of the Mice routine cspice_subpnt for an explanation of
%   how it applies aberration corrections.
%
%   Planetocentric coordinates are defined by a distance from a
%   central reference point, an angle from a reference meridian,
%   and an angle above the equator of a sphere centered at the
%   central reference point. These are the radius, longitude,
%   and latitude, respectively.
%
%   The program makes use of the following fundamental Mice
%   interface routines:
%
%      cspice_furnsh --- makes kernel information available to
%                        the user's program.
%
%      cspice_str2et --- converts strings representing time to counts
%                        of seconds past the J2000 epoch.
%
%      cspice_spkezr --- computes states of one object relative to
%                        another at a user specified epoch.
%
%      cspice_et2utc --- converts an ephemeris time J200 to
%                        a formatted UTC string.
%
%      cspice_subpnt --- calculate the position of the sub-observer point
%                        of one body with respect to another
%
%   For the sake of brevity, this program does NO error checking
%   on its inputs. Mistakes will cause a return to the MATLAB
%   prompt.
%
%-Required Reading
%
%      KERNEL        The CSPICE Kernel Pool
%      ROTATIONS     Rotations
%      SPK           S- and P- Kernel (SPK) Specification
%      TIME          Time routines in CSPICE
%
%   For questions about a particular function, refer to its
%   header.
%
%-Version
%
%   -Mice Version 1.0.0, 19-FEB-2008  (NJB) (EDW)
%
%-Index_Entries
%
%   None.
%
%-&

function subpt()

   SPICETRUE   = logical(1);
   SPICEFALSE  = logical(0);
   abcorr      = 'LT+S';
   answer      = 'n';

   %
   % An intro banner.
   %
   disp ( ' '                                                        )
   disp ( '             Welcome to SUBPT'                            )
   disp ( ' '                                                        )
   disp ( 'This program demonstrates the use of CSPICE in computing' )
   disp ( 'the apparent sub-observer point on a target body. The'    )
   disp ( 'computations use light time and stellar aberration'       )
   disp ( 'corrections.'                                             )
   disp ( ' '                                                        )

   %
   % Get the various inputs using interactive prompts:
   %

   disp( ' ' )

   %
   % First load the leapseconds file into the kernel pool, so
   % we can convert the UTC time strings to ephemeris seconds
   % past J2000.
   %
   leap = input( 'Enter the name of a leapseconds kernel file: ', 's');
   cspice_furnsh( leap )

   disp( ' ' )

   %
   % Get and load the physical constants kernel.
   %
   pck = input('Enter the name of a planetary constants kernel: ', 's' );
   cspice_furnsh( pck )

   disp( ' ' )

   %
   % Get and load the SPK.
   %
   spk = input( 'Enter the name of a binary SPK ephemeris file: ', 's');
   cspice_furnsh( spk )

   disp( ' ' )
   disp( 'Working ... Please wait' )
   disp( ' ' )

   %
   % Set-up for the user response loop
   %
   cont = SPICETRUE;

   %
   % Loop till the user quits.
   %
   while ( cont == SPICETRUE )

      %
      % Observer and target...
      %
      obs = input( 'Enter the name of the observing body: ', 's');
      disp( ' ' )

      targ = input('Enter the name of a target body: ', 's');
      disp( ' ' )

      fixfrm = input('Enter the name of the target body-fixed frame: ', 's');
      disp( ' ' )

      maxpts = input( 'Enter the number of points to calculate: ' );
      disp( ' ' )

      if ( maxpts <= 0 )
         maxpts = 1;
      end

      %
      % Input strings for the UTC time interval, or single UTC
      % time for a single evaluation.
      % 
      % Convert the UTC time interval to ET. ET stands for Ephemeris
      % Time and is in units of ephemeris seconds past Julian year
      % 2000. ET is the time system that is used internally in SPK
      % ephemeris files and reader subroutines.
      % 
      % DELTA is the increment between consecutive times, if
      % needed.
      % 
      if ( maxpts == 1 )

         %
         % Request for a single evaluation. No steps - no delta.
         %
         utcbeg = input( 'Enter the UTC time: ', 's');
         disp(' ')

         etbeg = cspice_str2et( utcbeg );
         delta = 0.;

         epoch  = etbeg;

      else

         % 
         % Request for a time interval with maxpts evaluations.
         %
         utcbeg = input( 'Enter the beginning UTC time: ', 's');
         disp(' ')

         utcend = input( 'Enter the ending UTC time: ', 's');
         disp(' ')

         etbeg = cspice_str2et( utcbeg );
         etend = cspice_str2et( utcend );

         delta  = ( etend - etbeg ) / (maxpts - 1. );
         epoch = [0:(maxpts-1)]*delta + etbeg;

      end

      %
      % Write the headings for the table of values.
      %
      disp( 'Planetocentric coordinates for the nearest point' )
      disp( 'on the target body to the observing body (deg).'  )

      txt = sprintf( 'Target body: %s          Observing body: %s', ...
                                                              targ, ...
                                                              obs );
      disp( txt )

      disp( ' ' )
      disp( '       UTC Time            Lat         Lon')
      disp( '----------------------------------------------')

      % 
      % Now, everything is set up for output
      %
      npts   = 1;

      while ( npts <= maxpts )

         %
         % Note: cspice_subpnt can also calculate a "sub-observer point" via
         % the intercept of the observer-target vector with the target
         % body's surface. The computation "method" argument value for
         % that calculation is
         % 
         %     "Intercept: ellipsoid"
         % 
         % The output sub-observer point `spoint' is expressed in the
         % body-fixed reference frame `fixfrm' specified by the user,
         % where the orientation of the frame is evaluated at the time
         % `trgepc'. `trgepc' is expressed in seconds past J2000 TDB, and
         % is equal to et-lt, where `lt' is the light time from the
         % sub-observer point to the observer. The output `srfvec' is
         % the apparent position of the sub-observer point relative to
         % the observer. `srfvec' is also expressed in the reference
         % frame `fixfrm'.
         % 
         % Please see the cspice_subpnt source file header for further
         % information.
         %
         [spoint, trgepc, srfvec] = cspice_subpnt( 'Near point: ellipsoid', ... 
                                                  targ, epoch(npts),        ...
                                                  fixfrm, abcorr, obs );
   
         [ radius, lon, lat] = cspice_reclat( spoint );

         %
         % Multiply lat and lon by the number of degrees per radian.
         %
         lon = lon * cspice_dpr;
         lat = lat * cspice_dpr;

         %
         % Convert the current EPOCH to UTC time for display.
         %
         utcout = cspice_et2utc( epoch(npts), 'C', 3 );

         %
         % Display results in a table format:
         %
         txt = sprintf ('  %.20s  %9.5f    %9.5f', utcout, ...
                                                   lat,    ... 
                                                   lon  );
         disp( txt )

         npts  = npts + 1;

      end

      disp( ' ' )
      answer = input( 'Continue? (Enter Y or N): ', 's');

      %
      % Perform a logical test to see if the user wants to
      % continue.
      %
      if ( strcmp( 'N', answer) || strcmp( 'n', answer) )
         cont = SPICEFALSE;
      end

   end

   %
   % It's always good form to unload kernels after use,
   % particularly in MATLAB due to data persistence.
   %
   cspice_kclear
   
   