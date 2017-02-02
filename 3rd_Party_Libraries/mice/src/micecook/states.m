%-Abstract
%
%   This "cookbook" program demonstrates the use of NAIF S- and P-
%   Kernel (SPK) files and subroutines to calculate the state
%   (position and velocity) of one solar system body relative to
%   another solar system body.
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
%   The program prompts the user for the following input:
%
%      - The name of a NAIF leapseconds kernel file.
%      - The name of a NAIF binary SPK ephemeris file.
%      - The name for the observing body.
%      - The name for the target body.
%      - Number of states to calculate.
%      - A time interval of interest.
%
%   Output
%
%      - The light time and stellar aberration corrected state of
%        the target body relative to the observing body plus
%        the magnitude to the position and velocity vectors.
%
%-Examples
%
%   None.
%
%-Particulars
%
%   The user supplies a NAIF leapseconds kernel file, a NAIF binary
%   SPK ephemeris file, valid names for both the target and
%   observing bodies, and the time to calculate the body's state.
%
%   The program makes use of the following fundamental CSPICE
%   interface routines:
%
%      cspice_furnsh  ---   makes kernel information available to
%                           the user's program.
%
%      cspice_str2et  ---   converts strings representing time to counts
%                           of seconds past the J2000 epoch.
%
%      cspice_spkezr  ---   computes states of one object relative to
%                           another at a user specified epoch.
%
%      cspice_et2utc  ---   converts an ephemeris time J200 to
%                           a formatted UTC string.
%
%   For the sake of brevity, this program does NO error checking
%   on its inputs. Mistakes will cause a return to the MATLAB
%   prompt.
%
%-Required Reading
%
%
%   For additional information, see NAIF IDS Required Reading, and the
%   headers of the functions cspice_furnsh, cspice_spkezr, cspice_et2utc
%   and cspice_str2et.
%
%-Version
%
%   -Mice Version 1.0.0, 19-FEB-2008  (EDW)
%
%-Index_Entries
%
%   None.
%
%-&

function states()

   %
   % Local constants.
   %
   format      = 'c';
   prec        = 0;
   maxpts      = 0;
   SPICETRUE   = logical(1);
   SPICEFALSE  = logical(0);

   %
   % Introduction.
   %
   disp( ' '                                                  )
   disp( '                Welcome to STATES'                  )
   disp( ' '                                                  )
   disp( 'This program demonstrates the use of NAIF S- and P-')
   disp( 'Kernel (SPK) files and subroutines by computing the')
   disp( 'state of a target body as seen from an observing'   )
   disp( 'body at a number of epochs within a given time'     )
   disp( 'interval.'                                          )
   disp( ' '                                                  )

   %
   % Get the various inputs using interactive prompts:
   %
   disp( ' ' )
   leap = input( 'Enter the name of a leapseconds kernel file: ', 's');
   disp( ' ' )

   %
   % First load the leapseconds file into the kernel pool, so
   % we can convert the UTC time strings to ephemeris seconds
   % past J2000.
   %
   cspice_furnsh( leap )

   spk = input( 'Enter the name of a binary SPK ephemeris file: ', 's');
   disp( ' ' )

   % 
   % Load the binary SPK file containing the ephemeris data
   % that we need.
   % 
   cspice_furnsh( spk  )

   %
   % Observer and target...
   %
   obs = input( 'Enter the name of the observing body: ', 's');
   disp( ' ' )

   targ = input('Enter the name of a target body: ', 's');
   disp( ' ' )

   while ( maxpts <= 0 )

      maxpts = input( 'Enter the number of states to be calculated: ');

      %
      % Check for a nonsensical input for the number of
      % look ups to perform.
      %
      % Make sure that the number of points is >= 1, to avoid a
      % division by zero error.
      %
      if ( maxpts <= 0 )
         disp( 'The number of states must be greater than 0.')
      end

      disp( ' ' )

   end

   %  
   % Query for the time interval.
   %
   if ( maxpts == 1 )

      utcbeg = input('Enter the UTC time: ', 's');
      disp( ' ' )

   elseif ( maxpts > 1 )

      utcbeg = input( 'Enter the beginning UTC time: ', 's' );
      disp( ' ' )

      utcend = input( 'Enter the ending UTC time: ', 's' );
      disp( ' ' )

   end

   frame  = input( 'Enter the inertial reference frame (e.g.:J2000): ', 's' );
   disp( ' ' )

   %
   % Output a banner for the aberration correction prompt.
   %
   disp( 'Type of correction                              Type of state    ')
   disp( '-------------------------------------------------------------    ')
   disp( '''LT+S''    Light-time and stellar aberration     Apparent state ')
   disp( '''LT''      Light-time only                       True state     ')
   disp( '''NONE''    No correction                         Geometric state')

   disp( ' ' )
   abcorr = input( 'Enter ''LT+S'', ''LT'', or ''NONE'': ', 's');

   disp( ' ' )
   disp( 'Working ... Please wait' )
   disp( ' ' )


   %
   %   Convert the UTC time strings into DOUBLE PRECISION ETs.
   %
   if ( maxpts == 1 )

      etbeg = cspice_str2et( utcbeg );

   elseif ( maxpts > 1 )

      etbeg = cspice_str2et ( utcbeg );
      etend = cspice_str2et ( utcend );

   end

   %
   % At each time, compute and print the state of the target body
   % as seen by the observer.  The output time will be in calendar
   % format, rounded to the nearest seconds.
   %
   % 'delta' is the increment between consecutive times.
   %
   % Construct the et (Ephemeris Time) vector.
   %
   if ( maxpts > 1 )

      delta = ( etend - etbeg ) / ( maxpts - 1.);
      et    = [0:(maxpts-1)]*delta + etbeg;

   elseif( maxpts == 1 )

      et    = etbeg; 

   end

   %
   % Compute the state of 'targ' from 'obs' at 'et' in the 'frame'
   % reference frame and aberration correction 'abcorr'.
   %
   [state, lt] = cspice_spkezr( targ, et, frame, abcorr, obs );

   %
   % Convert the ET (ephemeris time) into a UTC time string
   % for screen display.
   %
   utc = cspice_et2utc( et, format, prec );

   %
   % Initialize control variable for the output loop.
   %
   cont = SPICETRUE;
   i    = 1;

   %
   %   Perform the state look ups for the number of requested 
   %   intervals. The loop continues so long as the expression:
   %
   %            i <= maxpts  &&  cont == SPICETRUE
   %
   %   evaluates to true.
   %
   while ( (i <= maxpts)  &&  (cont == SPICETRUE) )

      txt = sprintf( 'For time %d of %d, the state of:', i, maxpts );
      disp( txt)

      txt = sprintf( 'Body            : %s', targ );
      disp( txt)

      txt = sprintf( 'Relative to body: %s', obs );
      disp( txt)

      txt = sprintf( 'In Frame        : %s', frame );
      disp( txt)

      txt = sprintf( 'At UTC time     : %s', utc(i,:) );
      disp( txt)

      disp( ' ' )
      disp('                 Position (km)              Velocity (km/s)'    )
      disp('            -----------------------     -----------------------')

      %
      % cspice_spkezr returns a 'state' and 'lt'. Extract i'th state vector
      % for output.
      %
      state_ele = state(:,i);

      txt = sprintf( '          X: %23.16e     %26.16e', state_ele(1), ...
                                                         state_ele(4) );
      disp( txt)

      txt = sprintf( '          Y: %23.16e     %26.16e', state_ele(2), ...
                                                         state_ele(5) );
      disp( txt)

      txt = sprintf( '          Z: %23.16e     %26.16e', state_ele(3), ...
                                                         state_ele(6) );
      disp( txt)

      txt = sprintf( '  MAGNITUDE: %23.16e     %26.16e', ...
                                   norm(state_ele(1:3)), ...
                                    norm(state_ele(4:6)) );
      disp( txt)

      %
      % One output cycle finished. Continue?
      %
      disp( ' ' )

      if ( i < maxpts )
         disp( ' ' )
         answer = input( 'Continue? (Enter Y or N): ', 's');
      end

      %
      % Perform a logical test to see if the user wants to
      % continue.
      %
      if ( strcmp( 'N', answer) || strcmp( 'n', answer) )
         cont = SPICEFALSE;
      end

      %
      % Increment the loop counter to mark the next cycle.
      %
      i  = i + 1;

   end

   %
   % It's always good form to unload kernels after use,
   % particularly in MATLAB due to data persistence.
   %
   cspice_kclear
   
   