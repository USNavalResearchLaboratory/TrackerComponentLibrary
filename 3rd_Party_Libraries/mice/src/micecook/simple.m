%-Abstract
%
%   This "cookbook" program demonstrates the use of SPICE SPK ephemeris
%   files and software.
%
%   Although this program lacks sophistication, it can serve
%   as a starting point from which you could build your own program.
%
%   The Mice subroutine cspice_furnsh (Furnish a program with SPICE
%   kernels) "loads" kernel files into the SPICE system. The calling
%   program indicates which files to load by passing their names to
%   cspice_furnsh.  It is also possible to supply cspice_furnsh with 
%   the name of a "metakernel" containing a list of files to load;
%   see the header of the function cspice_furnsh for an example.
%
%   cspice_spkezr (S/P Kernel, easier reader) computes states by by 
%   accessing the data loaded with cspice_furnsh (cspice_spkezr does not
%   require the name of an SPK file as input).
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
%   The user is prompted for the following input:
%
%      - The name of a NAIF leapseconds kernel file.
%      - The name of one binary NAIF SPK ephemeris file.
%      - The name for the first target body.
%      - The name for the second target body.
%      - The name for the observing body.
%      - A UTC time interval at which to determine states.
%
%   Output
%
%   The program calculates the angular separation of the two
%   target bodies as seen from the observing body.
%
%-Examples
%
%   None.
%
%-Particulars
%
%   The user enters the names for two target bodies and an
%   observer (these may be any objects in the solar system for
%   which the user has data) and the UTC time of interest.
%
%   For the sake of brevity, this program does NO error checking
%   on its inputs. Mistakes will cause a return to the MATLAB
%   prompt.
%
%-Required Reading
%
%   None.
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

function simple()

   SPICETRUE   = logical(1);
   SPICEFALSE  = logical(0);
   abcorr      = 'LT';
   answer      = 'n';
   MAXPTS      = 10;

   disp(' ')
   disp(' ') 
   disp('                    Welcome to SIMPLE'                  )
   disp(' ') 
   disp('This program calculates the angular separation of two'  )
   disp('target bodies as seen from an observing body.'          )
   disp(' ')
   disp('The angular separations are calculated for each of 10'  )
   disp('equally spaced times in a given time interval. A table' )
   disp('of the results is presented.')
   disp(' ')

   %
   % Set the time output format, the precision of that output
   % and the reference frame.  Note:  The angular separation has the
   % same value in all reference frames.  Let's use our favorite, J2000.
   % We need an aberration correction.  "LT+S", light time plus stellar
   % aberration, satisfies the requirements for this program.
   %
   ref  = 'J2000';
   corr = 'LT+S';

   %
   % First load the leapseconds file into the kernel pool, so
   % we can convert the UTC time strings to ephemeris seconds
   % past J2000.
   %
   leap = input( 'Enter the name of a leapseconds kernel file: ', 's');
   cspice_furnsh( leap )

   disp( ' ' )

   %
   % Get and load the SPK.
   %
   spk = input( 'Enter the name of a binary SPK ephemeris file: ', 's');
   cspice_furnsh( spk )

   disp( ' ' )

   %
   % Set-up for the user response loop
   %
   cont = SPICETRUE;

   %
   % Loop till the user quits.
   %
   while ( cont == SPICETRUE )

      obs = input( 'Enter the name of the observing body: ', 's');
      disp( ' ' )

      targ1 = input('Enter the name of the first target body: ', 's');
      disp( ' ' )

      targ2 = input('Enter the name of the second target body: ', 's');
      disp( ' ' )

      % 
      % Request for a time interval with maxpts evaluations.
      %
      utcbeg = input( 'Enter the beginning UTC time: ', 's');
      disp(' ')

      utcend = input( 'Enter the ending UTC time: ', 's');
      disp(' ')

      disp( ' ' )
      disp( 'Working ... Please wait' )
      disp( ' ' )

      etbeg = cspice_str2et( utcbeg );
      etend = cspice_str2et( utcend );

      delta = ( etend - etbeg ) / (MAXPTS - 1. );

      et = [0:(MAXPTS-1)]*delta + etbeg;

      %
      % Compute the state of targ1 and targ2 from obs at et then
      % calculate the angular separation between targ1 and targ2 
      % as seen from obs. Convert that angular value from radians 
      % to degrees.
      %
      [ pos1, lt1] = cspice_spkpos( targ1, et, ref, corr, obs );
      [ pos2, lt2] = cspice_spkpos( targ2, et, ref, corr, obs );

      y = cspice_vsep( pos1, pos2);
      y = y * cspice_dpr;

      %
      % Display the time and angular separation of the desired
      % target bodies for the requested observer for each of the
      % equally spaced evaluation times in the given time interval.
      %
      % If you have a graphics package, you may wish to write the
      % time and angular separation data to a file, and then plot
      % them for added effect.
      %
      disp( ' ' )
      txt = sprintf( 'The angular separation between bodies %s and %s,', ...
                                                         targ1, targ2 );
      disp( txt )

      txt = sprintf( 'as seen from body %s.', obs );
      disp( txt )

      disp( ' ')

      utcbeg = cspice_et2utc( etbeg, 'C', 0 );
      txt = sprintf( 'From: %s', utcbeg ); 
      disp( txt )

      utcend = cspice_et2utc( etend, 'C', 0 );
      txt = sprintf( 'To  : %s', utcend );
      disp( txt )

      utctim = cspice_et2utc( et, 'C', 0 );

      disp( ' ' )
      disp( '       UTC Time                 Separation' )
      disp( '----------------------------------------------' )

      for i=1:MAXPTS

         txt = sprintf( '  %.20s  %15.8f deg', utctim(i,:), y(i) );
         disp( txt )

      end

      disp( ' ')

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
   
   