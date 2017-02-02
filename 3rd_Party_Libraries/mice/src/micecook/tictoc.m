%-Abstract
%
%   This 'cookbook' example program demonstrates use of the
%   following two time conversion routines:
%
%                  cspice_str2et
%                  cspice_et2utc
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
%   The user will be prompted for:
%
%      - The name of a leapseconds Kernel file.
%
%   Output
%
%   This program will output to the terminal, several examples
%   of valid UTC time strings and their corresponding ET
%   (ephemeris time) values.
%
%-Examples
%
%   None.
%
%-Particulars
%
%   This program uses the routines cspice_str2et and cspice_et2utc.
%   These routines convert between UTC and ET representations of
%   time:
%
%      UTC    is a character string representation of Universal
%             Time Coordinated which may be in calendar, day
%             of year, or Julian date format.  UTC time strings
%             are human-readable and thus suitable as user input.
%
%      ET     which stands for Ephemeris Time, is a double precision
%             number of ephemeris seconds past Julian year 2000,
%             also called Barycentric Dynamical Time.  ET time is
%             used internally in CSPICE routines for reading
%             ephemeris files.
%
%   For the sake of brevity, this program does NO error checking
%   on its inputs. Mistakes will cause a return to the MATLAB
%   prompt.
%
%-Required Reading
%
%   Refer to Time Required Reading and the cspice_str2et and cspice_et2utc
%   module headers for additional information.
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

function tictoc()

   %
   % Local constants.
   %
   SPICETRUE   = logical(1);
   SPICEFALSE  = logical(0);
   cont        = SPICETRUE;
   prec        = 3;

   %
   %  Assign the example data.
   %
   utc = strvcat( '9 JAN 1986 03:12:59.22451', ...
                  '1/9/86 3:12:59.22451',      ...
                  '86-365//12:00',             ...
                  'JD 2451545',                ...
                  '77 JUL 1',                  ...
                  '1 JUL ''29'                 ...
                  );

   %
   %  Information for the user. 
   %
   disp( ' '                                                         )
   disp( '                 Welcome to TICTOC'                        )
   disp( ' '                                                         )
   disp( 'This program demonstrates the use of the time conversion ' )
   disp( 'utility routines: cspice_str2et and cspice_et2utc.'        )
   disp( ' '                                                         )

   %
   % Get and load the leapsecond kernel.
   %
   disp( ' ' )
   leap = input( 'Enter the name of a leapseconds kernel file: ', 's');
   cspice_furnsh( leap )

   disp( ' ' )
   disp( 'Working ... Please wait' )
   disp( ' ' )

   %
   % Convert the time strings to ephemeris time J2000.
   %
   et = cspice_str2et( utc );

   %
   % Retrieve the size of the time vector 'et'. 'sizes' returns
   % a 1x2 vector for a vector input. The second element has the
   % vector size.
   %
   NCASES = size( et, 2 );
   i      = 1;

   %
   %   The loop continues so long as the expression:
   %
   %            i <= NCASES  &&  cont == SPICETRUE
   %
   %   evaluates to true.
   %
   while ( (i <= NCASES)  &&  (cont == SPICETRUE) )

      %
      % Begin output.
      %
      txt = sprintf( '      Example UTC time      :  %s', utc(i,:) );
      disp( txt)

      %
      % Output the corresponding ephemeris time J2000.
      %
      disp( ' ' )
      txt = sprintf( '      Corresponding ET      :  %f', et(i)    );
      disp( txt)

      %
      % Convert the ephemeris time to a calendar format.
      %
      format  = 'C';
      timestr = cspice_et2utc( et(i), format, prec );
      txt     = sprintf( '      UTC calendar format   :  %s', timestr );
      disp( txt)

      %
      % Convert the ephemeris time to a day-of-year format.
      %
      format  = 'D';
      timestr = cspice_et2utc( et(i), format, prec );
      txt     = sprintf( '      UTC day of year format:  %s', timestr );
      disp( txt)

      %
      % Convert the ephemeris time to a Julian Day format.
      %
      format  = 'J';
      timestr = cspice_et2utc( et(i), format, prec );
      txt     = sprintf( '      UTC day of year format:  %s', timestr );
      disp( txt)

      %
      % One output cycle finished. Continue?
      %
      if ( i < NCASES )
         disp( ' ' )
         answer = input( 'Continue? (Enter Y or N): ', 's');
         disp( ' ' )
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
   
   