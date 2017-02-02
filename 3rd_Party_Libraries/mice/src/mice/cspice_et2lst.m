%-Abstract
%
%   CSPICE_ET2LST computes the local solar time at a given ephemeris epoch,
%   for an object on the surface of a body at a specified longitude.
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
%      et     the double precision scalar or 1xN array of ephemeris
%             time expressed as ephemeris seconds past J2000 at which
%             a local time is desired
%
%      body   an integer scalar SPICE ID-code of the body on which local
%             time is to be measured
%
%      lon    the scalar double precision longitude (either planetocentric
%             or planetographic) in radians of the site on the surface
%             of body for which local time should be computed
%
%      type   is the form of longitude supplied by the variable
%             lon.   Allowed values are "PLANETOCENTRIC" and
%             "PLANETOGRAPHIC".  Note the case of the letters
%             in type is insignificant.  Both "PLANETOCENTRIC"
%             and "planetocentric" are recognized.  Leading and
%             trailing blanks in type are not significant.
%
%   the call:
%
%      [ hr, min, sec, time, ampm] = cspice_et2lst( et, body, lon, type)
%
%   returns:
%
%      hr     a double precision scalar or double precision 1xN array describing
%             the integral number of the local "hour" of the site specified
%             at epoch 'et'
%
%             Note that an "hour" of local time does not have the same duration
%             as an hour measured by conventional clocks. It is simply a
%             representation of an angle.
%
%      mn     a double precision scalar or double precision 1xN array describing
%             the integral number of "minutes" past the hour of the local time
%             of the site at the epoch 'et'
%
%             Again note that a "local minute" is not the same as a minute you
%             would measure with conventional clocks
%
%      sc     a double precision scalar or double precision 1xN array describing
%             the integral number of "seconds" past the minute of the local time
%             of the site at the epoch `et'
%
%             Again note that a "local second" is not the same as a second
%             you would measure with conventional clocks.
%
%      time   the scalar string or NXM character array of output local time
%             on a "24 hour" local clock
%
%      ampm   the scalar string or NXM character array of output local time
%             on a "12 hour" local clock together with the traditional AM/PM
%             label to indicate whether the sun has crossed the local zenith
%             meridian.
%
%             All output arguments return with the same measure of
%             vectorization, N, as 'et'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Load a leapseconds kernel.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Define two UTC time strings to 'utc'
%      %
%      utc                        = strvcat( '2002 SEP 02 00:00:00', ...
%                                            '2002 SEP 30 00:00:00' );
%
%      %
%      % Convert 'utc' the ephemeris time, 'et'
%      %
%      et                         = cspice_str2et(utc);
%
%      %
%      % Define a planetographic longitude in degrees, convert the
%      % value to radians
%      %
%      dlon                       =  326.17;
%      rlon                       =  dlon * cspice_rpd;
%
%      %
%      % Convert inputs to Local Solar Time.
%      %
%      [hr, min, sec, time, ampm] = cspice_et2lst( et,   ...
%                                                  499,  ...
%                                                  rlon, ...
%                                                  'PLANETOGRAPHIC');
%
%      fprintf( ['The local time at Mars %6.2f degrees E ' ...
%               'planetographic longitude:\n'],            ...
%               dlon )
%      fprintf( '   at UTC %s, LST = %s\n', utc(1,:), ampm(1,:) )
%      fprintf( '   at UTC %s, LST = %s\n', utc(2,:), ampm(2,:) )
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      The local time at Mars 326.17 degrees E planetographic longitude:
%         at UTC 2002 SEP 02 00:00:00, LST = 03:25:35 A.M.
%         at UTC 2002 SEP 30 00:00:00, LST = 09:33:00 A.M.
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine et2lst_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 06-MAY-2009, EDW (JPL)
%
%      Added MICE.REQ reference to the Required Reading section.
%
%   -Mice Version 1.0.0, 07-MAR-2007, EDW (JPL)
%
%-Index_Entries
%
%   Compute the local time for a point on a body.
%
%-&

function [hr, min, sec, time, ampm] = cspice_et2lst( et, body, lon, type)

   switch nargin
      case 4

         et    = zzmice_dp(et);
         body  = zzmice_int(body);
         lon   = zzmice_dp(lon);
         type  = zzmice_str(type);

      otherwise

         error ( ['Usage: [ _hr_, _min_, _sec_, _`time`_, _`ampm`_] = ' ...
                 'cspice_et2lst( _et_, body, lon, `type`)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [hr, min, sec, time, ampm] = mice('et2lst_c', et, body, lon, type);

      %
      % Convert the integers returned from the interface to double precision
      % in case a user includes the return arguments in a calculation
      % with other doubles.
      %
      hr  = zzmice_dp(hr);
      min = zzmice_dp(min);
      sec = zzmice_dp(sec);

   catch
      rethrow(lasterror)
   end





