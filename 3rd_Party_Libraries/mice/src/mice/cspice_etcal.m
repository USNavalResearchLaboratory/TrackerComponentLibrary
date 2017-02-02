%-Abstract
%
%   CSPICE_ETCAL converts an ephemeris epoch measured in seconds past
%   the epoch of J2000 to a calendar string format using a
%   formal calendar free of leapseconds.
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
%      et   the double precision scalar or 1xN array of ephemeris
%           time expressed as ephemeris seconds past J2000
%
%   the call:
%
%      string = cspice_etcal(et)
%
%   returns:
%
%      string   the scalar string or NXM character array representing
%               the input ephemeris epoch 'et'. This string is based upon
%               extending the Gregorian Calendar backward and forward
%               indefinitely keeping the same rules for determining leap
%               years. Moreover, there is no accounting for leapseconds.
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
%      % Define a UTC time string.
%      %
%      TIMESTR = '2013 JUN 30 00:00:00.000';
%
%      %
%      % Convert the time string to ephemeris time.
%      %
%      et  = cspice_str2et( TIMESTR );
%
%      %
%      % Convert the ephemeris time to a time string, the conversion
%      % ignoring leapseconds. Note, this evaluation does not require
%      % loading a leapsecond kernel.
%      %
%      cal = cspice_etcal( et );
%
%      %
%      % Display the two time strings.
%      %
%      disp( ['Original times string: ' TIMESTR] )
%      disp( ['ETCAL time string    : ' cal    ] )
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Original times string: 2013 JUN 30 00:00:00.000
%      ETCAL time string    : 2013 JUN 30 00:01:05.184
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine etcal_c.
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
%   Convert ephemeris time to a formal calendar date
%
%-&

function [string] = cspice_etcal(et)

   switch nargin
      case 1

         et = zzmice_dp(et);

      otherwise

         error( 'Usage: [_`string`_] = cspice_etcal(_et_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [string] = mice('etcal_c', et );
   catch
      rethrow(lasterror)
   end



