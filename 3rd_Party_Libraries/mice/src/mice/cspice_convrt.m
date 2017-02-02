%-Abstract
%
%   CSPICE_CONVRT converts a measurement in a specified unit
%   to another specified unit.
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
%      x     a double precision scalar or 1xN array representing a
%            measurement in the units specified by 'in'
%
%      in    the scalar sting specifying the units associated with
%            measurement 'x'
%
%      out   represents the units desired for the measurement 'x'
%
%            Acceptable units for 'in' and 'out':
%
%              Angles:                 'RADIANS'
%                                      'DEGREES'
%                                      'ARCMINUTES'
%                                      'ARCSECONDS'
%                                      'HOURANGLE'
%                                      'MINUTEANGLE'
%                                      'SECONDANGLE'
%
%              Metric Distances:       'METERS'
%                                      'KM'
%                                      'CM'
%                                      'MM'
%
%              English Distances:      'FEET'
%                                      'INCHES'
%                                      'YARDS'
%                                      'STATUTE_MILES'
%                                      'NAUTICAL_MILES'
%
%              Astrometric Distances:  'AU'
%                                      'PARSECS'
%                                      'LIGHTSECS'
%                                      'LIGHTYEARS' julian lightyears
%
%              Time:                   'SECONDS'
%                                      'MINUTES'
%                                      'HOURS'
%                                      'DAYS'
%                                      'JULIAN_YEARS'
%                                      'TROPICAL_YEARS'
%                                      'YEARS' (same as julian years)
%
%      Neither 'in' nor 'out' are case sensitive.
%
%   the call:
%
%      y = cspice_convrt( x, in, out)
%
%   returns:
%
%      y   a double precision scalar or 1xN array representing the input
%          'x' measurement converted to the units defined by 'out'
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%    Example(1):
%
%      %
%      % Convert 300 miles (statute miles) to kilometers.
%      %
%      dist_sm = 300;
%      dist_km = cspice_convrt( dist_sm, 'statute_miles', 'km' )
%
%   MATLAB outputs:
%
%      dist_km =
%
%        482.8032
%
%   Example(2):
%
%      %
%      % Determine the number of lightyears in a vector of parsec values.
%      %
%      parsec     = [1, 3, 5];
%      lightyears = cspice_convrt( parsec, 'parsecs', 'lightyears' )
%
%   MATLAB outputs:
%
%      lightyears =
%
%          3.2616    9.7847   16.3078
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine convrt_c.
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
%   convert units
%
%-&

function [y] = cspice_convrt( x, in, out)

   switch nargin
      case 3

         x   = zzmice_dp(x);
         in  = zzmice_str(in);
         out = zzmice_str(out);

      otherwise

         error( 'Usage: [_y_] = cspice_convrt( _x_, `in`, `out`)' )

   end

   %
   % Call the MEX library.
   %
   try
      [y] = mice('convrt_c', x, in, out);
   catch
      rethrow(lasterror)
   end



