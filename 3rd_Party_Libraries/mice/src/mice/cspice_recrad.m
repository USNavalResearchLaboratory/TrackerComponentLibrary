%-Abstract
%
%   CSPICE_RECRAD converts rectangular (Cartesian) coordinates to
%   right ascension, declination coordinates.
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
%      rectan   a double precision 3x1 array or double precision
%               3xN array containing the rectangular coordinates of the
%               position or set of positions
%
%   the call:
%
%      [range, ra, dec] = cspice_recrad(rectan)
%
%   returns:
%
%      radius   a double precision scalar or 1XN-vector describing
%               the distance of the position from origin.
%
%      ra       a double precision scalar or 1XN-vector describing
%               the right ascension of the position as measured in
%               radians.
%
%      dec      a double precision scalar or 1XN-vector describing
%               the declination of the position as measured in radians.
%
%               'radius', 'ra', and 'dec' return with the same
%               vectorization measure (N) as 'rectan'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Output the right ascension and declination of the earth's pole
%      % in the J2000 frame approximately every month for the time
%      % interval January 1, 1990 to January 1, 2010 (UTC).
%      %
%      %
%      % Load a standard kernel set.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Define the time bounds for the time interval,
%      % 20 years,  convert to ephemeris time J2000.
%      %
%      utc_bounds = [ '1 Jan 1990'; '1 Jan 2010' ];
%      et_bounds = cspice_str2et( utc_bounds);
%
%      %
%      % Step in units of a month. 20 years ~ 240 months.
%      %
%      step = (et_bounds(2) - et_bounds(1)) / 240.;
%
%      %
%      % Create an array of 240 ephemeris times starting at
%      % et_bounds(1) in intervals of 'step'.
%      %
%      et = [0:239]*step + et_bounds(1);
%
%      %
%      % Set the conversion constant "radians to degrees."
%      %
%      r2d = cspice_dpr;
%
%      %
%      % Convert the 240-vector of 'et' to an array of corresponding
%      % transformation matrices (dimensions (3,3,240) ).
%      %
%      mat = cspice_pxform( 'IAU_EARTH', 'J2000', et);
%
%      %
%      % Extract the pole vector from the transformation matrix,
%      % convert to RA and DEC expressed in degrees.
%      %
%      % The last column in each matrix is the pole vector (z = (0,0,1))
%      % of the earth in IAU expressed in J2000. We need to copy the
%      % set of pole vectors to a 3xN array. Use reshape to do this.
%      %
%      pole = reshape( mat(:,3,:), 3,[] );
%
%      [radius, ra, dec] = cspice_recrad(pole);
%
%      ra  = ra * r2d;
%      dec = dec * r2d;
%
%      %
%      % Create an array of values for output.
%      %
%      output = [ et; ra; dec ];
%      txt = sprintf( '%17.8f %12.6f %12.6f\n' , output  );
%      disp(txt)
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      A partial output centered on et = 0:
%
%                         ...
%
%     -18408539.52023917   180.003739    89.996751
%     -15778739.49107254   180.003205    89.997215
%     -13148939.46190590   180.002671    89.997679
%     -10519139.43273926   180.002137    89.998143
%     -7889339.40357262   180.001602    89.998608
%     -5259539.37440598   180.001068    89.999072
%     -2629739.34523934   180.000534    89.999536
%           60.68392730   360.000000    90.000000
%      2629860.71309394   359.999466    89.999536
%      5259660.74226063   359.998932    89.999072
%      7889460.77142727   359.998397    89.998607
%     10519260.80059391   359.997863    89.998143
%     13149060.82976055   359.997329    89.997679
%     15778860.85892719   359.996795    89.997215
%     18408660.88809383   359.996261    89.996751
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine recrad_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   rectangular coordinates to ra and dec
%   rectangular to right_ascension and declination
%
%-&

function [range, ra, dec] = cspice_recrad(rectan)

   switch nargin
      case 1

         rectan = zzmice_dp(rectan);

      otherwise
         error ( 'Usage: [_range_, _ra_, _dec_] = cspice_recrad(_rectan(3)_)' )
   end

   %
   % Call the MEX library.
   %
   try
      [range, ra, dec] = mice('recrad_c',rectan);
   catch
      rethrow(lasterror)
   end

