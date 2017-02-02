%-Abstract
%
%   CSPICE_SXFORM returns the state transformation matrix from one
%   frame to another at a specified epoch.
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
%      from   the scalar string name of a reference frame in which
%             a position is known.
%
%      to     the scalar string name of a reference frame in which
%             it is desired to represent the position.
%
%      et     the double precision scalar or 1XN-vector of epochs in
%             ephemeris seconds past the epoch of J2000 (TDB) at which
%             the state transformation matrix should be evaluated.
%
%   the call:
%
%      xform = cspice_sxform( from, to, et )
%
%   returns:
%
%      xform   a double precision, 6x6 or 6x6xN array state
%              transformation matrix that transforms states from the
%              reference frame 'from' to frame 'to' at epoch 'et'
%
%              'xform' returns with the same vectorization measure (N)
%               as 'et'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Suppose you have geodetic coordinates of a station on the
%      % surface of Earth and that you need the inertial (J2000)
%      % state of this station.  The following code fragment
%      % illustrates how to transform the geodetic state of the
%      % station to a J2000 state.
%      %
%
%      %
%      % Load the SPK, PCK and LSK kernels.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Define a geodetic longitude, latitude, altitude
%      % coordinate set. These coordinates are defined in the
%      % non-inertial, earth fixed frame "IAU_EARTH".
%      %
%      lon = 118.25 * cspice_rpd;
%      lat = 34.05  * cspice_rpd;
%      alt = 0.;
%
%      %
%      % Define a UTC time of interest. Convert the 'utc' string
%      % to ephemeris time J2000.
%      %
%      utc = 'January 1, 1990';
%      et = cspice_str2et( utc );
%
%      %
%      % Retrieve the equatorial and polar axis of the earth (body 399).
%      %
%      abc = cspice_bodvrd( 'EARTH', 'RADII', 3 );
%      equatr =  abc(1);
%      polar  =  abc(3);
%
%      %
%      % Calculate the flattening factor for earth.
%      %
%      f =  ( equatr - polar  ) / equatr;
%
%      %
%      % Calculate the Cartesian coordinates on earth for the
%      % location at 'lon', 'lat', 'alt'.
%      %
%      estate = cspice_georec( lon, lat, alt, equatr, f);
%
%      %
%      % cspice_georec returned the position vector of the geodetic
%      % coordinates, but we want the state vector. Since it is a fixed
%      % location referenced in the "IAU_EARTH" frame, the location has
%      % no velocity. We need to extend estate to a 6-vector, the final
%      % three elements with value 0.d.
%      %
%      estate = [ estate; [0.; 0.; 0.] ];
%
%      %
%      % Retrieve the transformation matrix from "IAU_EARTH"
%      % to "J2000" at epoch 'et'.
%      %
%      xform = cspice_sxform( 'IAU_EARTH', 'J2000', et );
%
%      jstate = xform * estate;
%
%      disp( 'Scalar' )
%      txt = sprintf( 'Cartesian position in J2000 frame at epoch: %f ', et );
%      disp( txt )
%      txt = sprintf( '%16.8f %16.8f %16.8f ', jstate(1:3) );
%      disp( txt )
%
%      disp( 'Cartesian velocity in J2000 frame' )
%      txt = sprintf( '%16.8f %16.8f %16.8f ', jstate(4:6) );
%      disp( txt )
%
%      %
%      % Return the state transformation matrices from "IAU_EARTH"
%      % to "J2000" approximately every month for the time
%      % interval January 1, 1990 to January 1, 2010 (UTC).
%      %
%      %
%      % Define the time bounds for the time interval,
%      % 20 years,  convert to ephemeris time J2000.
%      %
%      utc_bounds = strvcat( '1 Jan 1990', '1 Jan 2010' );
%      et_bounds = cspice_str2et( utc_bounds );
%
%      %
%      % Step in units of a month. 20 years ~ 240 months.
%      %
%      step = (et_bounds(2) - et_bounds(1) ) / 240.;
%
%      %
%      % Create an array of 240 ephemeris times starting at
%      % et_bound(1) in intervals of 'step'.
%      %
%      et = [0:239]*step + et_bounds(1);
%
%      %
%      % Convert the 240-vector of 'et' to an array of corresponding
%      % transformation matrices (dimensions (6,6,240) ).
%      %
%      xform = cspice_sxform( 'IAU_EARTH', 'J2000', et );
%
%      %
%      % Show the dimensions of the 'xform'.
%      %
%      disp (' ' )
%      disp( 'Vector' )
%      disp( 'Dimension of xform:' )
%      disp( size(xform) )
%
%      %
%      % Apply the first and last of the transform matrices to the
%      % 'estate' vector.
%      %
%      % Transform the Cartesian state vector from "IAU_EARTH"
%      % to "J2000" at et(1) (initial epoch).
%      %
%      jstate = xform(:,:,1) * estate;
%
%      txt = sprintf( 'Cartesian position in J2000 frame at epoch: %f', et(1) );
%      disp( txt )
%      txt =  sprintf( '%24.8f %24.8f %24.8f', jstate(1:3) );
%      disp( txt )
%
%      disp( 'Cartesian velocity in J2000 frame ')
%      txt =  sprintf( '%24.12f %24.12f %24.12f', jstate(4:6) );
%      disp( txt )
%
%      disp (' ' )
%
%      %
%      % Same transformation, but at et(240) (final epoch).
%      %
%      jstate = xform(:,:,240) * estate;
%
%      txt = ...
%         sprintf( 'Cartesian position in J2000 frame at epoch: %f', et(240) );
%      disp( txt )
%
%      txt =  sprintf( '%24.8f %24.8f %24.8f', jstate(1:3) );
%      disp( txt )
%
%      disp( 'Cartesian velocity in J2000 frame ')
%      txt =  sprintf( '%24.12f %24.12f %24.12f', jstate(4:6) );
%      disp( txt )
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Scalar
%      Cartesian position in J2000 frame at epoch: -315575942.816070
%        -4131.46296088   -3308.37067191    3547.02152550
%      Cartesian velocity in J2000 frame
%            0.24124981      -0.30101944       0.00023422
%
%      Vector
%      Dimension of xform:
%           6     6   240
%
%      Cartesian position in J2000 frame at epoch: -315575942.816070
%                -4131.46296088        -3308.37067191       3547.02152550
%      Cartesian velocity in J2000 frame
%                0.241249810257       -0.301019439927       0.000234215852
%
%      Cartesian position in J2000 frame at epoch: 312946264.154758
%                 4533.62043540        2731.85929290        3546.67378733
%      Cartesian velocity in J2000 frame
%               -0.199210494903        0.330347334014       0.000192387677
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine sxform_c.
%
%   MICE.REQ
%   ROTATION.REQ
%   FRAMES.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   Find a state transformation matrix
%
%-&

function [xform] = cspice_sxform(from, to, et)

   switch nargin
      case 3

         from = zzmice_str(from);
         to   = zzmice_str(to);
         et   = zzmice_dp(et);

      otherwise

         error( 'Usage: [_xform(6,6)_] = cspice_sxform( `from`, `to`, _et_ )' )

   end

   %
   % Call the MEX library.
   %
   try
      [xform] = mice('sxform_c', from, to, et);
   catch
      rethrow(lasterror)
   end


