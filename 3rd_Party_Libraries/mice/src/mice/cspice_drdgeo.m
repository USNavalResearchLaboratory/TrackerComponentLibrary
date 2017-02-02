%-Abstract
%
%   CSPICE_DRDGEO computes the Jacobian of the transformation from
%   geodetic to rectangular coordinates.
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
%      lon   geodetic longitude of point (radians).
%
%            [1,n] = size(lon); double = class(lon)
%
%      lat   geodetic latitude of point (radians).
%
%            [1,n] = size(lat); double = class(lat)
%
%      alt   Altitude of point above the reference spheroid. Units of `alt'
%            must match those of `re'.
%
%            [1,n] = size(alt); double = class(alt)
%
%      re    equatorial radius of a reference spheroid. This spheroid is a
%            volume of revolution:  its horizontal cross sections are circular.
%             The shape of the spheroid is defined by an equatorial radius `re'
%            and a polar radius `rp'.  Units of 're' must match those of 'alt'.
%
%            [1,1] = size(re); double = class(re)
%
%      f     the flattening coefficient
%
%            [1,1] = size(f); double = class(f)
%
%               f = (re-rp) / re
%
%             where rp is the polar radius of the spheroid. (More importantly
%             rp = re*(1-f).) The units of `rp' match those of `re'.
%
%   the call:
%
%      jacobi = cspice_drdgeo( lon, lat, alt, re, f)
%
%   returns:
%
%      jacobi   the matrix of partial derivatives of the conversion between
%               geodetic and rectangular coordinates. It has the form
%
%               If [1,1] = size(lon) then [3,3]   = size(jacobi)
%               If [1,n] = size(lon) then [3,3,n] = size(jacobi)
%                                          double = class(jacobi)
%
%                  -                             -
%                 |  dx/dlon   dx/dlat  dx/dalt   |
%                 |                               |
%                 |  dy/dlon   dy/dlat  dy/dalt   |
%                 |                               |
%                 |  dz/dlon   dz/dlat  dz/dalt   |
%                  -                             -
%
%               evaluated at the input values of lon, lat and alt.
%
%               The formulae for computing x, y, and z from
%               geodetic coordinates are given below.
%
%                  x = [alt +        re/g(lat,f)]*cos(lon)*cos(lat)
%
%
%                  y = [alt +        re/g(lat,f)]*sin(lon)*cos(lat)
%
%                                    2
%                  z = [alt + re*(1-f) /g(lat,f)]*         sin(lat)
%
%               where
%
%                   re is the polar radius of the reference spheroid.
%
%                   f  is the flattening factor (the polar radius is
%                       obtained by multiplying the equatorial radius by 1-f).
%
%                   g( lat, f ) is given by
%
%                                2             2     2
%                      sqrt ( cos (lat) + (1-f) * sin (lat) )
%
%-Examples
%
%   None.
%
%-Particulars
%
%   It is often convenient to describe the motion of an object in
%   the geodetic coordinate system.  However, when performing
%   vector computations its hard to beat rectangular coordinates.
%
%   To transform states given with respect to geodetic coordinates
%   to states with respect to rectangular coordinates, one makes use
%   of the Jacobian of the transformation between the two systems.
%
%   Given a state in geodetic coordinates
%
%        ( lon, lat, alt, dlon, dlat, dalt )
%
%   the velocity in rectangular coordinates is given by the matrix
%   equation:
%
%                  t          |                                 t
%      (dx, dy, dz)   = jacobi|             * (dlon, dlat, dalt)
%                             |(lon,lat,alt)
%
%
%   This routine computes the matrix
%
%            |
%      jacobi|
%            |(lon,lat,alt)
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine drdgeo_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   Jacobian of rectangular w.r.t. geodetic coordinates
%
%-&

function [jacobi] = cspice_drdgeo( lon, lat, alt, re, f)

   switch nargin
      case 5

         lon = zzmice_dp(lon);
         lat = zzmice_dp(lat);
         alt = zzmice_dp(alt);
         re  = zzmice_dp(re);
         f   = zzmice_dp(f);

      otherwise

         error( ['Usage: [_jacobi(3,3)_] = '...
                 'cspice_drdgeo( _lon_, _lat_, _alt_, re, f)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [jacobi] = mice('drdgeo_c', lon, lat, alt, re, f);
   catch
      rethrow(lasterror)
   end




