%-Abstract
%
%   CSPICE_DRDLAT computes the Jacobian of the transformation from latitudinal
%   to rectangular coordinates.
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
%      radius   the distance of a point from the origin.
%
%               [1,n] = size(radius); double = class(radius)
%
%      lon      the angle of the point measured from the XZ plane in radians.
%               The angle increases in the counterclockwise sense about the 
%               +Z axis.
%
%               [1,n] = size(lon); double = class(lon)
%
%      lat      the angle of the point measured from the XY plane in radians.
%               The angle increases in the direction of the +Z axis.
%
%               [1,n] = size(lat); double = class(lat)
%
%   the call:
%
%      jacobi = cspice_drdlat( r, lon, lat)
%
%   returns:
%
%      jacobi   the matrix of partial derivatives of the conversion between
%               latitudinal and rectangular coordinates, evaluated at the input
%               coordinates. This matrix has the form
%
%               If [1,1] = size(radius) then [3,3]   = size(jacobi)
%               If [1,n] = size(radius) then [3,3,n] = size(jacobi)
%                                             double = class(jacobi)
%
%                   -                                -
%                  |  dx/dr     dx/dlon     dx/dlat   |
%                  |                                  |
%                  |  dy/dr     dy/dlon     dy/dlat   |
%                  |                                  |
%                  |  dz/dr     dz/dlon     dz/dlat   |
%                   -                                -
%
%               evaluated at the input values of r, lon and lat.
%               Here x, y, and z are given by the familiar formulae
%
%                  x = r * cos(lon) * cos(lat)
%                  y = r * sin(lon) * cos(lat)
%                  z = r *            sin(lat).
%
%-Examples
%
%   None.
%
%-Particulars
%
%   It is often convenient to describe the motion of an object
%   in latitudinal coordinates. It is also convenient to manipulate
%   vectors associated with the object in rectangular coordinates.
%
%   The transformation of a latitudinal state into an equivalent
%   rectangular state makes use of the Jacobian of the
%   transformation between the two systems.
%
%   Given a state in latitudinal coordinates,
%
%        ( r, lon, lat, dr, dlon, dlat )
%
%   the velocity in rectangular coordinates is given by the matrix
%   equation
%                  t          |                               t
%      (dx, dy, dz)   = jacobi|             * (dr, dlon, dlat)
%                             |(r,lon,lat)
%
%   This routine computes the matrix
%
%            |
%      jacobi|
%            |(r,lon,lat)
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine drdlat_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   Jacobian of rectangular w.r.t. latitudinal coordinates
%
%-&

function [jacobi] = cspice_drdlat( r, lon, lat)

   switch nargin
      case 3

         r   = zzmice_dp(r);
         lon = zzmice_dp(lon);
         lat = zzmice_dp(lat);

      otherwise

         error( 'Usage: [_jacobi(3,3)_] = cspice_drdlat( _r_, _lon_, _lat_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [jacobi] = mice('drdlat_c', r, lon, lat);
   catch
      rethrow(lasterror)
   end




