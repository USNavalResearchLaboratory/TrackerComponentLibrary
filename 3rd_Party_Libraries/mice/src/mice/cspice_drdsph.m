%-Abstract
%
%   CSPICE_DRDSPH computes the Jacobian of the transformation from
%   spherical to rectangular coordinates.
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
%      r       the distance of a point from the origin.
%
%              [1,n] = size(r); double = class(r)
%
%      colat   the angle between the point and the positive z-axis, in radians.
%
%              [1,n] = size(colat); double = class(colat)
%
%      lon     the angle of the point measured from the xz plane in radians.
%              The angle increases in the counterclockwise sense about
%              the +z axis.
%
%              [1,n] = size(lon); double = class(lon)
%
%   the call:
%
%      jacobi = cspice_drdsph( r, colat, lon)
%
%   returns:
%
%      jacobi   the matrix of partial derivatives of the conversion between
%               spherical and rectangular coordinates, evaluated at the input
%               coordinates. This matrix has the form
%
%               If [1,1] = size(r) then [3,3]   = size(jacobi)
%               If [1,n] = size(r) then [3,3,n] = size(jacobi)
%                                        double = class(jacobi)
%
%                   -                                 -
%                  |  dx/dr     dx/dcolat     dx/dlon  |
%                  |                                   |
%                  |  dy/dr     dy/dcolat     dy/dlon  |
%                  |                                   |
%                  |  dz/dr     dz/dcolat     dz/dlon  |
%                   -                                 -
%
%               evaluated at the input values of r, lon and lat.
%               Here x, y, and z are given by the familiar formulae
%
%                  x = r*cos(lon)*sin(colat)
%                  y = r*sin(lon)*sin(colat)
%                  z = r*cos(colat)
%
%-Examples
%
%   None.
%
%-Particulars
%
%   It is often convenient to describe the motion of an object in
%   the spherical coordinate system.  However, when performing
%   vector computations its hard to beat rectangular coordinates.
%
%   To transform states given with respect to spherical coordinates
%   to states with respect to rectangular coordinates, one uses
%   the Jacobian of the transformation between the two systems.
%
%   Given a state in spherical coordinates
%
%      ( r, colat, lon, dr, dcolat, dlon )
%
%   the velocity in rectangular coordinates is given by the matrix
%   equation:
%                  t          |                                   t
%      (dx, dy, dz)   = jacobi|              * (dr, dcolat, dlon )
%                             |(r,colat,lon)
%
%   This routine computes the matrix
%
%            |
%      jacobi|
%            |(r,colat,lon)
%
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine drdsph_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 09-NOV-2012, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   Jacobian of rectangular w.r.t. spherical coordinates
%
%-&

function [jacobi] = cspice_drdsph( r, colat, lon)

   switch nargin
      case 3

         r     = zzmice_dp(r);
         colat = zzmice_dp(colat);
         lon   = zzmice_dp(lon);

      otherwise

         error( 'Usage: [_jacobi(3,3)_] = cspice_drdsph( _r_, _colat_, _lon_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [jacobi] = mice('drdsph_c', r, colat, lon);
   catch
      rethrow(lasterror)
   end




