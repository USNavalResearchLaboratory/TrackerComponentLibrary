%-Abstract
%
%   CSPICE_DGEODR computes the Jacobian of the transformation from
%   rectangular to geodetic coordinates.
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
%      x,
%      y,
%      z    the rectangular coordinates of the point at which the Jacobian of
%           the map from rectangular to geodetic coordinates is desired.
%
%           [1,n] = size(z); double = class(z)
%
%      re   equatorial radius of a reference spheroid. This spheroid is a
%           volume of revolution: its horizontal cross sections are circular. 
%           The shape of the spheroid is defined by an equatorial radius `re'
%           and a polar radius `rp'.
%
%           [1,1] = size(re); double = class(re)
%
%      f    the flattening coefficient
%
%           [1,1] = size(f); double = class(f)
%
%               f = (re-rp) / re
%
%             where rp is the polar radius of the spheroid. (More importantly
%             rp = re*(1-f).) The units of `rp' match those of `re'.
%
%   the call:
%
%      jacobi = cspice_dgeodr( x, y, z, re, f)
%
%   returns:
%
%      jacobi   the matrix of partial derivatives of the conversion between
%               rectangular and geodetic coordinates, evaluated at the input
%               coordinates. This matrix has the form
%
%               [3,3] = size(jacobi); double = class(jacobi)
%
%                   -                            -
%                  |  dlon/dx   dlon/dy  dlon/dz  |
%                  |                              |
%                  |  dlat/dx   dlat/dy  dlat/dz  |
%                  |                              |
%                  |  dalt/dx   dalt/dy  dalt/dz  |
%                   -                            -
%
%               evaluated at the input values of x, y, and z.
%
%-Examples
%
%   None.
%
%-Particulars
%
%   When performing vector calculations with velocities it is
%   usually most convenient to work in rectangular coordinates.
%   However, once the vector manipulations have been performed,
%   it is often desirable to convert the rectangular representations
%   into geodetic coordinates to gain insights about phenomena
%   in this coordinate frame.
%
%   To transform rectangular velocities to derivatives of coordinates
%   in a geodetic system, one uses the Jacobian of the transformation
%   between the two systems.
%
%   Given a state in rectangular coordinates
%
%      ( x, y, z, dx, dy, dz )
%
%   the velocity in geodetic coordinates is given by the matrix
%   equation:
%                        t          |                     t
%      (dlon, dlat, dalt)   = jacobi|       * (dx, dy, dz)
%                                   |(x,y,z)
%
%   This routine computes the matrix
%
%            |
%      jacobi|
%            |(x, y, z)
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dgeodr_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   Jacobian of geodetic  w.r.t. rectangular coordinates
%
%-&

function [jacobi] = cspice_dgeodr( x, y, z, re, f)

   switch nargin
      case 5

         x = zzmice_dp(x);
         y = zzmice_dp(y);
         z = zzmice_dp(z);
         re= zzmice_dp(re);
         f = zzmice_dp(f);

      otherwise

         error( 'Usage: [_jacobi(3,3)_] = cspice_dgeodr( _x_, _y_, _z_, re, f)')

   end

   %
   % Call the MEX library.
   %
   try
      [jacobi] = mice('dgeodr_c', x, y, z, re, f);
   catch
      rethrow(lasterror)
   end




