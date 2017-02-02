%-Abstract
%
%   CSPICE_INEDPL calculates the intercept of a triaxial ellipsoid
%   and a plane.
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
%      a       [1,1] = size(a); double = class(a)
%
%      b       [1,1] = size(b); double = class(b)
%
%      c       [1,1] = size(c); double = class(c)
%
%              are the lengths of the semi-axes of a triaxial ellipsoid.
%              The ellipsoid is centered at the origin and oriented so that
%              its axes lie on the x, y and z axes. 'a', 'b', and 'c' are
%              the lengths of the semi-axes that respectively point in the
%              x, y, and z directions.
%
%      plane   a structure describing a SPICE plane. The intersection of
%              'plane' and the ellipsoid is sought.
%
%              [1,1] = size(plane); struct = class(plane)
%
%              The structure has the fields:
%
%                 normal:     [3,1] = size(normal); double = class(normal)
%                 constant:   [1,1] = size(constant); double = class(constant)
%
%   the call:
%
%      [ ellipse, found ] = cspice_inedpl( a, b, c, plane )
%
%   returns:
%
%      ellipse   a structure describing a SPICE ellipse that defines the
%                intersection of 'plane' and the ellipsoid.
%
%                [1,1] = size(ellipse); struct = class(ellipse)
%
%                The structure has the fields:
%
%                  center:    [3,1] = size(center); double = class(center)
%                  semiMajor: [3,1] = size(semiMajor); double = class(semiMajor)
%
%      found     the boolean indicating whether 'plane'
%                intersects the ellipsoid (true) or not (false).
%
%                [1,1] = size(found); logical = class(found)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%         KPL/MK
%
%         This meta-kernel is intended to support operation of SPICE
%         example programs. The kernels shown here should not be
%         assumed to contain adequate or correct versions of data
%         required by SPICE-based user applications.
%
%         In order for an application to use this meta-kernel, the
%         kernels referenced here must be present in the user's
%         current working directory.
%
%         The names and contents of the kernels referenced
%         by this meta-kernel are as follows:
%
%            File name                     Contents
%            ---------                     --------
%            de421.bsp                     Planetary ephemeris
%            pck00009.tpc                  Planet orientation and
%                                          radii
%            naif0009.tls                  Leapseconds
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( '/kernels/gen/lsk/naif0009.tls'
%                                '/kernels/gen/spk/de421.bsp'
%                                '/kernels/gen/pck/pck00009.tpc'
%                      )
%
%         \begintext
%
%      %
%      % Give a position relative to an ellipsoid, calculate
%      % the terminator on the ellipsoid as seen from the position.
%      % As an example, use the view of Earth from the sun.
%      %
%
%      %
%      % Standard SPK, LSK, PCK files.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Define the time to calculate the  terminator, the reference
%      % frame, and the light time correction.
%      %
%      TIME  = 'Oct 31 2002, 12:55:00 PST';
%      FRAME = 'J2000';
%      CORR  = 'LT+S';
%
%      %
%      % Convert the date string to ephemeris time.
%      %
%      et = cspice_str2et( TIME );
%
%      %
%      % calculate the position of Earth wrt the Sun.
%      %
%      [pos, ltime] = cspice_spkpos( 'EARTH', et, FRAME, CORR, 'SUN' );
%
%      %
%      % retrieve the triaxial radii of Earth.
%      %
%      radii = cspice_bodvrd( 'EARTH', 'RADII', 3 );
%
%      %
%      % Normalize the position to factors of the radii.
%      %
%      pos = [ pos(1)/radii(1)^2,
%              pos(2)/radii(2)^2,
%              pos(3)/radii(3)^2 ];
%
%      %
%      % Create the SPICE plane.
%      %
%      plane = cspice_nvc2pl( pos, 1. );
%
%      %
%      % Calculate the intercept.
%      %
%      [term, found] = cspice_inedpl( radii(1), radii(2), radii(3), plane );
%
%      %
%      % Show the ellipse.
%      %
%      center = term.center
%
%      smaj = term.semiMajor
%
%      smin = term.semiMinor
%
%      %
%      % What's the length measure of the semimajor axis.
%      %
%      smaj_norm = norm( term.semiMajor )
%
%      %
%      % What's the length measure of the semiminor axis?
%      %
%      smin_norm = norm( term.semiMinor )
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in IDL due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      center =
%
%           2.151203091605303e-01
%           1.554452696107440e-01
%           6.739164295284016e-02
%
%
%      smaj =
%
%          -3.735613602474674e+03
%           5.169706064912584e+03
%          -1.359546017505989e-11
%
%
%      smin =
%
%          -1.276335287715053e+03
%          -9.222759286967922e+02
%           6.159971549668506e+03
%
%
%      smaj_norm =
%
%           6.378139994119584e+03
%
%
%      smin_norm =
%
%           6.358055846565489e+03
%
%-Particulars
%
%   An ellipsoid and a plane can intersect in an ellipse, a single point, or
%   the empty set.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine inedpl_c.
%
%   MICE.REQ
%   ELLIPSES.REQ
%   PLANES.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 11-JUN-2013, EDW (JPL)
%
%       I/O descriptions edits to conform to Mice documentation format.
%
%   -Mice Version 1.0.0, 27-AUG-2012, EDW (JPL)
%
%-Index_Entries
%
%   intersection of ellipsoid and plane
%
%-&

function [ ellipse, found ] = cspice_inedpl( a, b, c, plane )

   switch nargin
      case 4

         a     = zzmice_dp(a);
         b     = zzmice_dp(b);
         c     = zzmice_dp(c);
         plane = zzmice_pln(plane);

      otherwise

         error ( 'Usage: [ ellipse, found ] = cspice_inedpl( a, b, c, plane )' )

   end

   try
      [ ellipse, found ] = mice( 'inedpl_c', a, b, c, plane );

      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch
      rethrow(lasterror)
   end


