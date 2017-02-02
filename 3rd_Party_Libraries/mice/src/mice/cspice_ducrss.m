%-Abstract
%
%   CSPICE_DUCRSS calculates the unit vector parallel to the cross product
%   of the position components of two state vectors and the time derivative
%   of this unit vector.
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
%      s1   a SPICE state(s);
%
%              s1 = (r1, dr1 ).
%                         --
%                         dt
%
%           [6,n] = size(s1); double = class(s1)
%
%      s2   a second SPICE state(s);
%
%              s2 = (r2, dr2 ).
%                        --
%                        dt
%
%           [6,n] = size(s2); double = class(s2)
%
%      An implicit assumption exists that 's1' and 's2' are specified
%      in the same reference frame. If this is not the case, the numerical
%      result has no meaning.
%
%   the call:
%
%      ducrss = cspice_ducrss ( s1, s2 )
%
%   returns:
%
%      ducrss   the unit vector(s) parallel to the cross product(s) of the
%               position components of 's1' and 's2' and the derivative
%               of the unit vector(s).
%
%               If the cross product of the position components is the zero
%               vector, then the position component of the output will be the
%               zero vector. The velocity component of the output will simply
%               be the derivative of the cross product of the position
%               components of 's1' and 's2'.
%
%               'ducrss' returns with the same measure of vectorization (N)
%               as 's1' and 's2'
%
%               Programmatically:
%
%                  crss   = cspice_dvcrss( s1, s2 )
%                  ducrss = cspice_dvhat( crss )
%
%               [6,n] = size(ducrss); double = class(ducrss)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      One can construct non-inertial coordinate frames from apparent
%      positions of objects or defined directions.  However, if one wants
%      to convert states in this non-inertial frame to states in an inertial
%      reference frame, the derivatives of the axes of the non-inertial
%      frame are required.
%
%      Define a reference frame with the apparent direction of the sun
%      as seen from earth as the primary axis (x). Use the earth pole vector
%      to define with the primary axis a primary plane of the frame.
%
%
%      %
%      % Load SPK, PCK, and LSK kernels, use a meta kernel for convenience.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Define the earth body-fixed pole vector (z). The pole
%      % has no velocity in the earth fixed frame "IAU_EARTH."
%      %
%      z_earth = [ 0, 0, 1, 0, 0, 0 ]';
%
%      %
%      % Calculate the state transformation between IAU_EARTH and J2000
%      % at an arbitrary epoch.
%      %
%      utc     = 'Jan 1, 2009';
%      et      = cspice_str2et( utc );
%      trans   = cspice_sxform( 'IAU_EARTH', 'J2000', et );
%
%      %
%      % Transform the earth pole vector from the IAU_EARTH frame to J2000.
%      %
%      z_j2000 = trans * z_earth;
%
%      %
%      % Calculate the apparent state of the sun from earth at the epoch
%      % 'et' in the J2000 frame.
%      %
%      target   = 'Sun';
%      observer = 'Earth';
%
%      [state, ltime] = cspice_spkezr( target, et, 'J2000', 'LT+S', observer);
%
%      %
%      % Define the z axis of the new frame as the cross product between
%      % the apparent direction of the sun and the earth pole. 'z_new' cross
%      % 'x_new' defines the y axis of the derived frame.
%      %
%      x_new = cspice_dvhat( state )
%      z_new = cspice_ducrss( state, z_j2000 )
%      y_new = cspice_ducrss( z_new, state   )
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in Matlab due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      x_new =
%
%           1.834466375917397e-01
%          -9.019196633169827e-01
%          -3.910092736476536e-01
%           2.024497675152527e-07
%           3.466010606102513e-08
%           1.503314202237741e-08

%
%      y_new =
%
%           7.884654015638601e-02
%          -3.829780802895584e-01
%           9.203863390571874e-01
%           8.238367850215384e-08
%           3.230941292533659e-08
%           6.386588623423665e-09
%
%
%      z_new =
%
%          -9.798625180410016e-01
%          -1.996715075815909e-01
%           8.572038510978363e-04
%           4.453114222022677e-08
%          -2.185310696303958e-07
%          -3.614002123088436e-11
%
%      These vectors define the transformation between the new frame and J2000.
%
%              -            -
%             |       :      |
%             |   R   :  0   |
%         M = | ......:......|
%             |       :      |
%             | dRdt  :  R   |
%             |       :      |
%              -            -
%
%      with
%
%         R    = [ x_new(1:3): y_new(1:3); z_new(1:3) ]
%
%         dRdt = [ x_new(4:6): y_new(4:6); z_new(4:6) ]
%
%-Particulars
%
%   The frame transformation described in the Example may also be implemented
%   using a dynamic frames kernel.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine ducrss_c.
%
%   FRAMES.REQ
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 09-NOV-2012, EDW (JPL)
%
%-Index_Entries
%
%   compute a unit cross product and its derivative
%
%-&

function [ducrss] = cspice_ducrss(s1, s2)

   switch nargin
      case 2

         s1 = zzmice_dp(s1);
         s2 = zzmice_dp(s2);

      otherwise

         error ( 'Usage: [_ducrss(6)_] = cspice_ducrss(_s1(6)_, _s2(6)_)' )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      [ducrss] = mice('ducrss_c', s1, s2);
   catch
      rethrow(lasterror)
   end



