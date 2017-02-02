%-Abstract
%
%   CSPICE_DVDOT returns the time derivative of the dot product of
%   two position vectors.
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
%      dvdot = cspice_dvdot( s1, s2 )
%
%   returns:
%
%      dvdot   the time derivative(s) of the dot product between the position
%              components of 's1' and 's2'.
%
%              'dvdot' returns with the same measure of vectorization (N)
%              as 's1' and 's2'.
%
%              [1,n] = size(dvdot); double = class(dvdot)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%     Suppose that given two state vectors (s1 and s2) whose position
%     components are unit vectors, and that we need to compute the
%     rate of change of the angle between the two vectors.
%
%     We know that the Cosine of the angle THETA between them is given
%     by
%
%        cos(theta) = dot(s1,s2)
%
%     Thus by the chain rule, the derivative of the angle is given
%     by:
%
%        sine(theta) dtheta/dt = cspice_dvdot(s1,s2)
%
%     Thus for values of theta away from zero we can compute
%
%     dtheta/dt as
%
%     dtheta = cspice_dvdot(s1,s2) / sqrt( 1 - dot(s1,s2)**2 )
%
%     Note that position components of s1 and s2 are parallel, the
%     derivative of the  angle between the positions does not
%     exist.  Any code that computes the derivative of the angle
%     between two position vectors should account for the case
%     when the position components are parallel.
%
%-Particulars
%
%   In this discussion, the notation
%
%      < V1, V2 >
%
%   indicates the dot product of vectors V1 and V2.
%
%   Given two state vectors s1 and s2 made up of position and velocity
%   components (r1,v1) and (r2,v2) respectively, cspice_dvdot calculates
%   the derivative of the dot product of p1 and p2, i.e. the time
%   derivative
%
%         d
%         -- < r1, r2 > = < v1, r2 > + < r1, v2 >
%         dt
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dvdot_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 20-APR-2010, EDW (JPL)
%
%-Index_Entries
%
%   time derivative of a dot product
%
%-&

function [dvdot] = cspice_dvdot(s1, s2)

   switch nargin
      case 2

         s1 = zzmice_dp(s1);
         s2 = zzmice_dp(s2);

      otherwise

         error ( 'Usage: [_dvdot_] = cspice_dvdot(_s1(6)_, _s2(6)_)' )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      [dvdot] = mice('dvdot_c', s1, s2);
   catch
      rethrow(lasterror)
   end

