%-Abstract
%
%   CSPICE_CGV2EL forms a SPICE ellipse from a center vector and two generating
%   vectors.
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
%      center   the location for an ellipse center.
%
%               [3,1] = size(center); double = class(center)
%
%      vec1 &
%      vec2     the two vectors defining the ellipse (the generating vectors)
%               with the 'center' in three-dimensional space. The ellipse is
%               the set of points
%
%                  center  +  cos(theta) vec1  +  sin(theta) vec2
%
%               where theta ranges over the interval (-pi, pi].
%
%               'vec1' and 'vec2' need not be linearly independent.
%
%               [3,1] = size(vec1); double = class(vec1)
%
%               [3,1] = size(vec2); double = class(vec2)
%
%   the call:
%
%      ellipse = cspice_cgv2el( center, vec1, vec2 )
%
%   returns:
%
%      ellipse   a structure describing a SPICE ellipse defined by the input
%                vectors. 
%
%                [1,1] = size(ellipse); struct = class(ellipse)
%
%                The structure has the fields:
%
%                center:    [3,1] = size(center); double = class(center)
%                semiMinor: [3,1] = size(semiMinor); double = class(semiMinor)
%                semiMajor: [3,1] = size(semiMajor); double = class(semiMajor)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Define the center and two linearly independent
%      % generating vectors of an ellipse (the vectors need not
%      % be linearly independent).
%      %
%      center = [ -1.;  1.; -1. ];
%      vec1   = [  1.;  1.;  1. ];
%      vec2   = [  1.; -1.;  1. ];
%
%      %
%      % Create the CSPICE_ELLIPSE structure.
%      %
%      ellipse = cspice_cgv2el( center, vec1, vec2 );
%
%      ellipse.semiMinor
%      ellipse.semiMajor
%      ellipse.center
%
%   MATLAB outputs for ellipse.semiMinor:
%
%         ans =
%
%             0.0000
%             1.4142
%             0.0000
%
%   MATLAB outputs for ellipse.semiMajor:
%
%         ans =
%
%             1.4142
%            -0.0000
%             1.4142
%
%   MATLAB outputs for ellipse.center:
%
%         ans =
%
%             -1
%              1
%             -1
%
%-Particulars
%
%   SPICE ellipses serve to simplify calling sequences and reduce
%   the chance for error in declaring and describing argument lists
%   involving ellipses.
%
%   The set of ellipse conversion routines is
%
%      cspice_cgv2el( Center and generating vectors to ellipse )
%      cspice_el2cgv( Ellipse to center and generating vectors )
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine cgv2el_c.
%
%   MICE.REQ
%   ELLIPSES.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 09-NOV-2012, EDW (JPL)
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.0, 30-DEC-2008, EDW (JPL)
%
%-Index_Entries
%
%   center and generating vectors to ellipse
%
%-&
function [ellipse] = cspice_cgv2el( center, vec1, vec2 )

   switch nargin

      case 3

         center = zzmice_dp(center);
         vec1   = zzmice_dp(vec1);
         vec2   = zzmice_dp(vec2);

      otherwise

         error ( ['Usage: [ellipse] = ' ...
                  'cspice_cgv2el( center(3), vec1(3), vec2(3) )'] )

   end

   %
   % Call the MEX library.
   %
   try
      [ellipse] = mice('cgv2el_c', center, vec1, vec2 );
   catch
      rethrow(lasterror)
   end
