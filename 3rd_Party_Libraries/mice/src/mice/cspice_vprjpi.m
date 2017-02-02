%-Abstract
%
%   CSPICE_VPRJPI calculates the vector in a specified plane that
%   maps under orthogonal projection to a specified vector in
%   another plane.
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
%      vin      a 3-vector.
%
%               [3,1] = size(vin); double = class(vin)
%
%      projpl   a SPICE plane that represents the geometric plane containing
%               'vin'. The structure has the fields:
%
%               SPICE plane STRUCTURE
%
%                  normal     constant  
%                             [3,1] = size(normal); double = class(normal)
%
%      invpl    a SPICE plane that represents the geometric plane containing
%               the inverse image of 'vin' under orthogonal projection onto
%               'projpl'. The structure has the fields:
%
%               SPICE plane STRUCTURE
%
%                  normal     constant 
%                             [3,1] = size(normal); double = class(normal)
%
%   the call:
%
%      [vout, found] = cspice_vprjpi( vin, projpl, invpl )
%
%   returns:
%
%      vout     inverse orthogonal projection of 'vin'. This is the vector
%               lying in the plane 'invpl' whose orthogonal projection onto the
%               plane 'projpl' is 'vin'. 'vout' is valid only when 'found' is
%               true. Otherwise, 'vout' is undefined.
%
%               [3,1] = size(vout); double = class(vout)
%
%      found    indicating whether the inverse orthogonal projection of 'vin'
%               could be computed. 'found' is true if so, false otherwise.
%
%               [1,1] = size(found); logical = class(found)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Construct 2 planes via cspice_nvc2pl. Define the normal
%      % vectors and constants for the planes.
%      %
%      norm1 = [ 0, 0, 1]';
%      norm2 = [ 1, 0, 1]';
%      con1  = 1.2;
%      con2  = 0.65;
%
%      %
%      % Use the normals and constants to create the plane structures,
%      % plane1 and plane2.
%      %
%      plane1 = cspice_nvc2pl( norm1, con1 );
%      plane2 = cspice_nvc2pl( norm2, con2 );
%
%      %
%      % Define a vector in plane1...
%      %
%      vec = [ 1, 1, 0]';
%
%      %
%      % Calculate the inverse projection to plane2.
%      %
%      [ vec_iproj, found] = cspice_vprjpi( vec, plane1, plane2);
%
%      if ( found )
%         disp( 'Found inverse vector:' )
%         vec_iproj
%      else
%         disp( 'Could not find the inverse vector.' )
%      end
%
%   MATLAB outputs:
%
%      Found inverse vector:
%
%      vec_iproj =
%
%           1.000000000000000e+00
%           1.000000000000000e+00
%          -3.500000000000000e-01
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine vprjpi_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   vector projection onto plane
%
%-&

function [vout, found] = cspice_vprjpi( vin, projpl, invpl )

   switch nargin

      case 3

         vin    = zzmice_dp( vin );
         projpl = zzmice_pln( projpl );
         invpl  = zzmice_pln( invpl );

      otherwise

         error ( ['Usage: [vout(3), found] = ' ...
                                  'cspice_vprjpi( vin(3), projpl, invpl )'] )

   end

   %
   % Call the MEX library.
   %
   % The developer decided to not complicate the interface call and so
   % use the individual fields of the 'plane' structure as arguments.
   %
   try
      [vout, found] = mice('vprjpi_c', vin, projpl, invpl );
   catch
      rethrow(lasterror)
   end

