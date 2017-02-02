%-Abstract
%
%   CSPICE_FRAME builds a right handed orthonormal frame x, y, z,
%   given an input vector 'vec' where the output vector 'x'
%   is parallel to the original input vector.
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
%      vec   the vector used to form the first vector of a right-handed
%            orthonormal triple.
%
%            [3,1] = size(vec); double = class(vec)
%
%   the call:
%
%      [ x, y, z ] = cspice_frame( vec )
%
%   returns:
%
%      x   
%      y   
%      z   [3,1] = size(x); double = class(x)
%          [3,1] = size(y); double = class(y)
%          [3,1] = size(z); double = class(z)
%
%          the vectors that form a right handed orthonormal frame, with
%          output 'x' parallel to the input 'vec'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Given some arbitrary vector
%      %
%      vec = [ 23., -3., 18. ]';
%
%      %
%      % Create an orthonormal frame with the
%      % x axis parallel to vec.
%      %
%      disp( 'Input vector parallel to output X vector' )
%
%      [ x, y, z ] = cspice_frame( vec );
%      fprintf( 'x  %16.8e   %16.8e   %16.8e\n', x )
%      fprintf( 'y  %16.8e   %16.8e   %16.8e\n', y )
%      fprintf( 'z  %16.8e   %16.8e   %16.8e\n\n', z )
%
%
%
%      %
%      % Alternative, make a frame with y parallel to
%      % vec...
%      %
%      disp( 'Input vector parallel to output Y vector' )
%
%      [ y, z, x ] = cspice_frame( vec );
%      fprintf( 'x  %16.8e   %16.8e   %16.8e\n', x )
%      fprintf( 'y  %16.8e   %16.8e   %16.8e\n', y )
%      fprintf( 'z  %16.8e   %16.8e   %16.8e\n\n', z )
%
%
%
%      %
%      % ...or a frame with z parallel to vec.
%      %
%      disp( 'Input vector parallel to output Z vector' )
%
%      [ z, x, y ] = cspice_frame( vec );
%      fprintf( 'x  %16.8e   %16.8e   %16.8e\n', x )
%      fprintf( 'y  %16.8e   %16.8e   %16.8e\n', y )
%      fprintf( 'z  %16.8e   %16.8e   %16.8e\n\n', z )
%
%   MATLAB outputs:
%
%      Input vector parallel to output X vector
%      x    7.83383109e-01    -1.02180405e-01     6.13082433e-01
%      y    6.16308262e-01     0.00000000e+00    -7.87505001e-01
%      z    8.04675803e-02     9.94765884e-01     6.29746281e-02
%
%      Input vector parallel to output Y vector
%      x    8.04675803e-02     9.94765884e-01     6.29746281e-02
%      y    7.83383109e-01    -1.02180405e-01     6.13082433e-01
%      z    6.16308262e-01     0.00000000e+00    -7.87505001e-01
%
%      Input vector parallel to output Z vector
%      x    6.16308262e-01     0.00000000e+00    -7.87505001e-01
%      y    8.04675803e-02     9.94765884e-01     6.29746281e-02
%      z    7.83383109e-01    -1.02180405e-01     6.13082433e-01
%
%-Particulars
%
%    Given an input vector x, this routine returns unit vectors x,
%    y, and z such that xyz forms a right-handed orthonormal frame
%    where the output x is parallel to the input x.
%
%    This routine is intended primarily to provide a basis for
%    the plane orthogonal to x.  There are no special properties
%    associated with y and z other than that the resulting xyz frame
%    is right handed and orthonormal. There are an infinite
%    collection of pairs (y,z) that could be used to this end.
%    Even though for a given x, y and z are uniquely determined, users
%    should regard the pair (y,z) as a random selection from this
%    infinite collection.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine frame_c.
%
%   MICE.REQ
%   FRAMES.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   build a right handed coordinate frame
%
%-&

function [x, y, z] = cspice_frame( vec )

   switch nargin
      case 1

         vec = zzmice_dp( vec );

      otherwise

         error ( 'Usage: [x(3), y(3), z(3)] = cspice_frame( vec(3) )' )

   end

   %
   % Call the MEX library.
   %
   try
      [x, y, z] = mice( 'frame_c', vec );
   catch
      rethrow(lasterror)
   end


