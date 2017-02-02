%-Abstract
%
%   CSPICE_STR2ET converts a string representing an epoch to a
%   double precision value representing the number of TDB seconds
%   past the J2000 epoch corresponding to the input epoch.
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
%      str   any scalar or NxM character array of strings recognized by
%            SPICE as an epoch.
%
%   the call:
%
%      et = cspice_str2et( str )
%
%   returns:
%
%      et   the scalar or 1XN-vector of double precision number of
%           TDB seconds past the J2000 epoch that corresponds to
%           the input 'str'.
%
%           'et' returns with the same vectorization measure (N) as 'str'.
%
%   Note: Reference the function cspice_tsetyr for information concerning
%   the translation of two digit representations of the century count.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Load a leapseconds kernel.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Define the epoch as a string.
%      %
%      date = 'Thu Mar 20 12:53:29 PST 1997';
%
%      %
%      % Convert a string to ephemeris time (ET).
%      %
%      et = cspice_str2et( date );
%
%      disp( 'Scalar:' )
%      txt = sprintf( '%20.8f', et );
%      disp( txt )
%
%      disp( ' ' )
%
%      %
%      % Define a vector of time strings:
%      %
%      time = strvcat( 'JD2454000.', ...
%                      'JD2464000.', ...
%                      'JD2474000.', ...
%                      'JD2484000.', ...
%                      'JD2494000.' );
%
%      %
%      % Convert the array of time strings 'time' to
%      % and array of ephemeris times 'et'.
%      %
%      et = cspice_str2et( time );
%
%      disp( 'Vector:' )
%      fprintf( '%20.8f\n', et' );
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%  MATLAB outputs:
%
%     Scalar:
%       -87836728.81438904
%
%     Vector:
%       212112064.18239054
%      1076112064.18491936
%      1940112064.18430591
%      2804112064.18263292
%      3668112064.18564129
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine str2et_c.
%
%   MICE.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%    Convert a string to TDB seconds past the J2000 epoch
%
%-&

function [et] = cspice_str2et(str)

   switch nargin
      case 1

         str = zzmice_str(str);

      otherwise

         error ( 'Usage: [_et_] = cspice_str2et(_`str`_)' )

   end

   try
      [et] = mice('str2et_c',str);
   catch
      rethrow(lasterror)
   end


