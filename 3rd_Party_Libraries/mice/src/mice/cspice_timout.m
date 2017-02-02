%-Abstract
%
%   CSPICE_TIMOUT converts an input epoch represented in TDB seconds
%   past the TDB epoch of J2000 to a character string formatted to
%   the specifications of a user's format picture.
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
%     et       a double precision scalar or 1xN array of time values
%              in seconds past the ephemeris epoch J2000.
%
%     pictur   a scalar string that specifies how the output should be
%              presented.  The string is made up of various markers
%              that stand for various components associated with
%              a time.
%
%   the call:
%
%      output = cspice_timout( et, pictur )
%
%   returns:
%
%      output   the scalar string or NxM character array of output time strings
%               equivalent to the input epoch 'et' in the format specified
%               by 'pictur'
%
%               'output' returns with the same vectorization measure (N)
%                as 'et'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Given a sample with the format of the UNIX date string
%      % local to California, create a SPICE time picture for use
%      % in cspice_timout.
%      %
%      sample = 'Thu Oct 1 11:11:11 PDT 1111';
%
%      %
%      % Load a leapseconds kernel file.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Create the pic string.
%      %
%      [ pic, ok, xerror ] = cspice_tpictr( sample );
%
%      %
%      % Check the error flag, 'ok', for problems.
%      %
%      if ( ~ok  )
%         error( xerror )
%      end
%
%      %
%      % Convert an ephemeris time to the 'pic' format.
%      %
%      % Using the ET representation for: Dec 25 2005, 1:15:00 AM UTC
%      %
%      et = 188745364.;
%
%      output = cspice_timout( et, pic );
%
%      disp( 'Scalar: ' )
%
%      txt = sprintf( 'ET              : %16.8f', et );
%      disp( txt )
%
%      disp( ['Converted output: ' output] )
%      disp( ' ' )
%
%      %
%      % Create an array of ephemeris times beginning
%      % at 188745364 with graduations of 10000.0
%      % ephemeris seconds.
%      %
%      et=[0:4] * 10000. + 188745364;
%
%      %
%      % Convert the array of ephemeris times 'et' to an
%      % array of time strings, 'output', in 'pic' format.
%      %
%      output = cspice_timout( et, pic );
%
%      disp( 'Vector:' )
%      for i=1:5
%         txt = sprintf( 'ET              : %16.8f', et(i) );
%         disp( txt)
%
%         disp( ['Converted output: ' output(i,:) ]  )
%         disp( ' ' )
%      end
%
%      %
%      %  It's always good form to unload kernels after use,
%      %  particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Scalar:
%      ET              : 188745364.00000000
%      Converted output: Sat Dec 24 18:14:59 PDT 2005
%
%      Vector:
%      ET              : 188745364.00000000
%      Converted output: Sat Dec 24 18:14:59 PDT 2005
%
%      ET              : 188755364.00000000
%      Converted output: Sat Dec 24 21:01:39 PDT 2005
%
%      ET              : 188765364.00000000
%      Converted output: Sat Dec 24 23:48:19 PDT 2005
%
%      ET              : 188775364.00000000
%      Converted output: Sun Dec 25 02:34:59 PDT 2005
%
%      ET              : 188785364.00000000
%      Converted output: Sun Dec 25 05:21:39 PDT 2005
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine timout_c.
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
%   Convert and format d.p. seconds past J2000 as a string
%
%-&

function [output] = cspice_timout( et, pictur )

   switch nargin
      case 2

         et     = zzmice_dp(et);
         pictur = zzmice_str(pictur);

      otherwise

         error( 'Usage: [_`output`_] = cspice_timout( _et_, `pictur` )' )

   end

   %
   % Call the MEX library.
   %
   try
      [output] = mice('timout_c', et, pictur);
   catch
      rethrow(lasterror)
   end




