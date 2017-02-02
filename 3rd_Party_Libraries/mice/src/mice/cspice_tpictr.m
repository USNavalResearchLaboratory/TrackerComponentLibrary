%-Abstract
%
%   CSPICE_TPICTR takes a time format sample then creates a
%   time format picture suitable for use by cspice_timout.
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
%      sample   a representative scalar time string to use
%               as a model to format time strings (the string need
%               not describe an actual date - only format matters)
%
%   the call:
%
%      [pictur, ok, errmsg] = cspice_tpictr(sample)
%
%   returns:
%
%      pictur   a scalar format picture string suitable for
%               use with the SPICE routine cspice_timout
%
%      ok       a scalar boolean indicating whether 'sample' parsed
%               without error (TRUE) or some parse error occurred
%               (FALSE)
%
%      errmsg   a scalar string containing the explanation of
%               the parse error
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Given a sample with the format of the UNIX date string,
%      % create a SPICE time picture for use in cspice_timout.
%      %
%      sample = 'Thu Oct 1 11:11:11 PDT 1111';
%
%      %
%      % Make the call. 'ok' returns false is an error occurred.
%      % The error description returns in the err variable.
%      %
%      [pictur, ok, errmsg] = cspice_tpictr( sample );
%
%      %
%      % If a false error flag, print the picture; if
%      % a true error flag, print the error message.
%      %
%      if ( ok )
%         disp( pictur )
%      else
%         disp( errmsg )
%      end
%
%   MATLAB outputs:
%
%      Wkd Mon DD HR:MN:SC PDT YYYY ::UTC-7
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine tpictr_c.
%
%   MICE.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.2.0, 10-MAY-2011, EDW (JPL)
%
%      "logical" call replaced with "zzmice_logical."
%
%   -Mice Version 1.0.1, 31-MAR-2010, EDW (JPL)
%
%      Renamed error message argument 'error' to 'errmsg'.
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   Use a sample time string to produce a time format picture
%
%-&

function [ pictur, ok, errmsg] = cspice_tpictr(sample)

   switch nargin
      case 1

         sample = zzmice_str(sample);

      otherwise

         error ( 'Usage: [ `pictur`, ok, `errmsg`] = cspice_tpictr( `sample`)' )

   end

   %
   % Call the MEX library.
   %
   try
      [pictur, ok, errmsg] =  mice('tpictr_c', sample );

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      ok = zzmice_logical(ok);
   catch
      rethrow(lasterror)
   end


