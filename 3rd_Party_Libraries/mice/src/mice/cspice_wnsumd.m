%-Abstract
%
%   CSPICE_WNSUMD summarizes the contents of a double precision window.
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
%      window   SPICE window containing zero or more intervals.
%
%               [2n,1] = size(window); double = class(window)
%
%   the call:
%
%      [ meas, avg, stddev, shortest, longest ] = cspice_wnsumd( window )
%
%   returns:
%
%      meas        total measure of the intervals in the input  window. This is
%                  just the sum of the measures of the individual intervals.
%
%                  [1,1] = size(meas); double = class(meas)
%
%      avg         average measure of the intervals in the input window.
%
%                  [1,1] = size(avg); double = class(avg)
%
%      stddev      standard deviation of the measures of the intervals in the
%                  input window.
%
%                  [1,1] = size(stddev); double = class(stddev)
%
%      shortest,
%      longest     indices of the left endpoint of, respectively, the shortest
%                  and longest intervals in the data contained in 'window'.
%
%                  [1,1] = size(longest); int32 = class(longest)
%
%                     The following algorithm describes the relation of
%                     'shortest' and 'longest' to the window data:
%
%                     The shortest interval:
%
%                        [ window(shortest), window(shortest+1) ]
%
%                     The longest interval:
%
%                        [ window(longest), window(longest+1) ]
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Define an array representing a window with six intervals.
%      % The values in 'window' have correct order for a
%      % SPICE window.
%      %
%      window = [ [  1.;  3.]; ...
%                 [  7.; 11.]; ...
%                 [ 18.; 18.]; ...
%                 [ 23.; 27.]; ...
%                 [ 30.; 69.]; ...
%                 [ 72.; 80.] ];
%
%      %
%      % Calculate the summary for 'window'.
%      %
%      [ meas, avg, stddev, shortest, longest] = cspice_wnsumd( window );
%
%      %
%      % 'shortest' and 'longest' refer to the indices of
%      % the 'cell' data array.
%      %
%
%      intrvl_short= (shortest+1)/2;
%      intrvl_long = (longest+1)/2;
%
%      fprintf( 'Measure           : %f\n', meas        )
%      fprintf( 'Average           : %f\n', avg         )
%      fprintf( 'Standard Dev      : %f\n', stddev      )
%      fprintf( 'Index shortest    : %i\n', shortest    )
%      fprintf( 'Index longest     : %i\n', longest     )
%      fprintf( 'Interval shortest : %i\n', intrvl_short)
%      fprintf( 'Interval longest  : %i\n', intrvl_long )
%
%      fprintf( 'Shortest interval : %f %f\n', window(shortest),  ...
%                                              window(shortest+1) )
%
%      fprintf( 'Longest interval  : %f %f\n', window(longest), ...
%                                              window(longest+1) )
%
%   MATLAB outputs:
%
%      Measure           : 57.000000
%      Average           : 9.500000
%      Standard Dev      : 13.413302
%      Index shortest    : 5
%      Index longest     : 9
%      Interval shortest : 3
%      Interval longest  : 5
%      Shortest interval : 18.000000 18.000000
%      Longest interval  : 30.000000 69.000000
%
%-Particulars
%
%   This routine provides a summary of the input window, consisting
%   of the following items:
%
%      - The measure of the window.
%
%      - The average and standard deviation of the measures
%        of the individual intervals in the window.
%
%      - The indices of the left endpoints of the shortest
%        and longest intervals in the window.
%
%   All of these quantities are zero if the window contains no
%   intervals.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine wnsumd_c.
%
%   MICE.REQ
%   WINDOWS.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.0, 15-DEC-2008, EDW (JPL)
%
%-Index_Entries
%
%   summary of a d.p. window
%
%-&

function  [ meas, avg, stddev, shortest, longest ] = cspice_wnsumd( window )

   switch nargin

      case 1

         window = zzmice_win( window);

      otherwise

         error( [ '[ meas, avg, stddev, shortest, longest ] ' ...
                  ' = cspice_wnsumd( window )' ] )

      end

   try
      [sumd] = mice( 'wnsumd_s',  [zeros(6,1);window] );

      meas     = reshape( [sumd.meas],     1, [] );
      avg      = reshape( [sumd.avg],      1, [] );
      stddev   = reshape( [sumd.stddev],   1, [] );
      shortest = reshape( [sumd.shortest], 1, [] );
      longest  = reshape( [sumd.longest],  1, [] );

   catch
      rethrow(lasterror)
   end
