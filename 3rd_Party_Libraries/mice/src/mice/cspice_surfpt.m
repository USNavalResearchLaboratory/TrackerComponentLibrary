%-Abstract
%
%     CSPICE_SURFPT determines the intersection of a line-of-sight vector with
%     the surface of an ellipsoid.
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
%    Given:
%
%       Arguments-
%
%       positn   the double precision position 3-vector position
%                of an observer with respect to the center of an
%                ellipsoid expressed in the body fixed coordinates of
%                the ellipsoid
%
%       u        the double precision direction 3-vector emanating from
%                'positn'
%
%       a,       the scalar double precision values of the ellipsoid's
%       b,       triaxial radii ellipsoid, where:
%       c
%                   'a' is length in kilometers of the semi-axis of the
%                    ellipsoid parallel to the x-axis of the body-fixed
%                    reference frame
%
%                   'b' is length in kilometers of the semi-axis of the
%                    ellipsoid parallel to the y-axis of the body-fixed
%                    reference frame
%
%                   'c' is length in kilometers of the semi-axis of the
%                    ellipsoid parallel to the z-axis of the body-fixed
%                    reference frame
%
%    the call:
%
%       [point, found] = cspice_surfpt ( positn, u, a, b, c )
%
%
%    returns:
%
%       point   a double precision 3-vector defining the location
%               on the ellipsoid at which the 'u' intercepts
%               the ellipsoid if the interception exists, 'point' returns
%               (0.d, 0.d, 0.d) if 'u' does not intersect the ellipsoid
%
%       found   a boolean indicating whether the intersection
%               between the ellipse and 'u' exists (TRUE) or
%               not (FALSE)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      Suppose that MGS has taken a picture of Mars at time 'etrec' with
%      the MOC narrow angle camera. We want to know the latitude and
%      longitude associated with two pixels projected to Mars'
%      surface:  the boresight and one along the boundary of the
%      field of view (FOV). Due to light time, the photons taken in
%      the picture left Mars at time 'etemit', when Mars was at a
%      different state than at time 'etrec'.
%
%      In order to solve this problem, we could use the 'cspice_sincpt'
%      routine for both pixels, but this would be slow.  Instead, we
%      will assume that the light time for each pixel is the same. We
%      will call 'cspice_sincpt' once to get the light time and surface point
%      associated with the boresight. Then, we will rotate the first
%      FOV boundary vector from the camera frame at 'etrec' to the
%      body-fixed Mars frame at 'etemit', and call the faster routine
%      'cspice_surfpt' to retrieve the surface point for the FOV boundary
%      vector.
%
%      This example problem could be extended to find the latitude
%      and longitude associated with every pixel in an instrument's
%      field of view, but this example is simplified to only solve
%      for two pixels:  the boresight and one along the boundary of
%      the field of view.
%
%      Assumptions:
%
%         1)  The light times from the surface points in the camera's
%             field of view to the camera are equal.
%
%         2)  The camera offset from the center of gravity of the
%             spacecraft is zero. If the data are more accurate
%             and precise, this assumption can be easily discarded.
%
%         3)  An ellipsoid shape model for the target body is
%             sufficient.
%
%         4)  The boundary field of view vector returned from 'cspice_getfov'
%             is associated with a boundary field of view pixel. If
%             this example were extended to include a geometric camera
%             model, this assumption would not be needed since the
%             direction vectors associated with each pixel would be
%             calculated from the geometric camera model.
%
%       Example program starts here.
%
%         % --------------------- Program Setup ---------------------
%
%         metakr = 'mgs_ex.tm';
%         camera = 'MGS_MOC_NA';
%         NCORNR = 4;
%         ABCORR = 'CN+S';
%
%         %
%         % Load kernels
%         %
%         cspice_furnsh( metakr );
%
%         %
%         % Convert the time the picture was taken from a
%         % UTC time string to seconds past J2000, TDB.
%         %
%         etrec = cspice_str2et( '2003 OCT 13 06:00:00 UTC' );
%
%         %
%         % Assume the one-way light times from different
%         % surface points on Mars to MGS within the camera's
%         % FOV are equal. This means the photons that make
%         % up different pixels were all emitted from Mars at
%         % 'etemit' and received by MGS at 'etrec'.  It would be
%         % slow to process images using 'cspice_sincpt' for every
%         % pixel.  Instead, we will use 'cspice_sincpt' on the
%         % boresight pixel and use 'cspice_surfpt' for the first FOV
%         % boundary pixel.  If this example program were extended
%         % to include all of the camera's pixels, 'cspice_surfpt' would
%         % be used for the remaining pixels.
%         %
%         % Get the MGS MOC Narrow angle camera (MGS_MOC_NA)
%         % ID code. Then look up the field of view (FOV)
%         % parameters by calling 'cspice_getfov'.
%         %
%         [camid, found] = cspice_bodn2c( camera );
%
%         if ( ~found )
%             txt = sprintf( ['SPICE(NOTRANSLATION)' ...
%                             'Could not find ID code for instrument %s.' ], ...
%                             camera );
%             error( txt )
%         end
%
%         %
%         % 'cspice_getfov' will return the name of the camera-fixed frame
%         % in the string OBSREF, the camera boresight vector in
%         % the array BSIGHT, and the FOV corner vectors in the
%         % array BOUNDS.
%         %
%         [shape, obsref, bsight, bounds] = cspice_getfov( camid, NCORNR);
%
%         fprintf( '\nObservation Reference Frame:  %s\n', obsref );
%
%         %
%         % ----------- Boresight Surface Intercept -----------
%         %
%         % Retrieve the time, surface intercept point, and vector
%         % from MGS to the boresight surface intercept point
%         % in IAU_MARS coordinates.
%         %
%         [ spoint, etemit, srfvec, found ] = ...
%                 cspice_sincpt( 'Ellipsoid', 'Mars', etrec,  'IAU_MARS', ...
%                                 ABCORR,     'MGS' , obsref,  bsight );
%
%         if ( ~found )
%             txt = sprintf(['SPICE(NOINTERCEPT)' ...
%                            'Intercept not found for boresight vector.']);
%             error( txt )
%         end
%
%         %
%         % Convert the intersection point of the boresight
%         % vector and Mars from rectangular into latitudinal
%         % coordinates. Convert radians to degrees.
%         %
%         [ radius, lon, lat ] = cspice_reclat( spoint );
%
%         lon = lon * cspice_dpr;
%         lat = lat * cspice_dpr;
%
%         fprintf( ['\n'                                         ...
%                   'Boresight surface intercept coordinates:\n' ...
%                   '    Radius    (km) :  %f\n'                 ...
%                   '    Latitude  (deg):  %f\n'                 ...
%                   '    Longitude (deg):  %f\n' ],              ...
%                    radius, lat, lon );
%
%         %
%         % ------ 1st Boundary FOV Surface Intercept (cspice_surfpt) -----
%         %
%         % Now we will transform the first FOV corner vector into the
%         % IAU_MARS frame so the surface intercept point can be
%         % calculated using cspice_surfpt, which is faster than
%         % cspice_subpnt.
%         %
%         % If this example program were extended to include all
%         % of the pixels in the camera's FOV, a few steps, such as
%         % finding the rotation matrix from the camera frame to the
%         % IAU_MARS frame, looking up the semi-axis values for Mars,
%         % and finding the position of MGS with respect to Mars could
%         % be done once and used for every pixel.
%         %
%         % Find the rotation matrix from the ray's reference
%         % frame at the time the photons were received (etrec)
%         % to IAU_MARS at the time the photons were emitted
%         % (etemit).
%         %
%         [rotate] = cspice_pxfrm2( obsref, 'IAU_MARS', etrec, etemit );
%
%         %
%         % Look up the semi-axis values for Mars.
%         %
%         radii = cspice_bodvrd( 'MARS', 'RADII', 3 );
%
%         %
%         % Find the position of the center of Mars with respect
%         % to MGS.  The position of the observer with respect
%         % to Mars is required for the call to 'cspice_surfpt'.  Note:
%         % the apparent position of MGS with respect to Mars is
%         % not the same as the negative of Mars with respect to MGS.
%         %
%         pos_mgs_wrt_mars = spoint - srfvec;
%
%         %
%         % The first boundary FOV pixel must be rotated into the
%         % IAU_MARS reference frame.
%         %
%         bndvec = rotate * bounds(:,1);
%
%         %
%         % Calculate the surface point of the boundary FOV
%         % vector.
%         %
%         [surface_point, found] = cspice_surfpt ( pos_mgs_wrt_mars, ...
%                                                  bndvec, radii(1), ...
%                                                  radii(2), radii(3) );s
%
%         if ( ~found )
%             txt = 'SPICE(NOTFOUND)Could not calculate surface point.';
%             error( txt )
%         end
%
%         surf_point_using_surfpt = surface_point;
%
%         %
%         % Convert the intersection point of the boundary
%         % FOV vector and Mars from rectangular into
%         % latitudinal coordinates. Convert radians
%         % to degrees.
%         %
%         [ radius, lon, lat ] = cspice_reclat( surface_point );
%
%         lon = lon * cspice_dpr;
%         lat = lat * cspice_dpr;
%
%         fprintf( ['\n'                                        ...
%                   'Boundary vector surface intercept\n'       ...
%                   'coordinates using cspice_surfpt:\n'        ...
%                   '    Radius    (km) :  %f\n'                ...
%                   '    Latitude  (deg):  %f\n'                ...
%                   '    Longitude (deg):  %f\n'                ...
%                   '    Emit time using boresight LT (s):  %10.9f\n'], ...
%                    radius, lat, lon, etemit );
%
%         %
%         % ------ 1st Boundary FOV Surface Intercept Verification ----
%         %
%         % For verification only, we will calculate the surface
%         % intercept coordinates for the first boundary vector
%         % using 'cspice_sincpt' and compare to the faster
%         % 'cspice_surfpt' method.
%         %
%         [ surface_point, etemit, srfvec, found ] = ...
%                  cspice_sincpt( 'Ellipsoid', 'Mars', etrec,  'IAU_MARS', ...
%                                  ABCORR,     'MGS' , obsref,  bounds(:,1) );
%
%         if ( ~found )
%             txt = sprintf( ['SPICE(NOTFOUND)' ...
%                             'Intercept not found for the boundary ' ...
%                             'FOV vector.'] );
%             error( txt )
%         end
%
%         %
%         % Convert the intersection point of the first boundary
%         % vector and Mars from rectangular into latitudinal
%         % coordinates. Convert radians to degrees.
%         %
%         [ radius, lon, lat ] = cspice_reclat( surface_point );
%
%         lon = lon * cspice_dpr;
%         lat = lat * cspice_dpr;
%
%         fprintf( ['\n'                                              ...
%                   'Boundary vector surface intercept\n'             ...
%                   'coordinates using cspice_sincpt:\n'              ...
%                   '    Radius    (km) :  %f\n'                      ...
%                   '    Latitude  (deg):  %f\n'                      ...
%                   '    Longitude (deg):  %f\n'                      ...
%                   '    Emit time using boresight LT (s):  %10.9f\n\n'], ...
%                    radius, lat, lon, etemit );
%
%         distance = cspice_vdist ( surf_point_using_surfpt, surface_point );
%
%         fprintf( [ 'Distance between surface points of the first\n' ...
%                    'boundary vector using cspice_surfpt and\n'      ...
%                    'cspice_sincpt:\n'                               ...
%                    '    Distance (mm):   %f\n' ],                   ...
%                    distance*10^6 );
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in MATLAB due to data persistence.
%         %
%         cspice_kclear
%
%
%   MATLAB outputs:
%
%             Observation Reference Frame:  MGS_MOC_NA
%
%             Boresight surface intercept coordinates:
%                 Radius    (km) :  3384.940410
%                 Latitude  (deg):  -48.479580
%                 Longitude (deg):  -123.436454
%
%             Boundary vector surface intercept
%             coordinates using cspice_surfpt:
%                 Radius    (km) :  3384.941136
%                 Latitude  (deg):  -48.477482
%                 Longitude (deg):  -123.474080
%                 Emit time using boresight LT (s):  119296864.181059480
%
%             Boundary vector surface intercept
%             coordinates using cspice_sincpt:
%                 Radius    (km) :  3384.941136
%                 Latitude  (deg):  -48.477482
%                 Longitude (deg):  -123.474079
%                 Emit time using boresight LT (s):  119296864.181059465
%
%             Distance between surface points of the first
%             boundary vector using cspice_surfpt and
%             cspice_sincpt:
%                 Distance (mm):   32.139880
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine surfpt_c.
%
%   MICE.REQ
%   ELLIPSES.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 24-OCT-2011, SCK (JPL)
%
%-Index_Entries
%
%   GF distance search
%
%-&

function [point, found] = cspice_surfpt ( positn, u, a, b, c )

   switch nargin

      case 5

         positn  = zzmice_dp(positn);
         u       = zzmice_dp(u);
         a       = zzmice_dp(a);
         b       = zzmice_dp(b);
         c       = zzmice_dp(c);

      otherwise

         error ( [ 'Usage: [point, found] = cspice_surfpt( `positn`, ' ...
                                                  '`u`, `a`, `b`, c )' ] )

   end

   %
   % Call surfpt_c.
   %
   try
      [surfpt] = mice('surfpt_s', positn, u, a, b, c );

      point  = reshape( [surfpt.spoint], 3, [] );
      found  = reshape( [surfpt.found ], 1, [] );
   catch
      rethrow(lasterror)
   end



















































