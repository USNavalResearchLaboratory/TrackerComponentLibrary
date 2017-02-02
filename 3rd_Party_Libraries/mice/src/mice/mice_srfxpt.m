%-Abstract
%
%   MICE_SRFXPT computes the surface intercept point of a specified ray
%   on a target body at a specified epoch, optionally corrected for light
%   time and stellar aberration, given an observer and a direction vector
%   defining a ray.
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
%      method   the scalar string providing the parameters to define the
%               computation method to use
%
%                  The only currently supported choice:
%
%                     "Ellipsoid"   The intercept computation uses
%                                   a triaxial ellipsoid to model
%                                   the surface of the target body.
%                                   The ellipsoid's radii must be
%                                   available in the kernel pool.
%
%      target   the scalar string name of the target body, 'target' being
%               case-insensitive, leading and trailing blanks are not
%               significant
%
%                  Optionally, you may supply a string containing the integer
%                  ID code for the object. For example both "MOON" and "301"
%                  are legitimate strings that indicate the moon is the
%                  target body.
%
%      et       the double precision, scalar or N-vector of epochs,
%               expressed as ephemeris seconds past J2000 TDB, at which to
%               compute the surface intercept point on the target body (this
%               epoch represents either the time of signal reception, or
%               transmission, depending on the selected 'abcorr')
%
%      abcorr   the scalar string name of the aberration correction to
%               apply when computing the observer-target state and the target
%               body orientation
%
%               For practical purposes, 'CN' (converged Newtonian)
%               represents the best correction choice.
%
%      obsrvr   the scalar string name of the observing body, 'obsrvr' being
%               case-insensitive, leading and trailing blanks are not
%               significant.
%
%                  Optionally, you may supply a string containing the integer
%                  ID code for the object. For example both "MOON" and "301"
%                  are legitimate strings that indicate the moon is the
%                  target body.
%
%      dref      the scalar string name of the reference frame containing the
%                'dvec' direction vector
%
%      dvec      the double precision 3-vector emanating from the observer
%
%                  'dvec' is specified relative to the reference frame
%                   designated by 'dref'.
%
%   the call:
%
%      [surf] = mice_srfxpt( method, target, et, abcorr, obsrvr, dref, dvec)
%
%   returns:
%
%      surf   the scalar or 1xN array of structures, each structure
%             consisting of five fields:
%
%                 spoint    the double precision 3-vector identifying the
%                           surface intercept point on 'target' of the ray
%                           'dvec' that emanates from 'obsrvr', with 'spoint'
%                           expressed in Cartesian coordinates relative to
%                           the body-fixed frame associated 'target'.
%
%                           The body-fixed target frame is evaluated at the
%                           epoch 'trgepc' NOT 'et'.
%
%                           The components of `spoint' are given in units of km.
%
%                  dist     the double precision, scalar distance in
%                           kilometers between the observer and surface
%                           intercept on the target body
%
%                 trgepc    the double precision, scalar "intercept epoch"
%                           expressed as ephemeris seconds past J2000 TDB where
%                          "intercept epoch" means the epoch at which the ray
%                           defined by 'obsrvr' and 'dvec' intercepts 'target'
%                           surface at 'spoint'
%
%                 obspos    the double precision 3-vector pointing from the
%                           center of 'target' at epoch 'trgepc' to 'obsrvr' at
%                           epoch 'et', with 'obspos' expressed in the target
%                           body-fixed reference frame
%
%                           The body-fixed target frame is evaluated at the
%                           epoch 'trgepc' NOT 'et'.
%
%                           The components of 'obspos' are given in units
%                           of km.
%
%                 found     a logical scalar indicating whether or not the ray
%                          'dvec' intersects 'target' (TRUE) or not (FALSE)
%
%             'surf' return with the same vectorization measure (N) as 'et'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Assign needed constants.
%      %
%      BUSID   = -94000;
%      MGS     = -94;
%      NCORNR  = 4;
%      ABCORR  = 'LT+S';
%      CAMERA  = 'MGS_MOC_NA';
%      DREF    = 'J2000';
%      METHOD  = 'ELLIPSOID';
%      OBSRVR  = 'MGS';
%      TARGET  = 'MARS';
%      UTC     = '2003 OCT 13 06:00:00 UTC';
%
%      %
%      %    Load kernel files:
%      %
%      %       - Leapseconds kernel
%      %       - MGS SCLK kernel
%      %       - Text PCK file
%      %       - Planetary SPK file
%      %       - MGS I-kernel
%      %       - MGS spacecraft bus C-kernel
%      %       - MGS SPK file
%      %
%      cspice_furnsh( '/kernels/gen/lsk/naif0008.tls' )
%      cspice_furnsh( '/kernels/MGS/sclk/MGS_SCLKSCET.00053.tsc' )
%      cspice_furnsh( '/kernels/MGS/pck/mars_iau2000_v0.tpc' )
%      cspice_furnsh( '/kernels/gen/spk/de405s.bsp' )
%      cspice_furnsh( '/kernels/MGS/ik/moc20.ti' )
%      cspice_furnsh( '/kernels/MGS/ck/mgs_sc_ext12.bc' )
%      cspice_furnsh( '/kernels/MGS/spk/mgs_ext12.bsp' )
%
%      %
%      % Convert the UTC request time to ET (seconds past
%      % J2000, TDB).
%      %
%      et = cspice_str2et( UTC );
%
%      %
%      % Get the MGS MOC Narrow angle camera (MGS_MOC_NA)
%      % ID code.  Then look up the field of view (FOV)
%      % parameters.
%      %
%      [camid, found] = cspice_bodn2c( CAMERA );
%
%      [shape, dref, bsight, bounds] = cspice_getfov( camid, NCORNR);
%
%      disp( ' ' )
%      disp( 'Surface Intercept Locations for Camera' )
%      disp( 'FOV Boundary and Boresight Vectors'     )
%      disp( ' ' )
%
%      txt = sprintf( '   Instrument:             %s', CAMERA);
%      disp( txt )
%
%      txt = sprintf( '   Epoch:                  %s', UTC);
%      disp( txt )
%
%      txt = sprintf( '   Aberration correction:  %s', ABCORR);
%      disp( txt )
%      disp( ' ' )
%
%      %
%      % Now compute and display the surface intercepts for the
%      % boresight and all of the FOV boundary vectors.
%      %
%      for i=1:NCORNR+1
%
%         if( i <= NCORNR )
%
%            %
%            % 'bounds' represents a 3 X NCORNR array with each row a bounds
%            % vector. Extract the vectors from 'bounds' using as a vector
%            % segment.
%            %
%            %    corner vector 0: bounds(:,1)
%            %    corner vector 1: bounds(:,2)
%            %    corner vector 2: bounds(:,3)
%            %    corner vector 3: bounds(:,4)
%            %
%            %
%            title = sprintf( 'Corner vector %d', i );
%            dvec = bounds(:,i);
%
%         else
%
%            title = sprintf( 'Boresight vector' );
%            dvec = bsight;
%
%         end
%
%         %
%         % Compute the surface intercept point using
%         % the specified aberration corrections.
%         %
%         [surf] = mice_srfxpt( METHOD, TARGET, et, ABCORR, OBSRVR, ...
%                               dref, dvec );
%
%         if( surf.found )
%
%            %
%            % Convert rectangular coordinates to planetocentric
%            % latitude and longitude.  Convert radians to degrees.
%            %
%            [ radius, lon, lat ] = cspice_reclat( surf.spoint );
%
%            lon = lon * cspice_dpr;
%            lat = lat * cspice_dpr;
%
%            %
%            % Display the results.
%            %
%            disp( title )
%            disp( ' ' )
%
%            txt = ...
%               sprintf('   Radius                   (km)  = %18.10e', radius );
%            disp( txt )
%
%            txt = ...
%               sprintf('   Planetocentric Latitude  (deg) = %18.10e', lat );
%            disp( txt )
%
%            txt = ...
%               sprintf('   Planetocentric Longitude (deg) = %18.10e', lon );
%            disp( txt )
%
%            txt = ...
%               sprintf('   Range                    (km)  = %18.10e', ...
%                                                               surf.dist );
%            disp( txt )
%            disp( ' ' )
%
%         else
%
%            disp( 'Intercept not found.' )
%
%         end
%
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Surface Intercept Locations for Camera
%      FOV Boundary and Boresight Vectors
%
%         Instrument:             MGS_MOC_NA
%         Epoch:                  2003 OCT 13 06:00:00 UTC
%         Aberration correction:  LT+S
%
%      Corner vector 1
%
%        Vector in MGS_MOC_NA frame =
%        1.8571383810e-06 -3.8015622659e-03  9.9999277403e-01
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849412615e+03
%           Planetocentric Latitude  (deg) =  -4.8477118861e+01
%           Planetocentric Longitude (deg) =  -1.2347365507e+02
%           Range                    (km)  =   3.8898362744e+02
%
%      Corner vector 2
%
%        Vector in MGS_MOC_NA frame =
%        1.8571383810e-06  3.8015622659e-03  9.9999277403e-01
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849398244e+03
%           Planetocentric Latitude  (deg) =  -4.8481272936e+01
%           Planetocentric Longitude (deg) =  -1.2339839939e+02
%           Range                    (km)  =   3.8897565851e+02
%
%      Corner vector 3
%
%        Vector in MGS_MOC_NA frame =
%       -1.8571383810e-06  3.8015622659e-03  9.9999277403e-01
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849398156e+03
%           Planetocentric Latitude  (deg) =  -4.8481298506e+01
%           Planetocentric Longitude (deg) =  -1.2339840260e+02
%           Range                    (km)  =   3.8897519958e+02
%
%      Corner vector 4
%
%        Vector in MGS_MOC_NA frame =
%       -1.8571383810e-06 -3.8015622659e-03  9.9999277403e-01
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849412527e+03
%           Planetocentric Latitude  (deg) =  -4.8477144435e+01
%           Planetocentric Longitude (deg) =  -1.2347365823e+02
%           Range                    (km)  =   3.8898316850e+02
%
%      Boresight vector
%
%        Vector in MGS_MOC_NA frame =
%        0.0000000000e+00  0.0000000000e+00  1.0000000000e+00
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849405358e+03
%           Planetocentric Latitude  (deg) =  -4.8479216591e+01
%           Planetocentric Longitude (deg) =  -1.2343603019e+02
%           Range                    (km)  =   3.8897626607e+02
%
%-Particulars
%
%   A sister version of this routine exists named cspice_srfxpt that returns
%   the structure field data as separate arguments.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine srfxpt_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 03-JAN-2006, EDW (JPL)
%
%-Index_Entries
%
%   surface intercept point
%
%-&

function [surf] = ...
          mice_srfxpt( method, target, et, abcorr, obsrvr, dref, dvec)

   switch nargin

      case 7

         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
         dref   = zzmice_str(dref);
         dvec   = zzmice_dp(dvec);

      otherwise

         error ( ['Usage: [_surf_ ] = '                ...
                  'mice_srfxpt( `method`, `target`,  ' ...
                  '_et_, `abcorr`, `obsrvr`, '         ...
                  '`dref`, dvec(6))']  )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [surf] = mice('srfxpt_s', method, target, et, abcorr, obsrvr, dref, dvec);
   catch
      rethrow(lasterror)
   end





