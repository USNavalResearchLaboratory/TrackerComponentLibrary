%-Abstract
%
%   CSPICE_XFMSTA transforms a state between coordinate systems.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
%   TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
%   WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
%   PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
%   BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
%   LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
%   INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
%   REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
%   REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
%   THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
%   CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
%   ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%
%      Arguments-
%
%      input_state       is a state vector in the input ('input_coord_sys')
%                        coordinate system representing position and velocity.
%
%                        [6,n] = size(input_state); double = class(input_state)
%
%                        All angular measurements must be in radians.
%
%                        Note: body radii values taken from the kernel
%                        pool are used when converting to or from geodetic or
%                        planetographic coordinates. It is the user's
%                        responsibility to verify the distance inputs are in
%                        the same units as the radii in the kernel pool,
%                        typically kilometers.
%
%      input_coord_sys   is the name of the coordinate system that the input
%                        state vector ('input_state') is currently in.
%
%                        [1,c1] = size(input_coord_sys)
%                        char = class(input_coord_sys)
%
%                          or
%
%                        [1,1] = size(output_coord_sys)
%                        cell = class(output_coord_sys)
%
%                        'input_coord_sys' may be any of the following:
%
%                              'RECTANGULAR'
%                              'CYLINDRICAL'
%                              'LATITUDINAL'
%                              'SPHERICAL'
%                              'GEODETIC'
%                              'PLANETOGRAPHIC'
%
%                        Leading spaces, trailing spaces, and letter case
%                        are ignored. For example, ' cyLindRical  ' would
%                        be accepted.
%
%      output_coord_sys  is the name of the coordinate system that the state
%                        should be converted to. Please see the description of
%                        'input_coord_sys' for details.
%
%                        [1,c2] = size(output_coord_sys)
%                        char = class(output_coord_sys)
%
%                          or
%
%                        [1,1] = size(output_coord_sys)
%                        cell = class(output_coord_sys)
%
%      body              is the name or NAIF ID of the body associated with the
%                        planetographic or geodetic coordinate system.
%
%                        [1,c3] = size(body); char = class(body)
%
%                          or
%
%                        [1,1] = size(body); cell = class(body)
%
%                        If one of the coordinate system choices is not
%                        geodetic or planetographic, 'body' may be an empty
%                        string (' ').
%
%                        Examples of accepted body names or IDs are:
%                                 'Earth'
%                                 '399'
%
%                        Leading spaces, trailing spaces, and letter case are
%                        ignored.
%
%   the call:
%
%      output_state = cspice_xfmsta ( input_state,      input_coord_sys, ...
%                                     output_coord_sys, body )
%
%   returns:
%
%      output_state      is the state vector that has been converted to the
%                        output coordinate system ('output_coord_sys').
%
%                        [6,n] = size(output_state)
%                        double = class(output_state)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%
%   Example(1):
%
%      Find the apparent state of Phoebe as seen by CASSINI in the
%      J2000 frame at three times starting at 2004 Jun 11 19:32:00.
%      Transform the state from rectangular to latitudinal coordinates.
%      For verification, transform the state back from latitudinal to
%      rectangular coordinates.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%         KPL/MK
%
%         File name: xfmsta_ex1.tm
%
%         This meta-kernel is intended to support operation of SPICE
%         example programs. The kernels shown here should not be
%         assumed to contain adequate or correct versions of data
%         required by SPICE-based user applications.
%
%         In order for an application to use this meta-kernel, the
%         kernels referenced here must be present in the user's
%         current working directory.
%
%         The names and contents of the kernels referenced
%         by this meta-kernel are as follows:
%
%                File name                     Contents
%                ---------                     --------
%                cpck05Mar2004.tpc             Planet orientation and
%                                              radii
%                naif0009.tls                  Leapseconds
%                020514_SE_SAT105.bsp          Satellite ephemeris for
%                                              Saturn
%                030201AP_SK_SM546_T45.bsp     CASSINI ephemeris
%                981005_PLTEPH-DE405S.bsp      Planetary ephemeris
%
%
%         \begindata
%
%         KERNELS_TO_LOAD = ( 'naif0009.tls'  ,
%                             '020514_SE_SAT105.bsp'  ,
%                             '030201AP_SK_SM546_T45.bsp'  ,
%                             '981005_PLTEPH-DE405S.bsp',
%                             'cpck05Mar2004.tpc'   )
%
%         End of meta-kernel
%
%       Example program starts here.
%
%         metakr = 'xfmsta_ex1.tm';
%         form = 'YYYY-MM-DD HR:MN:SC (TDB) ::TDB';
%
%         %
%         %   Load the meta kernel.
%         %
%         cspice_furnsh ( metakr );
%
%         %
%         %   Calculate the state at 2004 Jun 11 19:32:00 UTC.
%         %
%         times = [ '2004-JUN-11-19:32:00';
%                   '2004-JUN-11-19:40:00';
%                   '2004-JUN-11-19:48:00'];
%         et = cspice_str2et ( times );
%
%         %
%         %   Calculate the apparent state of Phoebe as seen by
%         %   CASSINI in the J2000 frame.
%         %
%         [state_rec, lt] = cspice_spkezr ( 'phoebe', et, 'iau_phoebe', ...
%                                           'lt+s', 'cassini' );
%
%         %
%         %   Transform the state from rectangular to latitudinal.
%         %   Notice that since neither the input nor output
%         %   coordinate frames are 'geodetic' or 'planetographic',
%         %   the input for the body name is a blank string.
%         %
%         state_lat = cspice_xfmsta ( state_rec, 'rectangular', ...
%                                     'latitudinal', ' ' );
%
%         %
%         %   Transform the state back to rectangular from latitudinal. The
%         %   result should be very close to 'state_rec'.
%         %
%         state_rec2 = cspice_xfmsta ( state_lat, 'latitudinal', ...
%                                      'rectangular', ' ');
%
%         %
%         %   Report the results.
%         %
%         fprintf('\nPhoebe as seen by Cassini - rectangular\n')
%         for i = 1:length(et)
%             fprintf('Time: %s\n', cspice_timout ( et(i), form))
%             fprintf('  Position: %16.6f %16.6f %16.6f\n', state_rec(1:3,i))
%             fprintf('  Velocity: %16.6f %16.6f %16.6f\n', state_rec(4:6,i))
%         end
%
%         fprintf('\nPhoebe as seen by Cassini - latitudinal\n')
%         for i = 1:length(et)
%             fprintf('Time: %s\n', cspice_timout ( et(i), form))
%             fprintf('  Position: %16.6f %16.6f %16.6f\n', state_lat(1:3,i))
%             fprintf('  Velocity: %16.6f %16.6f %16.6f\n', state_lat(4:6,i))
%         end
%
%         fprintf('\nVerification: Phoebe as seen by Cassini - rectangular\n')
%         for i = 1:length(et)
%             fprintf('Time: %s\n', cspice_timout ( et(i), form))
%             fprintf('  Position: %16.6f %16.6f %16.6f\n', state_rec2(1:3,i))
%             fprintf('  Velocity: %16.6f %16.6f %16.6f\n', state_rec2(4:6,i))
%         end
%
%         %
%         %   Unload the kernels.
%         %
%         cspice_kclear
%
%   MATLAB outputs:
%
%         Phoebe as seen by Cassini - rectangular
%         Time: 2004-06-11 19:33:04 (TDB)
%           Position:     -1982.639762      -934.530471      -166.562595
%           Velocity:         3.970832        -3.812496        -2.371663
%         Time: 2004-06-11 19:41:04 (TDB)
%           Position:      -257.857469     -2932.271650     -1304.929678
%           Velocity:         3.200173        -4.493939        -2.371569
%         Time: 2004-06-11 19:49:04 (TDB)
%           Position:      1075.597532     -5231.009547     -2443.280908
%           Velocity:         2.342482        -5.064820        -2.371563
%
%         Phoebe as seen by Cassini - latitudinal
%         Time: 2004-06-11 19:33:04 (TDB)
%           Position:      2198.169858        -2.701121        -0.075846
%           Velocity:        -1.780939         0.002346        -0.001144
%         Time: 2004-06-11 19:41:04 (TDB)
%           Position:      3219.867849        -1.658508        -0.417279
%           Velocity:         4.797399         0.001217        -0.000145
%         Time: 2004-06-11 19:49:04 (TDB)
%           Position:      5872.818108        -1.368003        -0.429078
%           Velocity:         5.926982         0.000239         0.000018
%
%         Verification: Phoebe as seen by Cassini - rectangular
%         Time: 2004-06-11 19:33:04 (TDB)
%           Position:     -1982.639762      -934.530471      -166.562595
%           Velocity:         3.970832        -3.812496        -2.371663
%         Time: 2004-06-11 19:41:04 (TDB)
%           Position:      -257.857469     -2932.271650     -1304.929678
%           Velocity:         3.200173        -4.493939        -2.371569
%         Time: 2004-06-11 19:49:04 (TDB)
%           Position:      1075.597532     -5231.009547     -2443.280908
%           Velocity:         2.342482        -5.064820        -2.371563
%
%   Example(2):
%
%      Example program starts here.
%
%         %   Initialize the cylindrical state.
%         %
%         state_cyl = [1 0.5 0.5 0.2 0.1 -0.2]';
%
%         %
%         %   Load kernels.
%         %
%         cspice_furnsh( 'cpck05Mar2004.tpc' );
%
%         %
%         %   Transform the state from cylindrical to planetographic.
%         %   Note that since one of the coordinate systems is
%         %   planetographic, the body name must be input.
%         %
%         state_plan = cspice_xfmsta ( state_cyl, 'cylindrical',  ...
%                                      'planetographic', 'earth' );
%
%         %
%         %   Transform the state back to cylindrical from
%         %   planetographic. The result should be very close to 'state_cyl'.
%         %
%         state_cyl2 = cspice_xfmsta ( state_plan, 'planetographic',...
%                                      'cylindrical', 'earth' );
%
%         %
%         %   Report the results.
%         %
%         fprintf('\nCylindrical State\n')
%         fprintf('  Position: %7.3f %7.3f %7.3f\n', state_cyl(1:3))
%         fprintf('  Velocity: %7.3f %7.3f %7.3f\n', state_cyl(4:6))
%
%         fprintf('\nPlanetographic State\n')
%         fprintf('  Position: %7.3f %7.3f %7.3f\n', state_plan(1:3))
%         fprintf('  Velocity: %7.3f %7.3f %7.3f\n', state_plan(4:6))
%
%         fprintf('\nVerification: Cylindrical State\n')
%         fprintf('  Position: %7.3f %7.3f %7.3f\n', state_cyl2(1:3))
%         fprintf('  Velocity: %7.3f %7.3f %7.3f\n', state_cyl2(4:6))
%
%         %
%         %   Unload the kernels.
%         %
%         cspice_kclear
%
%   MATLAB outputs:
%
%         Cylindrical State
%           Position:   1.000   0.500   0.500
%           Velocity:   0.200   0.100  -0.200
%
%         Planetographic State
%           Position:   0.500   1.548 -6356.238
%           Velocity:   0.100  -0.005  -0.195
%
%         Verification: Cylindrical State
%           Position:   1.000   0.500   0.500
%           Velocity:   0.200   0.100  -0.200
%
%-Particulars
%
%   Input Order
%   -------------------------------------------
%
%      The input and output states will be structured by the
%      following descriptions.
%
%      For rectangular coordinates, the state vector is the following
%      in which X, Y, and Z are the rectangular position components and
%      DX, DY, and DZ are the time derivatives of each.
%              ISTATE = ( X, Y, Z, DX, DY, DZ )
%
%      For cylindrical coordinates, the state vector is the following
%      in which R is the radius, LONG is the longitude, Z is the
%      height, and DR, DLONG, and DZ are the time derivatives of each.
%              ISTATE = ( R, LONG, Z, DR, DLONG, DZ )
%
%      For latitudinal coordinates, the state vector is the following
%      in which R is the radius, LONG is the longitude, LAT is the
%      latitude, and DR, DLONG, and DLAT are the time derivatives of
%      each.
%              ISTATE = ( R, LONG, LAT, DR, DLONG, DLAT )
%
%      For spherical coordinates, the state vector is the following in
%      which R is the radius, COLAT is the colatitude, LONG is the
%      longitude, and DR, DCOLAT, and DLONG are the time derivatives of
%      each.
%              ISTATE = ( R, COLAT, LONG, DR, DCOLAT, DLONG )
%
%      For geodetic coordinates, the state vector is the following in
%      which LONG is the longitude, LAT is the latitude, ALT is the
%      altitude, and DLONG, DLAT, and DALT are the time derivatives of
%      each.
%              ISTATE = ( LONG, LAT, ALT, DLONG, DLAT, DALT )
%
%      For planetographic coordinates, the state vector is the
%      following in which LONG is the longitude, LAT is the latitude,
%      ALT is the altitude, and DLONG, DLAT, and DALT are the time
%      derivatives of each.
%              ISTATE = ( LONG, LAT, ALT, DLONG, DLAT, DALT )
%
%   Input Boundaries
%   -------------------------------------------
%
%      There are intervals the input angles must fall within if
%      the input coordinate system is not rectangular.  These
%      intervals are provided below.
%
%         Input variable    Input meaning   Input interval [rad]
%         --------------    -------------   ------------------------
%             LONG           Longitude        0     <= LONG  <  2*pi
%             LAT            Latitude        -pi/2  <= LAT   <= pi/2
%             COLAT          Colatitude       0     <= COLAT <= pi
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine xfmsta_c.
%
%-Version
%
%   -Mice Version 1.0.0, 28-AUG-2012, SCK (JPL), EDW (JPL)
%
%-Index_Entries
%
%   state transformation between coordinate systems
%   convert state
%
%-&

function output_state = cspice_xfmsta ( input_state,      input_coord_sys, ...
                                        output_coord_sys, body )

   switch nargin
      case 4

         input_state      = zzmice_dp (input_state);
         input_coord_sys  = zzmice_str(input_coord_sys);
         output_coord_sys = zzmice_str(output_coord_sys);
         body             = zzmice_str(body);

      otherwise

         error ( ['Usage: [_output_state(6)_] = ' ...
                  'cspice_xfmsta( `_input_state(6)_`, `input_coord_sys`, ' ...
                  '`output_coord_sys`, `body`)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [output_state] = mice('xfmsta_c', input_state, input_coord_sys,...
                                        output_coord_sys, body);
   catch
      rethrow(lasterror)
   end









