%-Abstract
%
%   CSPICE_GETFOV returns the field-of-view parameters for a user
%   specified instrument.
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
%      instid   NAIF ID for the instrument of interest.
%
%               [1,1] = size(instid); int32 = class(instid)
%
%      room     the max number of double precision 'bounds' vectors to return.
%
%               [1,1] = size(room); int32 = class(room)
%
%   the call:
%
%      [shape, frame, bsight, bounds] = cspice_getfov( instid, room)
%
%   returns:
%
%      shape    the FOV shape for instrument 'instid'. Possible values:
%
%               [1,m] = size(shape); char = class(shape)
%
%                    "POLYGON"
%                    "RECTANGLE"
%                    "CIRCLE"
%                    "ELLIPSE"
%
%      frame    the name of the frame in which the FOV is defined.
%
%               [1,m] = size(frame); char = class(frame)
%
%      bsight   the vector pointing in the direction of the FOV center
%               (boresight).
%
%               [3,1] = size(bsight); double = class(bsight)
%
%      bounds   set of vectors pointing to the "corners" of the instrument FOV,
%               i.e. 'bounds' returns as a N columns of 3-vectors.
%
%               [3,n] = size(bounds); double = class(bounds)
%
%               Note: do not consider 'bounds' as a matrix in the
%               conventional sense. Its 3XN form serves only as
%               a container for the bounds vectors.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   The example program in this section loads the IK file
%   'example.ti' with the following contents defining four FOVs of
%   various shapes and sizes:
%
%      KPL/IK
%
%      The keywords below define a circular, 10-degree wide FOV with
%      the boresight along the +Z axis of the 'SC999_INST001' frame
%      for an instrument with ID -999001 using the "angles"-class
%      specification.
%
%      \begindata
%         INS-999001_FOV_CLASS_SPEC       = 'ANGLES'
%         INS-999001_FOV_SHAPE            = 'CIRCLE'
%         INS-999001_FOV_FRAME            = 'SC999_INST001'
%         INS-999001_BORESIGHT            = ( 0.0, 0.0, 1.0 )
%         INS-999001_FOV_REF_VECTOR       = ( 1.0, 0.0, 0.0 )
%         INS-999001_FOV_REF_ANGLE        = ( 5.0 )
%         INS-999001_FOV_ANGLE_UNITS      = ( 'DEGREES' )
%      \begintext
%
%      The keywords below define an elliptical FOV with 2- and
%      4-degree angular extents in the XZ and XY planes and the
%      boresight along the +X axis of the 'SC999_INST002' frame for
%      an instrument with ID -999002 using the "corners"-class
%      specification.
%
%      \begindata
%         INS-999002_FOV_SHAPE            = 'ELLIPSE'
%         INS-999002_FOV_FRAME            = 'SC999_INST002'
%         INS-999002_BORESIGHT            = ( 1.0, 0.0, 0.0 )
%         INS-999002_FOV_BOUNDARY_CORNERS = ( 1.0, 0.0, 0.01745506,
%                                             1.0, 0.03492077, 0.0 )
%      \begintext
%
%      The keywords below define a rectangular FOV with 1.2- and
%      0.2-degree angular extents in the ZX and ZY planes and the
%      boresight along the +Z axis of the 'SC999_INST003' frame for
%      an instrument with ID -999003 using the "angles"-class
%      specification.
%
%      \begindata
%         INS-999003_FOV_CLASS_SPEC       = 'ANGLES'
%         INS-999003_FOV_SHAPE            = 'RECTANGLE'
%         INS-999003_FOV_FRAME            = 'SC999_INST003'
%         INS-999003_BORESIGHT            = ( 0.0, 0.0, 1.0 )
%         INS-999003_FOV_REF_VECTOR       = ( 1.0, 0.0, 0.0 )
%         INS-999003_FOV_REF_ANGLE        = ( 0.6 )
%         INS-999003_FOV_CROSS_ANGLE      = ( 0.1 )
%         INS-999003_FOV_ANGLE_UNITS      = ( 'DEGREES' )
%      \begintext
%
%      The keywords below define a triangular FOV with the boresight
%      along the +Y axis of the 'SC999_INST004' frame for an
%      instrument with ID -999004 using the "corners"-class
%      specification.
%
%      \begindata
%         INS-999004_FOV_SHAPE            = 'POLYGON'
%         INS-999004_FOV_FRAME            = 'SC999_INST004'
%         INS-999004_BORESIGHT            = (  0.0,  1.0,  0.0 )
%         INS-999004_FOV_BOUNDARY_CORNERS = (  0.0,  0.8,  0.5,
%                                              0.4,  0.8, -0.2,
%                                             -0.4,  0.8, -0.2 )
%      \begintext
%
%
%   The program shown below loads the IK, fetches parameters for each
%   of the four FOVs and prints these parameters to the screen.
%
%      function getfov_t()
%
%      %
%      % Set maximum number of boundary vectors, number of
%      % instruments and instrument IDs.
%      %
%      MAXBND = 4;
%      NUMINS = 4;
%      insids = [ -999001, -999002, -999003, -999004 ];
%
%      %
%      % Load the IK file.
%      %
%      cspice_furnsh ( 'example.ti' );
%
%      %
%      % For each instrument ...
%      %
%      fprintf ( '--------------------------------------\n' );
%      for i = 1:NUMINS
%
%         %
%         % ... fetch FOV parameters and ...
%         %
%         [shape, frame, bsight, bounds] = ...
%                  cspice_getfov( insids(i), MAXBND );
%
%         %
%         % ... print them to the screen.
%         %
%         fprintf ( 'Instrument ID: %i\n', insids(i) );
%         fprintf ( '    FOV shape: %s\n', shape );
%         fprintf ( '    FOV frame: %s\n', frame );
%         fprintf ( 'FOV boresight: %f %f %f\n', ...
%                   bsight(1), bsight(2), bsight(3) );
%
%         fprintf ( '  FOV corners: \n' );
%         [m,n] = size(bounds);
%         for j= 1:n
%            fprintf ( '               %f %f %f\n', ...
%                      bounds(1,j), bounds(2,j), bounds(3,j) );
%         end
%
%         fprintf ( '--------------------------------------\n' );
%
%      end
%
%   The program produces the following output:
%
%      --------------------------------------
%      Instrument ID: -999001
%          FOV shape: CIRCLE
%          FOV frame: SC999_INST001
%      FOV boresight: 0.000000 0.000000 1.000000
%        FOV corners:
%                     0.087156 0.000000 0.996195
%      --------------------------------------
%      Instrument ID: -999002
%          FOV shape: ELLIPSE
%          FOV frame: SC999_INST002
%      FOV boresight: 1.000000 0.000000 0.000000
%        FOV corners:
%                     1.000000 0.000000 0.017455
%                     1.000000 0.034921 0.000000
%      --------------------------------------
%      Instrument ID: -999003
%          FOV shape: RECTANGLE
%          FOV frame: SC999_INST003
%      FOV boresight: 0.000000 0.000000 1.000000
%        FOV corners:
%                     0.010472 0.001745 0.999944
%                     -0.010472 0.001745 0.999944
%                     -0.010472 -0.001745 0.999944
%                     0.010472 -0.001745 0.999944
%      --------------------------------------
%      Instrument ID: -999004
%          FOV shape: POLYGON
%          FOV frame: SC999_INST004
%      FOV boresight: 0.000000 1.000000 0.000000
%        FOV corners:
%                     0.000000 0.800000 0.500000
%                     0.400000 0.800000 -0.200000
%                     -0.400000 0.800000 -0.200000
%      --------------------------------------
%
%-Particulars
%
%   This routine provides a common interface to retrieving the geometric
%   characteristics of an instrument field of view for a wide variety of
%   remote sensing instruments across many different space missions.
%
%   Given the NAIF instrument ID, (and having "loaded" the
%   instrument field of view description via the routine cspice_furnsh)
%   this routine returns the bore-sight of the instrument, the
%   "shape" of the field of view, a collection of vectors
%   that point along the edges of the field of view, and the
%   name of the reference frame in which these vectors are defined.
%
%   Currently this routine supports two classes of specifications
%   for FOV definitions: "corners" and "angles".
%
%   The "corners" specification requires the following keywords
%   defining the shape, boresight, boundary vectors, and reference
%   frame of the FOV be provided in one of the text kernel files
%   (normally an IK file) loaded into the kernel pool (in the
%   keywords below <INSTID> is replaced with the instrument ID as
%   passed into the module):
%
%      INS<INSTID>_FOV_CLASS_SPEC         must be set to 'CORNERS' or
%                                         omitted to indicate the
%                                         "corners"-class
%                                         specification.
%
%
%      INS<INSTID>_FOV_SHAPE              must be set to one of these
%                                         values:
%
%                                            'CIRCLE'
%                                            'ELLIPSE'
%                                            'RECTANGLE'
%                                            'POLYGON'
%
%      INS<INSTID>_FOV_FRAME              must contain the name of
%                                         the frame in which the
%                                         boresight and boundary
%                                         corner vectors are defined.
%
%      INS<INSTID>_BORESIGHT              must be set to a 3D vector
%                                         defining the boresight in
%                                         the FOV frame specified in
%                                         the FOV_FRAME keyword.
%
%      INS<INSTID>_FOV_BOUNDARY   or
%      INS<INSTID>_FOV_BOUNDARY_CORNERS   must be set to one (for
%                                         FOV_SHAPE = 'CIRCLE'), two
%                                         (for FOV_SHAPE =
%                                         'ELLIPSE'), three (for
%                                         FOV_SHAPE = 'RECTANGLE'),
%                                         or three or more (for
%                                         'POLYGON') 3D vectors
%                                         defining the corners of the
%                                         FOV in the FOV frame
%                                         specified in the FOV_FRAME
%                                         keyword.
%
%   The "angles" specification requires the following keywords
%   defining the shape, boresight, reference vector, reference and
%   cross angular extents of the FOV be provided in one of the text
%   kernel files (normally an IK file) loaded into the kernel
%   pool (in the keywords below <INSTID> is replaced with the
%   instrument ID as passed into the module):
%
%      INS<INSTID>_FOV_CLASS_SPEC         must be set to  'ANGLES' to
%                                         indicate the "angles"-class
%                                         specification.
%
%      INS<INSTID>_FOV_SHAPE              must be set to one of these
%                                         values:
%
%                                            'CIRCLE'
%                                            'ELLIPSE'
%                                            'RECTANGLE'
%
%      INS<INSTID>_FOV_FRAME              must contain the name of
%                                         the frame in which the
%                                         boresight and the computed
%                                         boundary corner vectors are
%                                         defined.
%
%      INS<INSTID>_BORESIGHT              must be set to a 3D vector
%                                         defining the boresight in
%                                         the FOV frame specified in
%                                         the FOV_FRAME keyword.
%
%      INS<INSTID>_FOV_REF_VECTOR         must be set to a 3D vector
%                                         that together with the
%                                         boresight vector defines
%                                         the plane in which the
%                                         first angular extent of the
%                                         FOV specified in the
%                                         FOV_REF_ANGLE keyword is
%                                         measured.
%
%      INS<INSTID>_FOV_REF_ANGLE          must be set to the angle
%                                         that is 1/2 of the total
%                                         FOV angular extent in the
%                                         plane defined by the
%                                         boresight and the vector
%                                         specified in the
%                                         FOV_REF_VECTOR keyword.
%
%      INS<INSTID>_FOV_CROSS_ANGLE        must be set to the angle
%                                         that is 1/2 of the total
%                                         FOV angular extent in the
%                                         plane containing the
%                                         boresight and perpendicular
%                                         to the plane defined by the
%                                         boresight and the vector
%                                         specified in the
%                                         FOV_REF_VECTOR keyword.
%                                         This keyword is not
%                                         required for FOV_SHAPE =
%                                         'CIRCLE'.
%
%      INS<INSTID>_FOV_ANGLE_UNITS        must specify units for the
%                                         angles given in the
%                                         FOV_REF_ANGLE and
%                                         FOV_CROSS_ANGLE keywords.
%                                         Any angular units
%                                         recognized by cspice_convrt
%                                         are acceptable.
%
%   This routine is intended to be an intermediate level routine.
%   It is expected that users of this routine will be familiar
%   with the SPICE frames subsystem and will be comfortable writing
%   software to further manipulate the vectors retrieved by this
%   routine.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine getfov_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.3, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%       I/O descriptions edits to conform to Mice documentation format.
%
%   -Mice Version 1.0.2, 24-APR-2010, EDW (JPL)
%
%      Minor edit to code comments eliminating typo.
%
%   -Mice Version 1.0.1, 05-FEB-2009, BVS (JPL)
%
%      Header update: added information about required IK keywords;
%      replaced old example with a new one more focused on getfov_c and
%      IK keywords.
%
%   -Mice Version 1.0.0, 07-DEC-2005, EDW (JPL)
%
%-Index_Entries
%
%   return instrument's FOV parameters
%
%-&

function [shape, frame, bsight, bounds] = cspice_getfov( instid, room )

   switch nargin
      case 2

         instid = zzmice_int( instid );
         room   = zzmice_int( room   );

      otherwise

         error ( ['Usage: [`shape`, `frame`, bsight(3), bounds(3,N)] = ' ...
                          'cspice_getfov( instid, room)']  )

   end

   %
   % Call the MEX library.
   %
   try
      [shape, frame, bsight, bounds] = mice( 'getfov_c', instid, room );
   catch
      rethrow(lasterror)
   end
