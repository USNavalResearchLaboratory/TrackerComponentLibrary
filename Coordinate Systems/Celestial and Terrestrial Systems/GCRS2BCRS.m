function vec=GCRS2BCRS(x,Jul1,Jul2,deltaTTUT1)
%%GCRS2BCRS  Convert vectors of position and possibly velocity from the
%            Geocentric Celestrial Reference System (GCRS), a type of
%            Earth-Centered Inertial (ECI) coordinate system, to the
%            Barycentric Celestial Reference System (BCRS).
%
%INPUTS: x     The NXnumVec collection of vectors to convert. N can be 3,
%              or 6. If the vectors are 3D, then they are position. 6D
%              vectors are assumed to be position and velocity.
%Jul1, Jul2    Two parts of a Julian date given in terrestrial time (TT).
%              The units of the date are days. The full date is the sum of
%              both terms. The date is broken into two parts to provide
%              more bits of precision. It does not matter how the date is
%              split.
%deltaTTUT1    An optional parameter specifying the difference between TT
%              and UT1 in seconds. This information can be obtained from
%http://www.iers.org/nn_11474/IERS/EN/DataProducts/EarthOrientationData/eop.html?__nnn=true
%              or 
%http://www.usno.navy.mil/USNO/earth-orientation/eo-products
%              If this parameter is omitted or if an empty matrix is
%              passed, then the value provided by the function getEOP
%              will be used instead. This parameter is used in the function
%              for the transformation to TDB.
%
%OUTPUTS: vec A 3XN or 6XN matrix of vectors converted from GCRS
%             coordinates to BCRS coordinates. If a 6XN matrix, the
%             velocity components are unchanged, since the BCRS is just a
%             change in origin compared to the GCRS.
%
%The BCRS is just offset in origin from the GCRS. The National Aeronautics
%and Space Administration's (NASA) Navigation and Ancillary Information
%Facility's (NAIF's) SPICE toolkit's function cspice_spkezr is used to get
%the location of Earth with respect to the solar system barycenter to
%perform the conversion. However, the main sticking point in a conversion
%to the BCRS is the timescale used. The function cspice_spkezr requires a
%time in barycentric dynamical time (TDB). The terrestrial time here is
%converted to TDB using TT2TDB without topocentric corrections. That is,
%the clock is just assumed to be located at the center of the Earth.
%
%Different celestial coordinate systems are compared in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "An overview of major terrestrial, celestial, and
%    temporal coordinate systems for target tracking", Report, U. S. Naval
%    Research Laboratory, to appear, 2016.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4)
   deltaTTUT1=[]; 
end

vec=zeros(size(x));

switch(size(x,1))
    case 3
    case 6%The velocity is not transformed.
        vec(4:6,:)=x(4:6,:);
    otherwise
        error('The input vector x has the wrong dimensionality.')
end

[TDB1,TDB2]=TT2TDB(Jul1,Jul2,deltaTTUT1,[0;0;0]);
%Convert to seconds past J2000.0 as a single double-precision floating
%point number. The date of J2000.0 in TDB is 2451545.0. There are
%precisely 86400 seconds in a Julian day.
TDBSec=86400*((TDB1-2451545.0)+TDB2);

%Get the ephemerides.
ScriptPath=mfilename('fullpath');
ScriptFolder = fileparts(ScriptPath);

%Load the ephemerides.
cspice_furnsh([ScriptFolder, '/data/de430.bsp'])

%Find the distance from the Earth to the Barycenter. 
state = cspice_spkezr('SOLAR SYSTEM BARYCENTER', ...
                                    TDBSec, ...
                                   'J2000', ...
                                   'NONE', ...
                                   'EARTH');
                  
%Convert the position components to meters
J2000Earth=state(1:3)*1e3;

%Align the Earth position vector with the GCRS/ICRS/BCRS axes.
rEarthGCRS=J2000F2ICRS(J2000Earth,Jul1,Jul2);

%Add the Earth position offset to the position vectors to get the BCRS
%position.
vec(1:3,:)=bsxfun(@plus,rEarthGCRS,x(1:3,:));

end

%LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.
