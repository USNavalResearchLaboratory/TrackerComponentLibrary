function [xEst,PEst]=spher2CartOrbitCubature(z,SR,systemType,useHalfRange,deltaT,zTx,zRx,M,isECEF,xi,w,longWay,GM,EarthRotRate)
%%SPHER2CARTORBITCUBATURE Solve Lambert's problem. That is, given two position
%            measurements of a ballistic (orbital) object and the time
%            difference of the measurements, determine the velocities of
%            the object at each of the two points assuming a Keplerian
%            gravitational model with just two bodies (the Earth (or Sun or
%            other planet), and the object (which is not massive). This
%            function assumes that measurements are given in noisy
%            spherical coordinates and it uses cubature integration to
%            solve the problem providing a covariance matrix. the
%            covariance matrix can help when starting tracks. The function
%            can work in a generic inertial coordinate system, or it can be
%            used with a generic Earth-centered Earth-fixed (ECEF)
%            coordinate system. That is, rather than being rigorously-
%            defined in terms of precession, nutation, Earth orientation
%            parameters, etc, the inputs are given in a generic ECEF
%            coordinate system that rotates at a rate of EarthRotRate
%            radians per second about the z-axis with respect to a generic
%            Earth centered inertial (ECI) coordinate system.
%
%INPUTS: z A 3X2 set of two points in 3D given in terms of range, azimuth
%          and elevation, with the angles in radians Note that many math
%          texts use a polar angle (pi/2-elevation) in place of elevation.
%          A polar angle is also known as a colatitude, an inclination
%          angle, a zenith angle, and a normal angle.
%       SR A 3X3X2 set of two lower-triangular square root covariance
%          matrices associated with the points. If a single 3X3 matrix is
%          passed, then it is assumed that SR is the same for both
%          measurements.
% systemType An parameter specifying the axes from which the angles are
%          measured. Possible values are
%          0 (The default if an empty matrix is passed) Azimuth is
%            measured counterclockwise from the x-axis in the x-y plane.
%            Elevation is measured up from the x-y plane (towards the
%            z-axis).
%          1 Azimuth is measured counterclockwise from the z-axis in the
%            z-x plane. Elevation is measured up from the z-x plane
%            (towards the y-axis).
%          2 This is the same as 0 except instead of being given elevation,
%            one desires the angle away from the z-axis, which is
%            (pi/2-elevation).
% useHalfRange An optional boolean value specifying whether the bistatic
%           (round-trip) range value has been divided by two. This normally
%           comes up when operating in monostatic mode (the most common
%           type of spherical coordinate system), so that the range
%           reported is a one-way range (or just half a bistatic range).
%           The default if an empty matrix is passed is false.
%    deltaT The positive time in seconds between when the two measurements
%           are taken.
%       zTx The 3X2 [x;y;z] location vectors of the transmitter in
%           Cartesian coordinates for each measurement. If only a single
%           3X1 vector is passed, then the transmitter locations are
%           assumed the same for both measurements.
%       zRx The 3X2 [x;y;z] location vectors of the receiver in Cartesian
%           coordinates for each measurement.  If only a single 3X1 vector
%           is passed, then the receiver locations are assumed the same for
%           both measurements.
%         M A 3X3X2 hypermatrix of the rotation matrices to go from the
%           alignment of the global coordinate system to that at the
%           receiver. The z-axis of the local coordinate system of the
%           receiver is the pointing direction of the receiver. If omitted,
%           then it is assumed that the local coordinate system is aligned
%           with the global and M=eye(3) --the identity matrix is used. If
%           only a single 3X3 matrix is passed, then it is assumed to be
%           the same for both measurements.
%    isECEF An optional boolean value indicating whether the r1Val and
%           everything are given in a generic ECEF coordinate system and
%           thus the outputs should be  given in ECEF coordinates. The
%           default if this parameter is omitted or an empty matrix is
%           passed is false.
%        xi A 6XnumCubaturePoints matrix of cubature points for the numeric
%           integration. If this and the final parameter are omitted or
%           empty matrices are passed, then fifthOrderCubPoints is used to
%           generate cubature points.
%         w A numCubaturePointsX1 vector of the weights associated with the
%           cubature points.
%   longWay Thisis an optional boolean value indicating whether the orbital
%           path goes the long way around during the time between
%           measurements For example, if the two points are 1 degree apart,
%           the target could have moved 1 degree in angle during the period
%           deltaT (longWay=false) or it could have moved 359 degrees
%           (longWay=true). The default if omitted or an empty matrix is
%           passed is false.
%        GM An optional value of the universal gravitational constant times
%           the mass of the Earth. If omitted, the value
%           Constants.EGM2008GM is used. The units are m^3/sec^2.
% EarthRotRate The rotation rate of the Earth in radians per second. This
%           parameter is only used if isECEF=true. if this parameter is
%           omitted or an empty matrix is passed, then
%           Constants.EGM2008EarthRotationRate is used.
%
%OUTPUTS: xEst The approximate 6X1 mean of the PDF of the target position
%              and velocity in Cartesian coordinates
%              ([x;y;z;xDot;yDot;zDot].
%         PEst The approximate 6X6 covariance matrix associated with xEst.
%
%In [1], a method of performing cubature measurement conversion with
%covariance matrices is given. using the function orbVelDet2Pt as the
%measurement conversion function, the same thing is applied here. See the
%comments to orbVelDet2Pt for more details on the solution to the Lambert
%problem.
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(systemType))
    systemType=0;
end

if(isempty(useHalfRange))
    useHalfRange=false;
end

if(nargin<8||isempty(M))
    M=eye(3,3);
end

if(nargin<9||isempty(isECEF))
    isECEF=false;
end

if(nargin<10||isempty(xi))
    [xi,w]=fifthOrderCubPoints(6);
end

if(nargin<12||isempty(longWay))
    longWay=false;
end

if(nargin<13||isempty(GM))
    GM=Constants.EGM2008GM;
end

if(nargin<14||isempty(EarthRotRate))
    EarthRotRate=Constants.EGM2008EarthRotationRate;
end

if(size(SR,3)==1)
    SR=repmat(SR,[1,1,2]);
end

if(size(M,3)==1)
    M=repmat(M,[1,1,2]);
end

if(size(zTx,2)==1)
    zTx=[zTx,zTx];
end

if(size(zRx,2)==1)
    zRx=[zRx,zRx];
end

zStacked=z(:);
SStacked=blkdiag(SR(:,:,1),SR(:,:,2));
numCubPoints=length(w);

%Transform the cubature points to the joint distribution of the
%measurements.
xiTrans=transformCubPoints(xi,zStacked,SStacked);

%Convert the transformed cubature points into global ECEF Cartesian
%coordinates.
r1Val=spher2Cart(xiTrans(1:3,:),systemType,useHalfRange,zTx(:,1),zRx(:,1),M(:,:,1));
r2Val=spher2Cart(xiTrans(4:6,:),systemType,useHalfRange,zTx(:,2),zRx(:,2),M(:,:,2));

%Solve Lambert's problem --the results will be in cell arrays, but
%there will only be one solution per cubature point (assuming nothing
%failed), so we can convert them directly into matrices.
[~,V2Comp]=orbVelDet2Pt(r1Val,r2Val,deltaT,0,isECEF,longWay,GM,EarthRotRate);
V2Comp=reshape(cell2mat(V2Comp),[3,numCubPoints]);

%Get the moments of the points to return
[xEst, PEst]=calcMixtureMoments([r2Val;V2Comp],w);
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
