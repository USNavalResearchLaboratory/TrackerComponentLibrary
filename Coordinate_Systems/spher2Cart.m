function cartPoints=spher2Cart(points,systemType,useHalfRange,zTx,zRx,M,flipNegRange)
%%SPHER2CART Convert points from spherical coordinates to Cartesian
%            coordinates. If no range is specified, then the angles are
%            converted into Cartesian unit vectors. An option allows for
%            the angles to be specified from different axes. Optionally, a
%            bistatic range can be used when considering bistatic
%            measurements in a local spherical coordinate system. 
%
%INPUTS: points One or more points given in terms of range, azimuth and
%           elevation, with the angles in radians, or in terms of just
%           azimuth and elevation if Cartesian unit vectors are desired. To
%           convert N points, points is a 3XN matrix with each column
%           having the format [range;azimuth;elevation] or it is a 2XN
%           matrix with each column having format [azimuth; elevation] if
%           unit vectors are desired. Note that many math texts use a polar
%           angle (pi/2-elevation) in place of elevation. A polar angle is
%           also known as a colatitude, an inclination angle, a zenith
%           angle, and a normal angle. systemType=2 supports a polar angle.
% systemType An optional parameter specifying the axis from which the
%           angles are measured in radians. Possible values are:
%           0 (The default if omitted) Azimuth is measured 
%             counterclockwise from the x-axis in the x-y plane. Elevation
%             is measured up from the x-y plane (towards the z-axis). This
%             is consistent with common spherical coordinate systems for
%             specifying longitude (azimuth) and geocentric latitude
%             (elevation).
%           1 Azimuth is measured counterclockwise from the z-axis in the
%             z-x plane. Elevation is measured up from the z-x plane
%             (towards the y-axis). This is consistent with some spherical
%             coordinate systems that use the z-axis as the boresight
%             direction of the radar.
%           2 This is the same as 0 except instead of being given
%             elevation, one is given the angle away from the z-axis, which
%             is (pi/2-elevation).
%           3 This is the same as 0 except azimuth is measured clockwise
%             from the y-axis in the x-y plane instead of counterclockwise
%             from the x-axis. This coordinate system often arises when
%             given "bearings" in a local East-North-Up coordinate system,
%             where the bearing directions are measured East of North.
% useHalfRange An optional boolean value specifying whether the bistatic
%           (round-trip) range value has been divided by two. This normally
%           comes up when operating in monostatic mode (the most common
%           type of spherical coordinate system), so that the range
%           reported is a one-way range (or just half a bistatic range).
%           The default if this parameter is not provided is false if zTx
%           is provided and true if it is omitted (monostatic). If no
%           range values are provided, an empty matrix can be passed.
%       zTx The 3XN [x;y;z] location vectors of the transmitters in global
%           Cartesian coordinates. If this parameter is omitted, then the
%           transmitters are assumed to be at the origin. If only a single
%           vector is passed, then the transmitter location is assumed the
%           same for all of the target states being converted. zTx can have
%           more than 3 rows; additional rows are ignored. If monostatic or
%           no range values are provided, an empty matrix can be passed.
%       zRx The 3XN [x;y;z] location vectors of the receivers in Cartesian
%           coordinates. If this parameter is omitted, then the receivers
%           are assumed to be at the origin. If only a single vector is
%           passed, then the receiver location is assumed the same for all
%           of the target states being converted. zRx can have more than 3
%           rows; additional rows are ignored. If no range values are
%           provided, an empty matrix can be passed.
%         M A 3X3XN hypermatrix of the rotation matrices to go from the
%           alignment of the global coordinate system to that at the
%           receiver. If omitted, then it is assumed that the local
%           coordinate system is aligned with the global and M=eye(3) --the
%           identity matrix is used. If only a single 3X3 matrix is passed,
%           then it is assumed to be the same for all of the N conversions.
% flipNegRange An optional boolean parameter. If true, the absolute value
%           of all range components, if present, is taken before
%           conversion.
%
%OUTPUTS: cartPoints For N points, cartPoints is a 3XN matrix of the
%                    converted points with each column having the format
%                    [x;y;z]. If no range values were passed, then all of
%                    the vectors are unit direction vectors.
%
%The conversion from spherical to Cartesian coordinates is given in [1].
%However, when considering the bistatic scenario, concepts from the
%bistatic r-u-v to Cartesian conversion described in [2] are used.
%
%REFERENCES:
%[1] R. L. Duncombe, "Computational techniques," in Explanatory Supplement
%    to the Astronomical Almanac, 3rd ed., S. E. Urban and P. K.
%    Seidelmann, Eds. Mill Valley, CA: University Science Books, 2013,
%    ch. 14.4.4.1.
%[2] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(zTx))
    zTx=[];
end

if(nargin<5||isempty(zRx))
    zRx=[];
end

if(nargin<6||isempty(M))
    M=[]; 
end

if(nargin<7||isempty(flipNegRange))
    flipNegRange=false;
end

if((nargin<4||isempty(zTx))&&(nargin<3||isempty(useHalfRange)))
    useHalfRange=true;
elseif(nargin<3||isempty(useHalfRange))
    useHalfRange=false;
end

if(nargin<2||isempty(systemType))
    systemType=0;
end

%Extract the coordinates
hasRange=size(points,1)>2;

if(hasRange)
    if(flipNegRange)
        r=abs(points(1,:));
    else
        r=points(1,:);
    end
    azimuth=points(2,:);
    elevation=points(3,:);
else
    azimuth=points(1,:);
    elevation=points(2,:);
end

if(systemType==2)
    elevation=pi/2-elevation;
    systemType=0;
elseif(systemType==3)
    azimuth=pi/2-azimuth;
    systemType=0;
end

%Get unit vectors pointing in the direction of the targets from the
%receiver in the LOCAL coordinate system of the receiver.
N=size(points,2);
uVecL=zeros(3,N);
switch(systemType)
    case 0
        uVecL(1,:)=cos(azimuth).*cos(elevation);
        uVecL(2,:)=sin(azimuth).*cos(elevation);
        uVecL(3,:)=sin(elevation);
    case 1
        uVecL(1,:)=sin(azimuth).*cos(elevation);
        uVecL(2,:)=sin(elevation);
        uVecL(3,:)=cos(azimuth).*cos(elevation);
    otherwise
        error('Invalid system type specified.')
end

if(nargin<6||isempty(M))
    M=repmat(eye(3),[1,1,N]);
elseif(size(M,3)==1)
    M=repmat(M,[1,1,N]);
end

%If just a standard spherical coordinate conversion is being performed.
if((hasRange==true)&&isempty(zTx)&&useHalfRange==true)
    cartPoints=bsxfun(@times,r,uVecL);
    
    %Deal with any rotation
    for curPoint=1:N
        cartPoints(:,curPoint)=M(:,:,curPoint)\cartPoints(:,curPoint);
    end
    return;
elseif((hasRange==true)&&isempty(zTx)&&useHalfRange==false)
    cartPoints=bsxfun(@times,r/2,uVecL);
    
    %Deal with any rotation
    for curPoint=1:N
        cartPoints(:,curPoint)=M(:,:,curPoint)\cartPoints(:,curPoint);
    end
    return;
end

%Allocate space for the return variables.
cartPoints=zeros(3,N);

%If just unit direction vectors are being returned.
if(hasRange==false)
    for curPoint=1:N
        %Convert to global coordinates, if a rotation is needed.
        cartPoints(:,curPoint)=M(:,:,curPoint)'*uVecL(:,curPoint);
    end
    return;
end

%If any type of bistatic range thing is done.
if(nargin<5||isempty(zRx))
    zRx=zeros(3,N);
elseif(size(zRx,2)==1)
    zRx=repmat(zRx,[1,N]);
end

if(nargin<4||isempty(zTx))
    zTx=zeros(3,N);
elseif(size(zTx,2)==1)
    zTx=repmat(zTx,[1,N]);
end

%The bistatic range is used in the conversions below. If the range provided
%to the function was divided by two, multiple it by two to geth the
%round-trip range.
if(useHalfRange)
   r=2*r; 
end

for curPoint=1:N
    %Convert the vector into the global coordinate system.
    uVec=uVecL(:,curPoint);
    
    %The transmitter location in the receiver's local coordinate system.
    zTxL=M(:,:,curPoint)*(zTx(1:3,curPoint)-zRx(1:3,curPoint));

    %The range from the transmitter to the target.
    r1=(r(curPoint)^2-norm(zTxL)^2)/(2*(r(curPoint)-dot(uVec,zTxL)));

    %This is the Cartesian location in the local coordinate system of
    %the receiver.
    zL=r1*uVec;
    %Convert to global Cartesian coordinates. The transpose of a rotation
    %matrix is its inverse.
    cartPoints(:,curPoint)=M(:,:,curPoint)'*zL+zRx(1:3,curPoint);
end
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
