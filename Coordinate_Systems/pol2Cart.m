function cartPoints=pol2Cart(z,systemType,useHalfRange,zTx,zRx,M)
%%POL2CART Convert points from bistatic polar coordinates to 2D Cartesian
%          coordinates. The angle can measured either counterclockwise
%          from the x-axis, which is standard in mathematics, or clockwise
%          from the y axis, which is more common in navigation.
%
%INPUTS: z A 2XN matrix of points in polar coordinates. Each column of the
%          matrix has the format [r;azimuth], with azimuth given in
%          radians. Alternatively, if one just wishes to have unit vectors
%          returned, this can be a 1XN matrix of just azimuthal angles, in
%          which case the useHalfRange, zTx,and zRx inputs are ignored.
% systemType An optional parameter specifying the axis from which the
%          angles are measured. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) The
%            azimuth angle is counterclockwise from the x axis.
%          1 The azimuth angle is measured clockwise from the y axis.
% useHalfRange A boolean value specifying whether the bistatic range value
%           should be divided by two. This normally comes up when operating
%           in monostatic mode, so that the range reported is a one-way
%           range. The default if this parameter is not provided, or an
%           empty matrix is passed, is true.
%       zTx The 2XN [x;y] location vectors of the transmitters in global
%           Cartesian coordinates. If this parameter is omitted or an
%           empty matrix is passed, then the transmitters are assumed to
%           be at the origin. If only a single vector is passed, then the
%           transmitter location is assumed the same for all of the target
%           states being converted. zTx can have more than 2 rows;
%           additional rows are ignored.
%       zRx The 2XN [x;y] location vectors of the receivers in Cartesian
%           coordinates.  If this parameter is omitted or an empty matrix
%           is passed, then the receivers are assumed to be at the origin.
%           If only a single vector is passed, then the receiver location
%           is assumed the same for all of the target states being
%           converted. zRx can have more than 2 rows; additional rows are
%           ignored.
%         M A 2X2XN hypermatrix of the rotation matrices to go from the
%           alignment of the global coordinate system to that at the
%           receiver. If omitted or an empty matrix is passed, then it is
%           assumed that the local coordinate system is aligned with the
%           global and M=eye(2) --the identity matrix is used. If only a
%           single 2X2 matrix is passed, then it is assumed to be the same
%           for all of the N conversions.
%
%OUTPUTS: cartPoints A 2XN or matrix of the points transformed into
%                    Cartesian coordinates. Each columns of cartPoints is
%                    of the format [x;y].
%
%The conversion utilizing bistatic polar measurements in 2D is similar to
%that using bistatic r-u-v measurements in 3D, which is discussed in [1].
%In both instances, one turns the direction components into a unit vector
%and the multiplied it by a one-way monostatic range that has to be
%computed.
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(z))
    cartPoints=zeros(2,0);
    return
end

N=size(z,2);

if(nargin<6||isempty(M))
    M=repmat(eye(2),[1,1,N]);
elseif(size(M,3)==1)
    M=repmat(M,[1,1,N]);
end

if(nargin<5||isempty(zRx))
    zRx=zeros(2,N);
elseif(size(zRx,2)==1)
    zRx=repmat(zRx,[1,N]);
end

if(nargin<4||isempty(zTx))
    zTx=zeros(2,N);
elseif(size(zTx,2)==1)
    zTx=repmat(zTx,[1,N]);
end

if(nargin<3||isempty(useHalfRange))
    useHalfRange=true;
end

if(nargin<2||isempty(systemType))
    systemType=0;
end

%Extract the components.
if(size(z,1)==1)
    %No ranges, so just get unit vectors.
    cartPoints=zeros(2,N);
    switch(systemType)
        case 0
            cartPoints(1,:)=cos(z);
            cartPoints(2,:)=sin(z);
        case 1
            cartPoints(1,:)=sin(z);
            cartPoints(2,:)=cos(z);
        otherwise
            error('Invalid system type specified.')
    end
    
    %Rotate the unit vectors into the global coordinate system.
    for k=1:N
        %Convert to global Cartesian coordinates.
        cartPoints(:,k)=M(:,:,k)\cartPoints(:,k);
    end
    return;
else
    rB=z(1,:);
    azimuth=z(2,:);
end

%The bistatic range is used in the conversions below.
if(useHalfRange)
   rB=2*rB; 
end

%Get unit vectors from the azimuth.
u=zeros(2,N);
switch(systemType)
    case 0
        u(1,:)=cos(azimuth);
        u(2,:)=sin(azimuth);
    case 1
        u(1,:)=sin(azimuth);
        u(2,:)=cos(azimuth);
    otherwise
        error('Invalid system type specified.')
end

cartPoints=zeros(2,N);
for curPoint=1:N
    if(all(zRx(:,curPoint)==zTx(:,curPoint)))
        %If it is monostatic, set r1 directly. This avoids NaNs with
        %zero-ranges in the monostatic case.
        r1=rB(curPoint)/2;
        %This is the Cartesian location in the local coordinate system of
        %the receiver.
        zL=r1*u(:,curPoint);
        %Convert to global Cartesian coordinates.
        cartPoints(:,curPoint)=M(:,:,curPoint)\zL+zRx(1:2,curPoint);
    else
        %The transmitter location in the receiver's local coordinate
        %system.
        zTxL=M(:,:,curPoint)*(zTx(1:2,curPoint)-zRx(1:2,curPoint));

        r1=(rB(curPoint)^2-norm(zTxL)^2)/(2*(rB(curPoint)-dot(u(:,curPoint),zTxL)));
        %This is the Cartesian location in the local coordinate system of
        %the receiver.
        zL=r1*u(:,curPoint);
        %Convert to global Cartesian coordinates.
        cartPoints(:,curPoint)=M(:,:,curPoint)\zL+zRx(1:2,curPoint);
    end
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
