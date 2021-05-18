function theta=getPolAngle(zC,systemType,zRx,M)
%%GETPOLANGLE Convert positions into local polar angles in radians.
%
%INPUTS: z A 2XN matrix of Cartesian points whoe polar angles are desired.
%          Each column of cartPoints is of the format [x;y].
% systemType An optional parameter specifying the axes from which the
%          angles are measured. Possible values are
%          0 (The default if omitted) The azimuth angle is counterclockwise
%            from the x axis.
%          1 The azimuth angle is measured clockwise from the y axis.
%      zRx The 2XN [x;y] location vectors of the receivers in Cartesian
%          coordinates. If this parameter is omitted or an empty matrix is
%          passed, then the receivers are assumed to be at the origin. If
%          only a single vector is passed, then the receiver location is
%          assumed the same for all of the target states being converted.
%        M A 2X2XN hypermatrix of the rotation matrices to go from the
%          alignment of the global coordinate system to that at the
%          receiver. If omitted or an empty matrix is passed, then it is
%          assumed that the local coordinate system is aligned with the
%          global and M=eye(2) --the identity matrix is used. If only a
%          single 2X2 matrix is passed, then it is assumed to be the same
%          for all of the N conversions.
%
%OUTPUT: theta The 1XN matrix of azimuth angles in radians
%
%The conversion a shifted rotated coordinate system in 2D is similar to
%that using bistatic r-u-v measurements in 3D, which is discussed in [1].
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(zC,2);

if(nargin<4||isempty(M))
    M=repmat(eye(2),[1,1,N]);
elseif(size(M,3)==1)
    M=repmat(M,[1,1,N]);
end

if(nargin<3||isempty(zRx))
    zRx=zeros(2,N);
elseif(size(zRx,2)==1)
    zRx=repmat(zRx,[1,N]);
end

if(nargin<2||isempty(systemType))
    systemType=0;
end

theta=zeros(1,N);
for curPoint=1:N
    %The target location in the receiver's coordinate system.
    z1=M(:,:,curPoint)*(zC(:,curPoint)-zRx(1:2,curPoint));

    %Extract the coordinates
    x=z1(1);
    y=z1(2);

    switch(systemType)
        case 0
            azimuth=atan2(y,x);
        case 1
            azimuth=atan2(x,y);
        otherwise
            error('Invalid system type specified.')
    end

    theta(:,curPoint)=azimuth;
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
