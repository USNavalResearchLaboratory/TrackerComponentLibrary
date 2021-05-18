function z=state2PolRR(xTar,systemType,useHalfRange,zTx,zRx,M,stateType)
%%STATE2POLRR Convert a 2D target state consisting of at least position
%             and velocity components into monostatic, one-way polar range
%             and angle and (if velocity components are provided) range
%             rate. The velocity components can either be a 2D Cartesian
%             vector or a heading angle in radians and a speed. The angles
%             can be measured either counterclockwise from the x-axis,
%             which is standard in mathematics, or clockwise from the
%             y-axis, which is more common in navigation.
%
%INPUTS: xTar A xDimXN matrix of N target states consisting of 2D
%          position and, if a range rate measurement is desired, velocity
%          components in the order xTar=[xPosition;xVelocity], if
%          stateType=0, or in the format [xPosition;headingAngle;speed] is
%          stateType=1. If stateType=-1, only position components need to
%          be provided and range rate is not returned. Other components in
%          the state are ignored.
% systemType An optional parameter specifying the axes from which the
%          angles are measured. Possible values are
%          0 (The default if omitted) The azimuth angle is counterclockwise
%            from the x axis.
%          1 The azimuth angle is measured clockwise from the y axis.
% useHalfRange A boolean value specifying whether the bistatic range value
%          should be divided by two. This normally comes up when operating
%          in monostatic mode, so that the range reported is a one-way
%          range. The default if this parameter is not provided (or an
%          empty matrix is provided) is true.
%      zTx The 4XN [x;y;xDot;yDot] position and velocity vectors of the
%          transmitters in global
%          Cartesian coordinates. If this parameter is omitted or an
%          empty matrix is passed, then the transmitters are assumed to be
%          at the origin. If only a single vector is passed, then the
%          transmitter location is assumed the same for all of the target
%          states being converted.
%      zRx The 2XN [x;y;xDot;yDot] position and velocity vectors of the
%          receivers in Cartesian coordinates. If this parameter is omitted
%          or an empty matrix is passed, then the receivers are assumed to
%          be at the origin. If only a single vector is passed, then the
%          receiver location is assumed the same for all of the target
%          states being converted.
%        M A 2X2 rotation matrix to go from the alignment of the global
%          coordinate system to that at the receiver. If omitted
%          or an empty matrix is passed, then it is assumed that the local
%          coordinate system is aligned with the global and M=eye(2) --the
%          identity matrix is used.
%   stateType An optional parameter specifying the format of the target
%          state. Possible values are
%          -1 The target state consists of 2D Cartesian position and no
%             range rate estimate will be returned. Other target state
%             components are ignored.
%           0 (The default if omitted) The target state is given in terms
%             of 2D position, and 2D Cartesian velocity. Other components
%             are ignored.
%           1 The target state is given in terms of 2D position, a
%             heading angle (in radians) and a speed. The heading angle
%             is measured in the same system as the measurement, given by
%             the next parameter. Other state components are ignored.
%
%OUPUTS: z If stateType!=-1, z is a 3XN matrix of the target states in xTar
%          converted into 2D polar coordinates consisting of range,
%          azimuthal angle, and range rate. Otherwise, if stateType==-1,
%          then z is a 2XN matrix of the position components of the target
%          state converted into 2D polar coordinates.
%
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xTar,1);

if(nargin<7||isempty(stateType))
   stateType=0; 
end

if(nargin<6||isempty(M))
	M=eye(2);
end

if(nargin<5||isempty(zRx))
	zRx=zeros(xDim,1);
end

if(nargin<4||isempty(zTx))
    zTx=zeros(xDim,1);
end

if(nargin<3||isempty(useHalfRange))
    useHalfRange=true;
end

if(nargin<2||isempty(systemType))
    systemType=0; 
end

%The number of state vectors to convert.
N=size(xTar,2);

%Convert the positions.
zPos=Cart2Pol(xTar(1:2,:),systemType,useHalfRange,zTx(1:2,:),zRx(1:2,:),M);

if(stateType==-1)
    z=zPos;
    return;
end

%Allocate space for the converted states.
z=zeros(3,N);
z(1:2,:)=zPos;

%Compute the bistatic range rates.
switch(stateType)
    case 0%The target state contains 2D position and velocity.
        z(3,:)=getRangeRate(xTar,useHalfRange,zTx,zRx);
    case 1%The target state contains 2D position, heading angle and speed.
        x=zeros(4,N);
        
        x(1:2,:)=xTar(1:2,:);
        theta=xTar(3,:);%Heading
        v=xTar(4,:);%Speed
        
        %Make the heading angle in radians North of East, not East of
        %North, if necessary.
        switch(systemType)
            case 0
            case 1
                theta=pi/2-theta;
            otherwise
                error('Unknown angle type specified');
        end
        
        x(3:4,:)=M\[v.*cos(theta);
                    v.*sin(theta)];
        z(3,:)=getRangeRate(x,useHalfRange,zTx,zRx);
    otherwise
        error('Unknown state type specified');
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
