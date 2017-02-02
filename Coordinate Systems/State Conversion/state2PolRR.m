function z=state2PolRR(xTar,stateType,angleType)
%%STATE2POLRR  Convert a 2D target state consisting of at least position
%              and velocity components into monostatic, one-way polar range
%              and angle and (if velocity components are provided) range
%              rate. The velocity components can either be a 2D Cartesian
%              vector or a heading angle in radians and a speed. The angles
%              can be measured either counterclockwise from the x-axis,
%              which is standard in mathematics, or clockwise from the
%              y-axis, which is more common in navigation.
%
%INPUTS: xTar A xDimXN matrix of N target states consisting of 2D
%             position and, if a range rate measurement is desired,
%             velocity components in the order xTar=[xPosition;xVelocity],
%             if stateType=0, or in the format
%             [xPosition;headingAngle;speed] is stateType=1. If
%             StateType=-1, only position components need to be provided
%             and range rate is not returned. Other components in the state
%             are ignored.
%   stateType An optional parameter specifying the format of the target
%             state. Possible values are
%            -1 The target state consists of 2D Cartesian position and no
%               range rate estimate will be returned. Other target state
%               components are ignored.
%             0 (The default if omitted) The target state is given in terms
%               of 2D position, and 2D Cartesian velocity. Other components
%               are ignored.
%             1 The target state is given in terms of 2D position, a
%               heading angle (in radians) and a speed. The heading angle
%               is measured in the same system as the measurement, given by
%               the next parameter. Other state components are ignored.
%  angleType  An optional parameter specifying the axes from which the
%             angles are measured. Possible vaues are
%             0 (The default if omitted) The azimuth angle is 
%               counterclockwise from the x-axis.
%             1 The azimuth angle is measured clockwise from the y-axis.
%
%OUPUTS: z   If stateType!=-1, z is a 3XN matrix of the target states in
%            xTar converted into 2D polar coordinates consisting of
%            monostatic (one-way) range, azimuthal angle, and range rate.
%            If stateType==-1, then z is a 2XN matrix of the position
%            components of the target state converted into 2D polar
%            coordinates.
%
%The range and range rate are both one-way.
%
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Set default parameters if the second and third inputs are omitted.
if(nargin<3)
    angleType=0;
end

if(nargin<2)
   stateType=0; 
end

%The number of state vectors to convert.
N=size(xTar,2);

%Convert the positions.
zPos=Cart2Pol(xTar(1:2,:),angleType);

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
        z(3,:)=getRangeRate(xTar,true,zeros(4,N),zeros(4,N),2);
    case 1%The target state contains 2D position, heading angle and speed.
        x=zeros(4,N);
        x(1:2,:)=xTar(1:2,:);
        theta=xTar(3,:);%Heading
        v=xTar(4,:);%Speed
        
        %Make the heading angle in radians North of East, not East of
        %North, if necessary.
        switch(angleType)
            case 0
            case 1
                theta=pi/2-theta;
            otherwise
                error('Unknown angle type specified');
        end
        
        x(3:4,:)=[v.*cos(theta);
                  v.*sin(theta)];
        z(3,:)=getRangeRate(x,true,zeros(4,N),zeros(4,N),2);
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
