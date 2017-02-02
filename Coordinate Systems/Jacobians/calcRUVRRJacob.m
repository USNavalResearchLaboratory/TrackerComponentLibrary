function J=calcRUVRRJacob(components,xState,useHalfRange,sTx,sRx,M)
%%CALCRUVRRJACOB    Calculate the Jacobian for the selected components of
%                   a monostatic or bistatic RUV measurement, possibly
%                   with non-relativistic range rate, ignoring atmospheric
%                   effects. A target state consists of position and (if
%                   range rate measurements are desired) velocity
%                   components.
%
%INPUTS: components An number specifying the components desired
%                   in the Jacobian matrix.
%                   0 is range, u-v angle and range rate.
%                   1 is range-only.
%                   2 is u-v direction only.
%                   3 is range rate only.
%                   4 is range and u-v direction.
%                   5 is range and range rate.
%                   6 is u-v direction and range rate. 
%       xState     The target state vector in the global coordinate system
%                   with [x;y;z] components if no range rate components are
%                   desired or with [x;y;z;xDot;yDot;zDot] components if
%                   range rate components are desired, where xDot yDot and
%                   zDot are components of a velocity vector.
%   useHalfRange    A boolean value specifying whether the bistatic range
%                   value has been divided by two. This normally comes up
%                   when operating in monostatic mode, so that the range
%                   reported is a one-way range. The default if this
%                   parameter is not provided is false.
%        sTx        The transmitter state vector in the global coordinate
%                   system with [x;y;z] components if no range rate
%                   components are desired or a with [x;y;z;xDot;yDot;zDot]
%                   components if range rate components are desired. If
%                   omitted, then a vector of zeros is used.
%        sRx        The receiver state vector in the global coordinate
%                   system with [x;y;z] components if no range rate
%                   components are desired or a with [x;y;z;xDot;yDot;zDot]
%                   components if range rate components are desired. If
%                   omitted, then a vector of zeros is used.
%        M          A rotation matrix from the global Coordinate system to
%                   the orientation of the coordinate system at the
%                   receiver. This is only necessary if UV direction
%                   components are desired. If omitted, it is assumed to be
%                   the identity matrix.
%
%OUTPUTS: J     The Jacobian matrix with derivatives with respect to global
%               position and velocity components. Each row is a
%               component of range, u, v and range rate with derivatives
%               taken with respect to [x,y,z,xdot,ydot,zdot]. The rows
%               are ordered as [bistatic range;u;v;bistatic range rate],
%               where a subset of those rows is provided as specified by
%               the components parameter. The bistatic range and range rate
%               components are divided by 2 if useHalfRange=true.
%
%A derivation of the components of the Jacobian is given in [1]. 
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%October 2013 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<6)
       M=eye(3); 
    end
    
    if(nargin<5)
       sRx=zeros(6,1); 
    end
    
    if(nargin<4)
       sTx=zeros(6,1); 
    end

    if(nargin<3)
        useHalfRange=false;
    end
   
    hasRange=false;
    hasDirection=false;
    hasRangeRate=false;
    
    switch(components)
        case 1
            hasRange=true;
        case 2
            hasDirection=true;
        case 3
            hasRangeRate=true;
        case 4
            hasRange=true;
            hasDirection=true;
        case 5
            hasRange=true;
            hasRangeRate=true;
        case 6
            hasDirection=true;
            hasRangeRate=true;
        otherwise
            hasRange=true;
            hasDirection=true;
            hasRangeRate=true;
    end

    %These are still in global Cartesian coordinates.
    dtr=xState(1:3,1)-sRx(1:3,1);
    dtl=xState(1:3,1)-sTx(1:3,1);

    tL=M*dtr;%The target location in local Coordinates.
    tMag=norm(tL);
    x=tL(1);
    y=tL(2);
    z=tL(3);
    
    if(hasDirection==true)
        %The gradient of u in local coordinates.
        du=zeros(6,1);
        du(1)=(y^2+z^2)/tMag^3;
        du(2)=-x*y/tMag^3;
        du(3)=-x*z/tMag^3;

        %The gradient of v in local coordinates.
        dv=zeros(6,1);
        dv(1)=-x*y/tMag^3;
        dv(2)=(x^2+z^2)/tMag^3;
        dv(3)=-y*z/tMag^3;

        %Now, the gradient vectors for the angular components must be rotated back
        %into the global coordinate system.
        du(1:3)=M\du(1:3);
        du(4:6)=M\du(4:6);
        dv(1:3)=M\dv(1:3);
        dv(4:6)=M\dv(4:6);
    end
    
    %The gradient of rB in global Coordinates
    if(hasRange==true||hasRangeRate==true)
        drB=zeros(6,1);
        drB(1:3)=dtr/norm(dtr)+dtl/norm(dtl);
    end

    %The gradient of the range-rate in local coordinates.
    if(hasRangeRate)
        drd=zeros(6,1);
        dvr=xState(4:6)-sRx(4:6);
        dvl=xState(4:6)-sTx(4:6);

        %The position derivaties
        drd(1:3)=A(dtr)*dvr/norm(dtr)^3+A(dtl)*dvl/norm(dtl)^3;

        %The velocity derivatives.
        drd(4:6)=drB(1:3);
    end

    if(hasRange&& hasRangeRate &&hasDirection)
        if(useHalfRange)
            J=[drB/2,du,dv,drd/2]';
        else
            J=[drB,du,dv,drd]';
        end
    elseif(hasRange&&hasRangeRate&&~hasDirection)
        if(useHalfRange)
            J=[drB/2,drd/2]';
        else
            J=[drB,drd]';
        end
    elseif(hasRange&&~hasRangeRate&&hasDirection)
        if(useHalfRange)
            J=[drB/2,du,dv]';
        else
            J=[drB,du,dv]';
        end
    elseif(hasRange&&~hasRangeRate&&~hasDirection)
        if(useHalfRange)
            J=drB'/2;
        else
            J=drB';
        end
    elseif(~hasRange&&hasRangeRate && hasDirection)
        J=[du,dv,drd]';
    elseif(~hasRange&&hasRangeRate&&~hasDirection)
        J=drd';
    else %(~hasRange&&~hasRangeRate&&hasDirection)
        J=[du,dv]';
    end
end

function val=A(vec)
%This is a simple helper function as defined in "Basic Tracking
%Using Nonlinear 3D Monostatic and Bistatic Measurements" by David F.
%Crouse.

    x=vec(1);
    y=vec(2);
    z=vec(3);
    val=[y^2+z^2,-x*y,   -x*z;
        -x*y,   x^2+z^2,-y*z;
        -x*z,   -y*z,   x^2+y^2];
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

