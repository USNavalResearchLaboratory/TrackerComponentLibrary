function J=calcPolarRRJacob(components,xState,systemType)
%%CALCPOLARRRJACOB  Calculate the Jacobian for the selected components of
%                   a 2D polar measurement, possibly with non-relativistic
%                   range rate, ignoring atmospheric effects. A target
%                   state consisits of position and (if range rate
%                   measurements are desired) velocity components.
%
%INPUTS: components An number specifying the components desired. Angles are
%                   in the Jacobian matrix. Angles are in radians. Possible
%                   values are:
%                   0 is range, azimuth angle and range rate.
%                   1 is range-only.
%                   2 is azimuth direction only.
%                   3 is range rate only.
%                   4 is range and azimuth direction.
%                   5 is range and range rate.
%                   6 is azimuth direction and range rate.
%       xState      The target state vector in local 2D coordinate system
%                   with [x;y] components if no range rate components are
%                   desired or with [x;y;z;xDot;yDot;zDot] components if
%                   range rate components are desired, where xDot yDot and
%                   zDot are components of a velocity vector.
%      systemType   An optional parameter specifying the axes from which
%                   the angles are measured. Possible vaues are
%                   0 (The default if omitted) The azimuth angle is
%                     counterclockwise from the x axis.
%                   1 The azimuth angle is measured clockwise from the y
%                     axis.
%
%OUTPUTS: J     The Jacobian matrix with derivatives with respect to
%               position and velocity components. Each row is a component
%               of range, azimuth, and range rate with derivatives taken
%               with respect to [x,y,xdot,ydot]. The rows are ordered as
%               [range;azimuth;range rate], where a subset of those rows is
%               provided as specified by the components parameter.
%
%The non-relativistic range-rate approximation is derived in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%July 2014 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<3)
        systemType=0;
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

    x=xState(1);
    y=xState(2);
    r2=x^2+y^2;%Squared range
    r=sqrt(r2);%Range
    
    vx=xState(3);
    vy=xState(4);

    if(hasDirection==true)
        %The gradients of azimuth in local coordinates.
        dAz=zeros(4,1);
        switch(systemType)
            case 0
                %Derivatives with respect to x.
                dAz(1)=-y/r2;

                %Derivatives with respect to y.
                dAz(2)=x/r2;
            case 1
                %Derivatives with respect to x.
                dAz(1)=y/r2;

                %Derivatives with respect to y.
                dAz(2)=-x/r2;
            otherwise
                error('Invalid system type specified.')
        end
    end
    
    if(hasRange==true||hasRangeRate==true)
        dr=zeros(4,1);
        
        %The derivative with respect to x.
        dr(1)=x/r;
        %The derivative with respect to y.
        dr(2)=y/r;
    end
    
    %The gradient of the range-rate.
    if(hasRangeRate)
       drd=[y*(-x*vy+y*vx)/r^3;
            x*(x*vy-y*vx)/r^3;
            x/r;
            y/r];
    end
    
    if(hasRange&& hasRangeRate &&hasDirection)
        J=[dr,dAz,drd]';
    elseif(hasRange&&hasRangeRate&&~hasDirection)
        J=[dr,drd]';
    elseif(hasRange&&~hasRangeRate&&hasDirection)
        J=[dr,dAz]';
    elseif(hasRange&&~hasRangeRate&&~hasDirection)
        J=dr';
    elseif(~hasRange&&hasRangeRate && hasDirection)
        J=[dAz,drd]';
    elseif(~hasRange&&hasRangeRate&&~hasDirection)
        J=drd';
    else %(~hasRange&&~hasRangeRate&&hasDirection)
        J=[dAz]';
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
