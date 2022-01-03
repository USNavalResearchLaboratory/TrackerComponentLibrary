function spherStates=stateCart2Sphere(x,systemType)
%%STATECART2SPHERE Convert a target state in 3D space consiting of
%               position, velocity and possibly acceleration terms into
%               monostatic spherical coordinates. With acceleration, this
%               would be going from
%               [x;y;z;xDot;yDot;zDot;xDDot;yDDot;zDDot] to
%               [r;theta;phi;rDot;thetaDot;phiDot;rDDot;thetaDDot;phiDDot]
%               where two Ds indicates a second derivative with respect to
%               time and theta and phi are azimuth and elevation.
%
%INPUTS: x The 6XN or 9XN set of Cartesian state vectors consisting of
%          position, velocity and possibly acceleration. Note that a
%          singularity exists at the origin, so non-finite terms can arise
%          in the output for entry i where x(1:3,i)=[0;0;0].
% systemType An optional parameter specifying the axes from which the
%          angles are measured in radians. Possible values are
%          0 (The default if omitted) Azimuth is measured counterclockwise
%            from the x-axis in the x-y plane. Elevation is measured up
%            from the x-y plane (towards the z-axis). This is consistent
%            with common spherical coordinate systems for specifying
%            longitude (azimuth) and geocentric latitude (elevation).
%          1 Azimuth is measured counterclockwise from the z-axis in the
%            z-x plane. Elevation is measured up from the z-x plane
%            (towards the y-axis). This is consistent with some spherical
%            coordinate systems that use the z axis as the boresight
%            direction of the radar.
%          2 This is the same as 0 except instead of being given
%            elevation, one desires the angle away from the z-axis, which
%            is (pi/2-elevation).
%
%OUTPUTS: spherStates The 6XN or 9XN set of target states given in spherical
%                   coordinates. The range is a one-way range and the
%                   angles are given in radians.
%
%The derivation is given in [1]. The function aCVSphere is an
%implementation of a related linear dynamic model.
%
%EXAMPLE:
%Here we note that the results are consistent with the inverse function:
%stateSpher2Cart.
% systemType=0;
% x=[100;-60;200;-3;12;160;108;-116;2];
% xRet=stateSpher2Cart(stateCart2Sphere(x,systemType),systemType);
% max(abs(x(:)-xRet(:)))
%One will see that the error is less than 1e-13, indicating good agreement.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic linear Cartesian dynamic models in local
%    coordinates," Naval Research Laboratory, Washington, DC, Tech. Rep.
%    NRL/MR/5344-19-9882, 24 Aug. 2019.
%
%August 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(systemType))
    systemType=0; 
end

numDim=size(x,1);
N=size(x,2);

posSphere=Cart2Sphere(x(1:3,:),systemType,true);
spherStates=zeros(numDim,N);

spherStates(1:3,:)=posSphere;

r=posSphere(1,:);
ur=bsxfun(@rdivide,x(1:3,:),r);

theta=posSphere(2,:);
phi=posSphere(3,:);

sinTheta=sin(theta);
cosTheta=cos(theta);
sinPhi=sin(phi);
cosPhi=cos(phi);

if(systemType==0||systemType==1)
    if(systemType==0)
        uTheta=[-sinTheta;cosTheta;zeros(1,N)];
        uPhi=[-cosTheta.*sinPhi;-sinTheta.*sinPhi;cosPhi];
    else
        uTheta=[cosTheta;zeros(1,N);-sinTheta];
        uPhi=[-sinTheta.*sinPhi;cosPhi;-cosTheta.*sinPhi];
    end
    
    rVel=x(4:6,:);
    
    rDot=sum(rVel.*ur,1);
    thetaDot=bsxfun(@rdivide,sum(rVel.*uTheta,1),r.*cosPhi);
    phiDot=bsxfun(@rdivide,sum(rVel.*uPhi,1),r);
    
    spherStates(4:6,:)=[rDot;thetaDot;phiDot];
    
    if(numDim>6)
        rAccel=x(7:9,:);
        
        tanPhi=sinPhi./cosPhi;
        
        ar=sum(rAccel.*ur,1);
        aTheta=sum(rAccel.*uTheta,1);
        aPhi=sum(rAccel.*uPhi,1);
        
        rDDot=ar+r.*phiDot.^2+r.*thetaDot.^2.*cosPhi.^2;
        thetaDDot=(1./r).*(-2*rDot.*thetaDot+aTheta./cosPhi+2*r.*thetaDot.*phiDot.*tanPhi);
        phiDDot=(1./r).*(aPhi-2*rDot.*phiDot-r.*thetaDot.^2.*cosPhi.*sinPhi);

        spherStates(7:9,:)=[rDDot;thetaDDot;phiDDot];
    end
else%systemType==2
    uTheta=[-sinTheta;cosTheta;zeros(1,N)];
    uPhi=[cosTheta.*cosPhi;sinTheta.*cosPhi;-sinPhi];
    
    rVel=x(4:6,:);
    
    rDot=sum(rVel.*ur,1);
    thetaDot=bsxfun(@rdivide,sum(rVel.*uTheta,1),r.*sinPhi);
    phiDot=bsxfun(@rdivide,sum(rVel.*uPhi,1),r);
    
    spherStates(4:6,:)=[rDot;thetaDot;phiDot];
    
    if(numDim>6)
        rAccel=x(7:9,:);
        
        cotPhi=cosPhi./sinPhi;
        
        ar=sum(rAccel.*ur,1);
        aTheta=sum(rAccel.*uTheta,1);
        aPhi=sum(rAccel.*uPhi,1);
        
        rDDot=ar+r.*phiDot.^2+r.*thetaDot.^2.*sinPhi.^2;
        thetaDDot=(1./r).*(-2*rDot.*thetaDot+aTheta./sinPhi-2*r.*thetaDot.*phiDot.*cotPhi);
        phiDDot=(1./r).*(aPhi-2*rDot.*phiDot+r.*thetaDot.^2.*cosPhi.*sinPhi);
        
        spherStates(7:9,:)=[rDDot;thetaDDot;phiDDot];
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
