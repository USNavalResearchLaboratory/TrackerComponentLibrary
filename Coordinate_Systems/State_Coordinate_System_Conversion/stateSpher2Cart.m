function cartStates=stateSpher2Cart(x,systemType)
%%STATESPHER2CART Convert a target state in 3D space in monostatic
%               spherical coordinates with first and possibly second time
%               derivatives of the position components into Cartesian
%               coordinates. The state has the format
%               [r;theta;phi;rDot;thetaDot;phiDot;rDDot;thetaDDot;phiDDot],
%               where two Ds indicate a second derivative with respect to
%               time and theta and phi are azimuth and elevation. In
%               Cartesian coordinates, the converted state has the
%               form [x;yz;;xDot;yDot;zDot;xDDot;yDDot;zDDot].
%
%INPUTS: x The 6XN or 9XN set of spherical state vectors consisting of
%          position, velocity and possibly acceleration. The angles are
%          given in radians. The range is a one-way range.
% systemType An optional parameter specifying the axes from which the
%           angles are measured in radians. Possible values are
%           0 (The default if omitted) Azimuth is measured
%             counterclockwise from the x-axis in the x-y plane. Elevation
%             is measured up from the x-y plane (towards the z-axis). This
%             is consistent with common spherical coordinate systems for
%             specifying longitude (azimuth) and geocentric latitude
%             (elevation).
%           1 Azimuth is measured counterclockwise from the z-axis in the
%             z-x plane. Elevation is measured up from the z-x plane
%             (towards the y-axis). This is consistent with some spherical
%             coordinate systems that use the z axis as the boresight
%             direction of the radar.
%           2 This is the same as 0 except instead of being given
%             elevation, one desires the angle away from the z-axis, which
%             is (pi/2-elevation).
%
%OUTPUTS: CartStates The 6XN or 9XN set of target states given in Cartesian
%                    coordinates.
%
%The derivation is given in [1]. The function aCVSphere is an
%implementation of a related linear dynamic model.
%
%EXAMPLE:
%Here we note that the results are consistent with the inverse function:
%stateCart2Sphere.
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
cartStates=zeros(numDim,N);

r=x(1,:);
theta=x(2,:);
phi=x(3,:);
rDot=x(4,:);
thetaDot=x(5,:);
phiDot=x(6,:);

sinTheta=sin(theta);
cosTheta=cos(theta);
sinPhi=sin(phi);
cosPhi=cos(phi);

if(systemType==0||systemType==1)
    if(systemType==0)
        ur=[cosTheta.*cosPhi;sinTheta.*cosPhi;sinPhi];
        uTheta=[-sinTheta;cosTheta;zeros(1,N)];
        uPhi=[-cosTheta.*sinPhi;-sinTheta.*sinPhi;cosPhi];
    else
        ur=[sinTheta.*cosPhi;sinPhi;cosTheta.*cosPhi];
        uTheta=[cosTheta;zeros(1,N);-sinTheta];
        uPhi=[-sinTheta.*sinPhi;cosPhi;-cosTheta.*sinPhi];
    end
    
    %Position
    cartStates(1:3,:)=bsxfun(@times,r,ur);
    
    %Velocity
    cartStates(4:6,:)=bsxfun(@times,rDot,ur)+bsxfun(@times,r.*thetaDot.*cosPhi,uTheta)+bsxfun(@times,r.*phiDot,uPhi);
    
    if(numDim>6)
        rDDot=x(7,:);
        thetaDDot=x(8,:);
        phiDDot=x(9,:);
        
        %Acceleration
        ar=rDDot-r.*thetaDot.^2.*cosPhi.^2-r.*phiDot.^2;
        aTheta=(2*rDot.*thetaDot+r.*thetaDDot).*cosPhi-2*r.*thetaDot.*phiDot.*sinPhi;
        aPhi=2*rDot.*phiDot+r.*phiDDot+r.*thetaDot.^2.*cosPhi.*sinPhi;
        
        cartStates(7:9,:)=bsxfun(@times,ar,ur)+bsxfun(@times,aTheta,uTheta)+bsxfun(@times,aPhi,uPhi);
    end
elseif(systemType==2)
    ur=[cosTheta.*sinPhi;sinTheta.*sinPhi;cosPhi];
    uTheta=[-sinTheta;cosTheta;zeros(1,N)];
    uPhi=[cosTheta.*cosPhi;sinTheta.*cosPhi;-sinPhi];
    
    %Position
    cartStates(1:3,:)=bsxfun(@times,r,ur);
    
    %Velocity
    cartStates(4:6,:)=bsxfun(@times,rDot,ur)+bsxfun(@times,r.*thetaDot.*sinPhi,uTheta)+bsxfun(@times,r.*phiDot,uPhi);
    
    if(numDim>6)
        rDDot=x(7,:);
        thetaDDot=x(8,:);
        phiDDot=x(9,:);
        
        %Acceleration
        ar=rDDot-r.*thetaDot.^2.*sinPhi.^2-r.*phiDot.^2;
        aTheta=(2*rDot.*thetaDot+r.*thetaDDot).*sinPhi+2*r.*thetaDot.*phiDot.*cosPhi;
        aPhi=2*rDot.*phiDot+r.*phiDDot-r.*thetaDot.^2.*cosPhi.*sinPhi;

        cartStates(7:9,:)=bsxfun(@times,ar,ur)+bsxfun(@times,aTheta,uTheta)+bsxfun(@times,aPhi,uPhi);
    end
else
   error('Invalid systemType specified.') 
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
