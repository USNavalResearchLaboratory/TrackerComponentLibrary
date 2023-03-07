function JTotal=calcSpherInvJacob(z,systemType)
%%SPHERINVJACOB Calculate the Jacobian for a 3D Cartesian position with
%          respect to monostatic spherical range, azimuth, and elevation
%          components. This produces derivatives of (x,y,z) with respect to
%          (range,Az,El). The function calcSpherJacob produces derivatives
%          of (range,Az,El) with respect to (x,y,z) in the more general
%          bistatic case.
%
%INPUTS: z The 3XN position vectors in the global spherical coordinate
%          system, each with [range;Az;El] components.
% systemType An optional parameter specifying the axis from which the
%          angles are measured in radians. Possible values are
%          0 (The default if omitted) Azimuth is measured 
%            counterclockwise from the x-axis in the x-y plane. Elevation
%            is measured up from the x-y plane (towards the z-axis). This
%            is consistent with common spherical coordinate systems for
%            specifying longitude (azimuth) and geocentric latitude
%            (elevation).
%          1 Azimuth is measured counterclockwise from the z-axis in the
%            z-x plane. Elevation is measured up from the z-x plane
%            (towards the y-axis). This is consistent with some spherical
%            coordinate systems that use the z axis as the boresight
%            direction of the radar.
%          2 This is the same as 0 except instead of being given
%            elevation, one desires the angle away from the z-axis, which
%            is (pi/2-elevation).
%          3 This is the same as 0 except azimuth is measured clockwise
%            from the y-axis in the x-y plane instead of counterclockwise
%            from the x-axis. This coordinate system often arises when
%            given "bearings" in a local East-North-Up coordinate system,
%            where the bearing directions are measured East of North.
%
%OUTPUTS: JTotal A 3X3XN Jacobian matrices where the rows in each matrix
%           are [x;y;z] in that order and the columns take the derivative
%           of the row component with respect to [r,azimuth,elevation] in
%           that order.
%
%This function evaluates analytic expressions for the Jacobian matrix that
%were derived by differentiating standard equations for the spherical
%coordinate systems.
%
%EXAMPLE:
%Here, we verify that the inverse of calcSpherJacob is equal to
%calcSpherInvJacob within reasonable finite precision limits.
% z=[1e3;0.3;0.25];
% systemType=2;
% J=calcSpherInvJacob(z,systemType)
% zC=spher2Cart(z,systemType);
% J1Inv=calcSpherJacob(zC,systemType);
% J1=inv(J1Inv)
%The difference between J and J1 values is small (and the values themselves
%are reasonable in magnitude), so we see that the results of this function
%are consistent with inverting the Jacobian from calcSpherJacob.
%
%June 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

if(nargin<2||isempty(systemType))
   systemType=0; 
end

N=size(z,2);

JTotal=zeros(3,3,N);

for curPoint=1:N
    r=z(1,curPoint);
    Az=z(2,curPoint);
    El=z(3,curPoint);
    
    sinAz=sin(Az);
    cosAz=cos(Az);
    sinEl=sin(El);
    cosEl=cos(El);
    
    J=zeros(3,3);
    switch(systemType)
        case 0
            %dx/dr
            J(1,1)=cosAz*cosEl;
            %dy/dr
            J(2,1)=cosEl*sinAz;
            %dz/dr
            J(3,1)=sinEl;
            %dx/dAz
            J(1,2)=-r*cosEl*sinAz;
            %dy/dAz
            J(2,2)=r*cosAz*cosEl;
            %dz/dAz
            J(3,2)=0;
            %dx/dEl
            J(1,3)=-r*cosAz*sinEl;
            %dy/dEl
            J(2,3)=-r*sinAz*sinEl;
            %dz/dEl
            J(3,3)=r*cosEl;
        case 1
            %dx/dr
            J(1,1)=cosEl*sinAz;
            %dy/dr
            J(2,1)=sinEl;
            %dz/dr
            J(3,1)=cosAz*cosEl;
            %dx/dAz
            J(1,2)=r*cosAz*cosEl;
            %dy/dAz
            J(2,2)=0;
            %dz/dAz
            J(3,2)=-r*cosEl*sinAz;
            %dx/dEl
            J(1,3)=-r*sinAz*sinEl;
            %dy/dEl
            J(2,3)=r*cosEl;
            %dz/dEl
            J(3,3)=-r*cosAz*sinEl;
        case 2
            %dx/dr
            J(1,1)=cosAz*sinEl;
            %dy/dr
            J(2,1)=sinAz*sinEl;
            %dz/dr
            J(3,1)=cosEl;
            %dx/dAz
            J(1,2)=-r*sinAz*sinEl;
            %dy/dAz
            J(2,2)=r*cosAz*sinEl;
            %dz/dAz
            J(3,2)=0;
            %dx/dEl
            J(1,3)=r*cosAz*cosEl;
            %dy/dEl
            J(2,3)=r*cosEl*sinAz;
            %dz/dEl
            J(3,3)=-r*sinEl;
        case 3
            %dx/dr
            J(1,1)=sinAz*cosEl;
            %dy/dr
            J(2,1)=cosAz*cosEl;
            %dz/dr
            J(3,1)=sinEl;
            %dx/dAz
            J(1,2)=r*cosEl*cosAz;
            %dy/dAz
            J(2,2)=-r*sinAz*cosEl;
            %dz/dAz
            J(3,2)=0;
            %dx/dEl
            J(1,3)=-r*sinAz*sinEl;
            %dy/dEl
            J(2,3)=-r*cosAz*sinEl;
            %dz/dEl
            J(3,3)=r*cosEl;
        otherwise
            error('Invalid system type specified.')
    end
    JTotal(:,:,curPoint)=J;
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
