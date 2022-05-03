function HTotal=spherAngHessian(xG,systemType,lRx,M)
%%SPHERANGHESSIAN Determine the Hessian matrix (a matrix of second partial
%          derivatives) of a 3D spherical azimuth and elevation measurement
%          with respect to 3D position. Relativity and atmospheric effects
%          are not taken into account.
%
%INPUTS: xG The 3XN position vectors in the global coordinate system, each
%          with [x;y;z] components.
% systemType An optional parameter specifying the axes from which the
%          angles are measured in radians. Possible vaues are
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
%          3 This is the same as 0 except azimuth is measured clockwise
%            from the y-axis in the x-y plane instead of counterclockwise
%            from the x-axis. This coordinate system often arises when
%            given "bearings" in a local East-North-Up coordinate system,
%            where the bearing directions are measured East of North.
%      lRx The 3X1 position vector of the receiver. If omitted, the
%          receiver is placed at the origin. All vectors in x are assumed
%          to be from the same receiver.
%        M A 3X3 rotation matrix from the global Coordinate system to the
%          orientation of the coordinate system at the receiver. If
%          omitted, it is assumed to be the identity matrix. All vectors in
%          x are assumed to have the same rotation matrix.
%
%OUTPUTS: HTotal A 3X3X2XN matrix such that HTotal(:,:,1,i) is the Hessian
%          matrix with respect to the azimuth component and HTotal(:,:,2,i)
%          is the Hessian matrix with respect to the elevation component.
%          The elements in the matrices for each component/ point are
%          ordered [d^2/(dxdx), d^2/(dxdy), d^2/(dxdz);
%                   d^2/(dydx), d^2/(dydy), d^2/(dydz);
%                   d^2/(dzdx), d^2/(dzdy), d^2/(dzdz)];
%          note that each matrix is symmetric (i.e.
%           d^2/(dydx)=d^2/(dxdy) ).
%
%EXAMPLE:
%Here, we verify that a numerically differentiated Hessian is consistent
%with the analytic one produced by this function.
% x=[100;-1000;500];
% lRx=[500;20;-400];
% epsVal=1e-5;
% systemType=0;
% M=randRotMat(3);
% 
% H=spherAngHessian(x,systemType,lRx,M);
% 
% J=spherAngGradient(x,systemType,lRx,M);
% JdX=spherAngGradient(x+[epsVal;0;0],systemType,lRx,M);
% JdY=spherAngGradient(x+[0;epsVal;0],systemType,lRx,M);
% JdZ=spherAngGradient(x+[0;0;epsVal],systemType,lRx,M);
% HNumDiff=zeros(3,3,2);
% HNumDiff(1,1,1)=(JdX(1,1)-J(1,1))/epsVal;
% HNumDiff(1,2,1)=(JdX(1,2)-J(1,2))/epsVal;
% HNumDiff(2,1,1)=HNumDiff(1,2,1);
% HNumDiff(1,3,1)=(JdX(1,3)-J(1,3))/epsVal;
% HNumDiff(3,1,1)=HNumDiff(1,3,1);
% HNumDiff(2,2,1)=(JdY(1,2)-J(1,2))/epsVal;
% HNumDiff(2,3,1)=(JdY(1,3)-J(1,3))/epsVal;
% HNumDiff(3,2,1)=HNumDiff(2,3,1);
% HNumDiff(3,3,1)=(JdZ(1,3)-J(1,3))/epsVal;
% 
% HNumDiff(1,1,2)=(JdX(2,1)-J(2,1))/epsVal;
% HNumDiff(1,2,2)=(JdX(2,2)-J(2,2))/epsVal;
% HNumDiff(2,1,2)=HNumDiff(1,2,2);
% HNumDiff(1,3,2)=(JdX(2,3)-J(2,3))/epsVal;
% HNumDiff(3,1,2)=HNumDiff(1,3,2);
% HNumDiff(2,2,2)=(JdY(2,2)-J(2,2))/epsVal;
% HNumDiff(2,3,2)=(JdY(2,3)-J(2,3))/epsVal;
% HNumDiff(3,2,2)=HNumDiff(2,3,2);
% HNumDiff(3,3,2)=(JdZ(2,3)-J(2,3))/epsVal;
% 
% max(abs((H(:)-HNumDiff(:))./H(:)))
%The relative error will be on the order of 1e-6 or better, indicating good
%agreement between the numerical Hessian matrix and the actual Hessian
%matrix.
%
%June 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(M))
   M=eye(3,3);
end

if(nargin<3||isempty(lRx))
   lRx=zeros(3,1); 
end

if(nargin<2||isempty(systemType))
   systemType=0; 
end

N=size(xG,2);

xLocal=M*bsxfun(@minus,xG(1:3,:),lRx(1:3));

HTotal=zeros(3,3,2,N);
for curPoint=1:N
    x=xLocal(1,curPoint);
    y=xLocal(2,curPoint);
    z=xLocal(3,curPoint);
    
    H=zeros(3,3,2);
    switch(systemType)
        case 0
            r2xy=x^2+y^2;
            r4=(x^2+y^2+z^2)^2;
            rxy=sqrt(r2xy);

            %dAzdxdx
            H(1,1,1)=(2*x*y)/r2xy^2;
            %dAzdydy
            H(2,2,1)=-(2*x*y)/r2xy^2;
            %dAzdzdz
            H(3,3,1)=0;
            %dAzdxdy
            H(1,2,1)=-(x-y)*(x+y)/r2xy^2;
            %dAzdydx
            H(2,1,1)=H(1,2,1);
            %dAzdxdz
            H(1,3,1)=0;
            %dAzdzdx
            H(3,1,1)=H(1,3,1);
            %dAzdydz
            H(2,3,1)=0;
            %dAzdzdy
            H(3,2,1)=H(2,3,1);

            %dEldxdx
            H(1,1,2)=(z*(2*x^4+x^2*y^2-y^2*(y^2+z^2)))/(rxy^3*r4);
            %dEldydy
            H(2,2,2)=-((z*(x^4-2*y^4+x^2*(-y^2+z^2)))/(rxy^3*r4));
            %dEldzdz
            H(3,3,2)=-(2*rxy*z)/r4;
            %dEldxdy
            H(1,2,2)=(x*y*z*(3*r2xy+z^2))/(rxy^3*r4);
            %dEldydx
            H(2,1,2)=H(1,2,2);
            %dEldxdz
            H(1,3,2)=-((x*(x^2+y^2-z^2))/(rxy*r4));
            %dEldzdx
            H(3,1,2)=H(1,3,2);
            %dEldydz
            H(2,3,2)=-((y*(x^2+y^2-z^2))/(rxy*r4));
            %dEldzdy
            H(3,2,2)=H(2,3,2);
        case 1
            r2xz=x^2+z^2;
            r4=(x^2+y^2+z^2)^2;
            rxz=sqrt(r2xz);

            %dAzdxdx
            H(1,1,1)=-((2*x*z)/r2xz^2);
            %dAzdydy
            H(2,2,1)=0;
            %dAzdzdz
            H(3,3,1)=(2*x*z)/r2xz^2;
            %dAzdxdy
            H(1,2,1)=0;
            %dAzdydx
            H(2,1,1)=H(1,2,1);
            %dAzdxdz
            H(1,3,1)=((x-z)*(x+z))/r2xz^2;
            %dAzdzdx
            H(3,1,1)=H(1,3,1);
            %dAzdydz
            H(2,3,1)=0;
            %dAzdzdy
            H(3,2,1)=H(2,3,1);

            %dEldxdx
            H(1,1,2)=(y*(2*x^4+x^2*z^2-z^2*(y^2+z^2)))/(rxz^3*r4);
            %dEldydy
            H(2,2,2)=-((2*y*rxz)/(x^2+y^2+z^2)^2);
            %dEldzdz
            H(3,3,2)=-((y*(x^4-2*z^4+x^2*(y-z)*(y+z)))/(rxz^3*r4));
            %dEldxdy
            H(1,2,2)=-((x*(x^2-y^2+z^2))/(rxz*r4));
            %dEldydx
            H(2,1,2)=H(1,2,2);
            %dEldxdz
            H(1,3,2)=(x*y*z*(3*x^2+y^2+3*z^2))/(rxz^3*r4);
            %dEldzdx
            H(3,1,2)=H(1,3,2);
            %dEldydz
            H(2,3,2)=-((z*(x^2-y^2+z^2))/(rxz*r4));
            %dEldzdy
            H(3,2,2)=H(2,3,2);
        case 2
            r2xy=x^2+y^2;
            r4=(x^2+y^2+z^2)^2;
            rxy=sqrt(r2xy);

            %dAzdxdx
            H(1,1,1)=(2*x*y)/r2xy^2;
            %dAzdydy
            H(2,2,1)=-((2*x*y)/r2xy^2);
            %dAzdzdz
            H(3,3,1)=0;
            %dAzdxdy
            H(1,2,1)=-(x-y)*(x+y)/r2xy^2;
            %dAzdydx
            H(2,1,1)=H(1,2,1);
            %dAzdxdz
            H(1,3,1)=0;
            %dAzdzdx
            H(3,1,1)=H(1,3,1);
            %dAzdydz
            H(2,3,1)=0;
            %dAzdzdy
            H(3,2,1)=H(2,3,1);

            %dEldxdx
            H(1,1,2)=(z*(-2*x^4-x^2*y^2+y^4+y^2*z^2))/(rxy^3*r4);
            %dEldydy
            H(2,2,2)=(z*(x^4-2*y^4+x^2*(-y^2+z^2)))/(rxy^3*r4);
            %dEldzdz
            H(3,3,2)=(2*rxy*z)/r4;
            %dEldxdy
            H(1,2,2)=-((x*y*z*(3*r2xy+z^2))/(rxy^3*r4));
            %dEldydx
            H(2,1,2)=H(1,2,2);
            %dEldxdz
            H(1,3,2)=(x*(x^2+y^2-z^2))/(rxy*r4);
            %dEldzdx
            H(3,1,2)=H(1,3,2);
            %dEldydz
            H(2,3,2)=(y*(x^2+y^2-z^2))/(rxy*r4);
            %dEldzdy
            H(3,2,2)=H(2,3,2);
        case 3
            r2xy=x^2+y^2;
            r4=(x^2+y^2+z^2)^2;
            rxy=sqrt(r2xy);

            %dAzdxdx
            H(1,1,1)=-(2*x*y)/r2xy^2;
            %dAzdydy
            H(2,2,1)=(2*x*y)/r2xy^2;
            %dAzdzdz
            H(3,3,1)=0;
            %dAzdxdy
            H(1,2,1)=(x-y)*(x+y)/r2xy^2;
            %dAzdydx
            H(2,1,1)=H(1,2,1);
            %dAzdxdz
            H(1,3,1)=0;
            %dAzdzdx
            H(3,1,1)=H(1,3,1);
            %dAzdydz
            H(2,3,1)=0;
            %dAzdzdy
            H(3,2,1)=H(2,3,1);

            %dEldxdx
            H(1,1,2)=(z*(2*x^4+x^2*y^2-y^2*(y^2+z^2)))/(rxy^3*r4);
            %dEldydy
            H(2,2,2)=-((z*(x^4-2*y^4+x^2*(-y^2+z^2)))/(rxy^3*r4));
            %dEldzdz
            H(3,3,2)=-(2*rxy*z)/r4;
            %dEldxdy
            H(1,2,2)=(x*y*z*(3*r2xy+z^2))/(rxy^3*r4);
            %dEldydx
            H(2,1,2)=H(1,2,2);
            %dEldxdz
            H(1,3,2)=-((x*(x^2+y^2-z^2))/(rxy*r4));
            %dEldzdx
            H(3,1,2)=H(1,3,2);
            %dEldydz
            H(2,3,2)=-((y*(x^2+y^2-z^2))/(rxy*r4));
            %dEldzdy
            H(3,2,2)=H(2,3,2);
        otherwise
            error('Invalid system type specified.')
    end
    %Adjust for the rotated cordinate system.
    HTotal(:,:,1,curPoint)=M'*H(:,:,1)*M;
    HTotal(:,:,2,curPoint)=M'*H(:,:,2)*M;
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
