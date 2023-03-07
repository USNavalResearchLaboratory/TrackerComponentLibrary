function HTotal=calcSpherInvHessian(z,systemType)
%%SPHERANGINVHESSIAN Determine the Hessian matrix (a matrix of second
%          partial derivatives) of a 3D Cartesian point with respect to
%          monostatic spherical range, azimuth, and elevation components.
%          This produces second derivatives of (x,y,z) with respect to
%          (r,Az,El). The function calcRuvHessian produces second
%          derivatives of (r,Az,El) with respect to (x,y,z) in a more
%          general context.
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
%OUTPUTS: HTotal The 3X3X3XN set of Hessian matrices, where H(:,:,1,i) is
%           the Hessian for the x component for the ith point,
%           H(:,:,2) is the Hessian for the y component, and H(:,:,3) is
%           the Hessian for the z component. The ordering of the
%           derivatives in each matrix is:
%                  [d^2/(drdr),  d^2/(drdAz),  d^2/(drdEl);
%                   d^2/(dAzdr), d^2/(dAzdAz), d^2/(dAzdEl);
%                   d^2/(dEldr), d^2/(dEldAz), d^2/(dEldEl)];
%           note that each matrix is symmetric (i.e. 
%           d^2/(dAzdr)=d^2/(drdAz) ).
%
%This function evaluates analytic expressions for the Hessian matrix that
%were derived by differentiating standard equations for the spherical
%coordinate systems.
%
%EXAMPLE:
%Here, we verify that a numerically differentiated Hessian is consistent
%with the analytic one produced by this function.
% z=[1e3;0.3;0.25];
% systemType=2;
% epsVal=1e-7;
% H=calcSpherInvHessian(z,systemType);
% J=calcSpherInvJacob(z,systemType);
% 
% Jdr=calcSpherInvJacob(z+[epsVal;0;0],systemType);
% JdAz=calcSpherInvJacob(z+[0;epsVal;0],systemType);
% JdEl=calcSpherInvJacob(z+[0;0;epsVal],systemType);
% 
% HNumeric=zeros(3,3,3);
% %d2xdrdr
% HNumeric(1,1,1)=(Jdr(1,1)-J(1,1))/epsVal;
% %d2ydrdr
% HNumeric(1,1,2)=(Jdr(2,1)-J(2,1))/epsVal;
% %d2zdrdr
% HNumeric(1,1,3)=(Jdr(3,1)-J(3,1))/epsVal;
% 
% %d2xdAzdAz
% HNumeric(2,2,1)=(JdAz(1,2)-J(1,2))/epsVal;
% %d2ydAzdAz
% HNumeric(2,2,2)=(JdAz(2,2)-J(2,2))/epsVal;
% %d2zdAzdAz
% HNumeric(2,2,3)=(JdAz(3,2)-J(3,2))/epsVal;
% 
% %d2xdEldEl
% HNumeric(3,3,1)=(JdEl(1,3)-J(1,3))/epsVal;
% %d2ydEldEl
% HNumeric(3,3,2)=(JdEl(2,3)-J(2,3))/epsVal;
% %d2zdEldEl
% HNumeric(3,3,3)=(JdEl(3,3)-J(3,3))/epsVal;
% 
% %d2xdrdAz
% HNumeric(1,2,1)=(JdAz(1,1)-J(1,1))/epsVal;
% HNumeric(2,1,1)=HNumeric(1,2,1);
% %d2ydrdAz
% HNumeric(1,2,2)=(JdAz(2,1)-J(2,1))/epsVal;
% HNumeric(2,1,2)=HNumeric(1,2,2);
% %d2zdrdAz
% HNumeric(1,2,3)=(JdAz(3,1)-J(3,1))/epsVal;
% HNumeric(2,1,3)=HNumeric(1,2,3);
% 
% %d2xdrdEl
% HNumeric(1,3,1)=(JdEl(1,1)-J(1,1))/epsVal;
% HNumeric(3,1,1)=HNumeric(1,3,1);
% %d2ydrdEl
% HNumeric(1,3,2)=(JdEl(2,1)-J(2,1))/epsVal;
% HNumeric(3,1,2)=HNumeric(1,3,2);
% %d2zdrdEl
% HNumeric(1,3,3)=(JdEl(3,1)-J(3,1))/epsVal;
% HNumeric(3,1,3)=HNumeric(1,3,3);
% 
% %d2xdAzdEl
% HNumeric(2,3,1)=(JdEl(1,2)-J(1,2))/epsVal;
% HNumeric(3,2,1)=HNumeric(2,3,1);
% %d2ydAzdEl
% HNumeric(2,3,2)=(JdEl(2,2)-J(2,2))/epsVal;
% HNumeric(3,2,2)=HNumeric(2,3,2);
% %d2zdAzdEl
% HNumeric(2,3,3)=(JdEl(3,2)-J(3,2))/epsVal;
% HNumeric(3,2,3)=HNumeric(2,3,3);
% 
% diffMat=HNumeric-H;
% max(abs(diffMat(:)./H(:)))
%The relative error will be on the order of 2e-7, indicating good agreement
%between the numerical Hessian matrix and the actual Hessian matrix.
%
%June 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

if(nargin<2||isempty(systemType))
   systemType=0; 
end

N=size(z,2);

HTotal=zeros(3,3,3,N);

for curPoint=1:N
    r=z(1,curPoint);
    Az=z(2,curPoint);
    El=z(3,curPoint);
    
    sinAz=sin(Az);
    cosAz=cos(Az);
    sinEl=sin(El);
    cosEl=cos(El);
    
    H=zeros(3,3,3);
    switch(systemType)
        case 0
            %d2xdrdr
            H(1,1,1)=0;
            %d2xdAzdAz
            H(2,2,1)=-r*cosAz*cosEl;
            %d2xdEldEl
            H(3,3,1)=-r*cosAz*cosEl;
            %d2xdrdAz
            H(1,2,1)=-cosEl*sinAz;
            %d2xdAzdr
            H(2,1,1)=H(1,2,1);
            %d2xdrdEl
            H(1,3,1)=-cosAz*sinEl;
            %d2xdEldr
            H(3,1,1)=H(1,3,1);
            %d2xdAzdEl
            H(2,3,1)=r*sinAz*sinEl;
            %d2xdEldAz
            H(3,2,1)=H(2,3,1);

            %d2ydrdr
            H(1,1,2)=0;
            %d2ydAzdAz
            H(2,2,2)=-r*cosEl*sinAz;
            %d2ydEldEl
            H(3,3,2)=-r*cosEl*sinAz;
            %d2ydrdAz
            H(1,2,2)=cosAz*cosEl;
            %d2ydAzdr
            H(2,1,2)=H(1,2,2);
            %d2ydrdEl
            H(1,3,2)=-sinAz*sinEl;
            %d2ydEldr
            H(3,1,2)=H(1,3,2);
            %d2ydAzdEl
            H(2,3,2)=-r*cosAz*sinEl;
            %d2ydEldAz
            H(3,2,2)=H(2,3,2);

            %d2zdrdr
            H(1,1,3)=0;
            %d2zdAzdAz
            H(2,2,3)=0;
            %d2zdEldEl
            H(3,3,3)=-r*sinEl;
            %d2zdrdAz
            H(1,2,3)=0;
            %d2zdAzdr
            H(2,1,3)=H(1,2,3);
            %d2zdrdEl
            H(1,3,3)=cosEl;
            %d2zdEldr
            H(3,1,3)=H(1,3,3);
            %d2zdAzdEl
            H(2,3,3)=0;
            %d2zdEldAz
            H(3,2,3)=H(2,3,3);
        case 1
            %d2xdrdr
            H(1,1,1)=0;
            %d2xdAzdAz
            H(2,2,1)=-r*cosEl*sinAz;
            %d2xdEldEl
            H(3,3,1)=-r*cosEl*sinAz;
            %d2xdrdAz
            H(1,2,1)=cosAz*cosEl;
            %d2xdAzdr
            H(2,1,1)=H(1,2,1);
            %d2xdrdEl
            H(1,3,1)=-sinAz*sinEl;
            %d2xdEldr
            H(3,1,1)=H(1,3,1);
            %d2xdAzdEl
            H(2,3,1)=-r*cosAz*sinEl;
            %d2xdEldAz
            H(3,2,1)=H(2,3,1);

            %d2ydrdr
            H(1,1,2)=0;
            %d2ydAzdAz
            H(2,2,2)=0;
            %d2ydEldEl
            H(3,3,2)=-r*sinEl;
            %d2ydrdAz
            H(1,2,2)=0;
            %d2ydAzdr
            H(2,1,2)=H(1,2,2);
            %d2ydrdEl
            H(1,3,2)=cosEl;
            %d2ydEldr
            H(3,1,2)=H(1,3,2);
            %d2ydAzdEl
            H(2,3,2)=0;
            %d2ydEldAz
            H(3,2,2)=H(2,3,2);

            %d2zdrdr
            H(1,1,3)=0;
            %d2zdAzdAz
            H(2,2,3)=-r*cosAz*cosEl;
            %d2zdEldEl
            H(3,3,3)=-r*cosAz*cosEl;
            %d2zdrdAz
            H(1,2,3)=-cosEl*sinAz;
            %d2zdAzdr
            H(2,1,3)=H(1,2,3);
            %d2zdrdEl
            H(1,3,3)=-cosAz*sinEl;
            %d2zdEldr
            H(3,1,3)=H(1,3,3);
            %d2zdAzdEl
            H(2,3,3)=r*sinAz*sinEl;
            %d2zdEldAz
            H(3,2,3)=H(2,3,3);
        case 2
            %d2xdrdr
            H(1,1,1)=0;
            %d2xdAzdAz
            H(2,2,1)=-r*cosAz*sinEl;
            %d2xdEldEl
            H(3,3,1)=-r*cosAz*sinEl;
            %d2xdrdAz
            H(1,2,1)=-sinAz*sinEl;
            %d2xdAzdr
            H(2,1,1)=H(1,2,1);
            %d2xdrdEl
            H(1,3,1)=cosAz*cosEl;
            %d2xdEldr
            H(3,1,1)=H(1,3,1);
            %d2xdAzdEl
            H(2,3,1)=-r*cosEl*sinAz;
            %d2xdEldAz
            H(3,2,1)=H(2,3,1);
        
            %d2ydrdr
            H(1,1,2)=0;
            %d2ydAzdAz
            H(2,2,2)=-r*sinAz*sinEl;
            %d2ydEldEl
            H(3,3,2)=-r*sinAz*sinEl;
            %d2ydrdAz
            H(1,2,2)=cosAz*sinEl;
            %d2ydAzdr
            H(2,1,2)=H(1,2,2);
            %d2ydrdEl
            H(1,3,2)=cosEl*sinAz;
            %d2ydEldr
            H(3,1,2)=H(1,3,2);
            %d2ydAzdEl
            H(2,3,2)=r*cosAz*cosEl;
            %d2ydEldAz
            H(3,2,2)=H(2,3,2);

            %d2zdrdr
            H(1,1,3)=0;
            %d2zdAzdAz
            H(2,2,3)=0;
            %d2zdEldEl
            H(3,3,3)=-r*cosEl;
            %d2zdrdAz
            H(1,2,3)=0;
            %d2zdAzdr
            H(2,1,3)=H(1,2,3);
            %d2zdrdEl
            H(1,3,3)=-sinEl;
            %d2zdEldr
            H(3,1,3)=H(1,3,3);
            %d2zdAzdEl
            H(2,3,3)=0;
            %d2zdEldAz
            H(3,2,3)=H(2,3,3);
        case 3
            %d2xdrdr
            H(1,1,1)=0;
            %d2xdAzdAz
            H(2,2,1)=-r*sinAz*cosEl;
            %d2xdEldEl
            H(3,3,1)=-r*sinAz*cosEl;
            %d2xdrdAz
            H(1,2,1)=cosEl*cosAz;
            %d2xdAzdr
            H(2,1,1)=H(1,2,1);
            %d2xdrdEl
            H(1,3,1)=-sinAz*sinEl;
            %d2xdEldr
            H(3,1,1)=H(1,3,1);
            %d2xdAzdEl
            H(2,3,1)=-r*cosAz*sinEl;
            %d2xdEldAz
            H(3,2,1)=H(2,3,1);

            %d2ydrdr
            H(1,1,2)=0;
            %d2ydAzdAz
            H(2,2,2)=-r*cosEl*cosAz;
            %d2ydEldEl
            H(3,3,2)=-r*cosEl*cosAz;
            %d2ydrdAz
            H(1,2,2)=-sinAz*cosEl;
            %d2ydAzdr
            H(2,1,2)=H(1,2,2);
            %d2ydrdEl
            H(1,3,2)=-cosAz*sinEl;
            %d2ydEldr
            H(3,1,2)=H(1,3,2);
            %d2ydAzdEl
            H(2,3,2)=r*sinAz*sinEl;
            %d2ydEldAz
            H(3,2,2)=H(2,3,2);

            %d2zdrdr
            H(1,1,3)=0;
            %d2zdAzdAz
            H(2,2,3)=0;
            %d2zdEldEl
            H(3,3,3)=-r*sinEl;
            %d2zdrdAz
            H(1,2,3)=0;
            %d2zdAzdr
            H(2,1,3)=H(1,2,3);
            %d2zdrdEl
            H(1,3,3)=cosEl;
            %d2zdEldr
            H(3,1,3)=H(1,3,3);
            %d2zdAzdEl
            H(2,3,3)=0;
            %d2zdEldAz
            H(3,2,3)=H(2,3,3);
        otherwise
            error('Invalid system type specified.')
    end
    HTotal(:,:,:,curPoint)=H;
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
