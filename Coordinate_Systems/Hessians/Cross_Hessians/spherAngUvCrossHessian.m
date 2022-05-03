function HTotal=spherAngUvCrossHessian(uv,systemType,Ms,Muv)
%%SPHERANGUVCROSSHESSIAN Determine second partial derivative matrices of 3D
%               spherical angular components with respect to u-v direction
%               cosines.
%
%INPUTS: uv A 2XN or 3XN (if the third components of the unit vector) set
%          of direction [u;v;w] cosines values in 3D. If the third
%          component of the unit vector is omitted, it is assumed to be
%          positive.
% systemType An optional parameter specifying the axes from which the
%          angles for the spherical coordinate system are measured in
%          radians. Possible vaues are
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
%   Ms,Muv If either the spherical coordinate system or the u-v coordinate
%          system is rotated compared to the global Cartesian coordinate
%          system, these optional 3X3 matrices provide the rotations. Ms
%          is a 3X3 matrix to go from the alignment of a global
%          Cartesian coordinate system to that in which the spherical
%          coordinates are computed. Similarly, Muv is a rotation matrix
%          to go from the alignment of a global Cartesian cordinate system
%          to that in which the u-v(-w) coordinates are computed. If
%          either of these in omitted or an empty matrix is passed, then
%          the missing one is replaced with the identity matrix.
%
%OUTPUTS: HTotal A 2X2X2XN matrix of second derivatives where
%                HTotal(:,:,i,j) is the Hessian matrix of the ith
%                component of [azimuth;elevation] evaluated at the jth
%                point in uv. The ordering of the second derivatives in the
%                i,jth Hessian matrix is [d/du^2, d/(dudv);
%                                         d/(dudv), d/dv^2]
%
%EXAMPLE:
%Here, we verify that the Hessian matrix computed by this function is close
%to the computed using forward differencing of the gradient.
% systemType=0;
% uv=[0.1;-0.2];
% 
% H=spherAngUvCrossHessian(uv,systemType);
% J=spherAngUvCrossGrad(uv,systemType);
% epsVal=1e-8;
% uv1=uv+[epsVal;0];
% J1=spherAngUvCrossGrad(uv1,systemType);
% dAz=(J1-J)/epsVal;
% 
% uv1=uv+[0;epsVal];
% J1=spherAngUvCrossGrad(uv1,systemType);
% dEl=(J1-J)/epsVal;
% 
% HNumDiff=zeros(2,2,2,1);
% 
% %Derivatives of azimuth
% HNumDiff(1,1,1)=dAz(1,1);%dAz/du^2
% HNumDiff(1,2,1)=dAz(1,2);%dAz/(dudv)
% HNumDiff(2,1,1)=HNumDiff(1,2,1);
% HNumDiff(2,2,1)=dEl(1,2);%dAz/dv^2
% 
% %Derivatives of elevation
% HNumDiff(1,1,2)=dAz(2,1);%dEl/du^2
% HNumDiff(1,2,2)=dAz(2,2);%dEl/(dudv)
% HNumDiff(2,1,2)=HNumDiff(1,2,2);
% HNumDiff(2,2,2)=dEl(2,2);%dEl/dv^2
% 
% max(abs(HNumDiff(:)-H(:)))
%One will see that the difference between the true Hessian and the Hessian
%from numeric differentiation is on the order of 8e-7, which indicates good
%agreement.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(Muv))
    Muv=eye(3,3);
end

if(nargin<3||isempty(Ms))
    Ms=eye(3,3);
end

if(nargin<2||isempty(systemType))
   systemType=0; 
end

hasW=size(uv,1)>2;
N=size(uv,2);

M=Ms*Muv';
m11=M(1,1);
m12=M(1,2);
m13=M(1,3);
m21=M(2,1);
m22=M(2,2);
m23=M(2,3);
m31=M(3,1);
m32=M(3,2);
m33=M(3,3);

HTotal=zeros(2,2,2,N);
for curPoint=1:N
    if(hasW==false)
        uvwCur=[uv(:,curPoint);sqrt(1-uv(1,curPoint)^2-uv(2,curPoint)^2)];
    else
        uvwCur=uv(:,curPoint);
    end
    u=uvwCur(1);
    v=uvwCur(2);
    w=uvwCur(3);

    uvwCur=M*uvwCur;
    u1=uvwCur(1);
    v1=uvwCur(2);
    w1=uvwCur(3);

    %First derivatives
    du1du=m11-m13*u/w;
    du1dv=m12-m13*v/w;
    dv1du=m21-m23*u/w;
    dv1dv=m22-m23*v/w;
    dw1du=m31-m33*u/w;
    dw1dv=m32-m33*v/w;

    w3=w^3;
    %Second derivatives
    d2u1dudu=m13*(v^2-1)/w3;
    d2u1dudv=-m13*u*v/w3;
    d2u1dvdv=m13*(u^2-1)/w3;
    d2v1dudu=m23*(v^2-1)/w3;
    d2v1dudv=-m23*u*v/w3;
    d2v1dvdv=m23*(u^2-1)/w3;
    d2w1dudu=m33*(v^2-1)/w3;
    d2w1dudv=-m33*u*v/w3;
    d2w1dvdv=m33*(u^2-1)/w3;

    switch(systemType)
        case 0
            denom1=(u1^2+v1^2)^2;
            denom2=sqrt(1-w1^2)^3;

            dazdu2=(-2*u1*(2*du1du*dv1du*u1-du1du^2*v1+dv1du^2*v1)+(2*du1du*dv1du+d2v1dudu*u1-d2u1dudu*v1)*(u1^2+v1^2))/denom1;
            dazdudv=(d2v1dudv*u1^3+v1^2*(du1dv*dv1du+du1du*dv1dv-d2u1dudv*v1)-u1^2*(du1dv*dv1du+du1du*dv1dv+d2u1dudv*v1)+u1*v1*(2*du1du*du1dv-2*dv1du*dv1dv+d2v1dudv*v1))/denom1;
            dazdv2=(-2*u1*(2*du1dv*dv1dv*u1-du1dv^2*v1+dv1dv^2*v1)+(2*du1dv*dv1dv+d2v1dvdv*u1-d2u1dvdv*v1)*(u1^2+v1^2))/denom1;
            deldu2=(d2w1dudu+dw1du^2*w1-d2w1dudu*w1^2)/denom2;
            deldudv=(d2w1dudv+dw1du*dw1dv*w1-d2w1dudv*w1^2)/denom2;
            deldv2=(d2w1dvdv+dw1dv^2*w1-d2w1dvdv*w1^2)/denom2;

%If there were no rotations, it would be:
%             u2v2=u^2+v^2;
%             dazdu2=(2*u*v)/u2v2^2;
%             dazdudv=(-u^2+v^2)/u2v2^2;
%             dazdv2=-((2*u*v)/u2v2^2);
%             deldu2=-((u^4+v^2-v^4)/(w*sqrt(u2v2))^3);
%             deldudv=-((u*v*(-1+2*u^2+2*v^2))/(w*sqrt(u2v2))^3);
%             deldv2=-((u^2-u^4+v^4)/(w*sqrt(u2v2))^3);
        case 1
            denom1=(u1^2+w1^2)^2;
            denom2=sqrt(1-v1^2)^3;

            dazdu2=(2*u1*(2*du1du*dw1du*u1-du1du^2*w1+dw1du^2*w1)+(-2*du1du*dw1du-d2w1dudu*u1+d2u1dudu*w1)*(u1^2+w1^2))/denom1;
            dazdudv=(u1^2*(du1dv*dw1du+du1du*dw1dv-d2w1dudv*u1)+u1*(-2*du1du*du1dv+2*dw1du*dw1dv+d2u1dudv*u1)*w1-(du1dv*dw1du+du1du*dw1dv+d2w1dudv*u1)*w1^2+d2u1dudv*w1^3)/denom1;
            dazdv2=(2*u1*(2*du1dv*dw1dv*u1-du1dv^2*w1+dw1dv^2*w1)+(-2*du1dv*dw1dv-d2w1dvdv*u1+d2u1dvdv*w1)*(u1^2+w1^2))/denom1;
            deldu2=(d2v1dudu+dv1du^2*v1-d2v1dudu*v1^2)/denom2;
            deldudv=(d2v1dudv+dv1du*dv1dv*v1-d2v1dudv*v1^2)/denom2;
            deldv2=(d2v1dvdv+dv1dv^2*v1-d2v1dvdv*v1^2)/denom2;

%If there were no rotations, it would be:
%             dazdu2=u/w^3;
%             dazdudv=v/w^3;
%             dazdv2=-((u*(-1+u^2+(-1+u^2)*v^2+2*v^4))/(w^3*(-1+v^2)^2));
%             deldu2=0;
%             deldudv=0;
%             deldv2=v/(1-v^2)^(3/2);
        case 2
            denom1=(u1^2+v1^2)^2;
            denom2=sqrt(1-w1^2)^3;

            dazdu2=(-2*u1*(2*du1du*dv1du*u1-du1du^2*v1+dv1du^2*v1)+(2*du1du*dv1du+d2v1dudu*u1-d2u1dudu*v1)*(u1^2+v1^2))/denom1;
            dazdudv=(d2v1dudv*u1^3+v1^2*(du1dv*dv1du+du1du*dv1dv-d2u1dudv*v1)-u1^2*(du1dv*dv1du+du1du*dv1dv+d2u1dudv*v1)+u1*v1*(2*du1du*du1dv-2*dv1du*dv1dv+d2v1dudv*v1))/denom1;
            dazdv2=(-2*u1*(2*du1dv*dv1dv*u1-du1dv^2*v1+dv1dv^2*v1)+(2*du1dv*dv1dv+d2v1dvdv*u1-d2u1dvdv*v1)*(u1^2+v1^2))/denom1;
            deldu2=(-dw1du^2*w1+d2w1dudu*(-1+w1^2))/denom2;
            deldudv=(-dw1du*dw1dv*w1+d2w1dudv*(-1+w1^2))/denom2;
            deldv2=(-dw1dv^2*w1+d2w1dvdv*(-1+w1^2))/denom2;

%If there were no rotations, it would be:
%             u2v2=u^2+v^2;
%             dazdu2=(2*u*v)/u2v2^2;
%             dazdudv=(-u^2+v^2)/u2v2^2;
%             dazdv2=-((2*u*v)/u2v2^2);
%             deldu2=(u^4+v^2-v^4)/(w*sqrt(u2v2))^3;
%             deldudv=(u*v*(-1+2*u^2+2*v^2))/(w*sqrt(u2v2))^3;
%             deldv2=(u^2-u^4+v^4)/(w*sqrt(u2v2))^3;
        case 3
            denom1=(u1^2+v1^2)^2;
            denom2=sqrt(1-w1^2)^3;

            dazdu2=(2*u1*(2*du1du*dv1du*u1-du1du^2*v1+dv1du^2*v1)+(-2*du1du*dv1du-d2v1dudu*u1+d2u1dudu*v1)*(u1^2+v1^2))/denom1;
            dazdudv=(u1^2*(du1dv*dv1du+du1du*dv1dv-d2v1dudv*u1)+u1*(-2*du1du*du1dv+2*dv1du*dv1dv+d2u1dudv*u1)*v1-(du1dv*dv1du+du1du*dv1dv+d2v1dudv*u1)*v1^2+d2u1dudv*v1^3)/denom1;
            dazdv2=(2*u1*(2*du1dv*dv1dv*u1-du1dv^2*v1+dv1dv^2*v1)+(-2*du1dv*dv1dv-d2v1dvdv*u1+d2u1dvdv*v1)*(u1^2+v1^2))/denom1;
            deldu2=(d2w1dudu+dw1du^2*w1-d2w1dudu*w1^2)/denom2;
            deldudv=(d2w1dudv+dw1du*dw1dv*w1-d2w1dudv*w1^2)/denom2;
            deldv2=(d2w1dvdv+dw1dv^2*w1-d2w1dvdv*w1^2)/denom2;

%If there were no rotations, it would be:
%             u2v2=u^2+v^2;
%             dazdu2=-(2*u*v)/u2v2^2;
%             dazdudv=(u-v)*(u+v)/u2v2^2;
%             dazdv2=(2*u*v)/u2v2^2;
%             deldu2=-((u^4+v^2-v^4)/(w*sqrt(u2v2))^3);
%             deldudv=-((u*v*(-1+2*u^2+2*v^2))/(w*sqrt(u2v2))^3);
%             deldv2=-((u^2-u^4+v^4)/(w*sqrt(u2v2))^3);
        otherwise
            error('Invalid system type specified.')
    end

    H=zeros(2,2,2);
    %Derivatives of azimuth
    H(1,1,1)=dazdu2;
    H(1,2,1)=dazdudv;
    H(2,1,1)=H(1,2,1);
    H(2,2,1)=dazdv2;
    %Derivatives of elevation
    H(1,1,2)=deldu2;
    H(1,2,2)=deldudv;
    H(2,1,2)=H(1,2,2);
    H(2,2,2)=deldv2;
    
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
