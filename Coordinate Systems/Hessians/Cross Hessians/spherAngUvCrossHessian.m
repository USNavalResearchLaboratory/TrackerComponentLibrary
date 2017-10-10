function HTotal=spherAngUvCrossHessian(uv,systemType)
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

if(nargin<2||isempty(systemType))
   systemType=0; 
end

hasW=size(uv,1)>2;

N=size(uv,2);

HTotal=zeros(2,2,2,N);

for curPoint=1:N
    u=uv(1,curPoint);
    v=uv(2,curPoint);
    
    if(hasW)
        w=uv(3,curPoint);
    else
        w=sqrt(1-u^2-v^2);
    end

    switch(systemType)
        case 0
            u2v2=u^2+v^2;

            dazdu2=(2*u*v)/u2v2^2;
            dazdudv=(-u^2+v^2)/u2v2^2;
            dazdv2=-((2*u*v)/u2v2^2);
            deldu2=-((u^4+v^2-v^4)/(w*sqrt(u2v2))^3);
            deldudv=-((u*v*(-1+2*u^2+2*v^2))/(w*sqrt(u2v2))^3);
            deldv2=-((u^2-u^4+v^4)/(w*sqrt(u2v2))^3);
        case 1
            dazdu2=u/w^3;
            dazdudv=v/w^3;
            dazdv2=-((u*(-1+u^2+(-1+u^2)*v^2+2*v^4))/(w^3*(-1+v^2)^2));
            deldu2=0;
            deldudv=0;
            deldv2=v/(1-v^2)^(3/2);
        case 2
            u2v2=u^2+v^2;

            dazdu2=(2*u*v)/u2v2^2;
            dazdudv=(-u^2+v^2)/u2v2^2;
            dazdv2=-((2*u*v)/u2v2^2);
            deldu2=(u^4+v^2-v^4)/(w*sqrt(u2v2))^3;
            deldudv=(u*v*(-1+2*u^2+2*v^2))/(w*sqrt(u2v2))^3;
            deldv2=(u^2-u^4+v^4)/(w*sqrt(u2v2))^3;
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
