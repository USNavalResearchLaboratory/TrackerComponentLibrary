function JTotal=spherAngUvCrossGrad(uv,systemType,Ms,Muv)
%%SPHERANGUVCROSSGRAD Determine the partial derivatives of u-v direction
%                  cosine components with respect to 3D spherical angles.
%
%INPUTS: uv A 2XnumPoints or 3XnumPoints (if the third components of the
%          unit vector) set of direction [u;v;w] cosines values in 3D. If
%          the third component of the unit vector is omitted, it is assumed
%          to be positive.
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
%          3 This is the same as 0 except azimuth is measured clockwise
%            from the y-axis in the x-y plane instead of counterclockwise
%            from the x-axis. This coordinate system often arises when
%            given "bearings" in a local East-North-Up coordinate system,
%            where the bearing directions are measured East of North.
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
%OUTPUTS: JTotal A 2X2XnumPoints set of Jacobian matrices where the rows
%            are azimuth and elevation the columns take the derivative of
%            the row components with respect to u and v in that order. Note
%            that singularities exist, in which case NaNs can be returned.
%
%The Jacobian is a matrix of partial derivatives. The Jacobian returned by
%uvSpherAngCrossGrad is the matrix inverse of that returned by this
%function.
%
%EXAMPLE:
%Here, we verify that the partial derivatives returned by this function are
%about equal to those returned via numeric differentiation (forward
%differencing).
% points=[[0.1;0.2],[0.1;0.1],[0;0]];%AzEl points
% epsVal=1e-9;
% systemType=3;
% Ms=randRotMat(3);
% Muv=randRotMat(3);
% 
% azEl=uv2SpherAng(points,systemType,Ms,Muv);
% azEl1=uv2SpherAng(bsxfun(@plus,points,[epsVal;0]),systemType,Ms,Muv);
% numDiffU=(azEl1-azEl)/epsVal;
% azEl1=uv2SpherAng(bsxfun(@plus,points,[0;epsVal]),systemType,Ms,Muv);
% numDiffV=(azEl1-azEl)/epsVal;
% gradVal=spherAngUvCrossGrad(points,systemType,Ms,Muv);
% 
% max(abs(numDiffU(:)-vec(gradVal(:,1,:))))
% max(abs(numDiffV(:)-vec(gradVal(:,2,:))))
%One will see that both numeric differences are on the order of 2e-7, which
%is a good amount of agreement.
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

N=size(uv,2);
JTotal=zeros(2,2,N);
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

    J=zeros(2,2);
    switch(systemType)
        case 0
            denom1=u1^2+v1^2;
            denom2=sqrt(1-w1^2);
            %Derivatives with respect to u.
            J(1,1)=(-v1*du1du+u1*dv1du)/denom1;
            J(2,1)=dw1du/denom2;
            %Derivatives with respect to v.
            J(1,2)=(-v1*du1dv+u1*dv1dv)/denom1;
            J(2,2)=dw1dv/denom2;

%If there were no rotations, it would be:
%             %Derivatives with respect to u.
%             J(1,1)=-v/(u^2+v^2);
%             J(2,1)=-u/(w*sqrt(u^2+v^2));
%             %Derivatives with respect to v.
%             J(1,2)=u/(u^2+v^2);
%             J(2,2)=-v/(w*sqrt(u^2+v^2));
        case 1
            denom1=u1^2+w1^2;
            denom2=sqrt(1-v1^2);
            
            %Derivatives with respect to u.
            J(1,1)=(w1*du1du-u1*dw1du)/denom1;
            J(2,1)=dv1du/denom2;
            %Derivatives with respect to v.
            J(1,2)=(w1*du1dv-u1*dw1dv)/denom1;
            J(2,2)=dv1dv/denom2;

%If there were no rotations, it would be:
%             %Derivatives with respect to u.
%             J(1,1)=1/w;
%             J(2,1)=0;
%             %Derivatives with respect to v.
%             J(1,2)=u*v/((1-v^2)*w);
%             J(2,2)=1/sqrt(1-v^2);
        case 2
            denom1=u1^2+v1^2;
            denom2=sqrt(1-w1^2);
            %Derivatives with respect to u.
            J(1,1)=(-v1*du1du+u1*dv1du)/denom1;
            J(2,1)=-dw1du/denom2;
            %Derivatives with respect to v.
            J(1,2)=(-v1*du1dv+u1*dv1dv)/denom1;
            J(2,2)=-dw1dv/denom2;

%If there were no rotations, it would be:
%             %Derivatives with respect to u.
%             J(1,1)=-v/(u^2+v^2);
%             J(2,1)=u/(w*sqrt(u^2+v^2));
%             %Derivatives with respect to v.
%             J(1,2)=u/(u^2+v^2);
%             J(2,2)=v/(w*sqrt(u^2+v^2));
        case 3
            denom1=u1^2+v1^2;
            denom2=sqrt(1-w1^2);
            %Derivatives with respect to u.
            J(1,1)=(v1*du1du-u1*dv1du)/denom1;
            J(2,1)=dw1du/denom2;
            %Derivatives with respect to v.
            J(1,2)=(v1*du1dv-u1*dv1dv)/denom1;
            J(2,2)=dw1dv/denom2;

%If there were no rotations, it would be:
%             %Derivatives with respect to u.
%             J(1,1)=v/(u^2+v^2);
%             J(2,1)=-u/(w*sqrt(u^2+v^2));
%             %Derivatives with respect to v.
%             J(1,2)=-u/(u^2+v^2);
%             J(2,2)=-v/(w*sqrt(u^2+v^2));
        otherwise
            error('Invalid system type specified')
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
