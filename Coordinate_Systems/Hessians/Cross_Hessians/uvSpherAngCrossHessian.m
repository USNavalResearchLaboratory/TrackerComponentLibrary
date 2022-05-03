function HTotal=uvSpherAngCrossHessian(azEl,systemType,includeW,Ms,Muv)
%%UVSPHERANGCROSSHESSIAN Determine second partial derivative matrices of 3D
%               direction cosines components u v and possibly w with
%               respect to 3D spherical angular coordinates.
%
%INPUTS: azEl A 2XN set of points in 3D [azimuth;elevation] in radians.
% systemType An optional parameter specifying the axes from which the
%          angles for the spherical coordinate system are measured in
%          radians. Possible vaues are:
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
% includeW If one wants to have second partial derivatives with respect to
%          the third component of the direction cosine unit vector (the z-
%          axis component), then set this to true. The default if omitted
%          or an empty matrix is passed is false.
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
%OUTPUTS: HTotal A 2X2X2XN (if includeW is false) or 2X2X3XN (if includeW
%                is true), matrix of second derivatives where
%                HTotal(:,:,i,j) is the Hessian matrix of the ith
%                component of [u;v;w] evaluated at the jth point in azEl.
%                The ordering of the second derivatives in the i,jth
%                Hessian matrix is [d/dAz^2, d/(dAzdEl);
%                                   d/(dAzdEl), d/dEl^2]
%
%EXAMPLE:
%Here, we verify that the Hessian matrix computed by this function is close
%to the computed using forward differencing of the gradient.
% systemType=2;
% includeW=true;
% azEl=[0.1;pi/4];
% 
% H=uvSpherAngCrossHessian(azEl,systemType,includeW);
% J=uvSpherAngCrossGrad(azEl,systemType,includeW);
% epsVal=1e-8;
% azEl1=azEl+[epsVal;0];
% J1=uvSpherAngCrossGrad(azEl1,systemType,includeW);
% dAz=(J1-J)/epsVal;
% 
% azEl1=azEl+[0;epsVal];
% J1=uvSpherAngCrossGrad(azEl1,systemType,includeW);
% dEl=(J1-J)/epsVal;
% 
% HNumDiff=zeros(2,2,3,1);
% 
% %Derivatives of u
% HNumDiff(1,1,1)=dAz(1,1);%du/dAz^2
% HNumDiff(1,2,1)=dAz(1,2);%du/(dAzdEl)
% HNumDiff(2,1,1)=HNumDiff(1,2,1);
% HNumDiff(2,2,1)=dEl(1,2);%du/dEl^2
% 
% %Derivatives of v
% HNumDiff(1,1,2)=dAz(2,1);%dv/dAz^2
% HNumDiff(1,2,2)=dAz(2,2);%dv/(dAzdEl)
% HNumDiff(2,1,2)=HNumDiff(1,2,2);
% HNumDiff(2,2,2)=dEl(2,2);%du/dEl^2
% 
% %Derivatives of w
% HNumDiff(1,1,3)=dAz(3,1);%dw/dAz^2
% HNumDiff(1,2,3)=dAz(3,2);%dw/(dAzdEl)
% HNumDiff(2,1,3)=HNumDiff(1,2,3);
% HNumDiff(2,2,3)=dEl(3,2);%dw/dEl^2
% 
% max(abs(HNumDiff(:)-H(:)))
%One will see that the difference between the true Hessian and the Hessian
%from numeric differentiation is on the order of 6e-9, which indicates good
%agreement.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(Muv))
    Muv=eye(3,3);
end

if(nargin<4||isempty(Ms))
    Ms=eye(3,3);
end

if(nargin<3||isempty(includeW))
    includeW=false;
end

if(nargin<2||isempty(systemType))
   systemType=0; 
end

N=size(azEl,2);

if(includeW)
    HDim=3;
else
    HDim=2;
end

HTotal=zeros(2,2,HDim,N);

M=Muv*Ms';
m11=M(1,1);
m12=M(1,2);
m13=M(1,3);
m21=M(2,1);
m22=M(2,2);
m23=M(2,3);
m31=M(3,1);
m32=M(3,2);
m33=M(3,3);

for curPoint=1:N
    az=azEl(1,curPoint);
    el=azEl(2,curPoint);
    
    sinAz=sin(az);
    cosAz=cos(az);
    sinEl=sin(el);
    cosEl=cos(el);
    switch(systemType)
        case 0
            dudaz2=-cosEl*(m11*cosAz+m12*sinAz);
            dudazdel=(-m12*cosAz+m11*sinAz)*sinEl;
            dudel2=-cosEl*(m11*cosAz+m12*sinAz)-m13*sinEl;
            dvdaz2=-cosEl*(m21*cosAz+m22*sinAz);
            dvdazdel=(-m22*cosAz+m21*sinAz)*sinEl;
            dvdel2=-cosEl*(m21*cosAz+m22*sinAz)-m23*sinEl;
            if(includeW)
                dwdaz2=-cosEl*(m31*cosAz+m32*sinAz);
                dwdazdel=(-m32*cosAz+m31*sinAz)*sinEl;
                dwdel2=-cosEl*(m31*cosAz+m32*sinAz)-m33*sinEl;
            end
%If there were no rotations, it would be:
%             dudaz2=-cosAz*cosEl;
%             dudazdel=sinAz*sinEl;
%             dudel2=-cosAz*cosEl;
%             dvdaz2=-cosEl*sinAz;
%             dvdazdel=-cosAz*sinEl;
%             dvdel2=-cosEl*sinAz;
%             if(includeW)
%                 dwdaz2=0;
%                 dwdazdel=0;
%                 dwdel2=-sinEl;
%             end
        case 1
            dudaz2=-cosEl*(m13*cosAz+m11*sinAz);
            dudazdel=(-m11*cosAz+m13*sinAz)*sinEl;
            dudel2=-cosEl*(m13*cosAz+m11*sinAz)-m12*sinEl;
            dvdaz2=-cosEl*(m23*cosAz+m21*sinAz);
            dvdazdel=(-m21*cosAz+m23*sinAz)*sinEl;
            dvdel2=-cosEl*(m23*cosAz+m21*sinAz)-m22*sinEl;
            if(includeW)
                dwdaz2=-cosEl*(m33*cosAz+m31*sinAz);
                dwdazdel=(-m31*cosAz+m33*sinAz)*sinEl;
                dwdel2=-cosEl*(m33*cosAz+m31*sinAz)-m32*sinEl;
            end
%If there were no rotations, it would be:
%             dudaz2=-cosEl*sinAz;
%             dudazdel=-cosAz*sinEl;
%             dudel2=-cosEl*sinAz;
%             dvdaz2=0;
%             dvdazdel=0;
%             dvdel2=-sinEl;
%             if(includeW)
%                 dwdaz2=-cosAz*cosEl;
%                 dwdazdel=sinAz*sinEl;
%                 dwdel2=-cosAz*cosEl;
%             end
        case 2
            dudaz2=-((m11*cosAz+m12*sinAz)*sinEl);
            dudazdel=cosEl*(m12*cosAz-m11*sinAz);
            dudel2=-m13*cosEl-(m11*cosAz+m12*sinAz)*sinEl;
            dvdaz2=-((m21*cosAz+m22*sinAz)*sinEl);
            dvdazdel=cosEl*(m22*cosAz-m21*sinAz);
            dvdel2=-m23*cosEl-(m21*cosAz+m22*sinAz)*sinEl;
            if(includeW)
                dwdaz2=-((m31*cosAz+m32*sinAz)*sinEl);
                dwdazdel=cosEl*(m32*cosAz-m31*sinAz);
                dwdel2=-m33*cosEl-(m31*cosAz+m32*sinAz)*sinEl;
            end

%If there were no rotations, it would be:
%             dudaz2=-cosAz*sinEl;
%             dudazdel=-cosEl*sinAz;
%             dudel2=-cosAz*sinEl;
%             dvdaz2=-sinAz*sinEl;
%             dvdazdel=cosAz*cosEl;
%             dvdel2=-sinAz*sinEl;
%             if(includeW)
%                 dwdaz2=0;
%                 dwdazdel=0;
%                 dwdel2=-cosEl;
%             end
        case 3
            dudaz2=-cosEl*(m12*cosAz+m11*sinAz);
            dudazdel=(-m11*cosAz+m12*sinAz)*sinEl;
            dudel2=-cosEl*(m12*cosAz+m11*sinAz)-m13*sinEl;
            dvdaz2=-cosEl*(m22*cosAz+m21*sinAz);
            dvdazdel=(-m21*cosAz+m22*sinAz)*sinEl;
            dvdel2=-cosEl*(m22*cosAz+m21*sinAz)-m23*sinEl;
            if(includeW)
                dwdaz2=-cosEl*(m32*cosAz+m31*sinAz);
                dwdazdel=(-m31*cosAz+m32*sinAz)*sinEl;
                dwdel2=-cosEl*(m32*cosAz+m31*sinAz)-m33*sinEl;
            end

%If there were no rotations, it would be:
%             dudaz2=-cosEl*sinAz;
%             dudazdel=-cosAz*sinEl;
%             dudel2=-cosEl*sinAz;
%             dvdaz2=-cosAz*cosEl;
%             dvdazdel=sinAz*sinEl;
%             dvdel2=-cosAz*cosEl;
%             if(includeW)
%                 dwdaz2=0;
%                 dwdazdel=0;
%                 dwdel2=-sinEl;
%             end
        otherwise
            error('Invalid system type specified.')
    end

    H=zeros(2,2,HDim); 

    %Derivatives of u
    H(1,1,1)=dudaz2;
    H(1,2,1)=dudazdel;
    H(2,1,1)=H(1,2,1);
    H(2,2,1)=dudel2;

    %Derivatives of v
    H(1,1,2)=dvdaz2;
    H(1,2,2)=dvdazdel;
    H(2,1,2)=H(1,2,2);
    H(2,2,2)=dvdel2;

    %Derivatives of w
    if(includeW)
        H(1,1,3)=dwdaz2;
        H(1,2,3)=dwdazdel;
        H(2,1,3)=H(1,2,3);
        H(2,2,3)=dwdel2;
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
