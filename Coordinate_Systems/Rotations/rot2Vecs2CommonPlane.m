function [vec1Rot,vec2Rot,R]=rot2Vecs2CommonPlane(vec1,vec2,zeroAxis,getAllRots)
%%ROT2VECS2COMMONPLANE Given two 3D vectors, vec1 and vec2, rotate the
%           vectors such that they fall on a common plane. That is, rotate
%           them such that their components for the selected axis are zero.
%           For example, if 'z' is the chosen axis (the default), then the
%           vectors will be rotated such that their components on the z-
%           axis are zero. This can be useful when dealing with coordinate
%           systems where a singularity exists in a particular direction.
%           For example, by zeroing the z-components of the vectors, one
%           avoids the singularities in the traditional spherical
%           coordinate system at elevations of +/-pi/2.
%
%INPUTS: vec1, vec2 The two 3X1 vectors that should be commonly aligned.
%          zeroAxis This specified which axis should be zeroed by the
%                   rotation. Possible values are 'x', 'y', and 'z'. If
%                   omitted or an empty matrix is passed, the default is
%                   'z'.
%        getAllRots The rotation matrix is formed by successive rotations
%                   about the two axes other than the one being zeroed.
%                   This means that 4 solutions are typically present. If
%                   this parameter is true, the vectors and rotation
%                   matrices for all 4 solutions are returned. Otherwise,
%                   just one solution is chosen. The default if omitted or
%                   an empty matrix is passed is false. Note that in
%                   some geometries, the 4 solutions will not always be
%                   unique.
%
%OUTPUTS: vec1Rot,vec2Rot The two 3X1 vectors rotated such that the desired
%                         axis is zero for both of them. If
%                         getAllRots=true, then these are 3X4 matrices of
%                         all solutions.
%                       R The 3X3 rotation matrix. If getAllRots is true,
%                         then this is a 3X3X4 set of all 4 rotation
%                         matrices.
%
%The algorithm is discussed with respect to a rotation about the z-axis.
%However, rotations about the other axes are derived in the same manner if
%one just changes the description.
%
%The goal is to find a rotation matrix the simultaneously zeros the z
%components of vec1 and vec2. An infinite number of solutions exist.
%However, we limit the number of possibilities by requiring the rotation
%matrix be the product of a rotation about the y axis followed by a
%rotation about the x axis. That is, the rotated and unrotated vectors are
%related as vRot=Rx*Ry*vec. This ends up limiting the number of rotation
%solutions to 4. The expressions for the z components of the rotated
%vectors form two equations, because we want them to be zero. Specifically,
%x3*cos(thetax)*cos(thetay)+x2*sin(thetax)-x1*cos(thetax)*sin(thetay)=0
%y3*cos(thetax)*cos(thetay)+y2*sin(thetax)-y1*cos(thetax)*sin(thetay)=0
%The four solutions (unique within factors of 2*pi) for thetax and thetay
%come from explcitly solving for thetax and thetay and use the atan2
%function. Note that one must specially address the case where vec1=vec2 or
%vec1=-vec2.
%
%EXAMPLE:
%Here, we rotate two vectors, verify that the z-components of the jointly
%rotated vectors are zero, and then we verify that the inverse rotation
%produces the original vectors.
% vec1=[-7.8;6.7;-6.9];
% vec2=[-0.1;0.1;-11.9];
% [vec1Rot,vec2Rot,R]=rot2Vecs2CommonPlane(vec1,vec2);
% 
% vec1Rot(3)
% vec2Rot(3)
% %One will see that the third element of the rotated vectors is zero
% %within finite precision limits.
% max(abs(R'*vec1Rot-vec1))
% max(abs(R'*vec2Rot-vec2))
%One will see that the inverse rotation brings back the original vectors,
%within finite precision limitations.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(zeroAxis))
    zeroAxis='z';
end

if(nargin<4||isempty(getAllRots))
    getAllRots=false;
end

x=vec1/norm(vec1);
y=vec2/norm(vec2);

x1=x(1);
x2=x(2);
x3=x(3);
y1=y(1);
y2=y(2);
y3=y(3);

switch(zeroAxis)
    case 'z'%Rotation xy to zero z.
        %The max operation in num1 and denom1 is just in case finite
        %precision errors make the argument of the square root negative.
        num1=sqrt(max(0,(x1^2+x3^2)*y2^2-2*x2*y2*(x1*y1+x3*y3)+x2^2*(y1^2+y3^2)));
        %Note that denom1 can be replaced by 1, because it appears in
        %both the numerator and denominator of the arctangent function and
        %nowhere else. On the other hand, num1 cannot be replaced by 1.
        %denom1=1;%sqrt(max(0,(x2^2+x3^2)*y1^2+x1^2*(y2^2+y3^2)+x2^2*y3^2+x3^2*y2^2-2*x2*x3*y2*y3-2*x1*y1*(x2*y2+x3*y3)));
        num2=x1*y3-x3*y1;
        num3=x1*y2-x2*y1;
        num4=x3*y2-x2*y3;

        rotSeries='xy';
    case 'y'%Rotation xz to zero y.
        num1=sqrt(max(0,x3^2*(y1^2+y2^2)-2*x1*x3*y1*y3-2*x2*x3*y2*y3+(x1^2+x2^2)*y3^2));
        %denom1=1;%sqrt(max(0,x3^2*(y1^2+y2^2)-2*x1*x3*y1*y3-2*x2*y2*(x1*y1+x3*y3)+x2^2*(y1^2+y3^2)+x1^2*(y2^2+y3^2)));
        num2=x2*y1-x1*y2;
        num3=x1*y3-x3*y1;
        num4=x3*y2-x2*y3;
        rotSeries='xz';
    case 'x'%Rotation yz to zero x.
        num1=sqrt(max(0,x3^2*(y1^2+y2^2)-2*x3*(x1*y1+x2*y2)*y3+(x1^2+x2^2)*y3^2));
        %denom1=1;%sqrt(max(0,x3^2*(y1^2+y2^2)-2*x1*x3*y1*y3-2*x2*y2*(x1*y1+x3*y3)+x2^2*(y1^2+y3^2)+x1^2*(y2^2+y3^2)));
        num2=x2*y1-x1*y2;
        num3=x2*y3-x3*y2;
        num4=x1*y3-x3*y1;
        rotSeries='yz';
    otherwise
        error('Unknown axis specified.')
end

%The following commented out lines are if we used the full ratio. However,
%due to the common (positive) denominators in the atan2 function, the
%denominators can be omitted.
% val1=-num1/denom1;
% val2=num2/denom1;
% val3=num3/num1;
% val4=num4/num1;

val1=-num1;
val2=num2;
val3=num3;
val4=num4;

if(getAllRots)
    theta12=zeros(2,4);
    R=zeros(3,3,4);
    vec1Rot=zeros(3,4);
    vec2Rot=zeros(3,4);
else
    theta12=zeros(2,1);
    R=zeros(3,3,1);
    vec1Rot=zeros(3,1);
    vec2Rot=zeros(3,1);
end

if(abs(num1)<eps()&&abs(num2)<=eps()&&abs(num3)<=eps()&&abs(num4)<=eps())
    %Special case: If vec1==vec2 or vec1==-vec2.

    %Rotate the vectors to be along one axis.
    R(:,:,1)=rotAxis2Vec(vec1,rotSeries(1))';
    vec1Rot(:,1)=R(:,:,1)*vec1;
    vec2Rot(:,1)=R(:,:,1)*vec2;

    %There are infinitely many solutions, but we just choose the ones
    %aligning with the coordinate axes.
    if(getAllRots)
        R(:,:,2)=-R(:,:,1);
        vec1Rot(:,2)=-vec1Rot(:,1);
        vec2Rot(:,2)=-vec2Rot(:,1);
        
        %Rotate the vectors to be along the other axis.
        R(:,:,3)=rotAxis2Vec(vec1,rotSeries(2))';
        vec1Rot(:,3)=R(:,:,3)*vec1;
        vec2Rot(:,3)=R(:,:,3)*vec2;

        R(:,:,4)=-R(:,:,3);
        vec1Rot(:,4)=-vec1Rot(:,3);
        vec2Rot(:,4)=-vec2Rot(:,3);
    end

    return;
end

theta12(1,1)=atan2(val2,val1);
theta12(2,1)=atan2(val4,val3);

if(getAllRots)
    theta12(1,2)=atan2(-val2,val1);
    theta12(2,2)=atan2(-val4,-val3);

    theta12(1,3)=atan2(val2,-val1);
    theta12(2,3)=atan2(-val4,-val3);

    theta12(1,4)=atan2(-val2,-val1);
    theta12(2,4)=atan2(val4,val3);
end

%Deal with finite precision issues resuling in NaNs.
theta12(~isfinite(theta12))=0;

R(:,:,1)=Euler2Ang2RotMat(theta12(1,1),theta12(2,1),rotSeries);

vec1Rot(:,1)=R(:,:,1)*vec1;
vec2Rot(:,1)=R(:,:,1)*vec2;

if(getAllRots)
    for k=2:4
        R(:,:,k)=Euler2Ang2RotMat(theta12(1,k),theta12(2,k),rotSeries);
        vec1Rot(:,k)=R(:,:,k)*vec1;
        vec2Rot(:,k)=R(:,:,k)*vec2;
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
