function dirVecs=cameraCoords2UVCoords(zCam,A,M,includeW)
%%CAMERACOORDS2UVCOORDS Consider a perfect pinhole camera with the
%   origin at the focus and the z axis being the pointing direction. The
%   coordinates of an [x;y;z] point in that system projected on the imaging
%   plane of the camera are
%   [xImag;yImage]=-(f/z)*[x;y] 
%   (More calibrated perspective cameras will use the A matrix as described
%   in [1] and below) Considering a point zCam=[xImage;yImage] measured by
%   the camera (probably measured as pixels on a charged coupled device
%   [CCD] and converted into distances across th CCD), this function 
%   converts that point into a uv direction cosines or a full 3D unit
%   vector. It can be givn in local camera coordinates or rotated into
%   global coordinates.
%
%INPUTS: zCam A 2XnumPts set of [x;y] coordinates on the CCD of the camera
%             measuring the directions. This is in distance units, not
%             pixels.
%           A A 3X3 matrix, as described below. The third row must be
%             [0,0,1].
%           M A 3X3 rotation matrix to go from the alignment of the
%             global coordinate system to that at the receiver. The z-axis
%             of the local coordinate system of the receiver is the
%             pointing direction of the receiver. If omitted or an empty
%             matrix is passed, then it is assumed that the local
%             coordinate system is aligned with the global and M=eye(3)
%             --the identity matrix is used. 
%    includeW If only the first two elements of a unit direction vctor
%             should be returned, then this is false. Otherwise, a full 3X1
%             unit direction vector is retuned. The default if omitted or
%             an empty matrix is passed is true.
%
%OUTPUT: dirVecs A 3XnumPoints set of the directions converted to 3D unit
%           vectors. If M is omitted, they are in the local coordinate
%           system of the sensor. Otherwise, they have been rotated into
%           the global coordnate system. If includeW is false, then only
%           the first two elements are returned. In the local coordinate
%           systm of the receiver, the third component is always positive
%           (in front of the receiver).
%
%The A matrix is a calibration matrix that allows for more general
%transformations. This type of perspective camera model is described in
%[1]. The A matrix is often estimated being of the form
%A=[alpha,gamma,x0;
%       0, beta,y0;
%       0,    0, 1];
%whereas here, the more general solution assuming A has the form
%A=[a11,a12,a13;
%   a21,a22,a23;
%     0,  0,a33];
%Given a vector xVec=[x;y;z] in the local coordinate system of the camera
%(the z axis is is the pointing direction and the origin is the focus), the
%2D coordinates measured on the imaging plane are obtained as
%xInt=A*xVec;                (1)
%zCam=xInt(1:2)/xInt(3);     (2)
%In this instance, xInt is essentially a set of homogeneous coordinates and
%the final step converts to actual coordinates on the imaging plane.
%The A matrix corresponds to that of a perfect pinhole camera with a focus
%of f is:
%A=[-f,0,0;
%   0,-f,0;
%   0, 0,1];
%To solve the problem, Equations 1 and 2 and th equations for local u-v
%coordinates:
%u=x/sqrt(x^2+y^2+z^2)
%v=y/sqrt(x^2+y^2+z^2)
%Are used to eliminate x,y, and z to obtain the expressions implemented
%here for u and v in terms of the zCam. If a 3D vector is desired, then
%w=sqrt(1-u^2-v^2) and one uses the vector [u;v;w]. If M is provided, that
%vector is rotated into the global coordinate system.
%
%The inverse of this function is uvCoords2CameraCoords.
%
%EXAMPLE:
%This demonstrates that the uvCoords2CameraCoords function is consistent
%with the outputs of the cameraCoords2UVCoords function. A random A matrix
%having the correct format is generated and is used to create the
%conversion. After converting uv values into camera coordinates and back,
%the error of the back-converted values is consistent with what one would
%expect due to finite precision limitations.
% uv=[-0.2,0.65,0,0.3;
%      0.5,0.7,0,-0.1];
% A=[rand(2,3);
%    0,0,1];
% uvBack=cameraCoords2UVCoords(uvCoords2CameraCoords(uv,A),A);
% AbsErr=max(max(abs(uvBack-uv)))
%
%REFERENCES:
%[1] J. Kannala, J. Heikkilä, and S. S. Brandt, "Geometric camera
%    calibration," in Wiley Encyclopedia of Computer Science and
%    Engineering, B. W. Wah, Ed., 2007, vol. 1.
%
%November 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(~all(A(3,1:2)==[0,0]))
    error('The third row of A has the wrong format.')
end

if(nargin<4||isempty(includeW))
    includeW=true;
end

if(nargin<3||isempty(M))
    M=eye(3,3);
end

a33=A(3,3);
a11=A(1,1)/a33;
a12=A(1,2)/a33;
a13=A(1,3)/a33;
a21=A(2,1)/a33;
a22=A(2,2)/a33;
a23=A(2,3)/a33;

xC=zCam(1,:);
yC=zCam(2,:);

denom=sqrt(a13^2*(a21^2+a22^2)-2*a11*a13*a21*a23-2*a12*a22*(a11*a21+a13*a23)+a12^2*(a21^2+a23^2)+a11^2*(a22^2+a23^2)+(-2*a13*(a21^2+a22^2)+2*(a11*a21+a12*a22)*a23)*xC+(a21^2+a22^2)*xC.^2+(2*a13*(a11*a21+a12*a22)-2*(a11^2+a12^2)*a23)*yC-2*(a11*a21+a12*a22).*xC.*yC+(a11^2+a12^2)*yC.^2);
signVal=sign(a11*a22-a12*a21);

u=signVal*(-a13*a22+a22*xC+a12*(a23-yC))./denom;
v=signVal*(a13*a21-a21*xC-a11*(a23-yC))./denom;

%Rotate into global coordinates, if required.
if(nargin>2&&~isempty(M))
    w=sqrt(1-u.^2-v.^2);
    uvw=M'*[u;v;w];
    if(includeW)
        dirVecs=uvw;
    else
        dirVecs=uvw(1:2,:);
    end
else
    if(includeW)
        w=sqrt(1-u.^2-v.^2);
        dirVecs=[u;v;w];
    else
        dirVecs=[u;v];
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
