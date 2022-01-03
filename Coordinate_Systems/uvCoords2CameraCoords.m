function [zCam,isBehind]=uvCoords2CameraCoords(uv,A,M)
%%UVCOORDS2CAMERACOORDS Convert a set of u-v direction cosines (or a full 
%          3D unit vector) direction vector into a measurement in the local
%          coordinate system of a perspective camera (distances on a CCD).
%          The model used is described in [1] and in the comments to
%          the inverse to this function, which is cameraCoords2UVCoords.
%
%INPUTS: uv A 2XnumPts set of uv direction cosines or a 3XnumPts set of
%           unit direction vectors. If M is provided, and only a 2XnumPts
%           set of vectors is given, the missing third component is assumed
%           to be positive. If after rotation into the local coordinate
%           system, the third coordinate is negative (which is not allowed
%           for targets in front of the camera), it is invalid and is
%           marked as being behind the sensor in isBehind.
%         A A 3X3 matrix, as described in the comments to
%           cameraCoords2UVCoords. The third row must be
%           [0,0,1]
%         M A 3X3XN rotation matrix to go from the alignment of the
%           global coordinate system to that at the receiver. The z-axis of
%           the local coordinate system of the receiver is the pointing
%           direction of the receiver. If omitted or an empty matrix is
%           passed, then it is assumed that the local coordinate system is
%           aligned with the global and M=eye(3) --the identity matrix is
%           used. 
%
%OUTPUTS: zCam The 2XnumPts set of direction given as distances on the CCD
%              of a camera.
%     isBehind A 1XnumPts set of boolean values indicating whether each
%              direction provided is behind the receiver and thus is not
%              valid.
%
%See the function cameraCoords2UVCoords for more details on the model and
%for an example conversion.
%
%REFERENCES:
%[1] J. Kannala, J. Heikkilä, and S. S. Brandt, "Geometric camera
%    calibration," in Wiley Encyclopedia of Computer Science and
%    Engineering, B. W. Wah, Ed., 2007, vol. 1.
%
%November 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(M))
    M=eye(3,3);
end

if(~all(A(3,1:2)==[0,0]))
    error('The third row of A has the wrong format.')
end

a33=A(3,3);
a11=A(1,1)/a33;
a12=A(1,2)/a33;
a13=A(1,3)/a33;
a21=A(2,1)/a33;
a22=A(2,2)/a33;
a23=A(2,3)/a33;

numMeas=size(uv,2);
isBehind=false(1,numMeas);
if(nargin>2&&~isempty(M))
    %Rotate into local coordinates.
    w=sqrt(1-uv(1,:).^2-uv(2,:).^2);
    uv=M*[uv;w];
    isBehind=uv(3,:)<0;
end

u=uv(1,:);
v=uv(2,:);

denom=sqrt(1-u.^2-v.^2);

zCam=[a13+(a11*u+a12*v)./denom;
      a23+(a21*u+a22*v)./denom];

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
