function [theta1,theta2,theta3]=quat2Euler3Angles(q,series,rotMatHanded,quatHanded)
%%QUAT2EULER3ANGLES Given a unit-magnitude quaternions representing a
%   rotation (a left handed or a right-handed quaternion), obtain the
%   equivalent set of Euler angles that represent the same rotation. The
%   Euler angle rotations are in order of theta3, theta2, and then theta1
%   about the specified axes. Note that the axes of rotation after theta3
%   are ROTATED axes. The rotations are either right or left-handed, with
%   right-handed rotations being the default. The components of vectors to
%   be rotated are assumed ordered [x;y;z].
%
%INPUTS: q A 4X1 unit-magnitude quaternion representing a rotation. The
%          handedness is specified by quatHanded.
%   series A character string specifying the series of axes about which the
%          three rotations for the Euler angles are taken. All possible
%          combinations of axes without repeating the same axis twice in a
%          row are valid (e.g. rotating about the x axis twice and then the
%          y axis). For example, 'xyz' means rotate theta3 about the z
%          axis, then rotate theta2 about the rotated y axis, then rotate
%          theta1 about the rotated x axis. All possible combinations of
%          values are: 'xzx', 'xyz', 'yxy', 'yzy', 'zyz', 'zxz', 'xzy',
%          'xyz', 'yxz', 'yzx', 'zyx', and 'zxy'.
%    handed The handedness of the Euler rotation angles. If omitted, it is
%           assumed that the rotation is right-handed (the standard).
%           Possible values are:
%           'right' The default if omitted. The rotation is right-handed.
%           'left'  The rotation is left-handed. The rotation angle is
%                   clockwise when one is looking into the rotation axis.
%    handed The handedness of the unit quaternion on the input. If omitted,
%           it is assumed that the quaternion is right-handed (the
%           standard). Possible values are:
%           'right' The default if omitted. The quaternion multiplication
%                   is assumed right-handed (standard).
%           'left'  The quaternion multiplication is assumed left-handed.
%
%OUTPUTS: theta,theta2,theta3 The three angles of rotation given in radians
%              about the three axes specified in series the rotations are
%              taken. The first rotation is theta3, then theta2, then
%              theta1.
%
%Algorithm 1 of [1] is implemented. However, the tests for theta2 being 0
%or pi/2 have been removed, because they produce incorrect results with
%some series; there is no need for those special conditions to be tested.
%
%EXAMPLE:
%Here, the maximum error between an original rotation matrix from Euler
%angles and that obtained by recreateating the Euler angles from a
%quaternion is computed over 1000 Monte Carlo runs. One will see that the
%absolute errors are on the order of what one might expect given finite
%precision errors.
% maxAbsErr=0;
% rotMatHanded='left';
% quatHanded='right';
% numRuns=1e3;
% for curRun=1:numRuns
%     theta1=2*pi*rand();
%     theta2=2*pi*rand();
%     theta3=2*pi*rand();
%     for series=["xzx", "xyz", "yxy", "yzy", "zyz", "zxz", "xzy", "xyz","yxz", "yzx", "zyx", "zxy"]
%         M=Euler3Ang2RotMat(theta1,theta2,theta3,series,rotMatHanded);
%         q=rotMat2Quat(M,quatHanded);
%         [theta1Back,theta2Back,theta3Back]=quat2Euler3Angles(q,series,rotMatHanded,quatHanded);
%         MBack=Euler3Ang2RotMat(theta1Back,theta2Back,theta3Back,series,rotMatHanded);
%         maxAbsErr=max(maxAbsErr,max(abs((M(:)-MBack(:)))));
%     end
% end
% maxAbsErr
%
%REFERENCES:
%[1] E. Bernandes and S. Viollet, "Quaternion to Euler angles conversion:
%    A direct, general and computationally efficient method," PLoS One,
%    vol. 17, no. 11, 10 Nov. 2022.
%
%May 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(quatHanded))
    quatHanded='right';
end

switch(quatHanded)
    case 'right'
    case 'left'
        %Switch to right-handed.
        q=[q(1);-q(2);-q(3);-q(4)];
    otherwise
        error('Invalid quaternions handedness provided.')
end

if(nargin<3||isempty(rotMatHanded))
    rotMatHanded='right';
end

series=char(series);%In case a straing object was passed.
switch(series(1))
    case 'x'
        k=1;
    case 'y'
        k=2;
    case 'z'
        k=3;
    otherwise
        error('Unknown series passed.')
end
switch(series(2))
    case 'x'
        j=1;
    case 'y'
        j=2;
    case 'z'
        j=3;
    otherwise
        error('Unknown series passed.')
end
switch(series(3))
    case 'x'
        i=1;
    case 'y'
        i=2;
    case 'z'
        i=3;
    otherwise
        error('Unknown series passed.')
end

%The paper assumes that q(1)>=0.
if(q(1)<0)
    q=-q;
end

if(i==k)
    notProper=false;
    k=6-i-j;
else
    notProper=true;
end
epsilon=(i-j)*(j-k)*(k-i)/2;

%Adjust i, j, and k for Matlab indexation (from 1, not 0).
i=i+1;
j=j+1;
k=k+1;

if(notProper)
    a=q(1)-q(j);
    b=q(i)+q(k)*epsilon;
    c=q(j)+q(1);
    d=q(k)*epsilon-q(i);
else
    a=q(1);
    b=q(i);
    c=q(j);
    d=q(k)*epsilon;
end
theta2=acos(2*(a^2+b^2)/(a^2+b^2+c^2+d^2)-1);
thetaPlus=atan2(b,a);
thetaMinus=atan2(d,c);
% if(theta2==0)
%     theta1=0;
%     theta3=2*thetaPlus-theta1;
% elseif(theta2==pi/2)
%     theta1=0;
%     theta3=2*thetaMinus+theta1;
% else
    theta1=thetaPlus-thetaMinus;
    theta3=thetaPlus+thetaMinus;
%end

if(notProper)
    theta3=epsilon*theta3;
    theta2=theta2-pi/2;
end

switch(rotMatHanded)
    case 'right'
    case 'left'
        theta1=-theta1;
        theta2=-theta2;
        theta3=-theta3;
    otherwise
        error('Invalid rotation matrix handedness provided.')
end

%In the paper, what is theta 1 and what is theta3 are flipped compared to
%the function Euler3Ang2RotMat.
temp=theta1;
theta1=theta3;
theta3=temp;

theta1=wrapRange(theta1,-pi,pi);
theta3=wrapRange(theta3,-pi,pi);
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
