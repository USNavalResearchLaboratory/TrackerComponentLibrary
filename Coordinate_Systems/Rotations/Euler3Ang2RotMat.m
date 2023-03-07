function M=Euler3Ang2RotMat(theta1,theta2,theta3,series,handed)
%%EULER3ANG2ROTMAT  Get a rotation matrix for a series of 3 Euler angles.
%              The axes of rotation are given by series. The rotations are
%              in order of theta3, theta2, and then theta1 about the
%              specified axes. Note that the axes of rotation after theta3
%              are ROTATED axes. The rotations are either right or
%              left-handed, with right-handed rotations being the default.
%              The components of vectors to be rotated are assumed ordered
%              [x;y;z].
%
%INPUTS: theta,theta2,theta3 The three angles of rotation given in radians
%              about the three axes specified in series the rotations are
%              taken. The first rotation is theta3, then theta2, then
%              theta1.
%       series A character string specifying the series of axes about which
%              the three rotations are taken. All possible combinations of
%              axes without repeating an axis are valid. For example, 'xyz'
%              means rotate theta3 about the z axis, then rotate theta2
%              about the rotated y axis, then rotate theta1 about the
%              rotated x axis. All possible combinations of values are:
%              'xzx', 'xyz', 'yxy', 'yzy', 'zyz', 'zxz', 'xzy', 'xyz',
%              'yxz', 'yzx', 'zyx', and 'zxy'.
%       handed The handedness of the rotation angle. If omitted, it is
%              assumed that the rotation is right-handed (the standard).
%              Possible values are:
%              'right' The default if omitted. The rotation is right-
%                      handed.
%              'left'  The rotation is left-handed. The rotation angle is
%                      clockwise when one is looking into the rotation
%                      axis.
%
%OUTPUTS: M The 3X3 rotation matrix such that M*v rotates the 3X1 vector v
%           according to the given rotation angles and axes.
%
%Euler angles are discussed in [1].
%
%EXAMPLE:
%Here, we rotate by some arbitrary angles and one can see that the
%rotations from this functions are equal to subsequently appending the
%rotations.
% theta1=2*pi*rand();
% theta2=2*pi*rand();
% theta3=2*pi*rand();
% max(max(abs(Euler3Ang2RotMat(theta1,theta2,theta3,'xzx')-Euler1Ang2RotMat(theta1,'x')*Euler1Ang2RotMat(theta2,'z')*Euler1Ang2RotMat(theta3,'x'))))
% max(max(abs(Euler3Ang2RotMat(theta1,theta2,theta3,'xyx')-Euler1Ang2RotMat(theta1,'x')*Euler1Ang2RotMat(theta2,'y')*Euler1Ang2RotMat(theta3,'x'))))
% max(max(abs(Euler3Ang2RotMat(theta1,theta2,theta3,'yxy')-Euler1Ang2RotMat(theta1,'y')*Euler1Ang2RotMat(theta2,'x')*Euler1Ang2RotMat(theta3,'y'))))
% max(max(abs(Euler3Ang2RotMat(theta1,theta2,theta3,'yzy')-Euler1Ang2RotMat(theta1,'y')*Euler1Ang2RotMat(theta2,'z')*Euler1Ang2RotMat(theta3,'y'))))
% max(max(abs(Euler3Ang2RotMat(theta1,theta2,theta3,'zyz')-Euler1Ang2RotMat(theta1,'z')*Euler1Ang2RotMat(theta2,'y')*Euler1Ang2RotMat(theta3,'z'))))
% max(max(abs(Euler3Ang2RotMat(theta1,theta2,theta3,'zxz')-Euler1Ang2RotMat(theta1,'z')*Euler1Ang2RotMat(theta2,'x')*Euler1Ang2RotMat(theta3,'z'))))
% max(max(abs(Euler3Ang2RotMat(theta1,theta2,theta3,'xzy')-Euler1Ang2RotMat(theta1,'x')*Euler1Ang2RotMat(theta2,'z')*Euler1Ang2RotMat(theta3,'y'))))
% max(max(abs(Euler3Ang2RotMat(theta1,theta2,theta3,'xyz')-Euler1Ang2RotMat(theta1,'x')*Euler1Ang2RotMat(theta2,'y')*Euler1Ang2RotMat(theta3,'z'))))
% max(max(abs(Euler3Ang2RotMat(theta1,theta2,theta3,'yxz')-Euler1Ang2RotMat(theta1,'y')*Euler1Ang2RotMat(theta2,'x')*Euler1Ang2RotMat(theta3,'z'))))
% max(max(abs(Euler3Ang2RotMat(theta1,theta2,theta3,'yzx')-Euler1Ang2RotMat(theta1,'y')*Euler1Ang2RotMat(theta2,'z')*Euler1Ang2RotMat(theta3,'x'))))
% max(max(abs(Euler3Ang2RotMat(theta1,theta2,theta3,'zyx')-Euler1Ang2RotMat(theta1,'z')*Euler1Ang2RotMat(theta2,'y')*Euler1Ang2RotMat(theta3,'x'))))
% max(max(abs(Euler3Ang2RotMat(theta1,theta2,theta3,'zxy')-Euler1Ang2RotMat(theta1,'z')*Euler1Ang2RotMat(theta2,'x')*Euler1Ang2RotMat(theta3,'y'))))
%All of the above differences will be zero or around eps(), indicating
%agreement within finite precision limits.
%
%REFERENCES:
%[1] M. D. Shuster, "A survey of attitude representations," The Journal of
%    Astronautical Sciences, vol. 41, no. 4, pp. 439-517, Oct. - Dec. 1993.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(handed))
    handed='right';
end

switch(handed)
    case 'right'
    case 'left'
        theta1=-theta1;
        theta2=-theta2;
        theta3=-theta3;
    otherwise
        error('Invalid handedness provided.')
end

c1=cos(theta1);
s1=sin(theta1);
c2=cos(theta2);
s2=sin(theta2);
c3=cos(theta3);
s3=sin(theta3);

%All of the formuae below are for right-handed rotations. Flipping the sign
%of the angles makes them left-handed rotations.
switch(series)
    %First the Euler angle combinations.
    case 'xzx'
        M=[   c2,         -c3*s2,           s2*s3;
           c1*s2, c1*c2*c3-s1*s3, -c3*s1-c1*c2*s3;
           s1*s2, c2*c3*s1+c1*s3,  c1*c3-c2*s1*s3];
    case 'xyx'
        M=[    c2,          s2*s3,           c3*s2;
            s1*s2, c1*c3-c2*s1*s3, -c2*c3*s1-c1*s3;
           -c1*s2, c3*s1+c1*c2*s3,  c1*c2*c3-s1*s3];
    case 'yxy'
        M=[ c1*c3-c2*s1*s3, s1*s2, c2*c3*s1+c1*s3;
                     s2*s3,    c2,         -c3*s2;
           -c3*s1-c1*c2*s3, c1*s2, c1*c2*c3-s1*s3];
    case 'yzy'
        M=[ c1*c2*c3-s1*s3, -c1*s2, c3*s1+c1*c2*s3;
                 c3*s2,         c2,          s2*s3;
           -c2*c3*s1-c1*s3,  s1*s2, c1*c3-c2*s1*s3];
    case 'zyz'
        M=[c1*c2*c3-s1*s3, -c3*s1-c1*c2*s3, c1*s2;
           c2*c3*s1+c1*s3,  c1*c3-c2*s1*s3, s1*s2;
                   -c3*s2,           s2*s3,    c2];
    case 'zxz'
        M=[c1*c3-c2*s1*s3, -c2*c3*s1-c1*s3, s1*s2;
           c3*s1+c1*c2*s3, c1*c2*c3-s1*s3, -c1*s2;
           s2*s3, c3*s2, c2];
    %Next, the Tait-Bryan-Cardan angle combinations.
    case 'xzy'
        M=[c2*c3,           -s2,    c2*s3;
           c1*c3*s2+s1*s3,	c1*c2,	-c3*s1+c1*s2*s3;
           c3*s1*s2-c1*s3,	c2*s1,	c1*c3+s1*s2*s3];
    case 'xyz'
        M=[c2*c3,           -c2*s3,         s2;
           c3*s1*s2+c1*s3,	c1*c3-s1*s2*s3,	-c2*s1;
           -c1*c3*s2+s1*s3,	c3*s1+c1*s2*s3,	c1*c2];
    case 'yxz'
        M=[c1*c3+s1*s2*s3,	c3*s1*s2-c1*s3,	c2*s1;
           c2*s3,           c2*c3,          -s2;
          -c3*s1+c1*s2*s3,	c1*c3*s2+s1*s3,	c1*c2];
    case 'yzx'
        M=[c1*c2,	-c1*c3*s2+s1*s3,    c3*s1+c1*s2*s3;
           s2,      c2*c3,              -c2*s3;
           -c2*s1,	c3*s1*s2+c1*s3,      c1*c3-s1*s2*s3];
    case 'zyx'
        M=[c1*c2,   -c3*s1+c1*s2*s3,    c1*c3*s2+s1*s3;
           c2*s1,   c1*c3+s1*s2*s3,     c3*s1*s2-c1*s3;
           -s2,     c2*s3,              c2*c3];
    case 'zxy'
        M=[c1*c3-s1*s2*s3,	-c2*s1, c3*s1*s2+c1*s3;
           c3*s1+c1*s2*s3,	c1*c2,  -c1*c3*s2+s1*s3;
           -c2*s3,          s2,     c2*c3];
    otherwise
        error('Invalid rotation series provided.')
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
