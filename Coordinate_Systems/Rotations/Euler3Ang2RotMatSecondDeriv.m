function MDDot=Euler3Ang2RotMatSecondDeriv(thetaState1,thetaState2,thetaState3,series,handed)
%%EULER3ANG2ROTMATSECONDDERIV This function provides the second time
%              derivative of a rotation matrix for a series of 3 Euler
%              angles, assuming that the angles vary with time. This is the
%              derivative of the output of the Euler3Ang2RotMatDeriv
%              function. The axes of rotation are given by series. The
%              rotations are in order of theta3, theta2, and then theta1
%              about the specified axes. Note that the axes of rotation
%              after theta3 are ROTATED axes. The rotations are either
%              right or left-handed, with right-handed rotations being the
%              default. The components of vectors to be rotated are assumed
%              ordered [x;y;z].
%
%INPUTS: thetaState1, thetaState2, thetaState3 These are the 3X1 states,
%              consisting of [angle; 1st angle derivative; second angle
%              derivative] with the angle given in radians for each of the
%              three rotations. The first rotation is theta3, then theta2,
%              then theta1.
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
%OUTPUTS: MDDot The 3X3 second derivative of the rotation matrix.
%
%This function hold analytic derivatives of the expressions given in
%Euler3Ang2RotMatDeriv.
%
%EXAMPLE:
%Angles, angle rates, and angle accelerations.
% thetaState1=[2*pi*rand(1);
%              randn(2,1)];
% thetaState2=[2*pi*rand(1);
%              randn(2,1)];
% thetaState3=[2*pi*rand(1);
%              randn(2,1)];
% deltaTEps=1e-8;%For finite differencing for comparison.
% FTheta=FPolyKal(deltaTEps,3,2);
% thetaState1Pred=FTheta*thetaState1;
% thetaState2Pred=FTheta*thetaState2;
% thetaState3Pred=FTheta*thetaState3;
% maxRelErr=0;
% for handed=["right","left"]
%     for series=["xzx", "xyz", "yxy", "yzy","zyz", "zxz", "xzy", "xyz","yxz", "yzx", "zyx", "zxy"]
%         MDDot=Euler3Ang2RotMatSecondDeriv(thetaState1,thetaState2,thetaState3,series,handed);
%         MDDotNumDiff=(Euler3Ang2RotMatDeriv(thetaState1Pred,thetaState2Pred,thetaState3Pred,series,handed)-Euler3Ang2RotMatDeriv(thetaState1,thetaState2,thetaState3,series,handed))/deltaTEps;
%         maxRelErr=max(maxRelErr,max(abs((MDDot(:)-MDDotNumDiff(:))./MDDotNumDiff(:))));
%     end
% end
% maxRelErr
%
%January 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(handed))
    handed='right';
end

theta1=thetaState1(1);
theta1Dot=thetaState1(2);
theta1DDot=thetaState1(3);
theta2=thetaState2(1);
theta2Dot=thetaState2(2);
theta2DDot=thetaState2(3);
theta3=thetaState3(1);
theta3Dot=thetaState3(2);
theta3DDot=thetaState3(3);

co1=cos(theta1);
si1=sin(theta1);
co2=cos(theta2);
si2=sin(theta2);
co3=cos(theta3);
si3=sin(theta3);

switch(series)
    %First the Euler angle combinations.
    case 'xzx'
        switch(handed)
            case 'right'
                MDDot=[-si2*theta2DDot-co2*theta2Dot^2,si2*si3*theta3DDot+co2*(-co3*theta2DDot+2*si3*theta2Dot*theta3Dot)+co3*si2*(theta2Dot^2+theta3Dot^2),co3*si2*theta3DDot+co2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-si2*si3*(theta2Dot^2+theta3Dot^2);
                        co1*co2*theta2DDot-si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-co1*si2*(theta1Dot^2+theta2Dot^2),-co1*(si3*theta1DDot+co3*si2*theta2DDot+2*co3*theta1Dot*theta3Dot-2*si2*si3*theta2Dot*theta3Dot)+si1*(si3*theta1Dot^2+2*co3*si2*theta1Dot*theta2Dot-co3*theta3DDot+si3*theta3Dot^2)-co2*(co3*si1*theta1DDot+si3*(co1*theta3DDot-2*si1*theta1Dot*theta3Dot)+co1*co3*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),-co1*co3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+si1*(co3*theta1Dot^2+si3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*co3*theta1Dot*theta3Dot+co3*theta3Dot^2)+co1*si3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2));
                        co2*si1*theta2DDot+co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-si1*si2*(theta1Dot^2+theta2Dot^2),co1*(co3*(-2*si2*theta1Dot*theta2Dot+theta3DDot)+co2*(co3*theta1DDot-2*si3*theta1Dot*theta3Dot)-si3*(theta1Dot^2+theta3Dot^2))-si1*(si3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),-co3*(si1*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*co2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)))];
            case 'left'
                MDDot=[-si2*theta2DDot-co2*theta2Dot^2,co2*(co3*theta2DDot-2*si3*theta2Dot*theta3Dot)-si2*(si3*theta3DDot+co3*(theta2Dot^2+theta3Dot^2)),co3*si2*theta3DDot+co2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-si2*si3*(theta2Dot^2+theta3Dot^2);
                        si1*si2*theta1DDot+co2*(-co1*theta2DDot+2*si1*theta1Dot*theta2Dot)+co1*si2*(theta1Dot^2+theta2Dot^2),-co1*(si3*theta1DDot+co3*si2*theta2DDot+2*co3*theta1Dot*theta3Dot-2*si2*si3*theta2Dot*theta3Dot)+si1*(si3*theta1Dot^2+2*co3*si2*theta1Dot*theta2Dot-co3*theta3DDot+si3*theta3Dot^2)-co2*(co3*si1*theta1DDot+si3*(co1*theta3DDot-2*si1*theta1Dot*theta3Dot)+co1*co3*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),co1*co3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)-si1*(co3*theta1Dot^2+si3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*co3*theta1Dot*theta3Dot+co3*theta3Dot^2)-co1*si3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2));
                        co2*si1*theta2DDot+co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-si1*si2*(theta1Dot^2+theta2Dot^2),si1*si3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(si3*theta1Dot^2-co3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*si3*theta1Dot*theta3Dot+si3*theta3Dot^2)+co3*si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),-co3*(si1*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*co2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)))];
           otherwise
                error('Invalid handedness provided.')
        end
    case 'xyx'
        switch(handed)
            case 'right'
                MDDot=[-si2*theta2DDot-co2*theta2Dot^2,co3*si2*theta3DDot+co2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-si2*si3*(theta2Dot^2+theta3Dot^2),co2*(co3*theta2DDot-2*si3*theta2Dot*theta3Dot)-si2*(si3*theta3DDot+co3*(theta2Dot^2+theta3Dot^2));
                        co2*si1*theta2DDot+co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-si1*si2*(theta1Dot^2+theta2Dot^2),-co3*(si1*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*co2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),si1*si3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(si3*theta1Dot^2-co3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*si3*theta1Dot*theta3Dot+si3*theta3Dot^2)+co3*si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2));
                        si1*si2*theta1DDot+co2*(-co1*theta2DDot+2*si1*theta1Dot*theta2Dot)+co1*si2*(theta1Dot^2+theta2Dot^2),co1*co3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)-si1*(co3*theta1Dot^2+si3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*co3*theta1Dot*theta3Dot+co3*theta3Dot^2)-co1*si3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),-co1*(si3*theta1DDot+co3*si2*theta2DDot+2*co3*theta1Dot*theta3Dot-2*si2*si3*theta2Dot*theta3Dot)+si1*(si3*theta1Dot^2+2*co3*si2*theta1Dot*theta2Dot-co3*theta3DDot+si3*theta3Dot^2)-co2*(co3*si1*theta1DDot+si3*(co1*theta3DDot-2*si1*theta1Dot*theta3Dot)+co1*co3*(theta1Dot^2+theta2Dot^2+theta3Dot^2))];
            case 'left'
                MDDot=[-si2*theta2DDot-co2*theta2Dot^2,co3*si2*theta3DDot+co2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-si2*si3*(theta2Dot^2+theta3Dot^2),si2*si3*theta3DDot+co2*(-co3*theta2DDot+2*si3*theta2Dot*theta3Dot)+co3*si2*(theta2Dot^2+theta3Dot^2);
                        co2*si1*theta2DDot+co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-si1*si2*(theta1Dot^2+theta2Dot^2),-co3*(si1*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*co2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),co1*(co3*(-2*si2*theta1Dot*theta2Dot+theta3DDot)+co2*(co3*theta1DDot-2*si3*theta1Dot*theta3Dot)-si3*(theta1Dot^2+theta3Dot^2))-si1*(si3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)));
                        co1*co2*theta2DDot-si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-co1*si2*(theta1Dot^2+theta2Dot^2),-co1*co3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+si1*(co3*theta1Dot^2+si3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*co3*theta1Dot*theta3Dot+co3*theta3Dot^2)+co1*si3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),-co1*(si3*theta1DDot+co3*si2*theta2DDot+2*co3*theta1Dot*theta3Dot-2*si2*si3*theta2Dot*theta3Dot)+si1*(si3*theta1Dot^2+2*co3*si2*theta1Dot*theta2Dot-co3*theta3DDot+si3*theta3Dot^2)-co2*(co3*si1*theta1DDot+si3*(co1*theta3DDot-2*si1*theta1Dot*theta3Dot)+co1*co3*(theta1Dot^2+theta2Dot^2+theta3Dot^2))];
            otherwise
                error('Invalid handedness provided.')
        end
    case 'yxy'
        switch(handed)
            case 'right'
                MDDot=[-co3*(si1*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*co2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),co2*si1*theta2DDot+co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-si1*si2*(theta1Dot^2+theta2Dot^2),co1*(co3*(-2*si2*theta1Dot*theta2Dot+theta3DDot)+co2*(co3*theta1DDot-2*si3*theta1Dot*theta3Dot)-si3*(theta1Dot^2+theta3Dot^2))-si1*(si3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)));
                        co3*si2*theta3DDot+co2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-si2*si3*(theta2Dot^2+theta3Dot^2),-si2*theta2DDot-co2*theta2Dot^2,si2*si3*theta3DDot+co2*(-co3*theta2DDot+2*si3*theta2Dot*theta3Dot)+co3*si2*(theta2Dot^2+theta3Dot^2);
                       -co1*co3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+si1*(co3*theta1Dot^2+si3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*co3*theta1Dot*theta3Dot+co3*theta3Dot^2)+co1*si3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),co1*co2*theta2DDot-si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-co1*si2*(theta1Dot^2+theta2Dot^2),-co1*(si3*theta1DDot+co3*si2*theta2DDot+2*co3*theta1Dot*theta3Dot-2*si2*si3*theta2Dot*theta3Dot)+si1*(si3*theta1Dot^2+2*co3*si2*theta1Dot*theta2Dot-co3*theta3DDot+si3*theta3Dot^2)-co2*(co3*si1*theta1DDot+si3*(co1*theta3DDot-2*si1*theta1Dot*theta3Dot)+co1*co3*(theta1Dot^2+theta2Dot^2+theta3Dot^2))];
            case 'left'
                MDDot=[-co3*(si1*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*co2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),co2*si1*theta2DDot+co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-si1*si2*(theta1Dot^2+theta2Dot^2),si1*si3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(si3*theta1Dot^2-co3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*si3*theta1Dot*theta3Dot+si3*theta3Dot^2)+co3*si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2));
                        co3*si2*theta3DDot+co2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-si2*si3*(theta2Dot^2+theta3Dot^2),-si2*theta2DDot-co2*theta2Dot^2,co2*(co3*theta2DDot-2*si3*theta2Dot*theta3Dot)-si2*(si3*theta3DDot+co3*(theta2Dot^2+theta3Dot^2));
                        co1*co3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)-si1*(co3*theta1Dot^2+si3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*co3*theta1Dot*theta3Dot+co3*theta3Dot^2)-co1*si3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),si1*si2*theta1DDot+co2*(-co1*theta2DDot+2*si1*theta1Dot*theta2Dot)+co1*si2*(theta1Dot^2+theta2Dot^2),-co1*(si3*theta1DDot+co3*si2*theta2DDot+2*co3*theta1Dot*theta3Dot-2*si2*si3*theta2Dot*theta3Dot)+si1*(si3*theta1Dot^2+2*co3*si2*theta1Dot*theta2Dot-co3*theta3DDot+si3*theta3Dot^2)-co2*(co3*si1*theta1DDot+si3*(co1*theta3DDot-2*si1*theta1Dot*theta3Dot)+co1*co3*(theta1Dot^2+theta2Dot^2+theta3Dot^2))];
           otherwise
                error('Invalid handedness provided.')
        end
    case 'yzy'
        switch(handed)
            case 'right'
                MDDot=[-co1*(si3*theta1DDot+co3*si2*theta2DDot+2*co3*theta1Dot*theta3Dot-2*si2*si3*theta2Dot*theta3Dot)+si1*(si3*theta1Dot^2+2*co3*si2*theta1Dot*theta2Dot-co3*theta3DDot+si3*theta3Dot^2)-co2*(co3*si1*theta1DDot+si3*(co1*theta3DDot-2*si1*theta1Dot*theta3Dot)+co1*co3*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),si1*si2*theta1DDot+co2*(-co1*theta2DDot+2*si1*theta1Dot*theta2Dot)+co1*si2*(theta1Dot^2+theta2Dot^2),co1*co3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)-si1*(co3*theta1Dot^2+si3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*co3*theta1Dot*theta3Dot+co3*theta3Dot^2)-co1*si3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2));
                        co2*(co3*theta2DDot-2*si3*theta2Dot*theta3Dot)-si2*(si3*theta3DDot+co3*(theta2Dot^2+theta3Dot^2)),-si2*theta2DDot-co2*theta2Dot^2,co3*si2*theta3DDot+co2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-si2*si3*(theta2Dot^2+theta3Dot^2);
                        si1*si3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(si3*theta1Dot^2-co3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*si3*theta1Dot*theta3Dot+si3*theta3Dot^2)+co3*si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),co2*si1*theta2DDot+co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-si1*si2*(theta1Dot^2+theta2Dot^2),-co3*(si1*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*co2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)))];
            case 'left'
                MDDot=[-co1*(si3*theta1DDot+co3*si2*theta2DDot+2*co3*theta1Dot*theta3Dot-2*si2*si3*theta2Dot*theta3Dot)+si1*(si3*theta1Dot^2+2*co3*si2*theta1Dot*theta2Dot-co3*theta3DDot+si3*theta3Dot^2)-co2*(co3*si1*theta1DDot+si3*(co1*theta3DDot-2*si1*theta1Dot*theta3Dot)+co1*co3*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),co1*co2*theta2DDot-si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-co1*si2*(theta1Dot^2+theta2Dot^2),-co1*co3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+si1*(co3*theta1Dot^2+si3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*co3*theta1Dot*theta3Dot+co3*theta3Dot^2)+co1*si3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2));
                         si2*si3*theta3DDot+co2*(-co3*theta2DDot+2*si3*theta2Dot*theta3Dot)+co3*si2*(theta2Dot^2+theta3Dot^2),-si2*theta2DDot-co2*theta2Dot^2,co3*si2*theta3DDot+co2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-si2*si3*(theta2Dot^2+theta3Dot^2);
                         co1*(co3*(-2*si2*theta1Dot*theta2Dot+theta3DDot)+co2*(co3*theta1DDot-2*si3*theta1Dot*theta3Dot)-si3*(theta1Dot^2+theta3Dot^2))-si1*(si3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),co2*si1*theta2DDot+co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-si1*si2*(theta1Dot^2+theta2Dot^2),-co3*(si1*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*co2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)))];
            otherwise
                error('Invalid handedness provided.')
        end
    case 'zyz'
        switch(handed)
            case 'right'
                MDDot=[-co1*(si3*theta1DDot+co3*si2*theta2DDot+2*co3*theta1Dot*theta3Dot-2*si2*si3*theta2Dot*theta3Dot)+si1*(si3*theta1Dot^2+2*co3*si2*theta1Dot*theta2Dot-co3*theta3DDot+si3*theta3Dot^2)-co2*(co3*si1*theta1DDot+si3*(co1*theta3DDot-2*si1*theta1Dot*theta3Dot)+co1*co3*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),-co1*co3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+si1*(co3*theta1Dot^2+si3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*co3*theta1Dot*theta3Dot+co3*theta3Dot^2)+co1*si3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),co1*co2*theta2DDot-si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-co1*si2*(theta1Dot^2+theta2Dot^2);
                       co1*(co3*(-2*si2*theta1Dot*theta2Dot+theta3DDot)+co2*(co3*theta1DDot-2*si3*theta1Dot*theta3Dot)-si3*(theta1Dot^2+theta3Dot^2))-si1*(si3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),-co3*(si1*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*co2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),co2*si1*theta2DDot+co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-si1*si2*(theta1Dot^2+theta2Dot^2);
                       si2*si3*theta3DDot+co2*(-co3*theta2DDot+2*si3*theta2Dot*theta3Dot)+co3*si2*(theta2Dot^2+theta3Dot^2),co3*si2*theta3DDot+co2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-si2*si3*(theta2Dot^2+theta3Dot^2),-si2*theta2DDot-co2*theta2Dot^2];
            case 'left'
                MDDot=[-co1*(si3*theta1DDot+co3*si2*theta2DDot+2*co3*theta1Dot*theta3Dot-2*si2*si3*theta2Dot*theta3Dot)+si1*(si3*theta1Dot^2+2*co3*si2*theta1Dot*theta2Dot-co3*theta3DDot+si3*theta3Dot^2)-co2*(co3*si1*theta1DDot+si3*(co1*theta3DDot-2*si1*theta1Dot*theta3Dot)+co1*co3*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),co1*co3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)-si1*(co3*theta1Dot^2+si3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*co3*theta1Dot*theta3Dot+co3*theta3Dot^2)-co1*si3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),si1*si2*theta1DDot+co2*(-co1*theta2DDot+2*si1*theta1Dot*theta2Dot)+co1*si2*(theta1Dot^2+theta2Dot^2);
                        si1*si3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(si3*theta1Dot^2-co3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*si3*theta1Dot*theta3Dot+si3*theta3Dot^2)+co3*si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),-co3*(si1*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*co2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),co2*si1*theta2DDot+co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-si1*si2*(theta1Dot^2+theta2Dot^2);
                        co2*(co3*theta2DDot-2*si3*theta2Dot*theta3Dot)-si2*(si3*theta3DDot+co3*(theta2Dot^2+theta3Dot^2)),co3*si2*theta3DDot+co2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-si2*si3*(theta2Dot^2+theta3Dot^2),-si2*theta2DDot-co2*theta2Dot^2];
            otherwise
                error('Invalid handedness provided.')
        end
    case 'zxz'
        switch(handed)
            case 'right'
                MDDot=[-co3*(si1*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*co2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),si1*si3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(si3*theta1Dot^2-co3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*si3*theta1Dot*theta3Dot+si3*theta3Dot^2)+co3*si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),co2*si1*theta2DDot+co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-si1*si2*(theta1Dot^2+theta2Dot^2);
                       co1*co3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)-si1*(co3*theta1Dot^2+si3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*co3*theta1Dot*theta3Dot+co3*theta3Dot^2)-co1*si3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),-co1*(si3*theta1DDot+co3*si2*theta2DDot+2*co3*theta1Dot*theta3Dot-2*si2*si3*theta2Dot*theta3Dot)+si1*(si3*theta1Dot^2+2*co3*si2*theta1Dot*theta2Dot-co3*theta3DDot+si3*theta3Dot^2)-co2*(co3*si1*theta1DDot+si3*(co1*theta3DDot-2*si1*theta1Dot*theta3Dot)+co1*co3*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),si1*si2*theta1DDot+co2*(-co1*theta2DDot+2*si1*theta1Dot*theta2Dot)+co1*si2*(theta1Dot^2+theta2Dot^2);
                       co3*si2*theta3DDot+co2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-si2*si3*(theta2Dot^2+theta3Dot^2),co2*(co3*theta2DDot-2*si3*theta2Dot*theta3Dot)-si2*(si3*theta3DDot+co3*(theta2Dot^2+theta3Dot^2)),-si2*theta2DDot-co2*theta2Dot^2];
            case 'left'
                MDDot=[-co3*(si1*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*co2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+si1*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),co1*(co3*(-2*si2*theta1Dot*theta2Dot+theta3DDot)+co2*(co3*theta1DDot-2*si3*theta1Dot*theta3Dot)-si3*(theta1Dot^2+theta3Dot^2))-si1*(si3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+co3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),co2*si1*theta2DDot+co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-si1*si2*(theta1Dot^2+theta2Dot^2);
                       -co1*co3*(theta1DDot+co2*theta3DDot-2*si2*theta2Dot*theta3Dot)+si1*(co3*theta1Dot^2+si3*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot+theta3DDot)+2*co2*co3*theta1Dot*theta3Dot+co3*theta3Dot^2)+co1*si3*(si2*theta2DDot+2*theta1Dot*theta3Dot+co2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),-co1*(si3*theta1DDot+co3*si2*theta2DDot+2*co3*theta1Dot*theta3Dot-2*si2*si3*theta2Dot*theta3Dot)+si1*(si3*theta1Dot^2+2*co3*si2*theta1Dot*theta2Dot-co3*theta3DDot+si3*theta3Dot^2)-co2*(co3*si1*theta1DDot+si3*(co1*theta3DDot-2*si1*theta1Dot*theta3Dot)+co1*co3*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),co1*co2*theta2DDot-si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot)-co1*si2*(theta1Dot^2+theta2Dot^2);
                       co3*si2*theta3DDot+co2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-si2*si3*(theta2Dot^2+theta3Dot^2),si2*si3*theta3DDot+co2*(-co3*theta2DDot+2*si3*theta2Dot*theta3Dot)+co3*si2*(theta2Dot^2+theta3Dot^2),-si2*theta2DDot-co2*theta2Dot^2];
            otherwise
                error('Invalid handedness provided.')
        end
    %Next, the Tait-Bryan-Cardan angle combinations.
    case 'xzy'
        switch(handed)
            case 'right'
                MDDot=[-co2*si3*theta3DDot+2*si2*si3*theta2Dot*theta3Dot-co3*(si2*theta2DDot+co2*(theta2Dot^2+theta3Dot^2)),-co2*theta2DDot+si2*theta2Dot^2,co2*co3*theta3DDot-si2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-co2*si3*(theta2Dot^2+theta3Dot^2);
                        si3*(co1*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)-si1*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2))-co3*(si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+co1*(-co2*theta2DDot-2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),-co1*si2*theta2DDot+2*si1*si2*theta1Dot*theta2Dot-co2*(si1*theta1DDot+co1*(theta1Dot^2+theta2Dot^2)),si1*(si3*(-2*co2*theta1Dot*theta2Dot+theta3DDot)-si2*(si3*theta1DDot+2*co3*theta1Dot*theta3Dot)+co3*(theta1Dot^2+theta3Dot^2))-co1*(co3*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+si3*(-co2*theta2DDot-2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)));
                        co1*co3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+si1*si3*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+co1*si3*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2)+co3*si1*(co2*theta2DDot+2*theta1Dot*theta3Dot-si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot)-si1*(si2*theta2DDot+co2*(theta1Dot^2+theta2Dot^2)),-co3*(si1*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+co1*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+si1*(co2*theta2DDot+2*theta1Dot*theta3Dot-si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)))];
            case 'left'
                MDDot=[-co2*si3*theta3DDot+2*si2*si3*theta2Dot*theta3Dot-co3*(si2*theta2DDot+co2*(theta2Dot^2+theta3Dot^2)),co2*theta2DDot-si2*theta2Dot^2,si2*si3*theta2DDot+co3*(-co2*theta3DDot+2*si2*theta2Dot*theta3Dot)+co2*si3*(theta2Dot^2+theta3Dot^2);
                        co3*si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+co1*si3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)-si1*si3*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2)+co1*co3*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),-co1*si2*theta2DDot+2*si1*si2*theta1Dot*theta2Dot-co2*(si1*theta1DDot+co1*(theta1Dot^2+theta2Dot^2)),-si1*(co3*theta1Dot^2+si3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+2*co3*si2*theta1Dot*theta3Dot+co3*theta3Dot^2)+co1*(co2*si3*theta2DDot-2*si3*theta1Dot*theta3Dot+co3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)-si2*si3*(theta1Dot^2+theta2Dot^2+theta3Dot^2));
                        co1*(co3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)-si3*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2))-si1*(si3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)+co3*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),si1*si2*theta2DDot+co1*(-co2*theta1DDot+2*si2*theta1Dot*theta2Dot)+co2*si1*(theta1Dot^2+theta2Dot^2),-co3*(si1*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+si1*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)))];
            otherwise
                error('Invalid handedness provided.')
        end
    case 'xyz'
        switch(handed)
            case 'right'
                MDDot=[-co2*si3*theta3DDot+2*si2*si3*theta2Dot*theta3Dot-co3*(si2*theta2DDot+co2*(theta2Dot^2+theta3Dot^2)),si2*si3*theta2DDot+co3*(-co2*theta3DDot+2*si2*theta2Dot*theta3Dot)+co2*si3*(theta2Dot^2+theta3Dot^2),co2*theta2DDot-si2*theta2Dot^2;
                       co1*(co3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)-si3*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2))-si1*(si3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)+co3*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),-co3*(si1*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+si1*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),si1*si2*theta2DDot+co1*(-co2*theta1DDot+2*si2*theta1Dot*theta2Dot)+co2*si1*(theta1Dot^2+theta2Dot^2);
                       co3*si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+co1*si3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)-si1*si3*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2)+co1*co3*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),-si1*(co3*theta1Dot^2+si3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+2*co3*si2*theta1Dot*theta3Dot+co3*theta3Dot^2)+co1*(co2*si3*theta2DDot-2*si3*theta1Dot*theta3Dot+co3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)-si2*si3*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),-co1*si2*theta2DDot+2*si1*si2*theta1Dot*theta2Dot-co2*(si1*theta1DDot+co1*(theta1Dot^2+theta2Dot^2))];
            case 'left'
                MDDot=[-co2*si3*theta3DDot+2*si2*si3*theta2Dot*theta3Dot-co3*(si2*theta2DDot+co2*(theta2Dot^2+theta3Dot^2)),co2*co3*theta3DDot-si2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-co2*si3*(theta2Dot^2+theta3Dot^2),-co2*theta2DDot+si2*theta2Dot^2;
                        co1*co3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+si1*si3*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+co1*si3*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2)+co3*si1*(co2*theta2DDot+2*theta1Dot*theta3Dot-si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),-co3*(si1*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+co1*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+si1*(co2*theta2DDot+2*theta1Dot*theta3Dot-si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot)-si1*(si2*theta2DDot+co2*(theta1Dot^2+theta2Dot^2));
                        si3*(co1*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)-si1*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2))-co3*(si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+co1*(-co2*theta2DDot-2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),si1*(si3*(-2*co2*theta1Dot*theta2Dot+theta3DDot)-si2*(si3*theta1DDot+2*co3*theta1Dot*theta3Dot)+co3*(theta1Dot^2+theta3Dot^2))-co1*(co3*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+si3*(-co2*theta2DDot-2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),-co1*si2*theta2DDot+2*si1*si2*theta1Dot*theta2Dot-co2*(si1*theta1DDot+co1*(theta1Dot^2+theta2Dot^2))];
            otherwise
                error('Invalid handedness provided.')
        end
    case 'yxz'
        switch(handed)
            case 'right'
                MDDot=[-co3*(si1*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+co1*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+si1*(co2*theta2DDot+2*theta1Dot*theta3Dot-si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),co1*co3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+si1*si3*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+co1*si3*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2)+co3*si1*(co2*theta2DDot+2*theta1Dot*theta3Dot-si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot)-si1*(si2*theta2DDot+co2*(theta1Dot^2+theta2Dot^2));
                        co2*co3*theta3DDot-si2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-co2*si3*(theta2Dot^2+theta3Dot^2),-co2*si3*theta3DDot+2*si2*si3*theta2Dot*theta3Dot-co3*(si2*theta2DDot+co2*(theta2Dot^2+theta3Dot^2)),-co2*theta2DDot+si2*theta2Dot^2;
                        si1*(si3*(-2*co2*theta1Dot*theta2Dot+theta3DDot)-si2*(si3*theta1DDot+2*co3*theta1Dot*theta3Dot)+co3*(theta1Dot^2+theta3Dot^2))-co1*(co3*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+si3*(-co2*theta2DDot-2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),si3*(co1*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)-si1*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2))-co3*(si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+co1*(-co2*theta2DDot-2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),-co1*si2*theta2DDot+2*si1*si2*theta1Dot*theta2Dot-co2*(si1*theta1DDot+co1*(theta1Dot^2+theta2Dot^2))];
            case 'left'
                MDDot=[-co3*(si1*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+si1*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),co1*(co3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)-si3*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2))-si1*(si3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)+co3*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),si1*si2*theta2DDot+co1*(-co2*theta1DDot+2*si2*theta1Dot*theta2Dot)+co2*si1*(theta1Dot^2+theta2Dot^2);
                        si2*si3*theta2DDot+co3*(-co2*theta3DDot+2*si2*theta2Dot*theta3Dot)+co2*si3*(theta2Dot^2+theta3Dot^2),-co2*si3*theta3DDot+2*si2*si3*theta2Dot*theta3Dot-co3*(si2*theta2DDot+co2*(theta2Dot^2+theta3Dot^2)),co2*theta2DDot-si2*theta2Dot^2;
                        -si1*(co3*theta1Dot^2+si3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+2*co3*si2*theta1Dot*theta3Dot+co3*theta3Dot^2)+co1*(co2*si3*theta2DDot-2*si3*theta1Dot*theta3Dot+co3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)-si2*si3*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),co3*si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+co1*si3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)-si1*si3*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2)+co1*co3*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),-co1*si2*theta2DDot+2*si1*si2*theta1Dot*theta2Dot-co2*(si1*theta1DDot+co1*(theta1Dot^2+theta2Dot^2))];
            otherwise
                error('Invalid handedness provided.')
        end
    case 'yzx'
        switch(handed)
            case 'right'
                MDDot=[-co1*si2*theta2DDot+2*si1*si2*theta1Dot*theta2Dot-co2*(si1*theta1DDot+co1*(theta1Dot^2+theta2Dot^2)),co3*si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+co1*si3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)-si1*si3*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2)+co1*co3*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),-si1*(co3*theta1Dot^2+si3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+2*co3*si2*theta1Dot*theta3Dot+co3*theta3Dot^2)+co1*(co2*si3*theta2DDot-2*si3*theta1Dot*theta3Dot+co3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)-si2*si3*(theta1Dot^2+theta2Dot^2+theta3Dot^2));
                        co2*theta2DDot-si2*theta2Dot^2,-co2*si3*theta3DDot+2*si2*si3*theta2Dot*theta3Dot-co3*(si2*theta2DDot+co2*(theta2Dot^2+theta3Dot^2)),si2*si3*theta2DDot+co3*(-co2*theta3DDot+2*si2*theta2Dot*theta3Dot)+co2*si3*(theta2Dot^2+theta3Dot^2);
                        si1*si2*theta2DDot+co1*(-co2*theta1DDot+2*si2*theta1Dot*theta2Dot)+co2*si1*(theta1Dot^2+theta2Dot^2),co1*(co3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)-si3*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2))-si1*(si3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)+co3*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),-co3*(si1*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+si1*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)))];
            case 'left'
                MDDot=[-co1*si2*theta2DDot+2*si1*si2*theta1Dot*theta2Dot-co2*(si1*theta1DDot+co1*(theta1Dot^2+theta2Dot^2)),si3*(co1*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)-si1*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2))-co3*(si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+co1*(-co2*theta2DDot-2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),si1*(si3*(-2*co2*theta1Dot*theta2Dot+theta3DDot)-si2*(si3*theta1DDot+2*co3*theta1Dot*theta3Dot)+co3*(theta1Dot^2+theta3Dot^2))-co1*(co3*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+si3*(-co2*theta2DDot-2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)));
                    -co2*theta2DDot+si2*theta2Dot^2,-co2*si3*theta3DDot+2*si2*si3*theta2Dot*theta3Dot-co3*(si2*theta2DDot+co2*(theta2Dot^2+theta3Dot^2)),co2*co3*theta3DDot-si2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-co2*si3*(theta2Dot^2+theta3Dot^2);
                    co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot)-si1*(si2*theta2DDot+co2*(theta1Dot^2+theta2Dot^2)),co1*co3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+si1*si3*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+co1*si3*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2)+co3*si1*(co2*theta2DDot+2*theta1Dot*theta3Dot-si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),-co3*(si1*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+co1*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+si1*(co2*theta2DDot+2*theta1Dot*theta3Dot-si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)))];
            otherwise
                error('Invalid handedness provided.')
        end
    case 'zyx'
        switch(handed)
            case 'right'
                MDDot=[-co1*si2*theta2DDot+2*si1*si2*theta1Dot*theta2Dot-co2*(si1*theta1DDot+co1*(theta1Dot^2+theta2Dot^2)),si1*(si3*(-2*co2*theta1Dot*theta2Dot+theta3DDot)-si2*(si3*theta1DDot+2*co3*theta1Dot*theta3Dot)+co3*(theta1Dot^2+theta3Dot^2))-co1*(co3*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+si3*(-co2*theta2DDot-2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),si3*(co1*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)-si1*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2))-co3*(si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+co1*(-co2*theta2DDot-2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)));
                        co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot)-si1*(si2*theta2DDot+co2*(theta1Dot^2+theta2Dot^2)),-co3*(si1*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+co1*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+si1*(co2*theta2DDot+2*theta1Dot*theta3Dot-si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),co1*co3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+si1*si3*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+co1*si3*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2)+co3*si1*(co2*theta2DDot+2*theta1Dot*theta3Dot-si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2));
                        -co2*theta2DDot+si2*theta2Dot^2,co2*co3*theta3DDot-si2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-co2*si3*(theta2Dot^2+theta3Dot^2),-co2*si3*theta3DDot+2*si2*si3*theta2Dot*theta3Dot-co3*(si2*theta2DDot+co2*(theta2Dot^2+theta3Dot^2))];
            case 'left'
                MDDot=[-co1*si2*theta2DDot+2*si1*si2*theta1Dot*theta2Dot-co2*(si1*theta1DDot+co1*(theta1Dot^2+theta2Dot^2)),-si1*(co3*theta1Dot^2+si3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+2*co3*si2*theta1Dot*theta3Dot+co3*theta3Dot^2)+co1*(co2*si3*theta2DDot-2*si3*theta1Dot*theta3Dot+co3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)-si2*si3*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),co3*si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+co1*si3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)-si1*si3*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2)+co1*co3*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2));
                    si1*si2*theta2DDot+co1*(-co2*theta1DDot+2*si2*theta1Dot*theta2Dot)+co2*si1*(theta1Dot^2+theta2Dot^2),-co3*(si1*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+si1*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),co1*(co3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)-si3*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2))-si1*(si3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)+co3*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)));
                    co2*theta2DDot-si2*theta2Dot^2,si2*si3*theta2DDot+co3*(-co2*theta3DDot+2*si2*theta2Dot*theta3Dot)+co2*si3*(theta2Dot^2+theta3Dot^2),-co2*si3*theta3DDot+2*si2*si3*theta2Dot*theta3Dot-co3*(si2*theta2DDot+co2*(theta2Dot^2+theta3Dot^2))];
            otherwise
                error('Invalid handedness provided.')
        end
    case 'zxy'
        switch(handed)
            case 'right'
                MDDot=[-co3*(si1*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)+co1*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(-co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+si1*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),si1*si2*theta2DDot+co1*(-co2*theta1DDot+2*si2*theta1Dot*theta2Dot)+co2*si1*(theta1Dot^2+theta2Dot^2),co1*(co3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)-si3*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2))-si1*(si3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)+co3*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)));
                       -si1*(co3*theta1Dot^2+si3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+2*co3*si2*theta1Dot*theta3Dot+co3*theta3Dot^2)+co1*(co2*si3*theta2DDot-2*si3*theta1Dot*theta3Dot+co3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)-si2*si3*(theta1Dot^2+theta2Dot^2+theta3Dot^2)),-co1*si2*theta2DDot+2*si1*si2*theta1Dot*theta2Dot-co2*(si1*theta1DDot+co1*(theta1Dot^2+theta2Dot^2)),co3*si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot+theta3DDot)+co1*si3*(theta1DDot+si2*theta3DDot+2*co2*theta2Dot*theta3Dot)-si1*si3*(theta1Dot^2+2*si2*theta1Dot*theta3Dot+theta3Dot^2)+co1*co3*(-co2*theta2DDot+2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2));
                       si2*si3*theta2DDot+co3*(-co2*theta3DDot+2*si2*theta2Dot*theta3Dot)+co2*si3*(theta2Dot^2+theta3Dot^2),co2*theta2DDot-si2*theta2Dot^2,-co2*si3*theta3DDot+2*si2*si3*theta2Dot*theta3Dot-co3*(si2*theta2DDot+co2*(theta2Dot^2+theta3Dot^2))];
            case 'left'
                MDDot=[-co3*(si1*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+co1*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2))+si3*(co1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+si1*(co2*theta2DDot+2*theta1Dot*theta3Dot-si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),co1*(co2*theta1DDot-2*si2*theta1Dot*theta2Dot)-si1*(si2*theta2DDot+co2*(theta1Dot^2+theta2Dot^2)),co1*co3*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+si1*si3*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+co1*si3*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2)+co3*si1*(co2*theta2DDot+2*theta1Dot*theta3Dot-si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2));
                        si1*(si3*(-2*co2*theta1Dot*theta2Dot+theta3DDot)-si2*(si3*theta1DDot+2*co3*theta1Dot*theta3Dot)+co3*(theta1Dot^2+theta3Dot^2))-co1*(co3*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)+si3*(-co2*theta2DDot-2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2))),-co1*si2*theta2DDot+2*si1*si2*theta1Dot*theta2Dot-co2*(si1*theta1DDot+co1*(theta1Dot^2+theta2Dot^2)),si3*(co1*(theta1DDot-si2*theta3DDot-2*co2*theta2Dot*theta3Dot)-si1*(theta1Dot^2-2*si2*theta1Dot*theta3Dot+theta3Dot^2))-co3*(si1*(si2*theta1DDot+2*co2*theta1Dot*theta2Dot-theta3DDot)+co1*(-co2*theta2DDot-2*theta1Dot*theta3Dot+si2*(theta1Dot^2+theta2Dot^2+theta3Dot^2)));
                        co2*co3*theta3DDot-si2*(si3*theta2DDot+2*co3*theta2Dot*theta3Dot)-co2*si3*(theta2Dot^2+theta3Dot^2),-co2*theta2DDot+si2*theta2Dot^2,-co2*si3*theta3DDot+2*si2*si3*theta2Dot*theta3Dot-co3*(si2*theta2DDot+co2*(theta2Dot^2+theta3Dot^2))];
            otherwise
                error('Invalid handedness provided.')
        end
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
