function MDot=rotMatDeriv2D(theta,thetaDot,isClockwise)
%%ROTMATDERIV2D Given an initial angle of rotation and the angle rate
%   (derivative of the angle with respect to time), find the derivative of
%   the 2D rotation matrix. This is the time derivative of rotMat2D
%   assuming that the rotation angle changes as a function of time.
%
%INPUTS: theta The initial rotation angle radians about which the
%              time-varying 2D rotation matrix would rotate a 2X1 vector at
%              the current time.
%     thetaDot The derivative of theta with respect to time.
%  isClockwise An optional boolean parameter indication whether the
%              rotation implied by theta is clockwise. The default if
%              omitted or an empty matrix is passed is false.
%
%OUTPUTS: MDot The 2X2 derivative of the rotation matrix (as one could
%              obtain using rotMat2D) with respect to time.
%
%EXAMPLE:
%This example shows that this function is consistent with numeric
%differentiation. The relative error of the finite differenced solution
%typically implies for 6-8 digits of precision.
% theta=2*pi*rand(1);
% thetaDot=randn(1);
% isClockwise=false;
% MDot=rotMatDeriv2D(theta,thetaDot,isClockwise);
% deltaT=1e-8;
% M=rotMat2D(theta,isClockwise);
% MDelta=rotMat2D(theta+deltaT*thetaDot,isClockwise);
% MDotNumDiff=(MDelta-M)/deltaT;
% RelErr=max(abs((MDotNumDiff(:)-MDot(:))./MDot(:)))
%
%December 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(isClockwise))
    isClockwise=false; 
end

sinTheta=sin(theta);
cosTheta=cos(theta);

if(isClockwise)
    MDot=thetaDot*[-sinTheta,cosTheta;
                   -cosTheta,-sinTheta];
else
    MDot=thetaDot*[-sinTheta,-cosTheta;
                    cosTheta,-sinTheta];
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
