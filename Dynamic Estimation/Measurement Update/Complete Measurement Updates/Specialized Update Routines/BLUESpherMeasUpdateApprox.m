function [xStateUpdate,PUpdate,innov,S,W]=BLUESpherMeasUpdateApprox(xStatePred,PPred,z,R)
%%BLUESPHERMEASUPDATEAPPROX Perform the measurement update step in the
%                   approximate best linear unbiased estimator (BLUE) for a
%                   Cartesian state consisting of 3D position and velocity
%                   components when given a measurement in monostatic
%                   spherical coordinates.
%
%INPUTS: xStatePred The 6X1 target state consisting of position components
%                   followed by velocity components:
%                   xStatePred=[x;y;z;xDot;yDot;zDot].
%             PPred The 6X6 covariance matrix associated with the predicted
%                   target state estimate.
%                 z A 3X1 monostatic spherical measurement (no refraction)
%                   with components ordered [range;azimuth;elevation] with
%                   azimuth and elevation angles in radians. Azimuth is 
%                   measured in the x-y plane, counterclockwise from the x
%                   axis and elevation is measured up from the x-y plane.
%                 R The 3X3 diagonal covariance matrix associated with the
%                   measurement z. Any cross terms present in R will be
%                   ignored.
%
%OUTPUTS: xStateUpdate The 6X1 updated state vector.
%              PUpdate The updated 6X6 state covariance matrix.
%             innov, S The 3X1 innovation and the 3X3 innovation covariance
%                      matrix are returned in case one wishes to analyze
%                      the consistency of the estimator or use those values
%                      in gating or likelihood evaluation.
%                    W The gain used in the update. This can be useful
%                      when gating and using the function
%                      calcMissedGateCov.
%
%The algorithm is the set of equations taken from Table I of [1].
%
%REFERENCES:
%[1] Z. Zhao, X. R. Li, and V. P. Jilkov, "Best linear unbiased filtering
%    with nonlinear measurements for target tracking," IEEE Transactions on
%    Aerospace and Electronic Systems, vol. 40, no. 4, pp. 1324-1336,
%    Oct. 2004.
%
%February 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Change the ordering of the elements in the state to match that used in the
%paper.
permuteIdx=[1;4;2;5;3;6];
undoPermuteIdx=[1;3;5;2;4;6];%The inverse permutation
xStateBar=xStatePred(permuteIdx);
PBar=PPred(permuteIdx,permuteIdx);

xBar=xStateBar(1);
yBar=xStateBar(3);
zBar=xStateBar(5);

sigmaR2=R(1,1);%Range standard deviation
sigmaTheta2=R(2,2);%Azimuth standard deviation
sigmaPhi2=R(3,3);%Elevation standard deviation

%The constants
lambda1=exp(-sigmaTheta2/2);
lambda2=(1/2)*(1+exp(-2*sigmaTheta2));
lambda3=(1/2)*(1-exp(-2*sigmaTheta2));
mu1=exp(-sigmaPhi2/2);
mu2=(1/2)*(1+exp(-2*sigmaPhi2));
mu3=(1/2)*(1-exp(-2*sigmaPhi2));

rBar=sqrt(xBar^2+yBar^2+zBar^2);
r1Bar=sqrt(xBar^2+yBar^2);
alpha=mu2*sigmaR2/rBar^2+mu3*zBar^2/r1Bar^2+mu3*sigmaR2*zBar^2/(r1Bar^2*rBar^2);
alpha1=(lambda2*mu2-lambda1^2*mu1^2)*xBar^2+lambda3*mu2*yBar^2;
alpha2=(lambda2*mu2-lambda1^2*mu1^2)*yBar^2+lambda3*mu2*xBar^2;
alpha3=(mu2-mu1^2)*zBar^2+mu3*(xBar^2+yBar^2);
alpha4=(mu2*(lambda2-lambda3)-lambda1^2*mu1^2)*xBar*yBar;
alpha5=(lambda1*(mu2-mu3)-lambda1*mu1^2)*zBar;

S(1,1)=lambda2*mu2*PBar(1,1)+lambda3*mu2*PBar(3,3)+alpha*(lambda2*xBar^2+lambda3*yBar^2)+alpha1;
S(2,2)=lambda2*mu2*PBar(3,3)+lambda3*mu2*PBar(1,1)+alpha*(lambda3*xBar^2+lambda2*yBar^2)+alpha2;
S(3,3)=mu2*PBar(5,5)+mu3*(PBar(1,1)+PBar(3,3)+mu2*sigmaR2*zBar^2/rBar^2+mu3*sigmaR2*r1Bar^2/rBar^2)+alpha3;
%The paper includes an open parenthesis with no closed parenthesis for the
%following equation. The closed-parenthesis has been added in the correct
%position.
S(1,2)=(lambda2-lambda3)*(mu2*PBar(1,3)+alpha*xBar*yBar)+alpha4;
S(2,1)=S(1,2);
S(1,3)=lambda1*(mu2-mu3)*(PBar(1,5)+sigmaR2*xBar*zBar/rBar^2)+alpha5*xBar;
S(3,1)=S(1,3);
S(2,3)=lambda1*(mu2-mu3)*(PBar(3,5)+sigmaR2*yBar*zBar/rBar^2)+alpha5*yBar;
S(3,2)=S(2,3);

W=mu1*[lambda1*PBar(:,1),lambda1*PBar(:,3),PBar(:,5)]/S;

zCart=spher2Cart(z);
innov=zCart-mu1*[lambda1*xBar;lambda1*yBar;zBar];
xStateUpdate=xStateBar+W*innov;
PUpdate=PBar-W*S*W';

%Make the ordering go back to the way it was.
xStateUpdate=xStateUpdate(undoPermuteIdx);
PUpdate=PUpdate(undoPermuteIdx,undoPermuteIdx);
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
