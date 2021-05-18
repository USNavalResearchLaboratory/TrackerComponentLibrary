function [xStateUpdate,PUpdate,innov,Pzz,W]=BLUEPolarMeasUpdateApprox(xStatePred,PPred,z,R)
%%BLUEPOLARMEASUPDATEAPPROX Perform the measurement update step in the
%                approximate best linear unbiased estimator (BLUE) for a
%                Cartesian state consisting of 2D position and velocity
%                when given a measurement in monostatic polar coordinates.
%
%INPUTS: xStatePred The 4X1 target state consisting of position components
%                   followed by velocity components:
%                   xStatePred=[x;y;xDot;yDot].
%             PPred The 3X3 covariance matrix associated with the predicted
%                   target state estimate.
%                 z A 2X1 one-way (monostatic) polar measurement (no
%                   refraction) with components ordered [range;azimuth]
%                   with azimuth in radians measured in radians,
%                   counterclockwise from the x-axis.
%                 R The 2X2 diagonal covariance matrix associated with the
%                   measurement z. Any cross terms present in R will be
%                   ignored.
%
%OUTPUTS: xStateUpdate The 4 X 1 updated state vector.
%              PUpdate The updated 4 X 4 state covariance matrix.
%           innov, Pzz The 2X1 innovation and the 2X2 innovation covariance
%                      matrix are returned in case one wishes to analyze
%                      the consistency of the estimator or use those values
%                      in gating or likelihood evaluation.
%                    W The gain used in the update. This can be useful
%                      when gating and using the function
%                      calcMissedGateCov.
%
%The algorithm is taken from Section IV of [1]. One equation that is not
%directly written in the above journal article is taken from the conference
%version [2].
%
%REFERENCES:
%[1] Z. Zhao, X. R. Li, and V. P. Jilkov, "Best linear unbiased filtering
%    with nonlinear measurements for target tracking," IEEE Transactions on
%    Aerospace and Electronic Systems, vol. 40, no. 4, pp. 1324-1336,
%    Oct. 2004.
%[2] Z. Zhao, X. R. Li, and V. P. Jilkov, "Optimal linear unbiased
%    filtering with polar measurements for target tracking," in Proceedings
%    of the 5th International Conference on Information Fusion, Annapolis,
%    MD, 8-11 Jul. 2002, pp. 1527-1534.
%
%February 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Change the ordering of the elements in the state to match that used in the
%paper.
permuteIdx=[1;3;2;4];
undoPermuteIdx=[1;3;2;4];%The inverse permutation
xStateBar=xStatePred(permuteIdx);
PBar=PPred(permuteIdx,permuteIdx);

xBar=xStateBar(1);
yBar=xStateBar(3);

covXTilde=PBar(1,1);
covYTilde=PBar(3,3);
covXTileYTilde=PBar(1,3);

sigmaR2=R(1,1);
sigmaTheta2=R(2,2);

%The Taylor series approximation to the expected value of
%(y^2-x^2)/(x^2+y^2) from Section VI of the journal article.
denom=(xBar^2+yBar^2)^3;
ExyRat1=(yBar^2-xBar^2)/(xBar^2+yBar^2)+...
        2*yBar^2*(yBar^2-3*xBar^2)*covXTilde/denom+...
        4*xBar*yBar*(xBar^2-yBar^2)*covXTileYTilde/denom-...
        2*xBar^2*(xBar^2-3*yBar^2)*covYTilde/denom;

%The Taylor series approximation to the expected value of x*y/(x^2+y^2)
%from Section VI is taken from the form written out explicitly in Equation
%18 of the conference paper.
ExyRat2=xBar*yBar/(xBar^2+yBar^2)+(1/2)*(...
        2*xBar*yBar*(xBar^2-3*yBar^2)*covXTilde/denom+...
        (6*xBar^2*yBar^2-xBar^4-yBar^4)*covXTileYTilde/denom+...
        2*xBar*yBar*(yBar^2-3*xBar^2)*covYTilde/denom);

lambda1=exp(-sigmaTheta2/2);
lambda2=(1/2)*(1+exp(-2*sigmaTheta2));
lambda3=(1/2)*(1-exp(-2*sigmaTheta2));

%The value of S11 is from Section IV of the journal article.
Pzz(1,1)=lambda2*covXTilde+lambda3*covYTilde+(1/2)*sigmaR2+lambda3*yBar^2+...
    (lambda2-lambda1^2)*xBar^2+(1/2)*sigmaR2*exp(-2*sigmaTheta2)*(-ExyRat1);
Pzz(2,2)=lambda2*covYTilde+lambda3*covXTilde+(1/2)*sigmaR2+lambda3*xBar^2+...
    (lambda2-lambda1^2)*yBar^2+(1/2)*sigmaR2*exp(-2*sigmaTheta2)*(ExyRat1);
Pzz(1,2)=exp(-2*sigmaTheta2)*covXTileYTilde+(exp(-2*sigmaTheta2)-...
    lambda1^2)*xBar*yBar-sigmaR2*exp(-2*sigmaTheta2)*ExyRat2;
Pzz(2,1)=Pzz(1,2);

%Equation 17
covXTildeZTilde=lambda1*[PBar(:,1),PBar(:,3)];
W=covXTildeZTilde/Pzz;
zCartPred=lambda1*[xBar;yBar];

zCart=pol2Cart(z);
innov=zCart-zCartPred;
xStateUpdate=xStateBar+W*(zCart-zCartPred);
PUpdate=PBar-W*Pzz*W';

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
