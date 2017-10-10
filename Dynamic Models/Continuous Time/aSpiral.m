function aVals=aSpiral(xPoints,t)
%%ASPIRAL The drift function for a spiraling target motion model in 3
%         dimensions. This formulation of the dynamic model is more
%         difficult to use to make a target go in a desired direction than
%         the function aSpiralSimp.
%
%INPUTS: xPoints The 10X1 target state at time t. It consists of position,
%                velocity, acceleration, and a scalar spiral rate. If x is
%                an xDim X numStates matrix, then the spiraling model is
%                evaluated for all of the state vectors.
%              t An unused time component so that aSpiral can be used with
%                Runge-Kutta methods that expect the function to take two
%                parameters.
%
%OUTPUTS: aVals The flat-Earth time-derivative of the state.
%
%A derivation of the flat-Earth spiraling dynamic model is given in [1],
%where various orthogonality assumptions are present to keep the speed of
%the target from increasing or decreasing.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear continuous-time dynamic
%    models," IEEE Aerospace and Electronic Systems Magazine, vol. 30, no.
%    2, Part II, pp. 4-41, Feb. 2015.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numPoints=size(xPoints,2);
    aVals=zeros(10,numPoints);
    
    aVals(1:6,:)=xPoints(4:9,:);
    
    rdot=xPoints(4:6,:);
    rddot=xPoints(7:9,:);
    omega=xPoints(10,:);
    
    normRDot = sqrt(xPoints(4,:).^2+xPoints(5,:).^2+xPoints(6,:).^2);
    normRDDotSqrd = xPoints(7,:).^2+xPoints(8,:).^2+xPoints(9,:).^2;

    aVals(7:9,:)=bsxfun(@rdivide,bsxfun(@times,bsxfun(@cross,rddot,rdot),omega),normRDot)-bsxfun(@times,rdot,normRDDotSqrd./(normRDot.^2));

    %Deal with errors that arise if the velocity is near zero, so 0/0
    %occurs.
    aVals(~isfinite(aVals))=0;
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
