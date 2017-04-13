function aVals=aSpiralSimp(xPoints,t)
%%ASPIRALSIMP The drift function for a non-ballistic spiraling target
%         motion model in 3 dimensions, formulated in a manner such that it
%         is simple to make a spiraling target follow a nominal trajectory.
%         This model is probably bad to design into a target tracking
%         algorithm, but might be good for designing spiraling trajectories
%         to use to test target tracking algorithms.
%
%INPUTS: xPoints The target state at time t. It consists of position 
%                (3 elements), instantaneous velocity (3 elements), the
%                velocity vector (3 elements) of the overall direction of
%                motion of the spiraling model (ground speed) and the
%                spiral rate of the target. Thus xDim=10. If x is an
%                xDim X numStates matrix, then the spiraling model is
%                evaluated for all of the state vectors.
%             t  An unused time component so that aSpiral can be used with
%                Runge-Kutta methods that expect the function to take two
%                parameters.
%
%OUTPUTS: aVals The flat-Earth time-derivative of the state.
%
%A derivation of the simplified flat-Earth spiraling dynamic model is given
%in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Simulating aerial targets in 3D accounting for the
%    Earth's curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, Jun. 2015.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numPoints=size(xPoints,2);
    aVals=zeros(6,numPoints);
    
    vl=xPoints(7:9,:);
    omega=xPoints(10,:);
    
    aVals(1:3,:)=xPoints(4:6,:);
    vs=xPoints(4:6,:)-vl;
    ul=bsxfun(@rdivide,vl,sqrt(sum(vl.*vl,1)));
    vsDot=omega.*bsxfun(@cross,ul,vs);
    aVals(4:6,:)=vsDot;
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
