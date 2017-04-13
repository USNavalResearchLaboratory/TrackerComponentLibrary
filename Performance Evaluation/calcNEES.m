function NEES=calcNEES(xTrue,xEst,PEst)
%CALCNEES Given a batch of estimates and associated true values, compute
%         the normalized estimation error squared (NEES). This version is
%         normalized such that the mean of a consistent estimator,
%         regardless of the dimensionality of the state, should be one.
%
%INPUTS: xTrue The truth data. If the truth is the same for all samples of
%              the estimate, then this is an xDimX1 vector. Otherwise, this
%              is an xDimXN matrix where N is the number of samples in the
%              batch.
%         xEst An xDim X N set of estimates.
%         PEst An xDim X xDim XN set of covariance matrices associated
%              with the estimates. If the covariance matrix is supposed to
%              be the same for all estimates, then PEst is just xDimXxDim.
%              These should not be singular.
%
%OUTPUTS: NEES The normalized estimation error squared of the estimates.
%
%The NEES is a measure of how well the provided by an estimator matches the
%actual accuracy of the estimates. The concept of the NEES is discussed in 
%[1].
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xEst,1);
numEst=size(xEst,2);

if(size(xTrue,2)==1)
    xTrue=repmat(xTrue,[1,numEst]);
end

if(size(PEst,3)==1)
    PEst=repmat(PEst,[1,1,numEst]);
end

NEES=0;
for curEst=1:numEst
    diff=xTrue(:,curEst)-xEst(:,curEst);

    NEES=NEES+invSymQuadForm(diff,PEst(:,:,curEst));%=NEES+diff'*inv(PEst(:,:,curEst))*diff;
end
NEES=NEES/(xDim*numEst);

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
