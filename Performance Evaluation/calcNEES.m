function NEES=calcNEES(xTrue,xEst,PEst,is3D)
%CALCNEES Given a batch of estimates and associated true values, compute
%         the normalized estimation error squared (NEES). This version is
%         normalized such that the mean of a consistent estimator,
%         regardless of the dimensionality of the state, should be one.
%
%INPUTS: xTrue The truth data. This is either an xDimXNumSamples matrix or
%              an xDimXNXnumSamples matrix. The latter formulation is
%              useful when the MSE over multiple Monte Carlo runs of an
%              entire length-N track is desired. In the first formulation,
%              we take N=1. Alternatively, if the same true value is used
%              for all numSamples, then xTrue can just be an xDimXN matrix.
%              Alternatively, if xTrue is the same for all numSamples and
%              all N, then just an xDimX1 matrix can be passed. N and
%              numSamples are inferred from xEst.
%         xEst An xDimXnumSamples set of estimates or an xDimXNXnumSamples
%              set of estimates (if values at N times are desired).
%         PEst An xDimXxDimXnumSamples or an xDimXxDimXNXnumSamples set of
%              covariance matrices associated with the estimates. If the
%              covariance matrix is supposed to be the same for all
%              estimates, then PEst is just xDimXxDim. These should not be
%              singular.
%         is3D An optional indicating that xEst is 3D. This is only used if
%              xEst is a matrix. In such an instance, there is an ambiguity
%              whether xEst is truly 2D, so N=1, or whether numSamples=1
%              and xEst is 3D with the third dimension being 1. If this
%              parameter is omitted or an empty matrix is passed, it is
%              assumed that N=1.
%
%OUTPUTS: NEES The normalized estimation error squared of the estimates.
%
%The NEES is a measure of how well the provided by an estimator matches the
%actual accuracy of the estimates. The concept of the NEES is discussed in 
%[1].
%
%EXAMPLE:
%Here, we demonstrate that the true covariance matrix is consistent with
%the samples of the distributon.
% R=[28,   4, 10;
%     4,  22, 16;
%     10, 16, 16];%The covariance matrix.
% xTrue=[10;-20;30];
% numRuns=100000;
% xEst=GaussianD.rand(numRuns,xTrue,R);
% NEESVal=calcNEES(xTrue,xEst,R)
%One will see that NEESVal is close to 1.
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(is3D))
    is3D=false; 
end

xDim=size(xEst,1);
if(ismatrix(xEst)&&is3D==false)
    N=1;
    numSamples=size(xEst,2);
    xEst=reshape(xEst,[xDim,1,numSamples]);
    
    if(size(xTrue,2)==1)
        %If the true values are the same for all samples.
        xTrue=repmat(xTrue,[1,1,numSamples]);
    else
        xTrue=reshape(xTrue,[xDim,1,numSamples]);
    end
    
    if(ismatrix(PEst))
        PEst=repmat(PEst,[1,1,N,numSamples]);
    elseif(ndims(PEst)==3) 
        PEst=reshape(PEst,[xDim,xDim,1,numSamples]);
    end
else
    N=size(xEst,2);
    numSamples=size(xEst,3);
    
    if(ismatrix(xTrue))
        if(size(xTrue,2)==1)
            %If the true values are the same for all samples and for all N.
            xTrue=repmat(xTrue,[1,N,numSamples]);
        else        
            %If the true values are the same for all samples.
            xTrue=repmat(xTrue,[1,1,numSamples]);
        end
    end
    
    if(ismatrix(PEst))
        PEst=repmat(PEst,[1,1,N,numSamples]);
    elseif(ndims(PEst)==3) 
        PEst=repmat(PEst,[1,1,1,numSamples]);
    end
end

NEES=zeros(1,N);
for k=1:N
    for curSamp=1:numSamples
        diff=xTrue(:,k,curSamp)-xEst(:,k,curSamp);

        NEES(k)=NEES(k)+invSymQuadForm(diff,PEst(:,:,k,curSamp));%=NEES+diff'*inv(PEst(:,:,k,curSamp))*diff;
    end
end
NEES=NEES/(xDim*numSamples);

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
