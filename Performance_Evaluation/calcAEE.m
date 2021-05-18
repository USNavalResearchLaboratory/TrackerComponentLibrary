function val=calcAEE(xTrue,xEst,is3D)
%%CALCAEE Compute the scalar average Euclidean error (AEE). In [1], it is
%         argued that the AEE is a better error metric to use than the
%         root-mean-squared error (RMSE). The AEE is essentially the
%         average magnitude of the error whereas the RMSE is the square
%         root of the average squared magnitude of the error.
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
%         is3D An optional indicating that xEst is 3D. This is only used if
%              xEst is a matrix. In such an instance, there is an ambiguity
%              whether xEst is truly 2D, so N=1, or whether numSamples=1
%              and xEst is 3D with the third dimension being 1. If this
%              parameter is omitted or an empty matrix is passed, it is
%              assumed that N=1.
%
%OUTPUTS: val The 1XN set of scalar AEE values.
%
%The AEE is given in Equation 2 of [1].
%
%EXAMPLE:
%Whereas the RMSE of a Gaussian random variables will be near the square
%root of the trace of the covariance matrix, the AEE is different, as is
%shown here.
% R=[28,   4, 10;
%     4,  22, 16;
%     10, 16, 16];%The covariance matrix.
% xTrue=[10;-20;30];
% numRuns=100000;
% xEst=GaussianD.rand(numRuns,xTrue,R);
% rootTrace=sqrt(trace(R))
% val=calcAEE(xTrue,xEst)
%One will see that the AEE is less than the square root of the trace of the
%covariance matrix.
%
%REFERENCES:
%[1] X. R. Li and Z. Zhao, "Measures of performance for evaluation of
%    estimators and filters," in Proceedings of SPIE: Conference on Signal
%    and Data processing of Small Targets, vol. 4473, San Diego, CA, 29
%    Jul. 2001, pp. 530-541.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(is3D))
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
end

val=zeros(1,N);
for k=1:N
    val(k)=sum(sqrt(sum((xEst(:,k,:)-xTrue(:,k,:)).^2,1)))/numSamples;
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
