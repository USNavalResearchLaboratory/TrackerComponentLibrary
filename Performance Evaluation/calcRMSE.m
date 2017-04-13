function val=calcRMSE(xTrue,xEst)
%%CALCRMSE Compute the scalar root-mean-squared error (RMSE) of estimates
%          compared to true values.
%
%INPUTS: xTrue The truth data. If the truth is the same for all samples of
%              the estimate, then this is an xDimX1 vector. Otherwise, this
%              is an xDimXN matrix where N is the number of samples in the
%              batch.
%         xEst An xDim X N set of estimates.
%
%OUTPUTS: val The scalar RMSE value.
%
%The RMSE is discussed with relation to alternative error measured in [1]
%and is given in Equation 1 of [1].
%
%EXAMPLE: 
%Here, we show that the RMSE of samples of a Gaussian random variable is
%related to the square root of the trace of a covariance matrix.
% R=[28,   4, 10;
%     4,  22, 16;
%     10, 16, 16];%The covariance matrix.
% xTrue=[10;-20;30];
% numRuns=100000;
% xEst=GaussianD.rand(numRuns,xTrue,R);
% rootTrace=sqrt(trace(R))
% val=calcRMSE(xTrue,xEst)
%One will see that the root-trace and the sample RMSE are close.
%
%REFERENCES:
%[1] X. R. Li and Z. Zhao, "Measures of performance for evaluation of
%    estimators and filters," in Proceedings of SPIE: Conference on Signal
%    and Data processing of Small Targets, vol. 4473, San Diego, CA, 29
%    Jul. 2001, pp. 530-541.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numEst=size(xEst,2);

if(size(xTrue,2)==1)
    xTrue=repmat(xTrue,[1,numEst]);
end

val=sqrt(sum(sum((xEst-xTrue).^2,2)/numEst));

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
