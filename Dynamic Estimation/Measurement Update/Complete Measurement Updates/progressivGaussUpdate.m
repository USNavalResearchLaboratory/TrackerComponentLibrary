function [xUpdate,SUpdate]=progressivGaussUpdate(xPred,SPred,zLike,xi,w,stepSizeBounds)
%%PROGRESSIVEGAUSSUPDATE Perform the measurement update state in a
%               progressive Gaussian filter. In this filter, the state is
%               always approximated by a Gaussian probability density
%               function, but the measurement has an arbitrary continuous
%               conditional density function that is provided. The
%               progressive update performs the update in a series of
%               steps that implement a type of homotopy transformation
%               (with renormalization each step) from the prior to the
%               posterior distribution.
%              
%INPUTS: xPred The xDim X 1 predicted target state.
%        SPred The xDim X xDim lower-triangular square root of the
%              predicted state covariance matrix.
%        zLike A function handle such that zLike(x) returns the conditional
%              likelikoods of the measurement given a matrix of N target
%              states. That is, an xDimXN matrix. The likelihoods can be
%              returned as a 1XN vector or an NX1 vector. It does not
%              matter if they are all off by a constant positive
%              multiplicative factor from their true values.
%           xi An xDim X numCubPoints matrix of cubature points.
%            w A numCubPoints X 1 vector of the weights associated with the
%              cubature points. These should be all greater than zero and
%              sum to one.
% stepSizeBounds An optional 2X1 or 1X2 vector giving bounds on the
%              adaptive step size used for the progressive measurement
%              update. These values are between 0 and 1, representing the
%              minimum/maximum change in the homotopy parameter each step.
%              If this parameter is omitted or an empty matrix is passed,
%              the default value of [1/100; 0.5] is used.
%
%OUTPUTS: xUpdate The xDim X 1 updated state vector. If zLike ever returned
%                 all zeros when evaluating a particular set of cubature
%                 points, then an empty matrix will be returned to indicate
%                 a failure of the algorithm.
%         SUpdate The updated xDim X xDim lower-triangular square root
%                 state covariance matrix, unless an algorithmic failure
%                 occurs, in which case an empty matrix is returned.
%
%The algorithm is essentially that of [1], which is a generalization and
%improvement of [2]. but it has been modified to allow the weights of the
%samples to be non-uniform. This just changes Equation 11 in [1] from an
%unweighted mixture into a weighted mixture. As a result, the ratio leading
%to the adaptive step size solution in Equation 14 changes as each f_l
%term in the ratio must be multiplied by the weight that is no in Equation
%11. This means that Equation 14 becomes
%Delta=log(R)/(log(f^L(x_s))+log(w_s)-log(f^L(x_l))-log(w_l))
%where w_s and w_l are the weights associated with x_s and x_l.
%
%REFERENCES:
%[1] J. Steinbring and U. D. Hanebeck, "Progressive Gaussian filtering
%    using explicit likelihoods," in Proceedings of the 17th International
%    Conference on Information Fusion, Salamanca, Spain, 7-10 Jun. 2014.
%[2] U. D. Hanebeck, "PGF 42: Progressive Gaussian filtering with a twist,"
%    in Proceedings of the 16th International Conference on Information
%    Fusion, Istanbul, Turkey, 9-12 Jul. 2013, pp. 1103-1110.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(stepSizeBounds))
    stepSizeBounds=[1/100;0.5];
end

%The number of cubature points.
M=size(xi,2);

xCur=xPred;
SCur=SPred;

logW=log(w);

gammaVal=0;
while(gammaVal<1)
    %Get the samples for the PDF at its current step in the progression.
    xiCur=transformCubPoints(xi,xCur,SCur);
    
    %The log-likelihoods of the transformed points.
    l=log(zLike(xiCur));
    
    %The logarithm of the unnormalized posterior w values (without Delta),
    %including the effects of the weights added to Equation 11.
    logWTilde=l(:)+logW;
    
    %The isfinite command is to get rid of zero likelihoods, which would
    %have led to -Inf as the log-likelihood, as the text suggests to avoid
    %numerical problems.
    logWMin=min(logWTilde(isfinite(logWTilde)));
    logWMax=max(logWTilde);
    
    %If zLike returned all zeros, then the update failed and we will return
    %empty matrices to indicate that.
    if(isempty(logWMin))
        xUpdate=[];
        SUpdate=[];
        return
    end
    
    %The step size
    Delta=-log(M)/(logWMin-logWMax);
    
    %Make sure that the step size is not too small or too large.
    Delta=max(stepSizeBounds(1),Delta);
    Delta=min(Delta,stepSizeBounds(2));
    
    %If we are about to step past the end (gammaVal=1), then 
    if(gammaVal+Delta>1)
        Delta=1-gammaVal;
    end
    
    %The exponential in Line 11 of Algorithm 1 in [1] would NOT apply to
    %the weight values. Thus, we cannot use the logWTilde term above when
    %getting the actual unnormalized posterior weights with the effects of
    %Delta included.
    fLDeltaw=exp(l(:)*Delta+logW);
    w=fLDeltaw/sum(fLDeltaw);
    
    %PCur is the covariance matrix, refferred to as C in the paper.
    [xCur, PCur]=calcMixtureMoments(xiCur,w);
    SCur=cholSemiDef(PCur,'lower');
    gammaVal=gammaVal+Delta;
end

xUpdate=xCur;
SUpdate=SCur;

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
