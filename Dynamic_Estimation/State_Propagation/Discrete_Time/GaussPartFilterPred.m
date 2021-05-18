function [xPred,PPred]=GaussPartFilterPred(xPrev,PPrev,transPDFSampFun,param4,wPart)
%%GAUSSPARTICLEFILTERPRED Perform the state prediction step in the Gaussian
%               particle filter. This can be used with an arbitrary dynamic
%               model, but it requires being able to draw random samples
%               from the propagated PDF conditioned on the prior state
%               value.
%
%INPUTS: xPrev The xDim X 1 state estimate at the previous time-step.
%        PPrev The xDim X xDim state covariance matrix at the previous
%              time-step.
% transPDFSampFun A function handle such that transPDFSampFun(xPrevVal)
%              draws a single sample from the propagated PDF conditioned on
%              the value of the prior PDF. This requires being able to
%              determine the conditional PDF of the propagated state. As
%              described below, in lieu of a better method, rejection
%              sampling can be used for the sampling function if the
%              transition PDF is known.
% param4,wPart These parameters determine the specifics of the algorithm
%              used. If param4 is a number and wPart is omitted or is an
%              empty matrix, then param4=numParticles is the number of
%              particles to use in the update step. If param4 and wPart are
%              both omitted or empty matrices are passed, then 1e4
%              particles is used. On the other hand, if param4=xiPart and
%              wPart is provided, then xiPart and wPart are the particles
%              returned by GaussPartFilterUpdate and the alternative time
%              update step is used.
%
%OUTPUTS: xPred The xDim X 1 predicted state estimate.
%         PPred The xDim X xDim predicted state covariance estimate.
%
%The algorithm is that given in Table I of [1]. There are two variants,
%depending on whether the particles from a previous measurement update step
%are available.
%
%Obtaining the transition PDF can be very difficult. On the other hand,
%obtaining a sampling function for the PDF need not be as difficult. Often,
%a discrete-time dynamic model is given in terms of a method that lends
%itself to sampling. For example, a model with additive noise is
%xPopagated= f(x)+w. Thus, just put x into f and add a random sample of w.
%For more general problems, one can use, for example, rejection sampling,
%which is described in Chapter 5.2 of [2]. In such an instance, one has a
%PDF g(x) that is easy to sample and one wishes to sample a PDF f(x) that
%is difficult to sample. g(x) and f(x) must have the same
%support. The rejection method is
%1) Determine a value c such that f(x)/g(x)<=c for all x.
%2) Generate a random variable y having density g.
%3) Generate a uniform(0,1) random variable u.
%4) If u<=f(y)/(c*g(y)) then accept y as the sample. otherwise, go to 2.
%
%REFERENCES:
%[1] J. H. Kotecha and P. M. Djuric, "Gaussian particle filtering,"
%    IEEE Transactions on Signal Processing, vol. 51, no. 10, pp.
%    2592-2601, Oct. 2003.
%[2] S. M. Ross, Simulation, 4th ed. Amsterdam: Academic Press, 2006.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xPrev,1);

if(nargin<4||isempty(param4))
%If no particles are given and no number of particles is specified, use 1e4
%particles by default.
    param4=1e4;
end

if(nargin<5||isempty(wPart))
    %%Time Update
    %1) Draw samples from the normal distributions.
    %In this case, param4 is the number of samples to draw.
    numParticles=param4;
    
    xiPart=bsxfun(@plus,xPrev,cholSemiDef(PPrev,'lower')*randn(xDim,numParticles));
    
    %2) Sample from the transition PDF conditioned on each sample state.
    for curP=1:numParticles
        xiPart(:,curP)=transPDFSampFun(xiPart(:,curP));
    end
 
    %3) Compute the mean and covariance.
    xPred=mean(xiPart,2);
    PPred=zeros(xDim,xDim);
    for curP=1:numParticles
        diff=xiPart(:,curP)-xPred;
        PPred=PPred+(diff*diff');
    end
    PPred=PPred/numParticles;
else
    %%Time Update Alternate
    %1) We are given the previous weighted samples.
    xiPart=param4;
    numParticles=size(xiPart,2);
    %2) Sample from the transition PDF conditioned on each previous sample.
    for curP=1:numParticles
        xiPart(:,curP)=transPDFSampFun(xiPart(:,curP));
    end

    %3) Compute the mean and covariance. This lacks the 1/M term given in
    %the paper as that appears to be a typo.
    [xPred, PPred]=calcMixtureMoments(xiPart,wPart);
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
