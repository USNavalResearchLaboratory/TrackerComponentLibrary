function [p,exitInfo]=fuseIndepProbEsts(pEsts,w,algorithm,maxIter,epsilon,alpha)
%%FUSEINDEPPROBESTS Given a set of independent discrete probability
%               estimates, fuse the the probability estimates. Each of the
%               estimates is assumed to be equally reliable unless
%               algorithm 2 is used (the default), in which case a
%               weighting vector w can be provided. The algorithm selected
%               determines the cost function used for fusion.
%
%INPUTS: pEsts A numProbsXnumEsts matrix such that pEsts(:,i) is a single
%              set of probability estimates. It is such that
%              sum(pEsts(:,i))=1 and all pEsts>=0.
%            w If algorithm=2, this is a numEstsX1 or 1XnumEsts set of
%              weights (w>0) that list the reliability of each of the
%              estimate sets. So, think of each estimate set being from a
%              particular expert and this is a relative weighting of how
%              accurate that expert is. If omitted or an empty matrix is
%              passed, all uniform weights are used.
%    algorithm The algorithm to use to determine the offsets. Possible
%              values are:
%              0 Use a quadratic programming solution as in Section 2.3 of
%                [1]. The cost function is thus the sum of squares of
%                offsets that need to be made to the individual pEst values
%                such that pEst(:,i)=pEst(:,j)
%                for all i and j. This doesn't use w. This is actually
%                equivalent to algorithm 2 with uniform weights, so it
%                makes more sense to just use algorithm 2.
%              1 As opposed to having the cost function be a sum of squared
%                offsets as in 0, this is the sum of absolute offsets. This
%                is a linear programming solution. This doesn't use w.
%              2 (The default if omitted or an empty matrix is passed) Use
%                the simple explicit solution derived in the comments
%                below. This is the only solution that supports using w. As
%                shown below, this is really just a weighted average of the
%                probabilities. It is optimal in terms of a quadratic cost
%                function, whereas algorithm 1 uses a different, more
%                difficult cost function.
%      maxIter The maximum number of iterations to use in the
%              convexQuadProg function for algorithm 0 or the
%              linProgInteriorPoint function for algorithm 1. If omitted or
%              an empty matrix is passed, the default of the respective
%              function is used.
%   epsilon, alpha These are parameters for the linProgInteriorPoint
%              function and are only used if algorithm=1. If omitted or
%              empty matrices are passed, then the defaults of that
%              function will be used.
%
%OUTPUTS: p An numProbsX1 vector holding the fused probabilty estimates.
%           Note that sum(p)=1.
%  exitInfo If algorithm=0, this is the exitCode output of the
%           convexQuadProg function. If algorithm=1, then this is the info
%           output of the linProgInteriorPoint function. This is not used
%           with algorthm=2 and will just be 0.
%
%Algorithms 0 and 1 that are implemented here are given in Section 2.3 of
%[1]. The work in the rest of [1] is more general than the solution used in
%this function. Algorithm 2 is more general than algorithm 0 and is
%actually explicit. It turns out that the quadratic programming in
%algorithm 0 is equivalent to just a weighted average. 
%
%ALgorithm 2 is derived as follows: Let fij be the probability that expert
%i chooses class j as being correct. w(i)>1 is a relative reliability
%ranking of the ith expert. There are NE experts and NC classes.
%To fuse into a consistent probability distribution covering all of the
%classes, we perform the optimization
%fHat=arg min_{fHat(j) for all j} sum_{i=1}^{N_E} w(i)*sum_{j=1}^{N_C}(fHat(j)-fij)^2
%such that sum(fHat)=1 and 0<=fHat(j)<=1
%This can be rewritten as a convex quadratic programming problem of the
%form
%fHat=arg min_{fHat} fHat'*Q*fHat+a'*fHat+b
%such that A1'*fHat=1 and A2'*fHat>=0 and A3'*fHat>=-1
%where I is an identity matrix and 
%Q=I*sum(w);
%a'=-2*sum_{i=1}^NE w(i)*f(i,:)
%b=sum_{j=1}^NC sum_{i=1}^NE w*fij^2
%A1=1
%A2=I
%A3=-1
%and f(i,:) selects the ith column of the matrix implied by fij.
%
%However, it is not necessary to impose the inequality constraints if all
%0<=fij<=1 for all (i,j). Any fHat(j) that is pushed over 1 necessitates
%one or more other fHat(j) to be pushed under 0 to satify the equality
%constraints. Because this is outside of the bounds of all fij, one would
%achieve a lower cost by clipping the fHat(j) to the region 0<=fHat(j)<=1.
%Thus, as long as we satisfy the equality constraints, we do not need to
%enforce the inequality constraints.
%
%The problem can be solved using Lagrangian relaxation. The Lagrangian is
%L=sum_{i=1}^NE w(i) sum_{j=1}^NC (fHat(j)-fij)^2+lambda*(sum(fHat)-1)
%The partial derivative with respect to the ith element is
%d/dfHat(j)=sum_{i=1}^NE 2*w(i)*(fHat(j)-fij)+lambda
%Setting the derivative equal to zero and solving for fHat(j) one gets
%fHat(j)=(sum_{i=1}^NE 2*w(i)*fij-lambda)/(2*sum(w))
%Substituting into the constraints, we get
%1=sum_{j=1}^NC (sum_{i=1}^NE 2*w(i)*fij-lambda)/(2*sum(w))
%Solving for lambda leads to 
%lambda=sum_{j=1}^NC sum_{i=1}^NE 2*w(i)*fij-2*sum(w)
%Substituting the solution for lambda back into the solution for fHat(j)
%and simplifying we get
%fHat(j)=sum_{i=1}^NE wTilde(i)*fij-sum_{k=1}^NC sum_i=1}^NE wTilde(i)*fik+1
%where
%wTilde=w/sum(w)
%We note that sum_{k=1}^NC sum_i=1}^NE wTilde(i)*fik=1 so that the result
%is nothing but a weighted average of the classifiers.
%fHat(j)=sum_{i=1}^NE wTilde(i)*fij
%That is how we get algorithm 2 without iteration. It also means that
%algorithm 0, which uses the same cost function, is nothing but a very
%inefficient method of taking a weighted average with uniform weights.
%
%EXAMPLE 1:
%Here, we use one of the many variants of an example of Zadeh's paradox in
%Dempster-Schafer theory. As in [2], consider one expert's opinion that a
%disease has a 99% chance of being meningitis and 1% that it is a
%concussion. A second doctor says that it has a 99% chance of being a tumor
%and a 1% chance of being a concussion. We fuse these estimates as:
% pEsts=[0.99, 0.00;
%        0.01, 0.01;
%        0.00, 0.99];      
% pFus1=fuseIndepProbEsts(pEsts,[],1)
% pFus2=fuseIndepProbEsts(pEsts,[],2)
%The result is pFus=[0.495;0.01;0.495], which is reasonable. The
%controversial Dempster-Schafer result would assing 100% probability to the
%concussion hypothesis.
%
%EXAMPLE 2:
%This is the same as example 1, but this time weightings are provided. One
%can see that if the two experts are ranked equally, then one gets the same
%result as before. However, if the first expert gets twice the rating, then
%the result is closer to his result.
% pEsts=[0.99, 0.00;
%        0.01, 0.01;
%        0.00, 0.99];
% w=[1,1];
% pFus=fuseIndepProbEsts(pEsts,w)
% w=[2,1];
% pFus=fuseIndepProbEsts(pEsts,w)
%
%REFERENCES:
%[1] S. Tsavachidis, "Eliminating incoherence from subjective estimates of
%    chance," Master's thesis, Rice University, Houston, TX, Feb. 2003.
%[2] R. Haenni, "Shedding light on Zadeh's criticism of Dempster's rule of
%    combination," in Proceedings of the 7th International Conference on
%    Information Fusion, Philadelphia, PA, 25-28 Jul. 2005.
%
%November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%Updated with Algorithm 2 in January 2025, David. F. Crouse
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6)
    alpha=[];
end

if(nargin<5)
    epsilon=[];
end

if(nargin<4)
    maxIter=[];
end

if(nargin<3||isempty(algorithm))
    algorithm=2;
end

numX=size(pEsts,1);
N=size(pEsts,2);

if(nargin<2||isempty(w))
    if(algorithm==2)
        %Uniform weighting.
        w=ones(N,1);
    end
end

switch(algorithm)
    case 0%Use quadratic programming.
        G=2*N*eye(numX,numX);
        a=-2*sum(pEsts,2);

        numEqConst=1;
        C=[ones(numX,1),eye(numX)];
        b=[1;zeros(numX,1)];

        epsVal=eps();%No need to use anything else.
        [p,~,exitInfo]=convexQuadProg(G,a,C,b,numEqConst,epsVal,maxIter);
    case 1%Use linear programming.
        %For the linear programing solution we have numX x variables that
        %are constrained to sum to 1 and to be >=0. We also have numX*N
        %slack variables. For each slack variable we have two inequality
        %constraints. We want to minimize the sum of the slack variables,
        %which are also constrained to be >=0.

        %The first numX variables are the X variables. Thus, the first
        %constraint is that they all sum to 1.
        numSlack=numX*N;
        %The cost is the sum of the slack variables.
        c=[zeros(numX,1);ones(numSlack,1)];
        %There is only one equality constraint: That the x values sum to 1.
        A=[ones(1,numX),zeros(1,numSlack)];
        b=1;

        ALeq=zeros(numX*2*N,numX+numSlack);
        bLeq=zeros(numX*2*N,1);

        selRow=1:numX;
        selColStart=1:numX;
        selColEnd=(numX+1):(2*numX);
        for curMeas=1:N
            ALeq(selRow,selColStart)=-eye(numX);
            ALeq(selRow,selColEnd)=-eye(numX);
            bLeq(selRow)=pEsts(:,curMeas);

            selRow=selRow+numX;
            ALeq(selRow,selColStart)=eye(numX);
            ALeq(selRow,selColEnd)=-eye(numX);
            bLeq(selRow)=pEsts(:,curMeas);

            selRow=selRow+numX;
            selColEnd=selColEnd+numX;
        end

        maximize=false;
        [~,xOpt,exitInfo]=linProgInteriorPoint(A,b,ALeq,bLeq,c,maximize,maxIter,epsilon,alpha);
        p=xOpt(1:numX);
    case 2%Use the simple explicit solution.
        fij=pEsts';
        wTilde=w/sum(w);
        fijWeighted=bsxfun(@times,wTilde(:),fij);
        p=sum(fijWeighted,1).';

        %p will be a valid set of probabilities unless finite precision
        %limitations affect it. Thus, these lines are just to make sure it
        %isn't slightly invalid and generally these lines won't do
        %anything.
        p=min(1,max(0,p));
        p=p/sum(p);
        
        exitInfo=0;%Not used.
    otherwise
        error('Unknown algorithm specified.')
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
