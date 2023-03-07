function [p,exitInfo]=fuseIndepProbEsts(pEsts,algorithm,maxIter,epsilon,alpha)
%%FUSEINDEPPROBESTS Given a set of independent discrete probability
%               estimates, fuse the the probability estimates. Each of the
%               estimates is assumed to be equally reliable. The algorithm
%               selected determines the cost function used for fusion.
%
%INPUTS: pEsts A numProbsXnumEsts matrix such that pEsts(:,i) is a single
%              set of probability estimates. It is such that
%              sum(pEsts(:,i))=1 and all pEsts>=0.
%    algorithm The algorithm to use to determine the offsets. Possible
%              values are:
%              0 (The default if omitted or an empty matrix is passed) Use
%                a quadratic programming solution. The cost function is
%                thus the sum of squares of offsets that need to be made to
%                the individual pEst values such that pEst(:,i)=pEst(:,j)
%                for all i and j.
%              1 As opposed to having the cost function be a sum of squared
%                offsets as in 0, this is the sum of absolute offsets. This
%                is a linear programming solution.
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
%OUTPUTS: p An xDimX1 vector holding the fused probabilty estimates. Note
%           that sum(p)=1.
%   exitInfo If algorithm=0, this is the exitCode output of the
%           convexQuadProg function. If algorithm=1, then this is the info
%           output of the linProgInteriorPoint function.
%
%Both algorithms implemented here are given in Section 2.3 of [1]. The work
%in the rest of [1] is more general than the solution used in this
%function.
%
%EXAMPLE:
%Here, we use one of the many variants of an example of Zadeh's paradox in
%Dempster-Schafer theory. As in [2], consider one expert's opinion that a
%disease has a 99% chance of being meningitis and 1% that it is a
%concussion. A second doctor says that it has a 99% chance of being a tumor
%and a 1% chance of being a concussion. We fuse these estimates as:
% pEsts=[0.99, 0.00;
%        0.01, 0.01;
%        0.00, 0.99];      
% pFus=fuseIndepProbEsts(pEsts)
%The result is pFus=[0.495;0.01;0.495], which is reasonable. The
%controversial Dempster-Schafer result would assing 100% probability to the
%concussion hypothesis.
%
%REFERENCES:
%[1] S. Tsavachidis, "Eliminating incoherence from subjective estimates of
%    chance," Master's thesis, Rice University, Houston, TX, Feb. 2003.
%[2] R. Haenni, "Shedding light on Zadeh's criticism of Dempster's rule of
%    combination," in Proceedings of the 7th International Conference on
%    Information Fusion, Philadelphia, PA, 25-28 Jul. 2005.
%
%November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5)
    alpha=[];
end

if(nargin<4)
    epsilon=[];
end

if(nargin<3)
    maxIter=[];
end

if(nargin<2||isempty(algorithm))
    algorithm=0;
end

numX=size(pEsts,1);
N=size(pEsts,2);

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
