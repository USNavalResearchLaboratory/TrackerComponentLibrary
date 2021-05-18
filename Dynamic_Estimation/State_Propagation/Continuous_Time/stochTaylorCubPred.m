function [xPred,PPred]=stochTaylorCubPred(xPrev,PPrev,aFun,BFun,t,deltaT,numSteps,isStrong,algorithm,simplified,xi,w)
%%STOCHTAYLORCUBPRED Given a (multivariate) Gaussian prior distribution
%           with mean xPrev and covariance matrix PPrev under a continuous
%           time dynamic model under It� calculus of the form:
%           dx=aFun(x,t)dt+BFun(x,t)*dw
%           where dw is the differential of a Wiener process and aFun and
%           BFun are drift and diffusion functions, approximate the mean
%           and covariance matrix of the distribution predicted forward
%           deltaT. This is done by computing the moments of a strong or
%           weak stochastic Ito-Taylor expansion, possibly taking more than
%           one step (using the previous expansion as the prior to the next
%           step).
%
%INPUTS: xPrev The dX1 prior target state.
%        PPrev The dXd prior state covariance matrix.
%         aFun A function handle for the drift function such that
%              [aCur,papy,p2apypy,papt]=aFun(y,t) where aCur is the value
%              of the drift function, papx(:,k) is the partial derivative
%              of aCur with respect to y(k), p2apxpx(:,k1,k2) is the second
%              partial derivative of aCur with respect to y(k1) and y(k2),
%              and papt is the partial derivative of aCur with respect to
%              t. Not all of the algorithms require all of the possible
%              outputs of a. See the descriptions of the algorithms to see
%              which outputs need to be implemented.
%         BFun A function handle for the diffusion function ssuch that
%              [BCur,pBpy,p2Bpypy,pBpt]=BFun(y,t). The outputs are related
%              to derivatives of the diffusion function in the same way the
%              outputs of aFun are derivatives of the drift function. Not
%              all of the algorithms require all of the possible outputs of
%              a. See the descriptions of the algorithms to see which
%              outputs need to be implemented.
%            t The time prior to prediction.
%       deltaT The time increment over which the prediction is taken.
%     numSteps The number of steps to take for the prediction. The default
%              if omitted or an empty matrix is passed is 1.
%     isStrong A parameter indicating whether a strong It�-Taylor expansion
%              step should be take (true) or a weak one (false). The
%              default if omitted or an empty matrix is passed is true.
%    algorithm A parameter specifying the algorithm to use. The meaning of
%              this parameter varies depending on the value of isStrong.
%              The default in all instances if omitted or an empty matrix
%              is passed is 0. If isStrong is true, then the algorithms
%              are:
%              0 Use the Euler-Maruyama method from Equation 2.4 in Chapter
%                10.2 of [1]. This is an order 0.5 method.
%              1 Use the Milstein scheme (order 1.0) for scalar noise from
%                Equation 3.2 of Chapter 10.3 of [1]. This requires pBpy
%                and that m=1.
%              2 Use the Milstein scheme (order 1.0) for general
%                multivariate noise from Equation 3.3 of Chapter 10.3 of
%                [1]. This requires pBpy.
%              3 Use the Milstein scheme (order 1.0) for diagonal noise
%                from Equation 3.12 of Chapter 10.3 of [1]. This requires
%                pBpy and that d=m.
%              4 Use the Milstein scheme (order 1.0) for commutative noise.
%                Equation 3.16 of Chapter 10.3 of [1] provides the solution
%                under Stratonovich calculus. However, one can modify
%                Equation 3.3 in the same manner to obtain a solution under
%                It� calculus, which is what is done here. This requires
%                pBpy and that d=m. Given diagonal noise, the result should
%                be the same as algorithm 3.
%              5 Use the strong order 1.5 Taylor method for autonomous
%                scalar problems from Equation 4.1 of Chapter 10.4 of [1].
%                This requires pBpy,p2Bpypy,papy, p2apypy, that d=m=1, and
%                that a and B don't depend on t.
%              6 Use the strong order 1.5 Taylor method for additive noise
%                from Equation 4.10 of Chapter 10.4 of [1]. This requires
%                papy, p2apypy, papt, and pBpt and that B not depend on y.
%              On the other hand, the possible values of algorithm if
%              isStrong=false are:
%              0 Use the order 1.0 simplified weak Euler scheme from
%                Equation 1.2 in Chapter 14.1 of [1] for the simplified
%                scheme; the unsimplified scheme is just the Euler-Maruyama
%                method and both have the same first and second moments.
%              1 Use the order 2.0 weak Taylor scheme for autonomous scalar
%                problems from Equation 2.1 in Chapter 14.2 of [1] or
%                Equation 2.2 for the simplified scheme. This requires
%                pBpy, p2Bpypy, papy, p2apypy and that d=m=1.
%              2 Use the order 2.0 weak Taylor scheme for general problems
%                from Equation 2.6 in Chapter 14.2 of [1] if simplified=0
%                or Equation 2.7 for the simplified scheme if simplified=1.
%                This requires pBpy, p2Bpypy, papy, p2apypy, papt, and
%                pBpt.
%              3 Use the order 2.0 weak Taylor scheme for problems with
%                scalar noise from Equation 2.5 in Chapter 14.2 of [1].
%                This requires pBpy, p2Bpypy, papy, p2apypy, papt, pBpt and
%                that m=1.
%   simplified This parameter is only used if isStrong=false. This is a
%              boolean parameter indicating whether the simplified version
%              of the algorithm should be used. Except for algorithm=0,3,
%              this can affect the results.
%           xi An dXnumCubPoints matrix of cubature points. If this and
%              the next parameter are omitted or empty matrices are passed,
%              then fifthOrderCubPoints(d) is used. It is suggested that xi
%              and w be provided to avoid needless recomputation of the
%              cubature points.
%            w A numCubPoints X 1 vector of the weights associated with the
%              cubature points.
%
%OUTPUTS: xPred The dX1 predicted (mean) state.
%         PPred The covariance matrix associated with xPred.
%
%The notion of using cubature integration to solve this problem is
%described in [2], though only an expansion for the special case with
%autonomous additive noise is considered. The implementation here is a
%generalization to many different Ito-Taylor expansions given in [1].
%Specific derivations are given in [3].
%
%EXAMPLE:
%We demonstrate the consistency of mean and covariance produced by the
%algorithms with respect to a multivariate Black-Scholes model (geometric
%Brownian).
% isStrong=false;%Use a weak approximation.
% simplified=0;
% numSteps=10;
% deltaT=1;
% xPrev=[1;2;3;1];
% SPrev=diag([4;1;2;0.5]);
% PPrev=SPrev*SPrev';
% a=[0.9;
%    1.7;
%    1.3;
%    0.1];
% D=[1.6154,  0.0284;
%    0.1034,  0.4361;
%    0.9386,  0.0641;
%    1.1955,  0.4186];%d=4,m=2.
% [xPredT,PPredT]=BlackScholesPredGaussPrior(xPrev,PPrev,a,D,deltaT);
% 
% t=0;
% aFun=@(x,t)aGeoBrownian(x,a);
% BFun=@(x,t)DGeoBrownian(x,D);
% 
% algorithm=0;
% [xPred,PPred]=stochTaylorCubPred(xPrev,PPrev,aFun,BFun,t,deltaT,numSteps,isStrong,algorithm,simplified);
% max(max(abs((xPred-xPredT)./xPredT)))
% max(max(abs((PPred-PPredT)./PPredT)))
% 
% algorithm=2;
% [xPred,PPred]=stochTaylorCubPred(xPrev,PPrev,aFun,BFun,t,deltaT,numSteps,isStrong,algorithm,simplified);
% max(max(abs((xPred-xPredT)./xPredT)))
% max(max(abs((PPred-PPredT)./PPredT)))
%One will see that algorithm 2 (an order 2.0 weak Taylor scheme) performs
%better than algorithms 0 (Euler-Maruyama method). The maximum relative
%error in the mean and covariance matrix for Euler's method should about
%0.1219 and 0.5056 for algorithm 0. However, they are only 0.0072 and
%0.0819 for algorithm 2.
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%[2] I. Arasaratnam, S. Haykin, and T. R. Hurd, "Cubature Kalman filtering
%    for continous-discrete systems: Theory and simulations," IEEE
%    Transactions on Signal Processing, vol. 58, no. 10, pp. 4977-4993,
%    Oct. 2010.
%[3] D. F. Crouse, "Ito-Taylor expansion moments for continuous-time state
%    propagation," NRL Memo, 2019.
%
%December 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xPrev,1);
if(nargin<11||isempty(xi))
    [xi,w]=fifthOrderCubPoints(xDim);
end

if(nargin<10||isempty(simplified))
    simplified=0; 
end

if(nargin<9||isempty(algorithm))
    algorithm=0;
end

if(nargin<8||isempty(isStrong))
    isStrong=true;
end

if(nargin<7||isempty(numSteps))
    numSteps=1;
end

deltaT=deltaT/numSteps;
x=xPrev;
P=PPrev;
for curStep=1:numSteps
    [x,P]=singleStochTaylorCubPredStep(x,P,aFun,BFun,t,deltaT,isStrong,algorithm,simplified,xi,w);
    t=t+deltaT;
end
xPred=x;
PPred=P;

end

function [xPred,PPred]=singleStochTaylorCubPredStep(xPrev,PPrev,aFun,BFun,t,deltaT,isStrong,algorithm,simplified,xi,w)
xDim=size(xPrev,1);

numCubPoints=length(w);

SPrev=chol(PPrev,'lower');
xiCur=transformCubPoints(xi,xPrev,SPrev);

SigmaMean=zeros(xDim,xDim);
if(isStrong)
    if(algorithm==0)
        for k=1:numCubPoints
            xCur=xiCur(:,k);
            aCur=aFun(xCur,t);
            BCur=BFun(xCur,t);

            [xPredCur,Sigma]=strongTaylorStepMeanCov(xCur,aCur,BCur,deltaT,algorithm);
            xiCur(:,k)=xPredCur;
            SigmaMean=SigmaMean+Sigma*w(k);
        end
    elseif(algorithm<=4)
        for k=1:numCubPoints
            xCur=xiCur(:,k);
            aCur=aFun(xCur,t);
            [BCur,pBpy]=BFun(xCur,t);

            [xPredCur,Sigma]=strongTaylorStepMeanCov(xCur,aCur,BCur,deltaT,algorithm,pBpy);
            xiCur(:,k)=xPredCur;
            SigmaMean=SigmaMean+Sigma*w(k);
        end
    elseif(algorithm==5)
        for k=1:numCubPoints
            xCur=xiCur(:,k);
            [aCur,papy,p2apypy]=aFun(xCur,t);
            [BCur,pBpy,p2Bpypy]=BFun(xCur,t);
            
            [xPredCur,Sigma]=strongTaylorStepMeanCov(xCur,aCur,BCur,deltaT,algorithm,pBpy,p2Bpypy,papy,p2apypy);
            xiCur(:,k)=xPredCur;
            SigmaMean=SigmaMean+Sigma*w(k);
        end
    elseif(algorithm==6)%Additive noise, not necessarily time-invariant.
        for k=1:numCubPoints
            xCur=xiCur(:,k);
            [aCur,papy,p2apypy,papt]=aFun(xCur,t);
            [BCur,~,~,pBpt]=BFun(xCur,t);
            
            [xPredCur,Sigma]=strongTaylorStepMeanCov(xCur,aCur,BCur,deltaT,algorithm,[],[],papy,p2apypy,papt,pBpt);
            xiCur(:,k)=xPredCur;
            SigmaMean=SigmaMean+Sigma*w(k);
        end
    else
        error('Unknown algorithm specified.')
    end
else%Use a weak stochastic Taylor method.
    if(algorithm==0)
        for k=1:numCubPoints
            xCur=xiCur(:,k);
            aCur=aFun(xCur,t);
            BCur=BFun(xCur,t);

            [xPredCur,Sigma]=weakTaylorStepMeanCov(xCur,aCur,BCur,deltaT,algorithm,simplified);
            xiCur(:,k)=xPredCur;
            SigmaMean=SigmaMean+Sigma*w(k);
        end
    elseif(algorithm==1)
        for k=1:numCubPoints
            xCur=xiCur(:,k);
            [aCur,papy,p2apypy]=aFun(xCur,t);
            [BCur,pBpy,p2Bpypy]=BFun(xCur,t);

            [xPredCur,Sigma]=weakTaylorStepMeanCov(xCur,aCur,BCur,deltaT,algorithm,simplified,pBpy,p2Bpypy,papy,p2apypy);
            xiCur(:,k)=xPredCur;
            SigmaMean=SigmaMean+Sigma*w(k);
        end
    elseif(algorithm==2||algorithm==3)
        for k=1:numCubPoints
            xCur=xiCur(:,k);
            [aCur,papy,p2apypy,papt]=aFun(xCur,t);
            [BCur,pBpy,p2Bpypy,pBpt]=BFun(xCur,t);

            [xPredCur,Sigma]=weakTaylorStepMeanCov(xCur,aCur,BCur,deltaT,algorithm,simplified,pBpy,p2Bpypy,papy,p2apypy,papt,pBpt);
            xiCur(:,k)=xPredCur;
            SigmaMean=SigmaMean+Sigma*w(k);
        end
    else
        error('Unknown algorithm specified.')
    end
end

[xPred,meanCov]=calcMixtureMoments(xiCur,w);

%From the Law of Total Covariance.
PPred=SigmaMean+meanCov;
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
