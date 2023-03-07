function [Q,q,exitCode]=ellipsOfEllipsSum(Q1,Q2,q1,q2,algorithm,AbsErr,RelErr,maxIter)
%%ELLIPSEOFELLIPSESUM Consider ellipses or ellipsoids defined by equations
%               of the form: (x-q)'*inv(Q)*(x-q)<=1, where Q is a dXd
%               positive-definite matrix q is the dX1 center of the
%               ellipse, and x is a point on or in the ellipsoid. We want
%               to find the parameters of the minimum volume ellipsoid that
%               encompasses the Minkowski sum of the ellipsoids. The
%               Minkowski sum of the ellipsoids is the shape containing
%               points x such that x=x1+x2, where x1 is a point in the
%               first ellipsoids and x2 is a point in the second ellipsoid.
%               A name for the minimum volume ellipsoid that encompasses
%               the Minkowsi sum of the two ellipsoids is the Löwner-John
%               ellipsoid. This function can alternatively return the
%               parameters for the ellipsoid with the minimum squared semi-
%               axis lengths that encloses the Minkowski sum.
%
%INPUTS: Q1, Q2 The dXd positive definite matrices defining the shapes of
%               the first and second ellipses being added.
%         q1,q2 The dX1 centers of the ellipses if one wishes to get q, the
%               center of the sum ellipse, on the output.If q is not
%               desired, these can be omitted or empty matrices can be
%               passed and an empty matrix will be returned for q.
%     algorithm An optional parameter indicating which algorithm to use.
%               Possible values are:
%               0 (The default if omitted or an empty matrix is passed) Use
%                 the iterative algorithm in [1] to find the minimum volume
%                 enclosing ellipsoid.
%               1 Use the algorithm in [1] for the 2D case, where a
%                 parameter beta is bounded. The desired beta value is
%                 found using a line search. An error is raised if d~=2.
%               2 Use the solution in [1] for the enclosing ellipsoid with
%                 the minimum squares semi-axis lengths. This is a
%                 non-iterative algorithm.
%        AbsErr A convergence parameter for algorithms 0 and 1. The sum
%               ellipsoid shape Q is parameterized by a scalar beta. This
%               is the absolute difference between iterative steps where
%               convergence is declared. The default if omitted or an empty
%               matrix is passed is 2^4*eps().
%        RelErr The relative change in the beta value for convergence used
%               in algorithm 0. The default if omitted or an empty
%               matrix is passed is 2^4*eps().
%       maxIter The maximum number of iterations for algorithm 0 or 1. If
%               this parameter is omitted or an empty matrix passed, the
%               default is 100.
%
%OUTPUTS: Q The dXd shape matrix of the ellipsoid fitting the sum.
%         q The dX1 center of the ellipsoid of the sum.
%  exitCode A parameter idicating how the algorithm exited. If algorithm=2,
%           then this is always 0. If algorithm=1, then this is the exit
%           code of the goldenSectionSearch function. If algorithm=0, then
%           possible values are:
%           0 The stepsize become too small accoridng to AbsTol or RelTol.
%           1 The maximum number of iterations elapsed.
%
%Extended target racking algorithms, such as in [2], represent
%extended/group targets as an ellipsoid. The center of the ellipsoid is
%distributed Gaussian and the shape matrix as inverse Wishart. Confidence
%intervals for the shape and the position can be each expressed as
%ellipsoids. An ad-hoc way to combine the confidence intervals is just to
%take the Minkowski sum of the ellipsoids.
%
%EXAMPLE:
%Here, two ellipses sharing the same origin are added (the result has the
%same origin too).
% Q1=[4,0;
%     0,1];
% Q2=[1,1.5;
%     1.5,4];
% figure(1)
% clf
% hold on
% drawEllipse([0;0],inv(Q1),1,'linewidth',2)
% drawEllipse([0;0],inv(Q2),1,'linewidth',2)
% Q=ellipsOfEllipsSum(Q1,Q2);
% drawEllipse([0;0],inv(Q),1,'--r','linewidth',2)
%The solid lines are the original ellipses and the dashed line is the sum.
%
%REFERENCES:
%[1] A. Halder, "On the parameterized computation of minimum volume outer
%    ellipsoid of Minkowski sum of ellipsoids," in Proceedings of the IEEE
%    Conference on Decision and Control, Miami Beach, FL, 17-19 Dec. 2018,
%    pp. 4040-4045.
%[2] M. Feldmann, D. Fränken, and W. Koch, "Tracking of extended objects
%    and group targets using random matrices," IEEE Transactions on Signal
%    Processing, vol. 59, no. 4, pp. 1409-1420, Apr. 2011.
%
%April 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<8||isempty(maxIter))
    maxIter=100;
end

if(nargin<7||isempty(RelErr))
    RelErr=2^4*eps();
end

if(nargin<6||isempty(AbsErr))
    AbsErr=2^4*eps();
end

if(nargin<5||isempty(algorithm))
    algorithm=0; 
end

if(nargin<3||isempty(q1))
    q=[];
else
    q=q1+q2;
end

switch(algorithm)
    case 0%Fixed point iteration.
        lambda=eig(Q1\Q2);
        
        beta=1;%Initial estimate.
        exitCode=1;
        for curIter=1:maxIter
            betaPrev=beta;
            
            lRat=1./(1+beta*lambda);
            
            beta=sqrt(sum(lRat)/sum(lambda.*lRat));
            
            err=abs(betaPrev-beta);
            if(err<=AbsErr||err<=RelErr*abs(beta))
                exitCode=0;
                break;
            end
        end
        Q=(1+1/beta)*Q1+(1+beta)*Q2;
    case 1%2D, bounded search.
        if(size(Q1,1)~=2)
            error('This algorithm is only for 2D matrices.')
        end
        
        lambda=eig(Q1\Q2);
        lambda1=lambda(1);
        lambda2=lambda(2);
        lambdaSum=lambda1+lambda2;
        lambdaProd=lambda1*lambda2;

        betaLB=(sqrt(lambdaSum*(lambdaSum+6*lambdaProd))-lambdaSum)/(6*lambdaProd);
        betaUB=min((lambda1+sqrt(lambda1*(lambda1+8)))/(2*lambda1),(lambda2+sqrt(lambda2*(lambda2+8)))/(2*lambda2));

        f=@(betaVal)costFun(betaVal,Q1,Q2);

        [betaMin,~,exitCode]=goldenSectionSearch(f,[betaLB;betaUB],AbsErr,maxIter);

        Q=(1+1/betaMin)*Q1+(1+betaMin)*Q2;
    case 2%Use the alternative cost function.
        beta=sqrt(trace(Q1)/trace(Q2));
        Q=(1+1/beta)*Q1+(1+beta)*Q2;
        exitCode=0;
    otherwise
        error('Unknown Algorithm specified.')
end
end

function val=costFun(betaVal,Q1,Q2)
%%COSTFUN The cost function for the 2D line-search algorithm in [1].
%
%REFERENCES:
%[1] A. Halder, "On the parameterized computation of minimum volume outer
%    ellipsoid of Minkowski sum of ellipsoids," in Proceedings of the IEEE
%    Conference on Decision and Control, Miami Beach, FL, 17-19 Dec. 2018,
%    pp. 4040-4045.
%
%April 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

Q=(1+1/betaVal)*Q1+(1+betaVal)*Q2;
val=log(det(Q));

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
