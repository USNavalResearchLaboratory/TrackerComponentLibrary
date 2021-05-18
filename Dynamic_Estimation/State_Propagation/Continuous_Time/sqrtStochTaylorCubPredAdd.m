function [xPred,SPred]=sqrtStochTaylorCubPredAdd(xPrev,SPrev,a,BVal,deltaT,numSteps,algorithm,xi,w)
%%SQRTSTOCHTAYLORCUBPREDADD Given a Gaussian prior distribution with mean
%           xPrev and lower-triangular square root covariance matrix SPrev
%           under a continuous time, autonomous dynamic model under Itô
%           calculus with additive noise of the form:
%           dx=a(x)dt+BVal*dw
%           where dw is the differential of a Wiener process and BVal is
%           constant matrix, approximate the mean and lower-triangular
%           square root covariance matrix of the distribution predicted
%           forward deltaT. This is done by computing the moments of a
%           stochastic Itô-Taylor expansion, possibly taking more than one
%           step (using the previous expansion as the prior to the next
%           step).
%
%INPUTS: xPrev The dX1 prior target state.
%        SPrev The dXd lower-triangular square root of the state covariance
%              matrix, such that the covariance matrix equals SPrev*SPrev'.
%            a A function handle such that [aCur,papx]=a(x) returns the
%              value of the dX1 drift coefficient (aCur) at x and the dXd
%              matrix of derivatives of the elements of the drift
%              coefficient with respect to the elements of x. papx(:,i) is
%              the derivtive of i with respect to the ith component of x.
%              papx does not need to be returned if algorithm=0.
%         BVal The constant dXm diffusion matrix. m is the dimensionality
%              of the Wiener process.
%       deltaT The time interval over which the prediction should take
%              place.
%     numSteps The number of steps to take for the prediction. The default
%              if omitted or an empty matrix is passed is 1.
%    algorithm A value specifying the type of stochastic Itô Taylor
%              expansion to use for the prediction. Possible values are:
%              0 Use the Euler-Maruyama method (order 0.5).
%              1 Use the order 1.5 Taylor expansion, as used in [1]. When
%                computing the cubature integrals, use the approximation of
%                [1] and assume that the terms involving higher papx do not
%                change rapidly, so papx does not need to be evaluated at
%                every cubature point.
%              2 (The default if omitted or an empty matrix is passed) Use
%                the order 1.5 Taylor expansion, as used in [1] but do
%                compute a new papx for every cubature point.
%           xi An xDim X numCubPoints matrix of cubature points. If this
%              and the next parameter are omitted or empty matrices are
%              passed, then fifthOrderCubPoints(d) is used. It is suggested
%              that xi and w be provided to avoid needless recomputation of
%              the cubature points.
%            w A numCubPoints X 1 vector of the weights associated with the
%              cubature points.
%
%OUTPUTS: xPred The dX1 predicted target state deltaT in the future.
%         SPred The dXd lower-triangular square root covariance matrix of
%               the predicted target state.
%
%The basic idea behind the prediction filter is described in [1]. The steps
%of the algorithm when holding papx constant are summarized in Appendix B
%of [1]. Specific derivations are given in [2].
%
%The option for the Euler-Maruyama method is a trivial modification
%obtained by truncating the expansion used in the paper.
%
%For algorithm 2, there is not a simple formulation of the square root of
%the covariance matrix, without putting a huge matrix into the tria
%function, so the mean predicted covariance matrix conditioned on the
%state (a cubature point) is computed and then a square root is taken prior
%to combining the result with the spears of the means term.
%
%EXAMPLE:
%Here, we use the same example as in stochTaylorCubPredAdd, but rather than
%verifying its accuracy with respect to Monte Carlo runs, we compare the
%results with stochTaylorCubPredAdd to demopnstrate that they produce
%equivalent answers within finite precision precision bounds.
% deltaT=1;
% omega=4*(pi/180);%Turn rate, radians per second.
% xInit=[1000;2650;0;150;omega];
% aFun=@(x,t)aCoordTurn2DOmega(x);
%           
% omegaRat=5e-4;
% BVal=[0;
%       0;
%       0;
%       0;
%       omegaRat];
% 
% SInit=diag([5;5;2;2;1e-3]);
% PInit=SInit*SInit';
% 
% numSteps=2;
% 
% disp('Euler-Maruyama.')
% [xPred0,PPred0]=stochTaylorCubPredAdd(xInit,PInit,aFun,BVal,deltaT,numSteps,0);
% [xPred00,SPred00]=sqrtStochTaylorCubPredAdd(xInit,SInit,aFun,BVal,deltaT,numSteps,0);
% norm(xPred00-xPred0)
% norm(SPred00*SPred00'-PPred0,'fro')/norm(PPred0,'fro')
% 
% disp('Order 1.5 Itô-Taylor with constant higher-order terms.')
% [xPred1,PPred1]=stochTaylorCubPredAdd(xInit,PInit,aFun,BVal,deltaT,numSteps,1);
% [xPred11,SPred11]=sqrtStochTaylorCubPredAdd(xInit,SInit,aFun,BVal,deltaT,numSteps,1);
% norm(xPred11-xPred1)
% norm(SPred11*SPred11'-PPred1,'fro')/norm(PPred1,'fro')
% 
% disp('Order 1.5 Itô-Taylor')
% [xPred2,PPred2]=stochTaylorCubPredAdd(xInit,PInit,aFun,BVal,deltaT,numSteps,2);
% [xPred22,SPred22]=sqrtStochTaylorCubPredAdd(xInit,SInit,aFun,BVal,deltaT,numSteps,2);
% norm(xPred22-xPred2)
% norm(SPred22*SPred22'-PPred2,'fro')/norm(PPred2,'fro')
%One will see that in all instances, the means are equal, and the
%equivalent covariance matrices are within a few eps() of each other.
%
%REFERENCES:
%[1] I. Arasaratnam, S. Haykin, and T. R. Hurd, "Cubature Kalman filtering
%    for continous-discrete systems: Theory and simulations," IEEE
%    Transactions on Signal Processing, vol. 58, no. 10, pp. 4977-4993,
%    Oct. 2010.
%[2] D. F. Crouse, "Itô-Taylor expansion moments for continuous-time state
%    propagation," NRL Memo, 2019.
%
%December 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xPrev,1);
if(nargin<8||isempty(xi))
    [xi,w]=fifthOrderCubPoints(xDim);
end

if(nargin<7||isempty(algorithm))
    algorithm=2;
end

if(nargin<6||isempty(numSteps))
    numSteps=1;
end

numCubPoints=length(w);
deltaT=deltaT/numSteps;
w=w(:).';
sqrtW=sqrt(w(:).');

x=xPrev;
S=SPrev;
for curStep=1:numSteps
    xiTrans=transformCubPoints(xi,x,S);
    if(algorithm==0)
        %Euler-Maruyama method for additive noise.
        for k=1:numCubPoints
            xiTrans(:,k)=xiTrans(:,k)+a(xiTrans(:,k))*deltaT;
        end
        x=sum(bsxfun(@times,xiTrans,w),2);
        %Center the cubature points.
        xiTrans=bsxfun(@minus,xiTrans,x);
        %Weight the points.
        xiTrans=bsxfun(@times,xiTrans,sqrtW);
        S=tria([xiTrans,sqrt(deltaT)*BVal]);
    elseif(algorithm==1)
        const1=(deltaT^2/2);
        for k=1:numCubPoints
            [aCur,papx,p2apxpx,papt]=a(xiTrans(:,k));
            L0a=L0Operator(aCur,BVal,papt,papx,p2apxpx);
            xiTrans(:,k)=xiTrans(:,k)+a(xiTrans(:,k))*deltaT+const1*L0a;
        end
        xPred=sum(bsxfun(@times,xiTrans,w),2);
        %Center the cubature points.
        xiTrans=bsxfun(@minus,xiTrans,xPred);
        %Weight the points.
        xiTrans=bsxfun(@times,xiTrans,sqrtW);
        
        [~,papx]=a(x);
        Lj=LjOperator(BVal,papx);
        
        x=xPred;
        S=tria([xiTrans,sqrt(deltaT)*(BVal+(deltaT/2)*Lj),sqrt(deltaT^3/12)*Lj]);
    elseif(algorithm==2)%Order 1.5 Taylor method for additive noise,
        %variable higher order terms.
        const1=deltaT^2/2;
        const2=(deltaT^3/3);
        Q=BVal*BVal';
        P0=zeros(xDim,xDim);
        for k=1:numCubPoints
            [aCur,papx,p2apxpx,papt]=a(xiTrans(:,k));
            L0=L0Operator(aCur,BVal,papt,papx,p2apxpx);
            Lj=LjOperator(BVal,papx);
        
            xiTrans(:,k)=xiTrans(:,k)+aCur*deltaT+L0*const1;
            BValLj=BVal*Lj';
            PCur=Q*deltaT+const2*(Lj*Lj')+const1*(BValLj+BValLj');
            P0=P0+w(k)*PCur;
        end
        x=sum(bsxfun(@times,xiTrans,w),2);
        %Center the cubature points.
        xiTrans=bsxfun(@minus,xiTrans,x);
        %Weight the points.
        xiTrans=bsxfun(@times,xiTrans,sqrtW);
        
        %Option 1 does not fully triangularize the square root.
        S0=cholSemiDef(P0,'lower',1);
        S=tria([xiTrans,S0]);
    else
        error('Unknown algorithm specified.')
    end
end
xPred=x;
SPred=S;

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
