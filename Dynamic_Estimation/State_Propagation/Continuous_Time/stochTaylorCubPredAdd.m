function [xPred,PPred]=stochTaylorCubPredAdd(xPrev,PPrev,a,BVal,deltaT,numSteps,algorithm,xi,w)
%%STOCHTAYLORCUBPREDALTADD Given a Gaussian prior distribution with mean
%           xPrev and covariance matrix PPrev under a continuous time,
%           autonomous dynamic model under Itô calculus with additive noise
%           of the form:
%           dx=a(x)dt+BVal*dw
%           where dw is the differential of a Wiener process and BVal is
%           constant matrix, approximate the mean and covariance matrix of
%           the distribution predicted forward deltaT. This is done by
%           computing the moments of a stochastic Itô-Taylor expansion,
%           possibly taking more than one step (using the previous
%           expansion as the prior to the next step).
%
%INPUTS: xPrev The dX1 prior target state.
%        PPrev The dXd prior state covariance matrix.
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
%           xi An dXnumCubPoints matrix of cubature points. If this
%              and the next parameter are omitted or empty matrices are
%              passed, then fifthOrderCubPoints(d) is used. It is suggested
%              that xi and w be provided to avoid needless recomputation of
%              the cubature points.
%            w A numCubPointsX1 vector of the weights associated with the
%              cubature points.
%
%OUTPUTS: xPred The dX1 predicted target state deltaT in the future.
%         PPred The dXd covariance matrix of the predicted target state.
%
%The basic idea behind the prediction filter is described in [1]. The steps
%of the algorithm when holding papx constant are summarized in Appendix A
%of [1]. Here, a typo in step 5 is corrected, since square roots were
%missing on two of the Q terms. Specific derivations are given in [2].
%
%The option for the Euler-Maruyama method is a trivial modification
%obtained by truncating the expansion used in the paper.
%
%EXAMPLE:
%Here, we consider the prediction of a target through a coordinated turn in
%2D. The target state is position, velocity and a turn rate. The dynamic
%model only adds noise to the turn rate. We compare the predicted mean and
%covariance matrix with that obtained via numMC Monte Carlo runs that are
%done on a much finder number of steps. The Monte Carlo runs of the random
%paths can take a while.
% numMC=1000;
% numStepsMC=300;
% deltaT=1;
% omega=4*(pi/180);%Turn rate, radians per second.
% xInit=[1000;2650;0;150;omega];
% xDim=length(xInit);
% 
% a=@(x,t)[x(3);
%        x(4);
%        -x(5)*x(4);
%        x(5)*x(3);
%        0];
% papx=@(x,t)[0,0,1,        0,0;
%             0,0,0,        1,0;
%             0,0,0,    -x(5),-x(4);
%             0,0,x(5),     0,x(3);
%             0,0,0,        0,0];     
% p2apxpx=zeros(xDim,xDim,xDim);
% p2apxpx(:,3,5)=[0;0;0;1;0];
% p2apxpx(:,5,3)=p2apxpx(:,3,5);
% p2apxpx(:,4,5)=[0;0;-1;0;0];
% p2apxpx(:,5,4)=p2apxpx(:,4,5);
% 
% papt=zeros(xDim,1);
% aFun=@(x,t)dealRobust(a(x),papx(x),p2apxpx,papt);
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
% xMC=zeros(xDim,numMC);
% for curRun=1:numMC
%     xCur=xInit+SInit*randn(xDim,1);
%     for curStep=2:numStepsMC
%         xPrev=xCur;
%         [aCur,papy,p2apypy]=aFun(xPrev);
%         xCur=strongStochTaylorStep(xPrev,aCur,BVal,deltaT/numStepsMC,6,[],[],papy,p2apypy);
%     end
%     xMC(:,curRun)=xCur;
% end
% [xPredMC,PPredMC]=calcMixtureMoments(xMC);
% numSteps=2;
% [xPred0,PPred0]=stochTaylorCubPredAdd(xInit,PInit,aFun,BVal,deltaT,numSteps,0);
% [xPred1,PPred1]=stochTaylorCubPredAdd(xInit,PInit,aFun,BVal,deltaT,numSteps,1);
% [xPred2,PPred2]=stochTaylorCubPredAdd(xInit,PInit,aFun,BVal,deltaT,numSteps,2);
% 
% disp('Euler-Maruyama.')
% mean(abs(xPred0-xPredMC)./abs(xPredMC))
% norm(PPred0-PPredMC,'fro')./norm(PPredMC,'fro')
% disp('Order 1.5 Itô-Taylor with constant higher-order terms.')
% mean(abs(xPred1-xPredMC)./abs(xPredMC))
% norm(PPred1-PPredMC,'fro')./norm(PPredMC,'fro')
% disp('Order 1.5 Itô-Taylor')
% mean(abs(xPred2-xPredMC)./abs(xPredMC))
% norm(PPred2-PPredMC)./norm(PPredMC,'fro')
%One will typically see that the order 1.5 prediction modifying the
%higher-order terms has a more accurate covariance matrix (in terms of
%relative Frobenius norm) than the other two methods. Also, the mean
%relative abolute error of the order 1.5 methods is typically better than
%The Euler-Maruyama method.
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

Q=BVal*BVal';

deltaT=deltaT/numSteps;
x=xPrev;
P=PPrev;
for curStep=1:numSteps
    S=chol(P,'lower');
    xiTrans=transformCubPoints(xi,x,S);
    if(algorithm==0)%Euler-Maruyama method for additive noise.
            for k=1:numCubPoints
                xiTrans(:,k)=xiTrans(:,k)+a(xiTrans(:,k))*deltaT;
            end
            [xPred,PPred0]=calcMixtureMoments(xiTrans,w);
            
            PPred=PPred0+deltaT*Q;
     elseif(algorithm==1)%Order 1.5 Taylor method for additive noise,
            %constant higher order terms.
            const1=(deltaT^2/2);
            for k=1:numCubPoints
                [aCur,papx,p2apxpx,papt]=a(xiTrans(:,k));
                L0a=L0Operator(aCur,BVal,papt,papx,p2apxpx);
                xiTrans(:,k)=xiTrans(:,k)+a(xiTrans(:,k))*deltaT+const1*L0a;
            end
            [xPred,PPred0]=calcMixtureMoments(xiTrans,w);
            
            [~,papx]=a(x);
            Lja=LjOperator(BVal,papx);

            %Modified from [1], since it has Lja*Q'+Q*Lja'. However, it
            %should be the square root of Q, which is BVal.
            PPred=PPred0+const1*(Lja*BVal'+BVal*Lja')+(deltaT^3/3)*(Lja*Lja')+deltaT*Q;
    elseif(algorithm==2)%Order 1.5 Taylor method for additive noise,
        %variable higher order terms.
        const1=(deltaT^2/2);
        const2=(deltaT^3/3);
        Sigma0=zeros(xDim,xDim);
        for k=1:numCubPoints
            [aCur,papx,p2apxpx,papt]=a(xiTrans(:,k));
            L0=L0Operator(aCur,BVal,papt,papx,p2apxpx);
            Lj=LjOperator(BVal,papx);
        
            xiTrans(:,k)=xiTrans(:,k)+aCur*deltaT+L0*const1;

            BValLj=BVal*Lj';
            SigmaCur=Q*deltaT+const2*(Lj*Lj')+const1*(BValLj+BValLj'); 
            
            Sigma0=Sigma0+w(k)*SigmaCur;
        end
        [xPred,PPred0]=calcMixtureMoments(xiTrans,w);
        PPred=PPred0+Sigma0;
    else
        error('Unknown algorithm specified.')
    end
    x=xPred;
    P=PPred;
end

xPred=x;
PPred=P;
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
