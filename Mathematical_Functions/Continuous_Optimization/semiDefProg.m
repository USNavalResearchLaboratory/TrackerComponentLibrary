function [X,y,Z,exitCode]=semiDefProg(C,AMats,a,algorithm,initVals,algParams)
%%SEMIDEFPROG Perform semidefinite programming in the real
%             domain with equality constraints using one of a number of
%             algorithms. This solves the optimization problem:
%                 maximize C(:)'*X(:)
%               subject to sum(sum(AMats(:,:,k).*X))=a(k) for all k
%                      and X is positive semidefinite.
%             The dual to this problem (i.e. a related problem derived from
%             the constraints that is also solved by this solver) is
%                 minimize a'*y
%               subject to sum_k AMats(:,:,k)*y(k)-C=Z
%                      and Z is positive semidefinite.
%
%INPUTS:    C A real nXn matrix.
%       AMats A real nXnXnumConst hypermatrix of numConst nXn matrices, one
%             for each constraint.
%           a A real numConstX1 vector.
%   algorithm An optional parameter specifying which algorithm to use.
%             Algorithms 0 is usually the fastest.
%             0 (The default if omitted or an empty matrix is passed) Use
%               the algorithm of [1] with the predictor-corrector
%               modification of [3]. The stepsize is computed as in [2].
%             1 Use the algorithm of [1]. The stepsize is computed as in
%               [2].
%             2 Use the XZ+ZX algorithm of [2] with Mehrota's predictor-
%               corrector rule using tau=0.9.
%             3 Use the XZ+ZX algorithm of [2] with the basic iteration.
%               The sigma parameter is computed as in [1] using tau=0.9.
%             4 Use the function splittingConicSolver to solve the problem.
%    initVals An optional structure providing initial values for the
%             algorithm. If this is omitted or an empty matrix is passed,
%             then the default values used in Section 5 of [3] are used.
%             When provided, the structure must contain members initVals.X,
%             initVals.y and initVals.Z. The initial estimate is not used
%             if algorithm=4.
%   algParams An optional structure containing parameters that affect the
%             convergence and evaluation of the algorithms. If algorithm=4,
%             then this is the params structre in the splittingConicSolver
%             function, except algParams.max_iters=50e3 as the default.
%             Possible members of the structure for the other algorithms
%             are
%             'epsVal1', 'epsVal2', 'epsVal3' These three entries are the
%                      relative error values for declaring convergence
%                      based on the duality gap and the satisfaction of
%                      the dual and primal constraints. Convergence is
%                      declared when
%                      abs(C(:)'*X(:)-a'*y)/(1+abs(a'*y))<epsVal1
%                      and
%                      sqrt(sum_k(sum(sum(AMats(:,:,k).*X))-a(k))^2)/(1+norm(a))<epsVal2
%                      and
%                      norm(sum_k AMats(:,:,k)*y(k)-C-Z,'fro')/(1+norm(C,'fro'))<epsVal3
%                      The default values for all three parameters if
%                      omitted are 1e-7.
%      primalInfeasVal The primal is declared infeasible if
%                      -a'*y/norm(sum_k AMats(:,:,k)*y(k)-Z,'fro')>primalInfeasVal
%                      The default value if omitted is 1e10.
%        dualInfeasVal The dual is declared infeasible if
%                      C(:)'*X(:)/norm(sum_k sum(sum(AMats(:,:,k).*X)))>dualInfeasVal
%                      The default value if omitted is 1e10.
%              maxIter The maximum number of iterations. The default value
%                      is 500.
%
%OUTPUTS: X The nXn solution to the primal or an empty matrix if the
%           problem is infeasible or an error occurred.
%         y The numConstX1 set of dual variables or an empty matrix is the
%           problem is infeasible or an error occurred.
%         Z The nXn set of dual matrix variables.
%  exitCode A parameter specifying how the algorithm terminated. If
%           algorithm=4, then this is the info structure from the
%           splittingConicSolver algorithm. Otherwise, possible
%           values are:
%           0 The algorithm converged.
%           1 The maximum number of iterations was reached.
%           2 The primal is infeasible.
%           3 The dual is infeasible.
%           4 A non-finite number or other numerical problems were
%             encountered.
%
%EXAMPLE 1:
%Here we have a very simple problem with two constraints:
% C=[1, 7, 0, 0, 0, 0, 0;
%    7, 2, 0, 0, 0,-2, 0;
%    0, 0, 3, 0, 6, 0, 0;
%    0, 0, 0, 4, 0, 0, 0;
%    0, 0, 6, 0, 5, 0, 0;
%    0, -2,0, 0, 0, 0, 0;
%    0, 0, 0, 0, 0, 0, 0];
% AMats=zeros(7,7,2);
% AMats(:,:,1)=[4, 6, 0, 0, 0, 0, 0;
%               6, 2, 0, 0, 0, 0, 0;
%               0, 0, 0, 0, 0, 0, 0;
%               0, 0, 0, 0, 0, 0, 0;
%               0, 0, 0, 0, 0, 0, 0;
%               0, 0, 0, 0, 0, 7, 0;
%               0, 0, 0, 0, 0, 0, 0];
% AMats(:,:,2)=[0, 0, 0, 0, 0, 0, 0;
%               0, 0, 0, 0, 0, 0, 0;
%               0, 0, 12,1, 24,0, 0;
%               0, 0, 1, 4, 0, 0, 0;
%               0, 0, 24,0, 48,0, 0;
%               0, 0, 0, 0, 0, 0, 0;
%               0, 0, 0, 0, 0, 0, 1];
% a=[1;2];
% [X0,y0,Z0,exitCode0]=semiDefProg(C,AMats,a)
%This problem is solved very well using methods 0 and 1. However, methods 2
%and 3 will generally not convergence, even though they work very well for
%the second example problem (albeit requiring a very large number of
%iterations). Though the stability of methods 2 and 3 is praised in [2],
%they often leave something to be desired, at least as implemented here.
%
%EXAMPLE 2:
%The second example is the example in Section 6.3 of [1]. The method is
%related to finding the maximum clique in a graph. Quite a lot of work goes
%into setting up the problem before it can be passed to the functions.
%Here, C is all ones and the constraints are that the trace of X equal 1
%and certain edges (symmetric values in X) be zero. The problem selected
%has a 50X50 X and 125 constraints.
% n=50;
% A=zeros(n,n);
% seed=3;
% prob=0.9;
% r=(4*seed+1)/16384^2;
% for i=1:(n-1)
%     for j=(i+1):n
%         r=mod(r*41475557,1);
%         A(i,j)=(r<(1-prob));
%     end
% end
% C=ones(n,n);
% 
% noEdgeIdx=find(A==1);
% k=length(noEdgeIdx)+1;
% AMats=zeros(n,n,k);
% 
% numEls=n*(n+1)/2;
% ASCS=zeros(k+numEls,numEls);
% for curIdx=1:(k-1)
%     ei=zeros(n,1);
%     ej=zeros(n,1);
%     [i,j]=ind2sub([n,n],noEdgeIdx(curIdx));
%     ei(i)=1;
%     ej(j)=1;
%     AMats(:,:,curIdx)=ei*ej'+ej*ei';
%     
%     idx=sub2VechInd(n,i,j);
%     ASCS(curIdx,idx)=1;
% end
% AMats(:,:,end)=eye(n,n);
% ASCS(k,:)=vech(eye(n));
% 
% %Must set semidefinite constraints, with appropriate scale factor.
% curVal=k+1;
% curEntry=1;
% for curCol=1:(n-1)
%     ASCS(curVal,curEntry)=1;
%     curEntry=curEntry+1;
%     curVal=curVal+1;
%     for curRow=(curCol+1):n
%         ASCS(curVal,curEntry)=sqrt(2);
%         curEntry=curEntry+1;
%         curVal=curVal+1;
%     end
% end
% 
% a=zeros(k,1);
% a(end)=1;
% 
% initVals.X=(1/n)*eye(n,n);
% lambda=1.1*n;
% initVals.Z=lambda*eye(n,n)-C;
% initVals.y=zeros(k,1);
% [X0,y0,Z0,exitCode0]=semiDefProg(C,AMats,a,0,initVals);
%One will find that C(:)'*X(:) is about 21.0910 regardless of the method
%used. However, while methods 0 and 1 take just a few seconds, methods 2
%and 3 are very slow. They are significantly faster if one sets usePinv to
%false.
%
%REFERENCES:
%[1] C. Helmberg, F. Rendl, R. J. Vanderbei, and H. Wolkowicz, "An interior
%    -point method for semidefinite programming," SIAM Journal on
%    Optimization, vol. 6, no. 2, pp. 342-361, May 1996.
%[2] F. Alizadeh, J.-P. A. Haeberly, and M. L. Overton, "Primal-dual
%    interior-point methods for semidefinite programming: Convergence
%    rates, stability and numerical results," SIAM Journal on Optimization,
%    vol. 8, no. 3, pp. 746-768, 1998.
%[3] B. Borchers, "CSDP, a C library for semidefinite programming,"
%    Optimization methods and Software, vol. 11, no. 1-4, pp. 613-623,
%    1999.
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(algorithm))
   algorithm=0;
end

%Initial feasible  solution.
if(nargin<5||isempty(initVals))
%If no initialization is provided, then use the default initialization
%described in Section 5 of [3].
    n=size(C,1);
    k=size(AMats,3);
    
    %Find alpha and beta from Equations 25 and 26 in [3].
    alpha=0;
    beta=0;
    for i=1:k
        AF=norm(AMats(:,:,i),'fro');
        alpha=max(alpha,n*(1+abs(a(k)))/(1+AF));
        beta=max(beta,AF);
    end
    CF=norm(C,'fro');
    beta=(1+max(beta,CF))/sqrt(n);

    initVals.X=alpha*eye(n,n);
    initVals.Z=beta*eye(n,n);
    initVals.y=zeros(k,1);
end

%Default parameters.
epsVal1=1e-7;%Duality gap constraint
epsVal2=1e-7;%Primal constraint gap
epsVal3=1e-7;%Dual constraint gap
primalInfeasVal=1e10;
dualInfeasVal=1e10;
maxIter=500;
usePinv=true;

%Extract any parameters that might override the options.
if(nargin>5&&~isempty(algParams)&&algorithm~=4)
    if(isfield(algParams,'epsVal1'))
        epsVal1=algParams.epsVal1;
    end
    if(isfield(algParams,'epsVal2'))
        epsVal2=algParams.epsVal2;
    end
    if(isfield(algParams,'epsVal3'))
        epsVal3=algParams.epsVal3;
    end
    if(isfield(algParams,'primalInfeasVal'))
        primalInfeasVal=algParams.primalInfeasVal;
    end
    if(isfield(algParams,'dualInfeasVal'))
        dualInfeasVal=algParams.dualInfeasVal;
    end
    if(isfield(algParams,'maxIter'))
        maxIter=algParams.maxIter;
    end
    if(isfield(algParams,'usePinv'))
        usePinv=algParams.usePinv;
    end
elseif(algorithm==4)
    if(nargin<5)
        algParams=[];%The default when using the splitting conic solver.
        algParams.max_iters=50e3;
    end
    if(~isfield(algParams,'max_iters'))
        algParams.max_iters=50e3;
    end
end

switch(algorithm)
    case 0%The XZ method of [1] with the predictor-corrector
          %modification of [3], using the approach to compute the stepsize
          %in [2].
        usePredictorCorrector=true;
        [X,y,Z,exitCode]=semiDefProgXZ(C,AMats,a,initVals,epsVal1,epsVal2,epsVal3,usePredictorCorrector,maxIter,primalInfeasVal,dualInfeasVal);
    case 1%The XZ method of [1] using the approach to compute the stepsize
          %in [2].
        usePredictorCorrector=false;
        [X,y,Z,exitCode]=semiDefProgXZ(C,AMats,a,initVals,epsVal1,epsVal2,epsVal3,usePredictorCorrector,maxIter,primalInfeasVal,dualInfeasVal);
    case 2%The XZ-ZX method of [2] with Mehrota's predictor-corrector rule.
        C=-C;
        [X,y,Z,exitCode]=semiDefProgZXXZMehrota(C,AMats,a,initVals,usePinv,epsVal1,epsVal2,epsVal3,maxIter,primalInfeasVal,dualInfeasVal);
    case 3%The XZ-ZX method of [2] as a single step.
        C=-C;
        [X,y,Z,exitCode]=semiDefProgZXXZ(C,AMats,a,initVals,usePinv,epsVal1,epsVal2,epsVal3,maxIter,primalInfeasVal,dualInfeasVal);
    case 4%Use the splittingConicSolver function.
        n=size(AMats,1);
        k=size(AMats,3);
        %The number of equality constraints
        cone.f=k;
        cone.s=n;
        
        numEls=n*(n+1)/2;
        ASCS=zeros(k+numEls,numEls);
        %Equality constraints
        for curConst=1:k
            ASCS(curConst,:)=vech(AMats(:,:,curConst));
        end
        
        %Now, the semidefinite constraint. The off-diagonal elements have
        %to be multiplied by 1/sqrt(2).
        ASCS((k+1):end,:)=diag(-vech(ones(n,n),1/sqrt(2)));
        b=[a;zeros(numEls,1)];
        c=-vech(C);%Negative for maximization.

        [x,y,~,info]=splittingConicSolver(ASCS,b,c,cone,algParams);
        X=vech2Mat(x,true,1/2);
        Z=vech2Mat(y((k+1):end),true,1/sqrt(2));
        y=y(1:k);
        exitCode=info;
    otherwise
        error('Unknown algorithm specified')
end

end

function [X,y,Z,exitCode]=semiDefProgZXXZ(C,AList,b,initVals,usePinv,epsVal1,epsVal2,epsVal3,maxIter,primalInfeasVal,dualInfeasVal)
%%SEMIDEFPROGZXXZ This function implements the ZX+XZ method of [1] using
%                 the basic iteration. An option is provided to use
%                 pseudoinverses when solving equations involving E and M.
%                 Such solutions are slower, but might be numerically
%                 stabler in some instances.
%
%REFERENCES:
%[1] F. Alizadeh, J.-P. A. Haeberly, and M. L. Overton, "Primal-dual
%    interior-point methods for semidefinite programming: Convergence
%    rates, stability and numerical results," SIAM Journal on Optimization,
%    vol. 8, no. 3, pp. 746-768, 1998.
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(AList,1);
m=size(AList,3);

tau=0.9;

X=initVals.X;
y=initVals.y;
Z=initVals.Z;

A=zeros(m,n*(n+1)/2);
for k=1:m
   A(k,:)=vech(AList(:,:,k),sqrt(2));
end

sigma=1/2;

CF=norm(C,'fro');

for curIter=1:maxIter
    x=vech(X,sqrt(2));
    XZ=X*Z;
    ZX=Z*X;
    %This is equal to sum(sum(X.*Z)), because X and Z are both symmetric.
    XDotZ=trace(XZ);

%%%Step 1 of the basic iteration
    mu=sigma*XDotZ/n;%Equation 2.17

%%%Step 2 of the basic iteration
    %Compute the initial deltas with mu=0.
    %Defined before Equation 2.8
    rp=b-A*x;

    Rd=C-Z-vech2Mat(A'*y,1/sqrt(2));
    Rc=mu*eye(n,n)-(1/2)*(XZ+ZX);%Equation 2.8
    
    %%%Test for convergence
    %Relative duality gap constraint.
    CXVal=sum(C(:).*X(:));
    
    diff1=abs(CXVal-b'*y)/(1+abs(b'*y));
    diff2=norm(rp)/(1+norm(b));
    diff3=norm(Rd,'fro')/(1+CF);
    
    if(diff1<epsVal1 && diff2<epsVal2&&diff3<epsVal3)
        exitCode=0;%Convergence achieved
        return;
    end

    %If the primal is infeasible.
    if(-b'*y/norm(Rd+C,'fro')>primalInfeasVal)
       X=[];
       Z=[];
       y=[];
       exitCode=2;
       return; 
    end
    
    %If the dual is infeasible.
    if(CXVal/norm(rp+b)>dualInfeasVal)
        X=[];
        Z=[];
        y=[];
        exitCode=3;
        return;
    end

    rd=vech(Rd,sqrt(2));
    rc=vech(Rc,sqrt(2));

    %These are defined after Equation 2.10
    E=kronSym(Z,eye(n,n));
    F=kronSym(X,eye(n,n));
    %Section 5 suggests using Equations 2.12 and 2.13  for computing deltaX
    %and deltaY for stability.
    if(usePinv)
        EInv=pinv(E);
        %Create the inverse of the Schur complement matrix of Equation 2.15.
        MInv=pinv(A*(EInv*F)*A');
        
        deltaY=MInv*(rp+A*EInv*(F*rd-rc));%Solving Equation 2.12
        deltaX=-EInv*(F*(rd-A'*deltaY)-rc);%Equation 2.13
    else
        M=A*(E\F)*A';
        deltaY=M\(rp+A*(E\(F*rd-rc)));%Solving Equation 2.12
        deltaX=-E\(F*(rd-A'*deltaY)-rc);%Equation 2.13
    end
    deltaZ=rd-A'*deltaY;%Equation 2.14
    
    DeltaX=vech2Mat(deltaX,1/sqrt(2));
    DeltaZ=vech2Mat(deltaZ,1/sqrt(2));
    
    if(any(~isfinite(deltaX(:)))||any(~isfinite(deltaY(:)))||any(~isfinite(deltaZ(:))))
        X=[];
        y=[];
        Z=[];
        exitCode=4;%Finite precision errors occured.
        return; 
    end
    
%%%Step 3 of the basic iteration does not apply to the XZ+ZX method.
%%%Step 4 of the basic iteration
    %Choose the steplengths
    alpha=findAlpha(X,DeltaX,tau);
    beta=findAlpha(Z,DeltaZ,tau);

    X=X+alpha*DeltaX;
    y=y+beta*deltaY;
    Z=Z+beta*DeltaZ;
    
    if(alpha+beta>1.8)
        sigma=1/4;
    else
        sigma=1/2;
    end
end

%Maximum number of iterations reached without convergence.
exitCode=1;

end

function [X,y,Z,exitCode]=semiDefProgZXXZMehrota(C,AList,b,initVals,usePinv,epsVal1,epsVal2,epsVal3,maxIter,primalInfeasVal,dualInfeasVal)
%%SEMIDEFPROGZXXZ This function implements the ZX+XZ method of [1] using
%                 Mehrota's predictor-corrector rule. An option is
%                 provided to use pseudoinverses when solving equations
%                 involving E and M. Such solutions are slower, but might
%                 be numerically stabler in some instances.
%
%REFERENCES:
%[1] F. Alizadeh, J.-P. A. Haeberly, and M. L. Overton, "Primal-dual
%    interior-point methods for semidefinite programming: Convergence
%    rates, stability and numerical results," SIAM Journal on Optimization,
%    vol. 8, no. 3, pp. 746-768, 1998.
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(AList,1);
m=size(AList,3);

X=initVals.X;
y=initVals.y;
Z=initVals.Z;

tau=0.9;

A=zeros(m,n*(n+1)/2);
for k=1:m
   A(k,:)=vech(AList(:,:,k),sqrt(2));
end

CF=norm(C,'fro');

for curIter=1:maxIter
    x=vech(X,sqrt(2));
    XZ=X*Z;
    ZX=Z*X;
    %This is equal to sum(sum(X.*Z)), because X and Z are both symmetric.
    XDotZ=trace(XZ);

%%%Step 1 of Mehrotra's predictor-corrector algorithm.
    %Compute the initial deltas with mu=0.
    %Defined before Equation 2.8
    rp=b-A*x;

    Rd=C-Z-vech2Mat(A'*y,1/sqrt(2));
    Rc=-(1/2)*(XZ+ZX);%Equation 2.8
    
    %%%Test for convergence
    %Relative duality gap constraint.
    CXVal=sum(C(:).*X(:));
    diff1=abs(CXVal-b'*y)/(1+abs(b'*y));
    diff2=norm(rp)/(1+norm(b));
    diff3=norm(Rd,'fro')/(1+CF);
    
    if(diff1<epsVal1 && diff2<epsVal2&&diff3<epsVal3)
        exitCode=0;%Convergence achieved
        return;
    end

    %If the primal is infeasible.
    if(-b'*y/norm(Rd+C,'fro')>primalInfeasVal)
       X=[];
       Z=[];
       y=[];
       exitCode=2;
       return; 
    end
    
    %If the dual is infeasible.
    if(CXVal/norm(rp+b)>dualInfeasVal)
        X=[];
        Z=[];
        y=[];
        exitCode=3;
        return;
    end

    rd=vech(Rd,sqrt(2));
    rc=vech(Rc,sqrt(2));

    %These are defined after Equation 2.10
    E=kronSym(Z,eye(n,n));
    F=kronSym(X,eye(n,n));
    %Section 5 suggests using Equations 2.12 and 2.13  for computing deltaX
    %and deltaY for stability.    
    if(usePinv)
        EInv=pinv(E);
        %Create the inverse of the Schur complement matrix of Equation 2.15.
        MInv=pinv(A*(EInv*F)*A');
        
        deltaY=MInv*(rp+A*EInv*(F*rd-rc));%Solving Equation 2.12
        deltaX=-EInv*(F*(rd-A'*deltaY)-rc);%Equation 2.13
    else
        M=A*(E\F)*A';
        
        deltaY=M\(rp+A*(E\(F*rd-rc)));%Solving Equation 2.12
        deltaX=-E\(F*(rd-A'*deltaY)-rc);%Equation 2.13
    end
    deltaZ=rd-A'*deltaY;%Equation 2.14
    
    DeltaX=vech2Mat(deltaX,1/sqrt(2));
    DeltaZ=vech2Mat(deltaZ,1/sqrt(2));
    
    if(any(~isfinite(deltaX(:)))||any(~isfinite(deltaY(:)))||any(~isfinite(deltaZ(:))))
        X=[];
        y=[];
        Z=[];
        exitCode=4;%Finite precision errors occured.
        return; 
    end
    
%%%Step 2 of Mehrotra's predictor-corrector algorithm.
    %Choose the steplengths
    alpha=findAlpha(X,DeltaX,tau);
    beta=findAlpha(Z,DeltaZ,tau);

    %Equation 7.1
    sigma=(sum(sum((X+alpha*DeltaX).*(Z+beta*DeltaZ)))/XDotZ)^3;
    mu=sigma*XDotZ/n;

%%%Step 3 of Mehrotra's predictor-corrector algorithm.
    %Recompute the delta with a new mu and a new value of Rc.
    Rc=mu*eye(n,n)-(1/2)*(XZ+ZX+DeltaX*DeltaZ+DeltaZ*DeltaX);
    rc=vech(Rc,sqrt(2));
    if(usePinv)
        deltaY=MInv*(rp+A*EInv*(F*rd-rc));%Solving Equation 2.12
        deltaX=-EInv*(F*(rd-A'*deltaY)-rc);%Equation 2.13
    else
        deltaY=M\(rp+A*(E\(F*rd-rc)));%Solving Equation 2.12
        deltaX=-E\(F*(rd-A'*deltaY)-rc);%Equation 2.13
    end
    deltaZ=rd-A'*deltaY;%Equation 2.14

    DeltaX=vech2Mat(deltaX,1/sqrt(2));
    DeltaZ=vech2Mat(deltaZ,1/sqrt(2));
    
    if(any(~isfinite(deltaX(:)))||any(~isfinite(deltaY(:)))||any(~isfinite(deltaZ(:))))
        X=[];
        y=[];
        Z=[];
        exitCode=4;%Finite precision errors occured.
        return; 
    end
    
    %Now, we must find new versions of alpha and beta
    alpha=findAlpha(X,DeltaX,tau);
    beta=findAlpha(Z,DeltaZ,tau);

    X=X+alpha*DeltaX;
    y=y+beta*deltaY;
    Z=Z+beta*DeltaZ;
end

%Maximum number of iterations reached without convergence.
exitCode=1;

end

function [X,y,Z,exitCode]=semiDefProgXZ(C,AMats,a,initVals,epsVal1,epsVal2,epsVal3,usePredictorCorrector,maxIter,primalInfeasVal,dualInfeasVal)
%%SEMIDEFPROGXZ This implements the algorithm of [1] with and without the
%               predictor-corrector modification of [2]. In both instances,
%               the stepsize choice of [3] is used. Only equality
%               constraints are considered.
%
%REFERENCES:
%[1] C. Helmberg, F. Rendl, R. J. Vanderbei, and H. Wolkowicz, "An interior
%    -point method for semidefinite programming," SIAM Journal on
%    Optimization, vol. 6, no. 2, pp. 342-361, May 1996.
%[2] B. Borchers, "CSDP, a C library for semidefinite programming,"
%    Optimization methods and Software, vol. 11, no. 1-4, pp. 613-623,
%    1999.
%[3] F. Alizadeh, J.-P. A. Haeberly, and M. L. Overton, "Primal-dual
%    interior-point methods for semidefinite programming: Convergence
%    rates, stability and numerical results," SIAM Journal on Optimization,
%    vol. 8, no. 3, pp. 746-768, 1998.
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(C,1);
k=size(AMats,3);

X=initVals.X;
y=initVals.y;
Z=initVals.Z;

CF=norm(C,'fro');

Aty=matFunT(AMats,y);
mu=sum(Z(:).*X(:))/(2*n);%Equation 9
for curIter=1:maxIter
    %We want to find the inverse of Z. inv works poorly with badly
    %conditioned matrices and pinv is slow. Thus, we will try to to a
    %square root using a Cholesky decomposition and forward/ backward
    %substitution (like inv). However, if it fails, we will just use pinv
    %anyway.
    [ZChol,cholDidFail]=chol(Z,'lower');
    if(cholDidFail)
        ZInv=pinv(Z);
    else
        opts.LT=true;
        opts.UT=false;
        %Solve using forward substitution.
        temp=linsolve(ZChol,eye(size(ZChol)),opts);
        opts.UT=true;
        opts.LT=false;
        %Solve using backward substitution.
        ZInv=linsolve(ZChol',temp,opts); 
    end
 
    Fd=-Aty+C+Z;%Equation 12

    %Allocate space
    O11=zeros(k,k);
    %Fill in values from Equation 16. 
    for curCol=1:k
        O11(:,curCol)=matFun(AMats,ZInv*AMats(:,:,curCol)*X);
    end
    
    [O11Chol,cholDidFail]=chol(O11,'lower');
    if(cholDidFail)
        O11Inv=pinv(O11);
    else
        opts.LT=true;
        opts.UT=false;
        %Solve using forward substitution.
        temp=linsolve(O11Chol,eye(size(O11Chol)),opts);
        opts.UT=true;
        opts.LT=false;
        %Solve using backward substitution.
        O11Inv=linsolve(O11Chol',temp,opts);
    end
    
    if(usePredictorCorrector)
        %Predictor step with mu=0
        deltaYHat=O11Inv*(-a+matFun(AMats,ZInv*Fd*X));
        AtDeltaYHat=matFunT(AMats,deltaYHat);
        deltaXHat=-X+ZInv*(Fd-AtDeltaYHat)*X;%Equation 17
        deltaXHat=(deltaXHat+deltaXHat')/2;%Force symmetry.
        deltaZHat=-Fd+AtDeltaYHat;%Equation 18

        %The corrector step
        deltaYBar=O11Inv*(matFun(AMats,mu*ZInv)-matFun(AMats,ZInv*deltaZHat*deltaXHat));
        AtDeltaYBar=matFunT(AMats,deltaYBar);
        deltaXBar=ZInv*(mu*eye(n,n)-deltaZHat*deltaXHat-AtDeltaYBar*X);
        deltaXBar=(deltaXBar+deltaXBar')/2;%Force symmetry.
        deltaZBar=AtDeltaYBar;

        deltaX=deltaXHat+deltaXBar;
        deltaY=deltaYHat+deltaYBar;
        deltaZ=deltaZHat+deltaZBar;
    else
        %Using v1 from Equation 4.20
        deltaY=O11Inv*(mu*matFun(AMats,ZInv)-a+matFun(AMats,ZInv*Fd*X));
        AtDeltaY=matFunT(AMats,deltaY);
        deltaX=mu*ZInv-X+ZInv*(Fd-AtDeltaY)*X;%Equation 4.16
        deltaX=(deltaX+deltaX')/2;%Force symmetry.
        deltaZ=-Fd+AtDeltaY;
    end

    if(any(~isfinite(deltaX(:)))||any(~isfinite(deltaY(:)))||any(~isfinite(deltaZ(:))))
        X=[];
        y=[];
        Z=[];
        exitCode=4;%Finite precision errors occured.
        return; 
    end
    
    alphaP=findAlpha(X,deltaX,0.9);
    X=X+alphaP*deltaX;%Take the step in the primal.
    alphaD=findAlpha(Z,deltaZ,0.9);

    %Take the step in the dual.
    Z=Z+alphaD*deltaZ;
    y=y+alphaD*deltaY;
    
    AX=matFun(AMats,X);
    
    Aty=matFunT(AMats,y);
    
    %The differences needed to test for convergence.
    diff1=abs(sum(C(:).*X(:))-a'*y)/(1+abs(a'*y));
    diff2=norm(AX-a)/(1+norm(a));
    diff3=norm(Aty-C-Z,'fro')/(1+CF);
    
    if(diff1<epsVal1&&diff2<epsVal2&&diff3<epsVal3)
        exitCode=0;
        return;
    end
    
    %The ad-hoc reduction when large steps are taken.
    if(alphaP+alphaD>1.8)
        mu=sum(Z(:).*X(:))/(4*n);
    else
        mu=sum(Z(:).*X(:))/(2*n);
    end

    %If the primal is infeasible.
    if(-a'*y/norm(Aty-Z,'fro')>primalInfeasVal)
       X=[];
       Z=[];
       y=[];
       exitCode=2;
       return; 
    end
    
    %If the dual is infeasible.
    if(sum(C(:).*X(:))/norm(AX)>dualInfeasVal)
        X=[];
        Z=[];
        y=[];
        exitCode=3;
        return;
    end
end

%If we get here, it got to the maximum number of iterations without
%converging.
exitCode=1;

end


function val=matFun(theMats,X)
%%MATFUN This function implements Equations 2 and 3 of [1].
%
%REFERENCES:
%[1] B. Borchers, "CSDP, a C library for semidefinite programming,"
%    Optimization methods and Software, vol. 11, no. 1-4, pp. 613-623,
%    1999.
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

k=size(theMats,3);
val=reshape(sum(sum(bsxfun(@times,theMats,X),1),2),k,1);
end

function val=matFunT(theMats,theVec)
%%MATFUNT This function implements Equations 5 and 6 of [1].
%
%REFERENCES:
%[1] B. Borchers, "CSDP, a C library for semidefinite programming,"
%    Optimization methods and Software, vol. 11, no. 1-4, pp. 613-623,
%    1999.
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
k=size(theMats,3);
val=sum(bsxfun(@times,reshape(theVec,1,1,k),theMats),3);
end

function alpha=findAlpha(X,DeltaX,tau)
%%FINDALPHA This function finds alpha using Equation 2.18 of [1] as well as
%           the identity with respect to the eigenvalues of a matrix after
%           that. Due to the symmetry of the equations, this function can
%           also be used to find beta from Equation 2.19 using Z and deltaZ
%           instead of X and deltaX.
%
%REFERENCES:
%[1] F. Alizadeh, J.-P. A. Haeberly, and M. L. Overton, "Primal-dual
%    interior-point methods for semidefinite programming: Convergence
%    rates, stability and numerical results," SIAM Journal on Optimization,
%    vol. 8, no. 3, pp. 746-768, 1998.
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

LInvX=pinv(cholSemiDef(X,'lower'));
alphaHat=1/max(eig(-LInvX*DeltaX*LInvX'));

%If deltaX is positive definite then it is possible that the maximum
%eigenvalue is negative. This just sets the step size to 1 in that case.
if(alphaHat<0)
    alphaHat=1;
end

alpha=min(1,tau*alphaHat);%Equation 2.18
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
