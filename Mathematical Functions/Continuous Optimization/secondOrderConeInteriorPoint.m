function [x,y,s,exitCode]=secondOrderConeInteriorPoint(A,b,c,blockLengths,params)
%%SECONDORDERCONEINTERIORPOINT Perform second order cone programming
%           (SOCP) using an interior point algorithm. This solves the
%           primal optimization problem:
%             minimize sum_{i=1}^n c_i'*x_i
%           subject to sum_{i=1}^n A_i*x_i=b
%             and x_i(1)-norm(x_i(2:end))>=0 (second order cone constraint)
%           which is equivalent to the dual optimization problem
%            maximize b'*y
%           subject to sum_{i=1}^n A_i'*y_i+s_i=c_i
%            and s_i(1)-norm(s_i(2:end))>=0
%           Many types of optimization problems can be formulated as SOCP
%           problems. This function only supports optimization in the real
%           domain.
%
%INPUTS: A A real mXsum(blockLengths) matrix containing coefficients for
%          all of the constraints on the primal. A must have full row rank.
%          With regard to the A_i in the problem formulation,
%          A=[A_1,a_2,...A_n]
%        b An mX1 real vector the coefficients of the dual objective
%          function to maximize.
%        c A sum(blockLengths)X1 vector containing the coefficients of the
%          primal objective function to minimize. With regard to the c_i in
%          the problem formulation, c=[c_1;c_2;...;c_n].
% blockLengths An nX1 vector containing the length of each of the cone
%          constraints. That is the first primal cone constraint is 
%          x(1)-norm(x(2:blockLengths(1)))>=0. The second is 
%          x(blockLengths(1)+1)-norm(x((blockLengths(1)+1):(blockLengths(1)+blockLengths(2)-1)))
%          and so on when considering the full stacked x vector (all of the
%          x_i stacked).
%   params An optional structure of parameters controlling the behaviour of
%          the algorithm. Possible members of the structure are:
%          gammaScal
%          'delta' A onstant parameter affecting the scale of a line search
%                  in the algorithm. 0<delta<1. The default if omitted is
%                  3/4.
%          'mu0'   The initial value of a positive real valued function
%                  affecting the smoothing function in the algorithm and
%                  thus the convergence. The default value if omitted is
%                  1/100.
%          'sigma' A constant parameter affecting the comvergence of the
%                  algorithm with delta. 0<sigma<1. The default if omitted
%                  is 1/4.
%          'maxIter' The maximum possible number of iterations to perform.
%                  The default if omitted is
%          'epsNormH' A tolerance value for declaring the convergence of
%                  the algorithm. This is the maximum value of norm(H) (H
%                  is defined in Equation 6 of [1]) below which convergence
%                  is declared. The default if omitted is 1e-6. Setting to
%                  zero will mean convergence is declared only when
%                  the minProgress contraint ocurs or the maximum number of
%                  iterations elapses.
%          'minProgress' A value such that the algorithm terminates if the
%                  change in normH is less than this value during multiple
%                  iterations. This allows termination when further
%                  iterations make little change. The default if omitted or
%                  an empty matrix is passed is epsNormH/(2*maxIter).
%                  Setting this to zero will make termination depend only
%                  on epsNormH and maxIter.
%
%OUTPUTS: x The sum(blockLengths)X1 set of stacked vectors solving the
%           primal optimization problem. x=[x_1;...;x_n] or an empty matrix
%           if an error occurred.
%         y The mX1 vector solving the dual optimization problem or an
%           empty matrix if an error occurred.
%         s The sum(blockLengths)X1 vector in the dual optimization
%           equality constraint or an empty matrix if an error occurred.
%  exitCode A value indicating how the algorithm terminated. Possible
%           values are
%           0 The accuracy constraint on norm(H) is satisfied (the
%             algorithm converged).
%           1 The maximum number of iterations elapsed without convergence.
%           2 Termination occurred because insufficient progress was made
%             (this depends on the params.minProgress value).
%           3 Termination occurred due to numerical problems.
%
%This function implements the algorithm described in [1]. The
%implementation slightly differs from the paper in that it can handle
%multiple cone constraints. The implementation does not take advantage of
%any sparsity in A.
%
%REFERENCES:
%[1] X. Chi and S. Liu, "A one-step smoothing Newton method for second-
%    order cone programming," Journal of Computational and Applied
%    Mathematics, vol. 223, no. 1, pp. 114-123, 1 Jan. 2009.
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%ni(i) is the size of the ith block.
%Implementing algorithm 3.1 in the paper....

numBlocks=length(blockLengths);%The number of constraints.
m=size(A,1);
k=size(A,2);%The length of x and b.

%Get a cell array where blockIndices{i} holds the indices for the elements
%of the ith block.
blockIndices=makeIndexList(blockLengths);

%Default parameters and initialization before step 1.
gammaScal=0.95;
delta=3/4;%A constant.
sigma=1/4;%A constant
mu0=1/100;%A positive value.
maxIter=100;
epsNormH=1e-6;%Convergence criterion on the norm of H.
%Termination criterion for taking a minimum stepsize
minProgress=epsNormH/(2*maxIter);

%If the user provided any parameters
if(nargin>4&&~isempty(params))
    if(isfield(params,'gammaScal'))
        gammaScal=params.gammaScal;
    end
    
    if(isfield(params,'delta'))
        delta=params.delta;
    end
    
    if(isfield(params,'sigma'))
        sigma=params.sigma;
    end
    
    if(isfield(params,'mu0'))
        mu0=params.mu0;
    end
    
    if(isfield(params,'maxIter'))
        maxIter=params.maxIter;
    end
    
    if(isfield(params,'epsNormH'))
        epsNormH=params.epsNormH;
    end
    
    if(isfield(params,'minProgress'))
        minProgress=params.minProgress;
    end
end

%Initialization before step 1
mu=mu0*ones(numBlocks,1);
zBar=[zeros(k+m,1);mu];%A constant

%We want initial x values that satisfy the cone constraints. We could set
%then all zero, but that would put it right on the boundary, so we will
%make all of the first elements of the cone constraints 1 and the rest
%zero, which puts it a bit away from the boundary.
x=zeros(k,1);
curIdx=1;
for curBlock=1:numBlocks
   x(curIdx)=1;
   curIdx=curIdx+blockLengths(curBlock);
end

%Initial dual variables.
y=zeros(m,1);%Dual variables initially just zero.
s=c-A'*y;%s must satisfy A'*y+s=c.
H=HFun(x,s,mu,A,b,blockIndices);
normH=norm(H);
normHPrev=Inf;

lkPrev=0;

%gammaVal satisfies both of the constraints for gamma: It is less than 1
%(as gammaScal is less than 1) and it is less than 1/normH.
gammaVal=gammaScal*min(1,1/normH);

if(any(~isreal(H)|~isfinite(H)))
    x=[];
    y=[];
    s=[];
    exitCode=3;%Numerical errors.
    return; 
end

for curIter=1:maxIter
%%Step 1
    if(normH<epsNormH)%If it converged.
        exitCode=0;
        return;
    elseif(normHPrev-normH<minProgress)
        %If normH did not sufficiently change between iterations, then
        %stop.
        exitCode=2;
        return;
    end
    
    %Rho is defined between equations 7 and 8.
    rho=gammaVal*normH*min(1,normH);
    
%%Step 2
    gradH=HFunGrad(x,s,mu,A,blockIndices);
    [HGradSqrt,didFail]=chol(gradH,'lower');
    if(didFail)%Using pinv for stability if poorly conditioned.
        deltaZ=lsqminnorm(gradH,(rho*zBar-H));
    else%Otherwise, solve using forward and backward substitution.
        opts.LT=true;
        opts.UT=false;
        %Solve using forward substitution.
        temp=linsolve(HGradSqrt,rho*zBar-H,opts);
        opts.UT=true;
        opts.LT=false;
        %Solve using backward substitution.
        deltaZ=linsolve(HGradSqrt',temp,opts);
    end
    deltaX=deltaZ(1:k);
    deltaY=deltaZ((k+1):(k+m));
    deltaMu=deltaZ((k+m+1):end);
    
%Steps 3 and 4
    %Search for the smallest allowable integer lk >=0. In this instance,
    %rather than just searching from 0 and going upwards, we start at the
    %last found value. In large problems, this can make a difference as lk
    %does not necessarily change much and is often not just 0.

    lk=lkPrev;
    lambda=delta^lk;
    xNew=x+lambda*deltaX;
    yNew=y+lambda*deltaY;
    muNew=mu+lambda*deltaMu;
    sNew=c-A'*yNew;
    HNew=HFun(xNew,sNew,muNew,A,b,blockIndices);
    normHNew=norm(HNew);
    boolValPrev=normHNew<=(1-sigma*(1-gammaVal*mu0)*lambda)*normH;
    
    while(1)
        if(boolValPrev==true)
            if(lk==0)
                break;
            end
            lk=lk-1;
        else
            lk=lk+1;
        end
        lambda=delta^lk;
        xNew=x+lambda*deltaX;
        yNew=y+lambda*deltaY;
        muNew=mu+lambda*deltaMu;
        sNew=c-A'*yNew;
        HNew=HFun(xNew,sNew,muNew,A,b,blockIndices);
        normHNew=norm(HNew);
        boolVal=normHNew<=(1-sigma*(1-gammaVal*mu0)*lambda)*normH;
        
        if(boolVal==true && boolValPrev==false)
            break;
        end
        boolValPrev=boolVal;
    end
    
    %We subtract one, because the loop above wants to try to bound the
    %correct lk value.
    lkPrev=max(lk-1,0);
    
    x=xNew;
    y=yNew;
    mu=muNew;
    s=sNew;
    H=HNew;
    normHPrev=normH;
    normH=normHNew;
end

%Maximum number of iterations completed.
exitCode=1;

end

function blockIndices=makeIndexList(blockLengths)
%%MAKEINDEXLIST This takes the vector of block lentghs and returns a cell
%               array where the ith element contains the indices of the
%               components in the ith block.

numBlocks=length(blockLengths);

blockIndices=cell(numBlocks,1);
curIdx=1;
for curBlock=1:numBlocks
    blockIndices{curBlock}=curIdx:(curIdx+blockLengths(curBlock)-1);
    curIdx=curIdx+blockLengths(curBlock);
end
end

function H=HFun(x,s,mu,A,b,blockIndices)
%%HFUN This implements the function H(z) from Equation 6. Note that s is
%%used instead of y, because c-A'*y is supposed to be passed as the second
%%argument to phiFun and that is equal to s.

H=[b-A*x;
   phiFun(x,s,mu,blockIndices);
   mu];
end

function gradG=HFunGrad(x,s,mu,A,blockIndices)
%%GRADH This function implements the gradient of H(z) in Equation 10. It
%       has been modified to handle more than one constraint.

m=size(A,1);
k=size(A,2);
numBlocks=length(blockIndices);

M=zeros(k,k);
N=zeros(k,k);
P=zeros(k,numBlocks);

%Allocate space and zero.
gradG=zeros(m+k+numBlocks,m+k+numBlocks);

gradG(1:m,1:k)=-A;%The first m rows.

for curBlock=1:numBlocks
    sel=blockIndices{curBlock};
    xCur=x(sel);
    sCur=s(sel);
    muCur=mu(curBlock);
    
    %e is defined toward the beginning of section 2.
    blockLength=length(sel);
    e=zeros(blockLength,1);
    e(1)=1;
    
    w1Bar=xCur+muCur*sCur;%s=c-A'*y
    w2Bar=muCur*xCur+sCur;
    wBar=vectorSqrt(vectorSquared(w1Bar)+vectorSquared(w2Bar)+2*muCur^2*e);
    
    %Lw is defined at the beginning of Section 2. Lw should be positive
    %definite, but is often pooly conditioned, so we use a pseudoinverse.
    Lw=Lx(wBar);
    [LwChol,cholDidFail]=chol(Lw,'lower');
    if(cholDidFail)
        LwInv=pinv(Lx(wBar));
    else
        opts.LT=true;
        opts.UT=false;
        %Solve using forward substitution.
        temp=linsolve(LwChol,eye(size(Lw)),opts);
        opts.UT=true;
        opts.LT=false;
        %Solve using backward substitution.
        LwInv=linsolve(LwChol',temp,opts);
    end

    Lw1=Lx(w1Bar);
    Lw2=Lx(w2Bar);
    
    M(sel,sel)=(1+muCur)*eye(blockLength,blockLength)-LwInv*(Lw1+muCur*Lw2);
    N(sel,sel)=(1+muCur)*eye(blockLength,blockLength)-LwInv*(muCur*Lw1+Lw2);
    P(sel,curBlock)=xCur+sCur-LwInv*(Lw1*sCur+Lw2*xCur+2*muCur*e);    
end

gradG((m+1):(m+k),1:k)=M;
gradG((m+1):(m+k),(k+1):(k+m))=-N*A';
gradG((m+1):(m+k),(k+m+1):(k+m+numBlocks))=P;
gradG((m+k+1):(m+k+numBlocks),(m+k+1):(m+k+numBlocks))=eye(numBlocks,numBlocks);

end

function val=Lx(x)
%%LX Create the arrow-shaped matrix Lx that is given toward the beginning
%    of section 2. x must be a single set of cone variables.

    x0=x(1);
    xBar=x(2:end);
    xBarLen=length(xBar);
    
    val=[x0,xBar';
         xBar,x0*eye(xBarLen,xBarLen)];
end


function phi=phiFun(x,s,mu,blockIndices)
%%PHIBAR This implements phi(x,s,mu) from Equation 4. It has been modified
%%to deal with multiple blocks. in this instance, the value for each block
%%is stacked.

numBlocks=length(blockIndices);

phi=zeros(length(x),1);

%jordan_sqrt(jordan_square
for curBlock=1:numBlocks
    sel=blockIndices{curBlock};
    
    %e is defined toward the beginning of section 2.
    numIdx=length(sel);
    e=zeros(numIdx,1);
    e(1)=1;
    
    muCur=mu(curBlock);
    xCur=x(sel);
    sCur=s(sel);

    phi(sel)=(1+muCur)*(xCur+sCur)-vectorSqrt(vectorSquared(xCur+muCur*sCur)+vectorSquared(muCur*xCur+sCur)+2*muCur^2*e);
end

end

function xSquared=vectorSquared(x)
%%VECTORSQUARED In Section 2, a definition of the square a second order
%cone vector in terms of a spectral factorization is given. It is
%%implemented in this function.

[lambda,u]=vectorSpectralFact(x);

xSquared=lambda(1)^2*u(:,1)+lambda(2)^2*u(:,2);

end

function sqrtX=vectorSqrt(x)
%%VECTORSQRT In Section 2, a definition of the square root of a second
%%order cone vector in terms of a spectral factorization is given. It is
%%implemented in this function.

[lambda,u]=vectorSpectralFact(x);

%sqrtX should be real.
sqrtX=sqrt(lambda(1))*u(:,1)+sqrt(lambda(2))*u(:,2);

end

function [lambda,u]=vectorSpectralFact(x)
%%SPECTRALFACT In section 2, a spectral factorization of a second order
%       cone vector is defined. The vector is factored into spectral values
%       lambda and spectral vectors u. This function computes the
%       factorization.

x0=x(1);
x1=x(2:end);
x1Mag=norm(x1);

%lambda should be positive. However, finite precision errors might push it
%negative.
lambda=abs([x0-x1Mag;x0+x1Mag]);

x1Norm=x1/x1Mag;
if(~all(isfinite(x1Norm)))
    x1Norm=zeros(size(x1));
    x1Norm(1)=1;
end

u=(1/2)*[[1;-x1Norm],[1;x1Norm]];
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
