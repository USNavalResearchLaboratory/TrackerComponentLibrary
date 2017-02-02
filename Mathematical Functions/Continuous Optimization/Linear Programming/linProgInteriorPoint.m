function [optCost,xOpt,exitFlag]=linProgInteriorPoint(A,b,ALeq,bLeq,c,maximize,maxIter,epsilon,alpha)
%%LINPROGINTERIORPOINT   Use the primal-dual interior point algorithm to
%                        solve a linear programming problem involving
%                        equality and/or inequality constraints. This
%                        solves the problem:
%                        minimize (maximize) c'*x
%                        given     A*x=b
%                              ALeq*x<=bLeq
%                                   x>=0
%                        The worst-case runtime of the interior point
%                        method is polynomial, though the algorithm is
%                        often slower than the simplex method, which has an
%                        exponential worst-case runtime.
%
%INPUTS: A An mXn matrix of equality constraints. If there are no equality
%          constraints, then use an empty matrix.
%        b A mX1 matrix of the right-hand of the equality constraints, or
%          an empty matrix if there are no equality cosntraints.
%     ALeq An mLeqXn matrix of the inequality constraints. If there  are no
%          inequality constraints, then use an empty matrix.
%     bLeq A mLeqX1 matrix of the right-hande of the inequality
%          constraints. or an empty matrix if there are no inequality
%          constraints.
%        c The nX1 cost vector.
% maximize A boolean variable specifying whether the problem is to maximize
%          or minimize the cost function. The default if omitted or an
%          empty matrix is passed is false.
%  maxIter An optional parameter specifying the maximum number of
%          iterations to use. The default if omitted or an empty matrix is
%          passed is 5000. Complicated problems might require additional
%          iterations. However, finite precision limitations can also cause
%          large problems to get stuck and not terminate.
%  epsilon An optional parameter specifying a tolerance for declaring the
%          dot product of a dual vector with x as zero, which defines when
%          the algorithm has converged. The default value if omitted or an
%          empty matrix is passed is 1e-12.
%   alpha  A scaling value for the stepsizes. This must be between 0 and 1.
%          This affects convergence. Normally, the default value sufficies,
%          but some specific problems might not work well with the default
%          value. If omitted or an empty matrix is passed, the default
%          value of 0.5 is used.
%
%OUTPUTS: optCost The optimal cost when the algorithm successfully
%                 terminates or when it terminates after having reached the
%                 maximum number of iterations. This is the minimum/maximum
%                 value of c'*x. If the algorithm did not terminate
%                 successfully, then an empty matrix is returned. If the
%                 cost is unbounded, then Inf or -Inf is returned.
%            xOpt The optimal vector x associated with the optimal cost. If
%                 the problem is unbounded or the algorithm did not
%                 terminate successfully, then an empty matrix is returned.
%        exitFlag A flag indicating the status upon termination. It can be
%                 0 Successful termination.
%                 1 The cost is unbounded.
%                 2 Maximum number of iterations reached.
%
%The c vector can contain infinite elements to force a variable x to be
%zero. This is handled by finding such elements, removing them, and running
%the optimization on the reduced problem. Note that it is not checked
%whether setting the forbidden assignments to zero is inconsistent with the
%constraints.
%
%The constraints can be redundant. Redundant constraints are not checked
%for consistency; rather they are simply removed. Thus, solutions might be
%obtained for infeasible problems.
%
%The algorithm is the infeasible primal-dual interior point method taken
%from Chapter 9.5 of [1] and has not been optimized for sparse matrices.
%Thus, on large, sparse problems, the algorithm can be slow.
%
%As an example, consider the 2D assignment problem. Given a rectangular
%cost matrix with more columns than rows, the goal is to choose one element
%per row, and at most one element per column so as to minimize the sum of
%the costs. Such a problem can be solved efficiently using the function
%[col4row, row4col, gain]=assign2D(C,false). However, it can also be solved
%using using a general linear programming solver, because it has been
%proven that the 0-1 constraint of assignments can be replaced with a
%non-negaitvity constraint without changing the optimal point.
%For example:
% 
% C=[Inf,  2,Inf,Inf,3;
%      7,Inf, 23,Inf,Inf;
%     17, 24,Inf,Inf,Inf;
%    Inf,  6, 13, 20,Inf];
% %Turn the 2D assignment formulation into a linear programming
% %formulation.
% c=C(:);
% %The inequality constraints specify that each column is assigned to at
% %most one row.
% ALeq=kron(eye(5),ones(1,4));
% bLeq=ones(5,1);
% %The equality constraints spacify that each row is assigned to precisely
% %one column.
% A=kron(ones(1,5),eye(4));
% b=ones(4,1);
% %The problem is solved using
% [optCost,xOpt,exitFlag]=linProgInteriorPoint(A,b,ALeq,bLeq,c)
% %To make the solution clearer, one can reshape the xOpt vector
% x=reshape(xOpt,4,5)
% %And the one's will be in the corresponding positions of the assigned
% elements of C, which can be compared to the output of the assign2D
% function.
%
%REFERENCES:
%[1] D. Bertsimas and J. N. Tsitsiklis, Introduction to Linear
%    Optimization. Belmont, MA: Athena Scientific, 1997.
%
%June 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<9||isempty(alpha))
    alpha=0.5;%An arbitrary value between 0 and 1.
end

if(nargin<8||isempty(epsilon))
    epsilon=1e-12;
end

if(nargin<7||isempty(maxIter))
    maxIter=10000;
end

if(nargin<6||isempty(maximize))
    maximize=false;
end

m=size(A,1);
nTotal=length(c);

%If inequality constraints exist, add slack variables to turn them into
%equality constraints. That is, Ax<=b becomes A*x+s=b.
if(~isempty(ALeq))
    numIneq=size(ALeq,1);
    
    %The extra slack variables have zero cost.
    c=[c;zeros(numIneq,1)];
    
    %The A matrix is enlarges for the slack variables. The slack variables
    %to dot change the existing equality constraints in A.
    if(~isempty(A))
        A=[A,   zeros(m,numIneq);
           ALeq, eye(numIneq)];
        b=[b;bLeq];
    else
        A=[ALeq,eye(numIneq)];
        b=bLeq;
    end
    nTotal=nTotal+numIneq;
    m=m+numIneq;
else
    numIneq=0;
end

if(maximize)
   c=-c; 
end

%If infinite costs are being passed in c, either the problem is unbounded,
%or those assignments are forbidden, in which case we can eliminate those
%elements from the problem and solve a reduced problem.
if(any(c==-Inf))%If the problem is unbounded.
    exitFlag=1;
    optCost=(-1)^(maximize)*(-Inf);
    xOpt=[];
    return;
elseif(any(c==Inf))%Reduce the problem to a subproblem without the
                   %forbidden assignments.
    removedxEls=(c==Inf);
    c=c(~removedxEls);
    A=A(:,~removedxEls);
    n=nTotal-sum(removedxEls);
else
    removedxEls=[];
    n=nTotal;
end

%Before the algorithm, we will get rid of redundant constraints. We will
%not check that the constraints are consistent.
indepRows=identifyIndepCols(A');
if(length(indepRows)~=m)
    m=length(indepRows);
    A=A(indepRows,:);
    b=b(indepRows);
end

%Step 1 (Initialization): The infeasible primal dual algorithm is used, so
%the initialization only has to satify a few non-negativity constraints for
%x and s.
%x>0 is required, so we set it to all ones. This is most likely infeasible.
x=ones(n,1);

%Initialize the dual to a solution that is usually not feasible but that
%satisfies the requirement that s be non-negative.
s=ones(n,1);
p=zeros(m,1);

%The dual and primal values are both feasible (but not optimal) for the
%Langrange variable mu=0.

rho=1;%An arbitrary initial value.

didConverge=false;
%If sDotx does not change from one iteration to the next, then the
%algorithm is likely stuck. In such an instance, we will reduce the rho
%parameter. As noted in Chapter 9.5, rho is generally only changed if the
%algorithm gets stuck.
sDotxOld=Inf;
for curIter=1:maxIter
    %Step 2 (Optimality Test)
    %Check for convergence
    sDotx=dot(s,x);
    if(sDotx<epsilon)
       didConverge=true;
       break;
    end
    
    %If sDotx has not changed enough, then the algorithm is probably stuck.
    %Thus, this shrinks rho if it gets stuck.
    if(abs(sDotxOld-sDotx)<=2*eps(sDotx))
        rho=rho/2;
    else
        sDotxOld=sDotx;
    end
    
    %Step 3 (Computation of Newton directions)
    
    %The barrier parameter
    mu=rho*sDotx/n;
    X=diag(x);
    S=diag(s);
    
    %Solving Equations 9.28-9.30 using the infeasible solution method on
    %page 436.
    LH=[A,          zeros(m,m), zeros(m,n);
     zeros(n,n),    A',         eye(n,n);
     S,             zeros(n,m), X];%left-hand side
    RH=-[A*x-b;
         A'*p+s-c;
         X*S*ones(n,1)-mu*ones(n,1)];%Right-hand side
     
    %Must use a pseudoinverse, because LH can be singular.
    dVec=pinv(LH)*RH;
    
    %The three Newton directions
    dx=dVec(1:n);
    dp=dVec((n+1):(n+m));
    ds=dVec((n+m+1):end);

    %Step 4(Find step lengths)

    %First, do the minimization step
    minxRat=Inf;
    minsRat=Inf;
    for i=1:n
        if(dx(i)<0)
            xRat=-x(i)/dx(i);
            if(xRat<minxRat)
                minxRat=xRat;
            end
        end
        if(ds(i)<0)
            sRat=-s(i)/ds(i);
            if(sRat<minsRat)
                minsRat=sRat;
            end
        end
    end
   
    
    %Reduce the scaling parameter if progress cannot be made, because all
    %of the ds and dx values are positive.
    if(~isfinite(minxRat)||~isfinite(minsRat))
        rho=rho/2;
        continue;
    end

    betaP=min(1,alpha*minxRat);
    betaD=min(1,alpha*minsRat);
    %Step 5 (Solution Update)
    x=max(0,x+betaP*dx);
    p=p+betaD*dp;
    s=max(0,s+betaD*ds);
end


%Compute the optimal cost, accounting for whether the sign of c was flipped
%in the beginning.
optCost=dot(c,x);
if(maximize)
    optCost=-optCost;
end

%Put back any removed elements
if(~isempty(removedxEls))
    xOpt=zeros(nTotal,1);
    
    xIdx=1;
    for curRetEl=1:nTotal
        if(removedxEls(curRetEl)==0)
            xOpt(curRetEl)=x(xIdx);
            xIdx=xIdx+1;
        end
    end
else
    xOpt=x;
end

%If slack variables were added for inequality constraints, remove them from
%x.
xOpt=xOpt(1:(nTotal-numIneq));

if(didConverge==false)
    exitFlag=2;
    return;
end
%The algorithm converged.
exitFlag=0;
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
