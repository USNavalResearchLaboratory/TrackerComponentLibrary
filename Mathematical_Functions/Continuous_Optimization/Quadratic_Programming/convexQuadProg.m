function [x,f,exitCode]=convexQuadProg(G,a,C,b,numEqConst,epsVal,maxIter)
%%CONVEXQUADPROG Perform quadratic programming on a convex problem with
%                equality and inequality constraints. Specifically, this
%                algorithm solves the optimization problem: 
%                minimize_x a'*x+(1/2)*x'*G*x
%                such that C(:,1:numEqConst)'*x=b(1:numEqConst)
%                     and C(:,(numEqConst+1):end)'*x>=b((numEqConst+1):end)
%                A dual active set algorithm for strictly convex problems
%                is used. Strict convexity means that G is positive
%                definite. The algorithm is robust to poorly conditioned
%                matrices. Thus, to handle the semidefinite case,
%                eigenvalue thresholding is applied prior to running the
%                algorithm.
%
%INPUTS: G An nXn real, positive (semi)definite symmetric matrix (n>0).
%        a An nX1 real vector.
%        C An nXm real matrix (m>=0). The first numEqConst constraints (the
%          equality constraints) must be linearly independent or the
%          algorithm will identify the problem as infeasible.
%        b An mX1 vector.
% numEqConst The number of constraints in C that are equality constraints
%          (numEqConst<=m). If this parameter is omitted or an empty matrix
%          is passed, the default numEqConst=0 is used.
%   epsVal A small value that is used to determine whether there is no step
%          in the primal space such that a constraint becomes feasible,
%          meaning that a step must be taken in the dual space and a
%          constraint dropped. This value should be >0 due to finite
%          precision errors. If this parameter is omitted or an empty
%          matrix is passed, the default value of eps() is used.
%  maxIter The maximum number of iterations to allow. The default value if
%          this parameter is omitted or an empty matrix is passed is 200.
%
%OUTPUTS: x The nX1 solution to the optimization problem or an empty matrix
%           if the problem is infeasible.
%         f The value of the objective function at x, or an empty ,matrix
%           if the problem is infeasible.
%  exitCode A value indicating how the algorithm terminated. Possible
%           values are:
%           0 The algorithm found the optimal value x.
%           1 The maximum number of iterations elapsed.
%          -1 The optimization problem is infeasible.
%           
%The algorithm of [1] is implemented. However, the algorithmic description
%in [1] omits a few minor details regarding deleting elements from uPlus
%and updating s removing a constraint. The changes are commented in the
%code. Also, the algorithm of [1] is only described for inequality
%constraints. The change to allow equality constraints is performed by
%forcing the first numEqConst columns of C to be chosen to be added to the
%active set and then making sure that they are never removed from the
%active set. If an attempt is made to remove them from the active set, then
%the corresponding dual variables are set to zero rather than removed.
%
%Another change from the algorithm of [1] deals with the qr formulation of
%certain terms. Though the qr formulation of the H and NStar matrices is
%used, the single column deletions and insertions are not performed. Though
%slower than the algorithm in the paper, one does not need to worry about
%cumulative finite precision errors due to many incremental updates.
%
%Also, algorithm [1] did not include an epsVal term to try to deal with
%finite precision effects limiting the ability to determine whether no
%primal step can activate a certain constraint.
%
%EXAMPLE 1:
%This is the example worked out step-by-step in [1]:
% G=[4, -2;-2, 4];
% b=[0;0;2];
% a=[6;0];
% C=[1,0,1;
%    0,1,1];  
% [x,f,exitCode]=convexQuadProg(G,a,C,b)
%The optimal solution can be seen to be x=[0.5;1.5] and f=6.5;
%
%EXAMPLE 2:
%In [1], the algorithm was tested on many random problems. Below is an
%implementation of the method described in the paper for n=81 dimensions
%and m=3*n=243 constraints. The first qStar=m/3=81 constraints end up being
%enforced.
% n=81;
% m=3*n;
% qStar=m/3;
% 
% G=rand(n,n);
% G=G+G';
% %G is not symmetric and all elements are in the range of -1 to 1
% SCur=sum(G(1,2:end));
% G(1,1)=1+SCur+rand(1);
% for curRow=2:n
%    SPrev=SCur;
%    SCur=sum(G(curRow,[1:(curRow-1),(curRow+1):end]));
%    G(curRow,curRow)=G(curRow-1,curRow-1)+SCur+SPrev+rand();
% end
% 
% xStar=-0.5+rand(n,1);
% C=-1+2*rand(n,m);
% C=bsxfun(@times,C,1./sqrt(sum(C.*C,1)));
% s=zeros(m,1);
% s((qStar+1):end)=rand(m-qStar,1);
% b=-(s-C'*xStar);%Sign changed from description in paper to be feasible.
% u=30*rand(m,1);
% u((qStar+1):end)=0;
% a=C*u-G*xStar;
% [x,f,exitCode]=convexQuadProg(G,a,C,b)
%
%EXAMPLE 3:
%This is an example of a problem having both inequality and equality
%constraints:
% G=[4, -2;-2, 4];
% b=[2;1.6];
% a=[6;0];
% C=[1,0;
%    1,1];
% numEqConst=1;%One equality constraint.
% [x,f,exitCode]=convexQuadProg(G,a,C,b,numEqConst)
%The optimal solution is x=[0.4;1.6] and f=6.56.
%
%REFERENCES:
%[1] D.Goldfarb and A. Idnani, "A numerically stable dual method for
%    solving strictly convex quadratic programs," Mathematical Programming,
%    vol. 27, no. 1, pp. 1-33, Sep. 1983.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(numEqConst))
   numEqConst=0; 
end

if(numEqConst>size(C,1))
   error('numEqConst is larger than the number of constraints in C')
end

if(nargin<6||isempty(epsVal))
    epsVal=eps();
end

if(nargin<7||isempty(maxIter))
    maxIter=200;
end

%This portion is to compensate for the possibility that G is positive
%semidefinite.
[V,D]=eig(G);
d=diag(D);
dNew=max(d,2*max(size(G))*eps(max(d)));
if(~all(dNew==d))
    %Only change G if some eigenvalues were too small.
    D=diag(dNew);
    G=V*D/V;
    G=(G+G')/2;%Ensure symmetry.
end

%If there are no constraints, the return the unconstrained solution.
if(isempty(b))
    x=-G\a;
    f=(1/2)*a'*x;
    exitCode=0;
    return
end

%Step 0: Find the unconstrained minimum and initialize parameters.
H=inv(G);
x=-G\a;
f=(1/2)*a'*x;
%Values are added to and deleted from A. Though this is inefficient, it is
%questionable whether any particular data structure could make it notably
%more efficient than what Matlab does without reprogramming this function
%in C.
A=[];
q=0;
u=[];
NStar=[];
L=chol(G,'lower');

skipStep1=false;
for curIter=1:maxIter
    %Step 1: Choose a violated constant, if any. We shall choose the most
    %violated constant.
    if(skipStep1==false)
        s=C'*x-b;
        %Force the equality constraints to be added first.
        if(curIter<=numEqConst)
            p=curIter;
        else
            [~,p]=min(s);

            if(any(p==A)||s(p)>=0)
                %If no constraints are violated, the current solution is
                %feasible and optimal.
                exitCode=0;
                return;
            end
        end
    
        uPlus=[u;0];%Dual variables with one for constraint p added.
        nPlus=C(:,p);
    end

    %Step 2: Check for feasibility and determine a new S-pair
    %Step 2a) Determine the step direction.
    z=H*nPlus;
    
    %r is empty when no constraints have been added.
    if(q==0)
        r=[];
    else
        r=NStar*nPlus;
    end
    
    %Step 2b) Determine the step length
    %Step 2bi) Find t1, the maximum step without violating dual
    %          feasibility.
    t1=Inf;
    l=[];
    for i=1:q
        %The added requirement that uPlus(i)>0 arises due to the equality
        %constraints. Otherwise, all of the dual variables enforced would
        %be greater than zero due to the nature of the dual problem.
        if(r(i)>0&&uPlus(i)>0)
            t1Cur=uPlus(i)/r(i);
            if(t1Cur<t1)
                t1=t1Cur;
                l=i;
            end
        end
    end

    %Step 2bii: Find the minimum step in the primal space such that the
    %pth constraint becomes feasible.
    zProd=z'*nPlus;
    
    %This comparison is not described in [1], but appears to be necessary
    %due to finite precision errors in some instances.
    if(zProd<epsVal)
        t2=Inf;
    else
        t2=-s(p)/zProd;
    end
    t=min(t1,t2);
    
    %Step 2c: Determine a new S-pair and take the step.
    %Step 2ci: If no step exists in the primal or dual spaces.
    if(~isfinite(t)||(t==0&&~all(uPlus==0)))
        %The problem is infeasible. The t=0 condition comes from adding
        %equality constraints to inequality constraints.
        x=[];
        f=[];
        exitCode=-1;
        return;
    end
    
    if(~isfinite(t2))
        %Step 2cii: If a step exists only in the dual space.
        uPlus=uPlus+t*[-r;1];
        if(l>numEqConst)%Do not delete equality constraints
            %Delete constraint A(l).
            A(l)=[];
            %Delete the dual variables associated with the constraint.
            %This deletion is not mentioned in the steps given in the paper,
            %but is necessary.
            uPlus(l)=[];
            q=q-1;

            %x did not change, so s(p) does not need to be recomputed.

            [NStar,H]=updateNStarH(q,A,C,L);
        else
            uPlus(l)=0;
        end
        skipStep1=true;
    else
        %Step 2ciii: Take a step in primal and dual space.
        x=x+t*z;
        f=f+t*z'*nPlus*((1/2)*t+uPlus(q+1));
        uPlus=uPlus+t*[-r;1];
        
        if(t==t2)%If a full step is taken, add constraint p
            u=uPlus;
            q=q+1;
            A(q)=p;

            [NStar,H]=updateNStarH(q,A,C,L);
            skipStep1=false;
        else%If a partial step is taken, then drop constraint k=A(l)
            if(l>numEqConst)%Do not delete equality constraints
                A(l)=[];
                %Delete the dual variables associated with the constraint.
                %This deletion is not mentioned in the steps given in the
                %paper, but is necessary.
                uPlus(l)=[];
                q=q-1;

                [NStar,H]=updateNStarH(q,A,C,L);
            else
                uPlus(l)=0;
            end
            %s must be recomputed for a partial step. This is not mentioned
            %in the paper.
            s=C'*x-b;
            skipStep1=true;
        end
    end
end

%Maximum number of iterations exceeded.
exitCode=1;
end

function [NStar,H]=updateNStarH(q,A,C,L)
%UPDATENSTARH This function is based on the discussion in Section 4 of [1].
%             It uses a qr decomposition, but does not use any of the qr
%             updating methods.
    N=C(:,A);

    opts.LT=true;
    opts.UT=false;
    %Solve B=L\N;
    B=linsolve(L,N,opts);
    [Q,R]=qr(B);
    R=R(1:q,1:q);
    
    opts.LT=false;
    opts.UT=true;
    %Solve J=(L')\Q;
    J=linsolve(L',Q,opts);
    J1=J(:,1:q);
    J2=J(:,(q+1):end);
    
    H=J2*J2';
    NStar=R\(J1');
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
