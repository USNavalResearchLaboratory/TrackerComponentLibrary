function [x,f,info]=convexQuadProgSOCP(G,a,C,b,numEqConst,params)
%%CONVEXQUADPROGSOCP Perform quadratic programming on a convex problem.
%                Specifically, this algorithm solves the optimization
%                problem 
%                minimize_x a'*x+(1/2)*x'*G*x
%                such that C(:,1:numEqConst)'*x=b
%                      and C(:,(numEqConst+1):end)'*x>=b
%                Convexity means that G is positive (semi)definite. This
%                reformulates the problem as a second order cone problem
%                and then uses splittingConicSolver to solve the problem.
%                Unlike the active set method of [2], this formulation
%                tends to converge slowly on large problems (large being,
%                for example, 81 dimensions).             
%
%INPUTS: G An nXn real, positive definite symmetric matrix (n>0).
%        a An nX1 real vector.
%        C An nXm real matrix (m>=0). The first numEqConst constraints (the
%          equality constraints) must be linearly independent or the
%          algorithm will identify the problem as infeasible.
%        b An mX1 vector.
% numEqConst The number of constraints in C that are equality constraints
%          (numEqConst<=m). If this parameter is omitted or an empty matrix
%          is passed, the default numEqConst=0 is used.
%   params Parameters that affect how the function splittingConicSolver,
%          which is used by this function, works. See the comments to
%          splittingConicSolver for more information.
%
%OUTPUTS: x The nX1 solution to the optimization problem or an empty matrix
%           if the problem is infeasible or unbounded.
%         f The value of the objective function at x, or an empty matrix
%           if the problem is infeasible or unbounded.
%      info The problem is formulated as a second-order cone optimization
%           problem. This is the information regarding the termination
%           state returned by the splittingConicSolver function.
%
%As shown in [1], the equivalent second order cone optimization problem is
%minimize t
%subject to norm(SG*x+SG\a)<=t
%       and C(:,1:numEqConst)'*x=b(1:numEqConst)
%       and C(:,(numEqConst+1):end)'*x>=b((numEqConst+1):end)
%       where SG=sqrtm(G)
%The function splittingConicSolver is used to solve the problem.
%
%EXAMPLE 1:
%This is the example worked out step-by-step in [2]:
% G=[4, -2;-2, 4];
% b=[0;0;2];
% a=[6;0];
% C=[1,0,1;
%    0,1,1];  
% [x,f,info]=convexQuadProgSOCP(G,a,C,b)
%The optimal solution can be seen to be x=[0.5;1.5] and f=6.5;
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
% [x,f,exitCode]=convexQuadProgSOCP(G,a,C,b,numEqConst)
%The optimal solution is x=[0.4;1.6] and f=6.56.
%
%REFERENCES:
%[1] S. Lobo, Miguel, L. Vandenberg, S. Boyd, and H. Lebret, "Applications
%    of second-order cone programming," Linear Algebra and its
%    Applications, vol. 284, no. 1, pp. 193-228, 1998.
%[2] D.Goldfarb and A. Idnani, "A numerically stable dual method for
%    solving strictly convex quadratic programs," Mathematical Programming,
%    vol. 27, no. 1, pp. 1-33, Sep. 1983.
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(numEqConst))
   numEqConst=0; 
end

if(nargin<6)
   params=[]; 
end

if(numEqConst>size(C,1))
   error('numEqConst is larger than the number of constraints in C')
end

P0=G;
%The real command can help deal with finite precision issues.
P0Root=real(sqrtm(P0));
q0=a;
Aeq=C(:,1:numEqConst)';
beq=b(1:numEqConst);
ALeq=-C(:,(numEqConst+1):end)';
bLeq=-b((numEqConst+1):end);

numIneqConst=length(bLeq);

%The extra variable t is appended to the end.
numDim=length(a);
c=zeros(numDim+1,1);
c(end)=1;

%Allocate space.
numVar=numEqConst+numIneqConst+numDim+1;
F=zeros(numVar,numDim+1);
g=zeros(numVar,1);

%First, the equality constraints.
F(1:numEqConst,1:numDim)=Aeq;
g(1:numEqConst)=beq;
curRow=numEqConst+1;
cone=[];
cone.f=numEqConst;

%Next, the inequality constraints
sel=curRow:(curRow+numIneqConst-1);
F(sel,1:numDim)=ALeq;
g(sel)=bLeq;
cone.l=numIneqConst;
curRow=curRow+numIneqConst;

%Finally, the quadratic constraint.
F(curRow,:)=-c';
g(curRow)=0;
curRow=curRow+1;

sel=curRow:(curRow+numDim-1);
F(sel,1:numDim)=-P0Root;
%g(sel)=P0Root\q0;
g(sel)=pinv(P0Root)*q0;
cone.q=numDim+1;

[x,~,~,info]=splittingConicSolver(F,g,c,cone,params);

if(~isempty(x))
    x=x(1:numDim);%Get rid of the auxiliary variable.
    f=a'*x+(1/2)*x'*G*x;
else
    f=[];
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
