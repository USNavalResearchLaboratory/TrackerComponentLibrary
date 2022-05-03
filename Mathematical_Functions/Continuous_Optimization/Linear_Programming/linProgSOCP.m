function [optCost,xOpt,info]=linProgSOCP(A,b,ALeq,bLeq,c,maximize,nonnegConst,params)
%%LINPRODSOCP Used a second-order cone programming algorithm to solve a
%             linear programming problem involving equality and/or
%             inequality constraints.  This solves the problem:
%                        minimize (maximize) c'*x
%                        given     A*x=b
%                              ALeq*x<=bLeq
%                                   x>=0 (if nonnegConst is true)
%            Second order cone programming is a generalization of linear
%            programming. The function splittingConicSolver is used to
%            perform the optimization.
%
%INPUTS: A An mXn matrix of equality constraints. If there are no equality
%          constraints, then use an empty matrix.
%        b A mX1 matrix of the right-hand of the equality constraints, or
%          an empty matrix if there are no equality constraints.
%     ALeq An mLeqXn matrix of the inequality constraints. If there  are no
%          inequality constraints, then use an empty matrix.
%     bLeq A mLeqX1 matrix of the right-hande of the inequality
%          constraints. or an empty matrix if there are no inequality
%          constraints.
%        c The nX1 cost vector. Unlike with linProgRevisedSimplex and
%          linProgInteriorPoint, this must contain all finite values.
% maximize A boolean variable specifying whether the problem is to maximize
%          or minimize the cost function. The default if omitted or an
%          empty matrix is passed is false.
% nonnegConst A parameter indicating whether the nonegtivity constraint
%          should be included. The default if omitted or an empty matrix is
%          passed is true.
%   params Parameters that affect how the function splittingConicSolver,
%          which is used by this function, works. See the comments to
%          splittingConicSolver for more information.
%
%OUTPUTS: optCost The optimal cost when the algorithm successfully
%                 terminates or when it terminates after having reached the
%                 maximum number of iterations. This is the minimum/maximum
%                 value of c'*x. If the algorithm did not terminate
%                 successfully, then an empty matrix is returned. If the
%                 cost is unbounded, then Inf or -Inf is returned. An INf
%                 or -Inf can also be returned if accuracy errors prevented
%                 the algorithm from getting a solution.
%            xOpt The optimal vector x associated with the optimal cost. If
%                 the problem is unbounded or the algorithm did not
%                 terminate successfully, then an empty matrix is returned.
%            info The problem is formulated as a second-order cone
%                 optimization problem. This is the information regarding
%                 the termination state returned by the
%                 splittingConicSolver function.
%
%This just formulates the inputs in an appropriate manner so that the
%function splittingConicSolver can be used to solve the problem.
%
%EXAMPLE:
%%Consider the problem
%maximize x4
%given 3*x1-2*x2-x3-x4=0
%       x1-2*x2+x3=5
%       -2*x1-4*x2+x3<=-1
%       x1+12*x2-5*x3+x4<=4
%       x1>=0
%       x2>=0
%
% c=[0;0;0;1];
% A=[3,-2,-1,-1;
%    1,-2,1,0];
% b=[0;5];
% ALeq=[-2,-4,1,0;
%        1, 12,-5,1];
% bLeq=[-1;4];
% maximize=true;
% [optCost,xOpt,info]=linProgSOCP(A,b,ALeq,bLeq,c,maximize)
%One gets the solution optCost=8.6 and xOpt=[3;0;1.6;8.6];
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(maximize))
    maximize=false;
end

if(nargin<7||isempty(nonnegConst))
   nonnegConst=true; 
end

if(nargin<8||isempty(params))
   params=[]; 
end

if(maximize)
   c=-c; 
end

numDim=length(c);
numEqConst=length(b);
numIneqConst=length(bLeq);

if(nonnegConst)
    F=zeros(numEqConst+numDim+numIneqConst,numDim);
    g=zeros(numEqConst+numDim+numIneqConst,1);
else
    F=zeros(numEqConst+numIneqConst,numDim);
    g=zeros(numEqConst+numIneqConst,1);
end

%First, add the equality constraints.
cone=[];
cone.f=numEqConst;
F(1:numEqConst,:)=A;
g(1:numEqConst)=b;
curRow=numEqConst+1;

%Next, we will add the non-negativity constraints on all of the x terms.
if(nonnegConst)
    sel=curRow:(curRow+numDim-1);
    F(sel,:)=-eye(numDim,numDim);
    g(sel)=0;
    curRow=curRow+numDim;
end

%Finally, we will add the inequality constraints.
sel=curRow:(curRow+numIneqConst-1);
F(sel,:)=ALeq;
g(sel)=bLeq;

%The total number of inequality constraints.
if(nonnegConst)
    cone.l=numDim+numIneqConst;
else
    cone.l=numIneqConst;
end

[xOpt,~,~,info]=splittingConicSolver(F,g,c,cone,params);

optCost=info.pobj;

if(maximize)
    optCost=-optCost;
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
