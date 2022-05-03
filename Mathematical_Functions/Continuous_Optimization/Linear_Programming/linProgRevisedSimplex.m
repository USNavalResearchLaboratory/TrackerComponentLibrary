function [optCost,xOpt,exitCode]=linProgRevisedSimplex(A,b,ALeq,bLeq,c,maximize,maxIter,epsilon)
%%LINPROGREVISEDSIMPLEX  Use the revised simplex algorithm to solve a
%                        linear programming problem involving equality
%                        and/or inequality constraints. This solves the
%                        problem:
%                        minimize (maximize) c'*x
%                        given     A*x=b
%                              ALeq*x<=bLeq
%                                   x>=0
%                        The worst-case runtime of the revised simplex
%                        algorithm is O(2^n), where n is the length of x.
%                        However, in practice, it is often faster than an
%                        interior point method, which has a polynomial
%                        worst-case runtime.
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
%        c The nX1 cost vector. As shown in the example below, this can
%          sometimes contain infinite values.
% maximize A boolean variable specifying whether the problem is to maximize
%          or minimize the cost function. The default if omitted or an
%          empty matrix is passed is false.
%  maxIter An optional parameter specifying the maximum number of
%          iterations to use. The default if omitted or an empty matrix is
%          passed is 5000. Complicated problems might require additional
%          iterations. However, finite precision problems can also cause
%          large problems to get stuck and not terminate.
%  epsilon An optional parameter specifying a tolerance for declaring
%          values zero. The default value if omitted or an empty matrix is
%          passed is 1e-9. Setting this parameter to zero can cause the
%          optimization problem to cycle or fail by incorrectly creating an
%          invalid basis.
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
%        exitCode A flag indicating the status upon termination. It can be
%                 0 Successful termination.
%                 1 The cost is unbounded.
%                 2 Maximum number of iterations reached.
%                 3 Finite precision problems caused an invalid basis to be
%                   formed (changing epsilon might help, or the problem
%                   might just be poorly suited for double-precision
%                   arithmetic.
%                 4 The problem is not feasible.
%                 5 The subroutine to find the initial feasible basis
%                   failed.
%
%The c vector can contain infinite elements to force a variable x to be
%zero. This is handled by finding such elements, removing them, and running
%the optimization on the reduced problem. Note that it is not checked
%whether setting the forbidden assignments to zero is inconsistent with the
%constraints.
%
%The constraints can be redundant. Redundant constraints are not checked
%for consistency.
%
%The algorithm is mostly taken from Chapter 3.3 of [1] and has not been
%optimized for sparse matrices. Thus, on large, sparse problems, the
%algorithm can be slow. Bland's rule from Chapter 3.4 is used to prevent
%cycling. The introduction of comparisons to epsilon rather than zero is
%not in the book, but was found to be necessary for even small problems.
%Preprocessing to handle inequality constrains by adding slack variables
%and turning them into equality constraints is not described in the text.
%The method of Chapter 3.5 for finding an intial feasible solution is used.
%Chapter 3.3 describes periodically recomputing the inverse basis matrix
%from the actual bases to avoid accumulating roundoff errors. This is done
%every 30th iteration. Finite precision problems can potentially lead to
%redundant basis vectors being added. If this occurs, then the algorithm
%reports as failing.
%
%As an example, consider the 2D assignment problem. Given a rectangular
%cost matrix with more columns than rows, the goal is to choose one element
%per row, and at most one element per column so as to minimize the sum of
%the costs. Such a problem can be solved efficiently using the function
%[col4row, row4col, gain]=assign2D(C,false). However, it can also be solved
%using using a general linear programming solver, because it has been
%proven that the 0-1 constraint of assignments can be replaced with a
%non-negaitvity constraint without changing the optimal point.
%
%EXAMPLE:
% C=[Inf,  2,   Inf,Inf,3;
%      7,  Inf, 23, Inf,Inf;
%     17,  24,  Inf,Inf,Inf;
%    Inf,  6,   13, 20, Inf];
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
% [optCost,xOpt,exitCode]=linProgRevisedSimplex(A,b,ALeq,bLeq,c)
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
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<8||isempty(epsilon))
    epsilon=1e-12;
end

if(nargin<7||isempty(maxIter))
    maxIter=5000;
end

if(nargin<6||isempty(maximize))
    maximize=false;
end

m=size(A,1);
nTotal=length(c);

%If inequality constraints exist, add slack variables to turn them into
%equality constraints.
if(~isempty(ALeq))
    numIneq=size(ALeq,1);
    
    %The extra slack variables have zero cost.
    c=[c;zeros(numIneq,1)];
    
    %The A matrix is enlarged for the slack variables. The slack variables
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
    exitCode=1;
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

%Before the algorithm, we could get rid of redundant constraints. However,
%doing so would also remove one constraint in a bounded region (unbounding
%the region). For example, if the constraints were of the form a<x<b, then
%getting rid of redundant constraints like
% indepRows=identifyIndepCols(A');
% if(length(indepRows)~=m)
%     m=length(indepRows);
%     A=A(indepRows,:);
%     b=b(indepRows);
% end
%would unbound the region. Consequently, we keep all constraints.

%To meet the requirement in determining the initial basis in Chapter 3.5
%that all b are positive, we flip the sign of b and the appropriate rows of
%A accordingly.
sel=b<0;
b(sel)=-b(sel);
A(sel,:)=-A(sel,:);

%Step 1 in the book.
%The above steps reformulated the problem into
%                   minimize c'*x
%                      given A*x=b
%                             x>=0
%Here, we have to determine an initial feasible basis to perform the actual
%optimization. To do this, we use the method of Section 3.5, which actually
%involves solving a different, augmented problem for which an easy
%initialization is known. The augmented problem is
%                   minimize y
%                      given A*x+y=b
%                             x>=0, y>=0
%A feasible initialization of that problem isx=0,y=b.
%Thus, first that problem is solved.
AInit=[A,   eye(m)];
cInit=[zeros(n,1);ones(m,1)];
xInit=[zeros(n,1);b];
basisIdx=(n+1):(n+m);
[xInit,basisIdx,exitCode,foundSol]=solveSimplexGivenBasis(AInit,b,cInit,xInit,basisIdx,epsilon,maxIter);
if(exitCode~=0||foundSol==false)
    optCost=[];
    xOpt=[];
    exitCode=5;%Could not get an initial feasible basis.
    return;
end

%If any of the extra variables are nonzero, then the problem was not
%feasible.
if(any(abs(xInit((n+1):end))>epsilon))
    optCost=[];
    xOpt=[];
    %The problem is infeasible, because nonzero y values ended up in the
    %basis.
    exitCode=4;
    return;
end

%We have to see whether any of the artifical variables are in the basis. If
%so, we must get rid of them and add any other vectors to the basis. This
%will not change the result as the corresponding entries in x will be zero.
basisIdx=sort(basisIdx,'ascend');

%Mark which basis vectors were not assigned.
sel=zeros(n+m,1);
sel(basisIdx)=1;
%This has the indices of unassigned basis vectors.
unassignedBases=find(~sel(1:n));

curIdx=m;
%While we are considering one of the extra vectors in the basis.
while(curIdx>0&&basisIdx(curIdx)>n)
    %The current basis index is one of the extra vectors. We will set it to
    %an unassigned basis vector. The unassigned vectors are used in order.
    newBasisIdx=(m-curIdx)+1;
    basisIdx(curIdx)=unassignedBases(newBasisIdx);
    
    curIdx=curIdx-1;
end

%We now have a feasible set of bases for the problem.
x=xInit(1:n);

[x,~,exitCode,foundSol]=solveSimplexGivenBasis(A,b,c,x,basisIdx,epsilon,maxIter);
switch(exitCode)
    case 1
        optCost=(-1)^(maximize)*(-Inf);
        xOpt=[];
        return
    case 3
       optCost=[];
       xOpt=[];
       return;
    otherwise%Assume 0, no error.
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

if(foundSol==false)
    exitCode=2;
else
    exitCode=0;
end

end

function [x,basisIdx,exitCode,foundSol]=solveSimplexGivenBasis(A,b,c,x,basisIdx,epsilon,maxIter)
    n=length(x);
    m=size(A,1);
    
    if(rcond(A(:,basisIdx))<eps())
        BInv=pinv(A(:,basisIdx));
    else
        BInv=inv(A(:,basisIdx));
    end
    
    exitCode=0;
    foundSol=false;
    for curIter=1:maxIter
        %Step 2 in the book.
        %The components of the cost function associated with the chosen
        %bases.
        cB=c(basisIdx); 
        p=cB'*BInv;
        %Compute the reduced costs.
        cRed=c-(p*A)';

        negRedCostIdx=[];
        for curRedCost=1:n
        %We are using Bland's rule to avoid cycling, as mentioned at the
        %end of Chapter 3.4. That means, choose the first index that
        %satisfies the criterion.
            if(cRed(curRedCost)<-epsilon)
                negRedCostIdx=curRedCost;
                break;
            end
        end

        %If none of the reduced costs are negative, then the optimal
        %solution in basis form has been found.
        if(isempty(negRedCostIdx))
            foundSol=true;
            break;
        end

        %Steps 3. Compute u and check for an unbounded solution (within
        %precision bounds). 
        u=BInv*A(:,negRedCostIdx);
        if(all(u<epsilon))
            exitCode=1;
            return;
        end

        %Step 4. Find the positive component of u that minimizes a ratio.
        theta=Inf;%To hold the optimal ratio value.
        for i=1:m
            if(u(i)>=epsilon)
                ratio=x(basisIdx(i))/u(i);
                %Since this is less than and not less than or equal to,
                %this satisfies the criterion for Bland's rule to avoid
                %cycling.
                if(ratio<theta)
                    theta=ratio;
                    minIdx=i;
                end
            end
        end

        %Steps 5 and 6: Every 30th iteration, the matrix BInv is rebuilt
        %from the bases and the corresponding values of x are directly
        %recomputed. Otherwise, use the method in the book for iteratively
        %updating BInv and x.
        if(mod(curIter,30)==0)
            basisIdx(minIdx)=negRedCostIdx;
            %We will examine whether duplicate elements exist in basisIdx.
            %If so, then a numerical error has occurred and the algorithm
            %should terminate.
            if(length(unique(basisIdx))~=m)
               exitCode=3;
               return 
            end
            if(rcond(A(:,basisIdx))<eps())
                BInv=pinv(A(:,basisIdx));
            else
                BInv=inv(A(:,basisIdx));
            end

            x=zeros(n,1);
            x(basisIdx)=BInv*b;
        else
            %Step 5. Update the feasible basic solution. Also, we will
            %update the basis index.
            x(basisIdx)=x(basisIdx)-theta*u;
            x(negRedCostIdx)=theta;
            basisIdx(minIdx)=negRedCostIdx;

            %Step 6. Update the BInv matrix. We have to perform parallel
            %row operations on Binv and u until u becomes a unit vector
            %with a 1 in the position of minIdx.

            %First, we will perform the operation that would turn the
            %minIdx element of u to 1.
            BInv(minIdx,:)=BInv(minIdx,:)/u(minIdx);
            %Next, we will perform the operations that would zero out all
            %of the other rows of u.
            for curRow=1:m
                if(curRow==minIdx) 
                    continue;
                end
                addVal=-u(curRow);
                BInv(curRow,:)=BInv(curRow,:)+addVal*BInv(minIdx,:);
            end
        end
    end

    %Examine whether duplicate elements exist in basisIdx.
    if(length(unique(basisIdx))~=m)
       exitCode=3;
       return 
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
