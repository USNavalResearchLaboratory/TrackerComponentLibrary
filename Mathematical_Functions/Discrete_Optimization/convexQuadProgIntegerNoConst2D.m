function [x,minCost]=convexQuadProgIntegerNoConst2D(Q,b,nonNeg)
%%CONVEXQUADPROGINTEGERNOCONST2D Perform quadratic programming on a 2D
%       convex problem with no constraints except that the solutions must
%       be integers that are either unconstrained or are non-negative.
%       This functions solves the optimization problem
%         minimize_x x'*Q*x+2*b'*x
%         subject to x is a 2X1 set of integers (optionally non-negative
%                    integers).
%
%INPUTS: Q An 2X2 real, positive definite symmetric matrix.
%        b An 2X1 real vector.
%   nonNeg This is true if all the x elements should be non-negative and
%          false otherwise. The default if omitted or an empty matrix is
%          passed is false.
%
%OUTPUTS: x The 2X1 optimal solution (or an empty matrix if the optimal
%           cost is infinite).
%   minCost The cost of the optimal constrained solution. If the cost isn't
%           finite, then this is an empty matrix.
%
%One could implement the branch-and bound algorithm for large scale
%problems as in [1]. However, in the 2D case, there are so few branches, it
%isn't worth it. Instead, we simply branch as in [1] (with a
%simple modification for non-negative constraints) and we traverse all 8
%hypotheses and then choose the lowest cost one. In [1], when branching,
%one creates smaller Q and b submatrices to optimize over rather than
%redoing the entire problem with an added eualaity constraint. In the 2D
%case, fixing an index means that the subproblem is just optimizing a scalar
%quadratic equation. Thus, we don't need to explicitely form submatrices;
%we just directly evaluate the scalar solution to the quadratic.
%
%EXAMPLE 1:
%Here is a simple example. We get the unconstrained solution. Then, we get
%the integer constrained solution and see that the cost is higher. Finally,
%we get the integer constrained solution with non-negativity constraints.
%That is NOT just clipping the integer constrained solution to 0. To show
%that, we compute the cost of the clipped solution and we see that it is
%much larger than the solution returned by the function.
% Q=[32,  20;
%    20, 22]; 
% b=[-178;-1200];
% [xOpt,minCost]=convexQuadProgEqConst(Q,b)
% nonNeg=false;
% [xOptInt,minCostInt]=convexQuadProgIntegerNoConst2D(Q,b,nonNeg)
% nonNeg=true;
% [xNonNegInt,minCostNonNegInt]=convexQuadProgIntegerNoConst2D(Q,b,nonNeg)
% %We show that the zero-clipped solution is much higher than the one with
% %non-negativity constraints.
% xTest=max(xOptInt,0);
% testCost=xTest'*Q*xTest+2*b'*xTest
%
%EXAMPLE 2:
%One might be tempted the think that one can always obtain the globally
%optimal solution by solving the unconstrained problem and then just
%rouning the result or by just trying the 4 points obtained by taking the
%floor and ceiling of each element of the unconstrained solution. However,
%that is not the case. In this example, we show the unconstrained solution
%and then the optimal integer constrained differ by more than 1 index.
% Q=[25, 0.3;
%   0.3, 0.01];
% b=[-156.9;
%    -3.31];
% xUnconst=convexQuadProgEqConst(Q,b)
% [x,minCost]=convexQuadProgIntegerNoConst2D(Q,b)
% %Now, we show that the cost of just rouning up or down the non-integer part
% %of xUnconst does not produce a globally optimal solution. Only the first
% %element of xUnConst is not an integer, so we only need to consider points
% %rounding that one up or down (as opposed to 4 points rounding each
% %element).
% xJustRoundUp=ceil(xUnconst);
% costJustRoundUp=xJustRoundUp'*Q*xJustRoundUp+2*b'*xJustRoundUp
% xJustRoundDown=floor(xUnconst);
% costJustRoundDown=xJustRoundDown'*Q*xJustRoundDown+2*b'*xJustRoundDown
%
%REFERENCES:
%[1] F. Koerner, "An efficient branch and bound algorithm to solve the
%    quadratic integer programming problem," Computing, vol. 30, pp. 253-
%    260, Sep. 1983.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(nonNeg))
    nonNeg=false;
end

xOpt=convexQuadProgEqConst(Q,b);

%Impose the non-negativity constraint. When branching on each value, if the
%non-negativity constraint was enforce,d both branches will be the same.
%However, we aren't going to bother checking to see whether a branch should
%be skipped, because there are so few branches.
if(nonNeg)
    xOpt=max(xOpt,0);
end

x=[];
minCost=Inf;

%Branch on the first index.
for x1=[floor(xOpt(1)),ceil(xOpt(1))]
    %The part of the cost function that depends only on x1 and not x2
    %(which is fixed).
    constCostTerm=2*b(1)*x1+Q(1,1)*x1^2;
    %The optimal solution for x2 given x1 and no constraints.
    x2Node=-(b(2)+Q(1,2)*x1)/(Q(2,2));
    if(nonNeg)
        %Apply a non-negativity constraint.
        x2Node=max(x2Node,0);
    end

    %Now, branch on the second index. With it fixed, we just evaluate the
    %cost (There is nothing else to optimize).
    for x2=[floor(x2Node),ceil(x2Node)]
        costCur=2*b(2)*x2+2*Q(1,2)*x1*x2+Q(2,2)*x2^2+constCostTerm;
        if(costCur<minCost)
            minCost=costCur;
            x=[x1;x2];
        end
    end
end

%Now branch on the second index.
for x2=[floor(xOpt(2)),ceil(xOpt(2))]
    %The part of the cost function that depends only on x2 and not x1
    %(which is fixed).
    constCostTerm=2*b(2)*x2+Q(2,2)*x2^2;
    %The optimal solution for x1 given x2 and no constraints.
    x1Node=-(b(1)+Q(1,2)*x2)/(Q(1,1));

    if(nonNeg)
        %Apply a non-negativity constraint.
        x1Node=max(x1Node,0);
    end

    %Now, branch on the first index.
    for x1=[floor(x1Node),ceil(x1Node)]
        costCur=2*b(1)*x1+Q(1,1)*x1^2+2*Q(1,2)*x1*x2+constCostTerm;
        if(costCur<minCost)
            minCost=costCur;
            x=[x1;x2];
        end
    end
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
