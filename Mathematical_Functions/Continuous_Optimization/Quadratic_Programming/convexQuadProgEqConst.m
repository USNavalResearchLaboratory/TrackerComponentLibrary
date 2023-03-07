function [xOpt,xCost]=convexQuadProgEqConst(Q,b,D,d,algorithm)
%%CONVEXQUADPROGEQCONST Solve a convex quadratic programming problem that
%       has only equality constraints, no inequality constraints.
%       Specifically, this algorithm solves the optimization problem: 
%        minimize_x x'*Q*x+2*b'*x
%        subject to D*x=d
%       Because there are only equality constraints, there is an explicit
%       solution to the problem, unlike in the convexQuadProg function,
%       which must iteratively solve the problem with inequality
%       constraints. 
%
%INPUTS: Q An nXn real positive definite, symmetric matrix.
%        b An nX1 real vector.
%        D A numConstXn real matrix. This can be an omitted or empty matrix
%          can be passed if there are no constraints. numConst<=n. There
%          shouldn't be any redundant constraints.
%        d A numConstX1 real vector.
% algorithm An optional parameter specifying how the problem is solved.
%          Possible values are:
%          -1 (The default if omitted or an empty matrix is passed) Just
%             use the straightforward solution as derived below.
%           0-2 Transform the problem into a format that can be solved by
%             constrainedLSEq and use the same-numbered algorithm in
%             constrainedLSEq.
%
%OUTPUTS: xOpt The nX1 optimal constrained solution to x.
%        xCost The scalar cost of the solution in xOpt.
%
%The solution is a straightforward application of Lagrangian relaxation
%with equality constraints as discussed in Chapter 4.1.3 of [1]. The
%Lagrangian of this problem is: 
%L=x'*Q*x+2*b'*x+2*Lambda'*(D*a-d)
%Taking the gradient with respect to x, setting it equal to 0 and solving
%for x, one obtains
%a=-inv(Q)*(b+D'*Lambda);
%One has to find the Lambda vector that satisfies the constraint
%D*a=d.
%Thus, one substitutes in a and solved for Lambda to get
%Lambda=-inv(D*inv(Q)*D')*(d+D*inv(Q)*b)
%This is the notin behin algorithm -1 of this function.
%
%EXAMPLE:
%On a simple example problem, we show that the algorithms in the function
%all get the same result and that it is also the same as the output of
%convexQuadProg when called with no equality constraints (all equality
%limited by finite precision errors).
% Q=[32,     7,    12,    17;
%     7,    22,    17,    22;
%    12,    17,    32,    27;
%    17,    22,    27,    40];
% b=[1;2;3;4];
% D=[2,0,0,1];
% d=4.5;
% xOpt=convexQuadProgEqConst(Q,b,D,d,-1)
% xOpt=convexQuadProgEqConst(Q,b,D,d,0)
% xOpt=convexQuadProgEqConst(Q,b,D,d,1)
% xOpt=convexQuadProgEqConst(Q,b,D,d,2)
% xOpt=convexQuadProg(Q,b,D',d,1)
%
%REFERENCES:
%[1] D. P. Bertsekas, Nonlinear Programming, 3rd ed. Belmont, MA: Athena
%    Science, 2016.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(algorithm))
    algorithm=-1;
end

if(nargin<3||isempty(D))
    %The unconstrained problem.
    xOpt=-Q\b;
else
    %The constrained problem.
    switch(algorithm)
        case -1
            DQInv=D/Q;
            lambda=-(DQInv*D')\(d+DQInv*b);
            xOpt=-Q\(b+D'*lambda);
        otherwise
            A=chol(Q,'upper');
            bLS=-linsolve(A',b);
            xOpt=constrainedLSEq(A,bLS,D,d,algorithm);
    end
end

if(nargout>1)
    xCost=xOpt'*Q*xOpt+2*b'*xOpt;
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
