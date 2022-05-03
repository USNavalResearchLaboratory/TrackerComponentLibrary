function [x,y,s,info]=splittingConicSolver(A,b,c,cone,params)
%%SPLITTINGCONICSOLVER This function calls a library to solve many types of
%      cone programming problems. The primal problem is
%            minimize c'*x
%        subject to s=b-A*x and s in a region K.
%      where both x in in R^n and s is in R^m (m-dimensional space of real
%      numbers). K is the product of cones in the order free cone, lp cone,
%      second order cone(s), semi-definite cone(s), primal exponential
%      cones, and dual exponential cones. The meaning of the cones is
%      defined below. The equivalent dual problem is
%               min -b'*y
%        subject to -A'*y=c and y is in the dual to the cone K.
%      The algorithm used is based on the alternating direction method of
%      multipliers (ADMM).
%
%INPUTS: A, b The real mXn matrix A and mX1 real vector b in the constraint
%             s=A*x-b. A can be a sparse matrix. Thr rows of A are ordered
%             to correspond with the different cones used.
%           c The real nX1 weighting vector for the cost value c'*x being
%             maximized.
%        cone A structure specifying the cones applied for every row of A.
%             The possible elements are listed below in the order that they
%             correspond to rows in A. The sample code below shows how each
%             of the cones is used. Possible values in the structure are:
%             'f' This is the number of free variables in the primal
%                 problem (equality constraints on x with respect to b).
%                 That is, s(1:cone.f) are zero.
%             'l' This is the number of non-negative components in
%                 the primal, meaning that s((cone.f+1):cone.l)>=0.
%             'q' This is a vector that lists the dimensions of the
%                 second order cone (quadratic, Lorentz) constraints. The
%                 lengths must be >=2. For example, if cone.q=[q1;...;qn]
%                 for n constraints, then the first constraint is
%                 s(cone.f+cone.l+1)>=norm(s((cone.f+cone.l+2):(cone.f+cone.l+cone.q(1)))
%                 The second constraint is
%                 s(cone.f+cone.l+cone.q(1)+1)>=norm(s((cone.f+cone.l+cone.q(1)+2):(cone.f+cone.l+cone.q(1)+cone.q(2)))
%                 and so on.
%             's' This is a vector that lists the dimensions of the
%                 positive semidefinite constraints. However, it is assumed
%                 that all of the input data has been modified such that
%                 the off-diagonal terms have been scaled by sqrt(2). For
%                 example, for cone.s=[s1;...;sn] the first constraint is
%                 vech2Mat(s((cone.f+cone.l+sum(cone.q)+1):(cone.f+cone.l+sum(cone.q)+k1),true,1/sqrt(2))
%                 be positive semidefinite, where k1=s1*(s1+1)/2 is the
%                 size of one side of the square matrix. The second such
%                 constraint would be that
%                 vech2Mat(s((cone.f+cone.l+sum(cone.q)+k1+1):(cone.f+cone.l+sum(cone.q)+k1+k2),true,1/sqrt(2))
%                 be positive semidefinite, where k2=s2*(s2+1)/2, and so
%                 on.
%            'ep' The number of exponential cones. Like all of the
%                 previous cones, theses only act on elements of s after
%                 those used by the previous cones. Each exponential cone
%                 takes three rows of A. Calling the s values taken in
%                 order as (a,b,c), the exponential cone is b*exp(a/b)<=c
%                 and b>0.
%            'ed' The number of dual exponential cones. These are
%                 similar to the exponential cones in that three rows of A
%                 are taken for each cone. For three values in s being
%                 (a,b,c), the dual exponential cone is
%                 -a*exp(b/a)<=exp(1)*c and a<0
%      params This is an optional structure of parameters that affect how
%             the algorithm runs. Possible elements are
%             'alpha' The over relaxation parameter to use. This can
%                  be between 0 and two. The default value is 1.5.
%             'rho_x' The momentum of the x term in the
%                  optimization. The default value is 1e-3.
%             'max_iters' The maximum number of alternating
%                  direction method of multipliers (ADMM) iterations. The
%                  default if omitted is 2500.
%             'eps' The convergence tolerance. The default value is
%                  1e-9.
%             'verbose' Indicates whether verbose mode should be
%                  used. The default value is false. If verbose mode is on,
%                  text with a lot of information will be displayed.
%             'normalize' Indicates whether a heuristic data
%                  rescaling method should be used. The default value is
%                  false.
%             'scale' The factor by which the data should be
%                  rescaled if params.normalize=true. The default value is
%                  1.
%             'cg_rate' The rate at which the conjugate gradient
%                  tolerance is tightened. The default value is 2.
%             'use_indirect' If true, then the indirect algorithm will be
%                  used. The indirect algorithm often requires more
%                  iterations, though it is available in a compiled
%                  version.
%
%OUTPUTS: x The nX1 solution to the primal problem. An empty matrix is
%           returned if the algorithm fails or the problem is infeasible or
%           unbounded.
%         y The mX1 solution to the dual problem. An empty matrix is
%           returned if the algorithm fails or the problem is infeasible or
%           unbounded.
%         s The s value on wehich the constraint s=A*x-b and the conic
%           constraints were applied. Note that the off-diagonal elements
%           of s corresponding to semidefinite constraints must be rescaled
%           by a factor of 1/sqrt(2) due to the requirements for the
%           scaling of the parameters going into this function. An empty
%           matrix is returned if the algorithm fails or the problem is
%           infeasible or unbounded.
%      info A structure providing information regarding the termination
%           state of the algorithm. The fields are
%           'iter' The number of iterations taken.
%           'status' A string describing how the algorithm terminated.
%             Ideally, this is 'Solved' and not 'Solved/Inaccurate',
%             'Infeasible', 'Unbounded', 'Infeasible/Inaccurate', or
%             'Unbounded/Inaccurate' when using the compiled indirect
%             algorithm. When using the non-compiled algorithms, this can
%             only be 'Solved', 'Undetermined', 'Infeasible' or
%             'Unbounded'.
%           'pobj' The value of the primal objective function.
%           'dobj' The value of the dual objective function.
%           'resPri' The primal equality residual.
%           'resDual' The dual equality residual.
%           'resInfeas' The residual of the infeasibility certification
%                   (Only returned if compiled).
%           'resUnbdd' The unbounded certification residual (only returned
%               if compiled)
%           'relGap' The relative duality gap.
%           'setupTime' The time to setup the problem in millisecond (only
%             returned if compiled).
%           'solveTime' The time to solve the problem in milliseconds (only
%             returned if compiled.
%
%This function does nothing but call the scs_matlab function from the
%splitting conic solver (SCS) library, except this function changes some of
%the default behaviours and has empty matrices returned instead of NaNs in
%the event of a failure.
%The library is described in [1] and the associated web page is
%https://github.com/cvxgrp/scs
%However, the Matlab implementation, which is used here is from
%https://github.com/bodono/scs-matlab
%
%EXAMPLE 1:
%This example shows how to use inequality constraints and how to impose
%multiple constraints on a single variable through the use of multiple
%variables.
% Here, we are to solve a linear programming problem with only inequality
% constraints. Such a problem has the form
% minimize a'*x
% subject to
% B*x<=d
% x>=0
% 
% The solution is
% numVar=length(a);
% c=a;
% A=[B;
%   -eye(numVar,numVar)];
% b=[d;zeros(numVar,1)];
% cone.l=2*numVar;
%
% %As an example, consider
% a=[1;-5];
% B=[2, 1;
%    -1, 2];
% d=[3;3];
% %Transforming into the values needed by the solver, we get
% c=a;
% A=[B;
%    -eye(2,2)];
% b=[d;zeros(2,1)];
% cone=[];
% cone.l=4;
% [x,y,s,info]=splittingConicSolver(A,b,c,cone)
%The solution is in the first two elements of x and is [0.6;1.8;0;0], where
%the final two variables can be discarded.
%
%EXAMPLE 2:
%This example shows how to use second order cone constraints. A second
%order cone (also known as a quadratic cone) optimization problem
%has the form
%minimize c'*x
%subject to
%norm(A1*x+b1)<=d1'*x+f1
%norm(A2*x+b2)<=d2'*x+f2
%...
%norm(An*x+bn)<=dn'*x+fn
%
%To be able to solve such a problem, use
% A=-[d1';
%    A1;
%    d2';
%    A2;
%    ...
%    dn';
%    An];
% c=c;
% b=[f1;b1;f2;b2...;fn;bn];
%and
% cone.q(1)=size(A1,1)+1;
% cone.q(2)=size(A2,1)+1;
% ...
% cone.q(n)=size(A2,1)+1;
%
%As an example, consider the problem where
% c=[-2;3;4];
% A1=[-13, 3,   5;
%     -12, 13, -6];
% A2=[-3,  6,  7;
%     2,   9,  1;
%     -1, -19, 3];
% b1=[-3;-2];
% b2=[0;4;-5];
% d1=[-13;-6;42];
% d2=[-2;6;-10];
% f1=-12;
% f2=27;
% %Transforming into the types of parameters needed for this solver, we
% %have
% c=c;
% A=-[d1';
%    A1;
%    d2';
%    A2];
% b=[f1;b1;f2;b2];
% cone=[];
% cone.q(1)=3;
% cone.q(2)=4;
% [x,y,s,info]=splittingConicSolver(A,b,c,cone)
%The answer for x is about x=[-0.9614; -1.2123;0.0227]
%
%EXAMPLE 3:
%In this example, we give an example of using semidefinite programming to
%deal with inequalities involving symmetric matrices. A matrix inequality
%refers to comparisons between the eigenvalues of the matrices.
%Here the optimization problem is
%minimize c'*x
%such that sum_{i=1}^N x(i)*B1i < D1
%          sum_{i=1}^N x(i)*B2i < D2 and
%          ....
%          sum_{i=1}^N x(i)*Bni < Dn
%where N=length(x) and the B and D  are all symmetric matrices. In this
%instance, we note that the jth constraint is equivalent to saying that
%-sum_{i=1}^N x(i)*Bji+S=Dj
%where S is a symmetric positive (semi)-definite matrix. Using this, we can
%see how to formulate the problem:
% c=-c;
%Note the required scaling by sqrt(2):
% A=[vech(B11,sqrt(2)),vech(B12,sqrt(2)),...,vech(B1N,sqrt(2));
%     vech(B21,sqrt(2)),vech(B22,sqrt(2)),...,vech(B2N,sqrt(2));
%     ...
%     vech(Bn1,sqrt(2)),vech(Bn2,sqrt(2)),...,vech(BnN,sqrt(2))];
% b=[vech(D1,sqrt(2));vech(D2,sqrt(2));...;vech(Dn,sqrt(2))];
% cone.s=[length(D1);length(D2);...length(Dn)];
%
%As an example, consider
% c=[1;-2;1];
% B11=[-8, -11;
%      -11, 3];
% B12=[7,  -42;
%      -42, 8];
% B13=[-2, -8;
%      -8,  1];
% B21=[-21, -12,  1;
%      -12,  10,  7;
%       1,   7,   5];
% B22=[0,   10,   16;
%      10, -16, -10;
%      16, -10,  3];
% B23=[-3,  2, -19;
%       2, -6,  8;
%      -19, 8,  7];
% D1=[22, -8;
%     -8,  37];
% D2=[14, 10,  40;
%     10,  99, 10;
%     40, 10, 15];
% %Now, transform the parameters into the form needed for the optimization
% %algorithm.
% c=c;
% A=[vech(B11,sqrt(2)),vech(B12,sqrt(2)),vech(B13,sqrt(2));
%     vech(B21,sqrt(2)),vech(B22,sqrt(2)),vech(B23,sqrt(2))];
% b=[vech(D1,sqrt(2));vech(D2,sqrt(2))];
% cone=[];
% cone.s=[2;3];
% [x,y,s,info]=splittingConicSolver(A,b,c,cone)
%One will see that the result is about x=[-0.34299;0.9879;-1.5830];
%
%REFERENCES:
%[1] B. O'Donoghue, E. Chu, N. Parikh, and S. Boyd, "Conic optimization via
%    operator splitting and homogeneous self-dual embedding," Journal of
%    Optimization Theory and Applications, 2016, to appear in print.
%    [Online]. Available: http://www.stanford.edu/%7Eboyd/papers/scs.html
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(params))
   params=[]; 
end

%Change the default behaviour to non-verbose.
if(~isfield(params,'verbose'))
   params.verbose=0.0;
elseif(params.verbose)
    %Make sure that it is a double as the compiled scs_indirect function
    %expects.
    params.verbose=1.0;
end

%This is for when the Matlab version is used.
params.extra_verbose=params.verbose;

%The function expects doubles.
if(isfield(params,'normalize'))
    params.normalize=double(params.normalize);
end

%Lower the default eps.
if(~isfield(params,'eps'))
   params.eps=1e-9; 
end

if(~isfield(params,'use_indirect'))
   params.use_indirect=true; 
end

data.A=sparse(A);
data.b=b;
data.c=c;

[x,y,s,info]=scs_matlab(data,cone,params);
isUnbounded=strcmp(info.status,'Unbounded');
info.pobj=c'*x;

if(any(~isfinite(x))||isUnbounded)
    x=[];
end

if(any(~isfinite(y))||isUnbounded)
    y=[];
end

if(any(~isfinite(s))||isUnbounded)
    s=[];
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
