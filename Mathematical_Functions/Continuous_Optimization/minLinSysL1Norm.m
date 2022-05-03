function [x,info]=minLinSysL1Norm(A,b,algorithm,params)
%%MINLINSYSL1NORM This function solves the optimization problem
%                 minimize norm(A*x-b,1) over x. The result differs from
%                 solving such a system using pinv(A)*b, because the l1,
%                 not the l2 norm is used. Complex valued arguments are
%                 allowed. This algorithm does not work well if A does not
%                 have full column rank.
%
%INPUTS: A An mXn matrix, which can be real or complex. A can be sparse.
%        b An mX1 vector, which can be real or complex.
%algorithm An optional parameter specifying which algorithm to use to solve
%          the problem. The problem is solved as a second order cone
%          problem. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            function secondOrderConeInteriorPoint.
%          1 Use the function splittingConicSolver.
%   params Parameters that affect how the second order cone programming
%          algorithm works. See the comments to either
%          secondOrderConeInteriorPoint for algorithm 0 or
%          splittingConicSolver for algorithm 1 for more details. The
%          default values are the same as for those solvers, except the
%          maximum number of iterations for algorithm 1 is increased to
%          30e3 due to the slow convergence of that algorithm with moderate
%          sized matrices.
%
%OUTPUTS: x The nX1 solution to the problem or an empty matrix is no
%           solution is possible.
%      info The problem is formulated as a second-order cone optimization
%           problem. This is the information regarding the termination
%           state. This is either the return value exitCode from
%           secondOrderConeInteriorPoint or info from splittingConicSolver,
%           depending on the selected algorithm.
%
%As described in [1], the problem can be formulated as a second order cone
%optimization problem. Thus, available second order cone programming
%algorithms are used to solve the problem.
%
%EXAMPLE 1:
% A=[magic(4);
%    1,2,3,50];
% b=[160;-12;64;-123;1];
% [x,info]=minLinSysL1Norm(A,b)
% %The solution of about x=[12.3;13.9;-24.5;0.7] is found.
% norm(A*x-b,1)
% x1=pinv(A)*b;
% norm(A*x1-b,1)
% %And one sees that the l1 error from this estimator is less than the l1
% %error from using a pseudoinverse (18.3333 versus 22).
%
%EXAMPLE 2:
%Larger dense systems can be a bit slow. For example,
% A=rand(50,33);b=rand(50,1);
% [x,info]=minLinSysL1Norm(A,b)
% norm(A*x-b,1)
% x1=pinv(A)*b;
% norm(A*x1-b,1)
%Will often take 2000 to 25000 iterations.
%
%REFERENCES:
%[1] S. Lobo, Miguel, L. Vandenberg, S. Boyd, and H. Lebret, "Applications
%    of second-order cone programming," Linear Algebra and its
%    Applications, vol. 284, no. 1, pp. 193-228, 1998.
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
   algorithm=0; 
end

if(nargin<4||isempty(params))
   params=[];
   %This default only matters if algorithm=1.
   params.max_iters=30e3;
end

p=size(A,1);
q=size(A,2);

AReal=real(A);
AImag=imag(A);
BReal=real(b);
BImag=imag(b);

%p auxiliary variables are used.
c=[zeros(2*q,1);ones(p,1)];

%Allocate space
F=zeros(3*p,2*q+p);
g=zeros(3*p,1);%Allocate space

e=zeros(2*q+p,1);

curStart=1;
for i=1:p
    span=curStart:(curStart+2);
    e(2*q+i)=1;
    ACur=[AReal(i,:), -AImag(i,:);
          AImag(i,:), AReal(i,:)];
    F(span,:)=-[e';
                ACur,zeros(2,p)];
    g(span)=[0;-BReal(i);-BImag(i)];
    e(2*q+i)=0;
    
    curStart=curStart+3;
end

cone=[];
cone.q=3*ones(p,1);

switch(algorithm)
    case 0
        [~,x,~,info]=secondOrderConeInteriorPoint(F',-c,g,cone.q,params);
    case 1
        [x,~,~,info]=splittingConicSolver(F,g,c,cone,params);
    otherwise
        error('Unknown algorithm specified')
end

%Extract the real and imaginary parts and get rid of the auxiliary
%variables.
if(~isempty(x))
    x=x(1:q)+1i*x((q+1):(2*q));
end

%If the original problem was all real, then 
if(all(isreal(A(:)))&&all(isreal(b)))
    x=real(x);
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
