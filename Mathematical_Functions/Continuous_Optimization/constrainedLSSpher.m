function [x,lambda,exitFlag]=constrainedLSSpher(A,b,alpha,varargin)
%%CONSTRAINEDLSSPHER Find x to minimize norm(A*x-b,2) under the constraint
%                    that norm(x,2)<=alpha. This is essentialy contraining
%                    x to a sphere of radius alpha. This only solves real
%                    systems.
%
%INPUTS: A A real mXn matrix with m>=n.
%        b A real mX1 vector.
%    alpha The inequality constraint value. If this parameter is omitted or
%          an empty matrix is passed, the default of 1 is used.
% varargin Any parameters that one might wish to pass to the fzero
%          function. These are generally comma-separated things, such as
%          'TolX',1e-12.
%
%OUTPUTS: x The optimal value of x subject to the spherical constraint.
%    lambda The Lagrangian multiplier used in the optimization. lambda>=0.
%           If lambda=0, then the constraint on x did not have to be
%           enforced.
%  exitFlag The exit flag from the fzero function, if used. Otherwise, this
%           is 0.
%
%This is an implementation of Algorithm 6.2.1 in Chapter 6.2.1 of [1]. The
%zero of the secular equation is found using the fzero function after
%bounding the equation.
%
%If after solving the unconstrained system norm(x)>alpha, then one must
%find the zero of a polynomial equation given in Algorithm 6.2.1 of [1].
%The equation is a sum of squared terms. Since this is an inequality
%constrained problem, we know that the Lagrangian multiplier must be >=0.
%The lower bound is thus 0. For an upper bound, we need to find a value of
%lambda that changes the sign of the objective function. To do this, we
%take the largest term in the sum and multiply it by r (r is the number of
%terms in the sum). Then, we find when the lambda for that to solve the
%optimization problem, multiply the result by 2 and we have an upper bound.
%Of course, we do not know which is the largest term, so we evaluate all of
%them and choose the largest.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4)
    varargin{1}=[];
end

if(nargin<3||isempty(alpha))
    alpha=1;
end


r=rank(A);

[U,Sigma,V]=svd(A,0);
sigma=diag(Sigma);

%Sums are only up to r, so get rid of the extra elements.
U=U(:,1:r);
V=V(:,1:r);
sigma=sigma(1:r);

bTilde=U'*b;
x=V*(bTilde./sigma);

if(norm(x)>alpha)
    %The abs is because there are two solutions and the correct one will be
    %the positive one.
    upperBound=2*max((abs(bTilde.*sigma)*sqrt(r)-alpha*sigma.^2)./alpha);
    f=@(lambda)objFun(lambda,alpha,bTilde,sigma);
    [lambda,~,exitFlag]=fzero(f,[0,upperBound],varargin{:});
    x=V*((sigma.*bTilde)./(sigma.^2+lambda));
else
    lambda=0;
    exitFlag=0;
end
end

function val=objFun(lambda,alpha,bTilde,sigma)
%This is the function in Algorithm 6.2.1 that must be zeroed to find the
%Lagrange multiplier.

    val=sum(((sigma.*bTilde)./(sigma.^2+lambda)).^2)-alpha^2;
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
