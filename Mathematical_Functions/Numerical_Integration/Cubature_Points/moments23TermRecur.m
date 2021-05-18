function [alpha,beta]=moments23TermRecur(mu)
%%MOMENTS23TERMRECUR Given 2*N+1 moments of a scalar weighting function
%          from 0 to 2*N, the kth moment is
%          mu(k+1)=integral_a^b w(x)*x^k dx where w(x) is a weighting
%          function and a and b need not be finite, obtain the three term
%          recusion for orthogonal polynomials that share the same moments.
%          The three term recursion coefficients cna be used with the
%          orthoPolyZerosFromRecur to obtain quadrature points and weights
%          when using the same weighting function. The three-term recursion
%          has the form:
%          x*p_{k-1}(x)=beta_{k-1}*p_{k-2}(x)+alpha_k*p_{k-1}+beta_kp_k(x)
%          for k=1,...,N starting with p_{-1}(x)=0. This can be useful for
%          obtaining cubature points for integration over continuous
%          probability distributions with easy-to-express noncentral
%          moments (e.g. from the moment generting function of the
%          distribution.
%
%INPUTS: mu A length (2*N+1) array holding the moments from orders 0 to 2*N
%           of the desired weighting function (which could be, but needn't
%           be, a continuous probability distribution).
%
%OUTPUTS: alpha A length N vector of coefficients.
%          beta A length N-1 vector of coefficients.S
%
%The algorithm is taken from Section 4 of [1]. To use these with the
%orthoPolyZerosFromRecur function, one must pass alpha, beta, and mu(1).
%See the example below.
%
%EXAMPLE:
%In this example, we get the values for a three-term recursion of the
%standard normal distribution. We then plus the values into
%orthoPolyZerosFromRecur and get cubature points for integration voer the
%distribution. We then show that these points are identical tot he ones
%obtained from the quadraturePoints1D function (differences are within
%finite precision limits).
% meanVal=0;
% variance=1;
% N=5;
% mu=zeros(2*N+1,1);
% for k=0:2*N
%     mu(k+1)=GaussianD.momentGenFun(meanVal,variance,k);
% end
% [alpha,beta]=moments23TermRecur(mu);
% [xi,w]=orthoPolyZerosFromRecur(N,alpha,beta,[],mu(1));
% [xiGauss,wGauss]=quadraturePoints1D(N,0);
% max(abs((xi-xiGauss)))
% max(abs(w-wGauss))
%
%REFERENCES:
%[1] G. H. Golub and J. H. Welsh, "Calculation of Gauss quadrature rules,"
%    Mathematics of Computation, vol. 23, pp. 221-230, 1969.
%
%January 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=(length(mu)-1)/2;

if(fix(N)~=N)
    error('This function requires an odd number of moments.')
end

%The Gram Matrix in Section 4.
M=zeros(N+1,N+1);
for i1=1:(N+1)
    for i2=1:(N+1)
        M(i1,i2)=mu(i1+i2-1);
    end
end

%Using cholSemiDef instead of chol, it is more robust to semidefinite
%matrices.
R=cholSemiDef(M,'upper');

%Equation 4.3 in [1] for both alpha and beta.
alpha=zeros(N,1);
alpha(1)=R(1,2)/R(1,1);
for i=2:N
    alpha(i)=R(i,i+1)/R(i,i)-R(i-1,i)/R(i-1,i-1);
end

beta=zeros(N-1,1);
for i=1:(N-1)
    beta(i)=R(i+1,i+1)/R(i,i);
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
