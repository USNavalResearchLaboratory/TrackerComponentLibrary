function [a,b,c]=moments2OrthoPolyRecurCoeffs(mu,format)
%%MOMENTS2ORTHOPOLYRECURCOEFFS Given n values of the real integral
%          mu(k+1)=integral_a^b x^k w(x) dx for integer k=0 to n-1 (the
%          moments), determine a set of coefficients for a three-term 
%          recursion for polynomials that are orthonormal to w(x). If
%          p_n(x) is an nth order polynomial, then the set of orthonormal
%          polynomials to w(x) have the form:
%          integral_a^b w(x)p_m(x)p_n(x) dx =1 if m==n and =0 when m~=n.
%          A three-term recursion for these polynomials can be expressed as
%          b(i)*p_i(x)=(x-a(i))*p_{(i-1)}(x)-b(i-1)*p_{i-2}(x)
%          where we start from p_{-1}(x)=0. Alternatively, another form is 
%          p_i(x)=(a(i)*x+b(i))*p_{(i-1)}(x)-c(i)*p_{i-2}(x)
%          for for a length-n mu vector, where n=2*N+1, this function can
%          output a of length N and b of length N-1, for the firxt type of
%          recursion, or a and b of length N-1 and c of length N-2 for the
%          second type of recursion. This function only works if all mu>=0
%          and mu(1)>0. The three-element expansion can be used as in [1]
%          to generate cubature points (see the function
%          orthoPolyZerosFromRecur).
%
%INPUTS: mu A length n vector of moments of w(x) as described above. m(1)>0
%           and mu(2:n)>=0.
%    format A parameter selecting the format of the output. Possible values
%           are:
%           0 (The default if omitted or an empty matrix is passed) The
%             output is a and b for the first type of 3-element recursion.
%             c is set to an empty matrix.
%           1 The output is a, b, and c for the second type of recursion.
%
%OUTPUTS: a, b, and c The parameters for the 3-term recursion, the format
%                     of which is selected using the format output.
%
%This function implements the transformation that is given in Section 4 of
%[1].
%
%EXAMPLE:
%Here, we generate the moments for w(x)=exp(-x) for an integral from 0 to
%Inf using monomialIntExpAlt. We then obtain the three-term recursion from
%this function and pass the result to orthoPolyZerosFromRecur. The
%resulting cubature points and weights are compared to those obtained for
%the same integral using quadraturePoints1D, which uses a different
%technique.
% N=5;
% n=2*N+1;
% 
% mu=zeros(n,1);
% for i=1:n
%     mu(i)=monomialIntExpAlt(i-1);
% end
% [alpha,beta]=moments2OrthoPolyRecurCoeffs(mu,0);
% [xi1,w1]=orthoPolyZerosFromRecur(N,alpha,beta,[],mu(1))
% [xi,w]=quadraturePoints1D(N,4)
% %One will see that within finite precision limits, xi1 and xi are the
% %same and w and w1 are the same. Also note that the moments can be
% %verified. For example, the fifth-order moment agreement can be
% relErr=(mu(6)-sum(xi1.^5.*w'))/mu(6)
%where the relative error will be on the order of 1e-15.
%
%REFERENCES:
%[1] G. H. Golub and J. H. Welsh, "Calculation of Gauss quadrature rules,"
%    Mathematics of Computation, vol. 23, pp. 221-230, 1969.
%
%May 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(format))
    format=0;
end

if(mu(1)==0||any(mu<0))
    error('It is required that all moments are >=0 and mu(1)>0.') 
end

n=length(mu);

N=fix((n-1)/2);
L=N+1;

%The "Gram matrix".
M=hankel(mu(1:L),mu(L:n));

try
    R=chol(M,'upper');
catch
    error('It appears that finite precision errors preclude the decomposition of a matrix. Try using  a shorted mu vector.')
end

alphaVals=zeros(N,1);
betaVals=zeros(N-1,1);

%This implements Equation 4.3 of [1].
alphaVals(1)=R(1,2)/R(1,1);
betaVals(1)=R(2,2)/R(1,1);

for j=2:(N-1)
    alphaVals(j)=R(j,j+1)/R(j,j)-R(j-1,j)/R(j-1,j-1);
    betaVals(j)=R(j+1,j+1)/R(j,j);
end
j=N;
alphaVals(j)=R(j,j+1)/R(j,j)-R(j-1,j)/R(j-1,j-1);

if(format==0)
    a=alphaVals;
    b=betaVals;
    c=[];
else
    a=1./betaVals;
    b=-alphaVals(1:(N-1))./betaVals;
    
    c=zeros(N-1,1);
    c(1)=0;%This term is unused.
    c(2:(N-1))=betaVals(1:(N-2))./betaVals(2:(N-1));
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
