function vals=Laguerre(x,n,alpha,findAll)
%%%LAGUERRE Evaluate a generalized Laguerre polynomial (also known as an
%           associated Laguerre polynomial) at one or more points. Laguerre
%           polynomials arise when computing the mean of the Chi
%           distribution. They also arise in various problems in quantum
%           mechanics.
%
%INPUTS: x The 1XN or NX1 set of points at which the Laguerre polynomial
%          should be evaluted.
%        n The integer value of the order of the polynomial. If a negative
%          parameter is passed, then an empty matrix will be returned.
%          Results for very large parameters can be slow as a recurrence
%          relation is used.
%    alpha The real scalar prameter in a generalized Laguerre polynomial.
%          In a standard Laguerre polynomial, this is zero. If this
%          parameter is omitted or an empty matrix is passed, then alpha=0
%          is used.
%  findAll An optional parameter indicating whether the values of all
%          Laguerre polynomials from order zero to n should be returned or
%          whether just the value of the polynomial of order n should be
%          returned. The default if this parameter is omitted or an empty
%          matrix is passed is false, which means that only the value of
%          the polynomial of order n is returned.
%
%OUTPUTS: vals If findAll is false, then this is a NX1 set of the values of
%              the Laguerre polynomial of order n evaluated at the ponts in
%              x. Otherwise, this is a NX(n+1) set of all Laguerre
%              polynomials of order 0 to n evaluated at x.
%
%This function implements the recursion given in  [1].
%
%REFERENCES:
%[1] Abramowitz, M. and Stegun, I. A. (Eds.). "Orthogonal Polynomials."
%    Table 22.7 in Ch. 22 in Handbook of Mathematical Functions with
%    Formulas, Graphs, and Mathematical Tables, 9th printing. New York:
%    Dover, pp. 782, 1972.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(findAll))
   findAll=false; 
end

if(nargin<3||isempty(alpha))
   alpha=0; 
end

if(n<0)
   vals=[]; 
   return;
end

x=x(:);
numPoints=length(x);

%Deal with the first two special cases.
if(n==0)
   vals=ones(numPoints,1);
   return;
end

if(findAll)
    vals=zeros(numPoints,n+1);
    vals(:,1)=1;
    vals(:,2)=1+alpha-x;
    for k=2:n
        %The two-value recursion in Ch. 22.7 of [1].
        vals(:,k+1)=(1/k)*(2*(k-1)+alpha+1-x).*vals(:,k-1+1)-((k-1+alpha)/k)*vals(:,k-2+1);
    end
else%If only the nth order polynomial value is supposed to be returned.
    if(n==1)
        vals=1+alpha-x;
        return;
    end

    %The 2 holds the last two values.
    vals=zeros(numPoints,2);

    %For n=0
    vals(:,2)=1;

    %For n=1
    vals(:,1)=1+alpha-x;

    for k=2:n
        temp=vals(:,1);
        %The two-value recursion in Ch. 22.7 of [1].
        vals(:,1)=(1/k)*(2*(k-1)+alpha+1-x).*vals(:,1)-((k-1+alpha)/k)*vals(:,2);
        vals(:,2)=temp;
    end

    %Just return the last set of values.
    vals=vals(:,1);
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
