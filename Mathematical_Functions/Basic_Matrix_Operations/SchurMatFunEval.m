function fX=SchurMatFunEval(f,X)
%%SCHURMATFUNEVAL Evaluate the matrix function f(X), where the function f
%         is a function that has a Taylor series representation using the
%         Schur-Parlett method. In this approach, f(x) only needs to be
%         evaluated on _scalar_ quantities. Thus, this function can be used
%         to turn scalar functions into matrix functions. Matrix functions
%         with a Taylor series expansions satisfy the identities
%         X*f(X)=f(X)*X and f(inv(B)*X*B)=B*f(X)*inv(B) and include
%         functions such as the matrix exponential and square root.
%
%INPUTS: f A function handle that takes a scalar input and returns a scalar
%          output.
%        X The real or complex nXn matrix that should be evaluated via the
%          function f.
%
%OUTPUTS: fX The value f(X), where the function f applied to scalar terms
%            is applied to the matrix X. This is not the same as applying f
%            to every scalar element in X.
%
%This function implements Algorithm 9.1.1 of Chapter 9.1.4 of [1].
%
%EXAMPLE:
%Here, we show that the result agrees with the mpow function to take a
%matrix square root.
% X=randn(6,6);
% f=@(x)sqrt(x);
% fX=SchurMatFunEval(f,X);
% diffVal=fX-mpower(X,1/2)
%One will see that the entries in diffVal are all on the order of 1e-14 or
%less, indicating that the two are probably equal within finite precision
%bounds.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(X);

%We must use the 'complex' option or else T is not guaranteed to be
%upper triangular.
[Q,T]=schur(X,'complex');

%The diagonal terms.
F=zeros(n,n);
for i=1:n
    F(i,i)=f(T(i,i));
end

for p=1:(n-1)
    for i=1:(n-p)
        j=i+p;
        s=T(i,j)*(F(j,j)-F(i,i));
        for k=(i+1):(j-1)
            s=s+T(i,k)*F(k,j)-F(i,k)*T(k,j);
        end
        F(i,j)=s/(T(j,j)-T(i,i));
    end
end

if(isreal(X))
    fX=real(Q*F*Q');
else
    fX=Q*F*Q';
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
