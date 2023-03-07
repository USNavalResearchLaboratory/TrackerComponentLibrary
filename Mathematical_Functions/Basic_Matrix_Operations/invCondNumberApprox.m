function rCond=invCondNumberApprox(X)
%%INVCONDNUMBERAPPROX Find the approximate inverse condition number of a
%       real, square matrix using the l1 norm. The output is very similar
%       to the rcond function built into Matlab, when given a real matrix.
%       It approximates the inverse condition number without directly
%       inverting the matrix X.
%
%INPUTS: X an nXn real matrix.
%
%OUTPUTS: rCond The approximate inverse condition number (a real scalar).
%               If any non-finite values are in X, then a NaN is returned.
%
%The inverse condition number is 1/(norm(X,1)*norm(inv(X),1)). This
%function uses algorithm 4.1 from [1] to approximate norm(inv(X),1) without
%performing the matrix inversion.
%
%REFERENCES:
%[1] J. Higham, Nicholas, "FORTRAN codes for estimating the one-norm
%    of a real or complex matrix with applications to condition
%    estimation," ACM Transactions on Mathematical Software, vol. 14, no.
%    4, pp. 381-396, Dec. 1988.
%
%June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(size(X,1)~=size(X,2))
    error('X must be a square matrix.')
end

if(isempty(X))
    rCond=Inf;
    return
end

if(any(~isfinite(X(:))))
    rCond=NaN;
    return;
end

%The maximum absolute column sum of the matrix.
normX=norm(X,1);
if(normX==0)
    rCond=0;
    return;
end

rCond=1/(normX*approxInvMat1Norm(X));

end

function [gammaVal,v]=approxInvMat1Norm(A)
%%APPROXINVMAT1NORM This implements the approximate matrix norm function of
%               in algorithm 4.1 of [1].
%
%REFERENCES:
%[1] J. Higham, Nicholas, "FORTRAN codes for estimating the one-norm
%    of a real or complex matrix with applications to condition
%    estimation," ACM Transactions on Mathematical Software, vol. 14, no.
%    4, pp. 381-396, Dec. 1988.
%
%June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[L,U,P]=lu(A);
%L must be lower triangular and U upper triangular The permutation in P can
%be ignored, because it does not change the result.

n=size(A,1);
v=backSubstitution(U,forwardSubstitution(L,ones(n,1)/n));
if(n==1)
    gammaVal=abs(v);
    return;
end
gammaVal=norm(v,1);
zeta=signP(v);
x=backSubstitution(L',forwardSubstitution(U',zeta));
I=eye(n,n);
k=2;
[~,j]=max(abs(x));
while(1)
    v=backSubstitution(U,forwardSubstitution(L,I(:,j)));
    gammaBar=gammaVal;
    gammaVal=norm(v,1);
    zetaNew=signP(v);
    if(all(zetaNew==zeta)||gammaVal<gammaBar)
        break;
    end
    zeta=zetaNew;
    x=backSubstitution(L',forwardSubstitution(U',zeta));
    k=k+1;
    
    [~,jNew]=max(abs(x));
    if(jNew==j||k>5)
        break;
    end
    j=jNew;
end

signVal=1;
for i=1:n
    x(i)=signVal*(1+(i-1)/(n-1)); 
    signVal=-signVal;
end
x=backSubstitution(U,forwardSubstitution(L,x));

testVal=2*norm(x,1)/(3*n);
if(testVal>gammaVal)
   v=x;
   gammaVal=testVal;
end

if(isnan(gammaVal))
    gammaVal=Inf;
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
