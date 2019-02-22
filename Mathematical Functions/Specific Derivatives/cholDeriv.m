function dS=cholDeriv(SQ,dQ,Q,upperOrLower,algorithm)
%%CHOLDERIV Let the symmetric, positive definite matrix Q be a function of
%           a scalar variable x. dQ is the derivative of Q with respect to
%           x and SQ is the Cholesky decomposition of Q with respect to x.
%           That means that Q=SQ*SQ' if SQ is lower triangular and Q=SQ'*SQ
%           if SQ is upper triangular. This function finds the derivative
%           of SQ with respect to x.
%
%INPUTS: SQ The numDimXnumDim Cholesky decomposition (lower triangular or
%           upper triangular) of Q. One can get this with chol or
%           cholSemiDef.
%        dQ The numDimXnumDim derivative fo the matrix Q with respect to x.
%         Q The numDimXnumDim symmetric positive definite matrix.
% upperOrLower A string indicating whether SQ is an upper-triangular or
%           lower-triangular Cholesky decomposition matrix. If omitted or
%           an empty matrix is passed, the default value of 'upper' is
%           used. This can take the values 'upper' and 'lower'.
% algorithm A parameter indicating the algorithm to use. Possible values
%           are:
%           0 (The default if omitted or an empty matrix is passed) Use
%             Theorem A.1 in Appendix A.2 of [1].
%           1 Use algorithm A.2 of Appendix A.2 of [1].
%
%EXAMPLE:
%Here we have an example of a matrix that is a nonlinear function of x and
%we find the derivative of its Cholesky decomposition three ways, including
%using numerical differentiation. All of the results will agree within
%resonable precision bounds.
% QFun=@(x)((magic(4)+magic(4)')*x.^2+50*sqrt(x)*eye(4));
% dQFun=@(x)(2*(magic(4)+magic(4)')*x+(50/(2*sqrt(x)))*eye(4));
% x=1;
% Q=QFun(x);
% dQ=dQFun(x);
% SQ=chol(Q,'lower');
% dS0=cholDeriv(SQ,dQ,Q,'lower',0)
% dS1=cholDeriv(SQ,dQ,Q,'lower',1)
% epsVal=1e-6;
% Q1=QFun(x+epsVal);
% SQ1=chol(Q1,'lower');
% dSNumDiff=(SQ1-SQ)/epsVal
%
%REFERENCES:
%[1] S. Särkkä, Bayesian Filtering and Smoothing. Cambridge: Cambridge
%    University Press, 2013.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(algorithm))
    algorithm=0;
end

if(nargin<4||isempty(upperOrLower))
    upperOrLower='upper';
end

numDim=size(SQ,1);

switch(upperOrLower)
    case 'lower'
        isTransposed=false;
    case 'upper'
        SQ=SQ';
        dQ=dQ';
        Q=Q';
        isTransposed=true;
    otherwise
        error('An invalid value was entered for the upperOrLower parameter.')
end

if(algorithm==0)
    SQInv=pinv(SQ);
    M=SQInv*dQ*SQInv';
    %Appy the Phi Function of A.8.
    for i=1:numDim
        M(i,i)=M(i,i)/2;
        idx=(i+1):numDim;
        M(i,idx)=0;
    end
    dS=SQ*M;
elseif(algorithm==1)
    dS=zeros(numDim,numDim);
    for i=1:numDim
        idx=1:(i-1);
        SQ(i,i)=sqrt(Q(i,i)-sum(SQ(i,idx).^2));
        dS(i,i)=(dQ(i,i)-2*sum(dS(i,idx).*SQ(i,idx)))/(SQ(i,i)*2);
        for j=(i+1):numDim
            SQ(j,i)=(Q(j,i)-sum(SQ(j,idx).*SQ(i,idx)))/SQ(i,i);
            tempVal=dQ(j,i)-sum(dS(j,idx).*SQ(i,idx)+SQ(j,idx).*dS(i,idx));
            dS(j,i)=tempVal/SQ(i,i)-(dS(i,i)/SQ(i,i))*SQ(j,i);
        end
    end
else
    error('Unknown algorithm specified.')
end

if(isTransposed)
    dS=dS';
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
