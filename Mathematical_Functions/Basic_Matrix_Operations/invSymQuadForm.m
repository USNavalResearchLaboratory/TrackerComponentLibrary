function dist=invSymQuadForm(x,M,matType)
%%INVSYMQUADFORM Compute the quadratic form x'*inv(R)*x in a manner that
%                should be more numerically stable than directly evaluating
%                the matrix inverse, where R is a symmetric, positive
%                definite matrix. The parameter M can either be the matrix
%                R directly (R=M) (the default), or the lower-triangular
%                square root of R (R=M*M'). The distance  can be computed
%                for a single matrix and multiple vectors x. If x consists
%                of vector differences and R is a covariance matrix, then
%                this can be viewed as a numerically robust method of
%                computing Mahalanobis distances.
%
%INPUTS: x An mXN matrix of N vectors whose quadratic forms with the
%          inverse of R are desired.
%        M The mXm real, symmetric, positive definite matrix R, or the
%          square root of R, as specified by the following parameter.
%  matType This optional parameter specified the type of matrix that M is.
%          Possible values are
%          0 (The default if omitted or an empty matrix is passed) M is the
%            matrix R in the form x'*inv(R)*x.
%          1 M is the lower-triangular square root of the matrix R.
%
%OUTPUTS: dist A 1XN vector of the values of x*inv(R)*x for every vector in
%              the matrix x.
%
%As one can find in many textbooks, solving A*x=b using matrix inversion is
%generally a bad idea. This relates to the problem at hand, because one can
%write
%x'*inv(R)*x=x'inv(C*C')*x
%           =x'*inv(C)'*inv(C)*x
%where C is the lower-triangular Cholesky decomposition of R. Next, say
%that
%y=inv(C)*x.
%Then, we are just computing y'*y.
%What is a stable way to find y? Well, we can rewrite the equation as
%C*y=x
%Since C is lower-triangular, we can find x using forward substitution.
%This should be the same as one of the many ways that Matlab solves the
%equation C*y=x when one uses the \ operator. One can explicitly tell
%Matlab that the matrix is lower triangular when using the linsolve
%function, thus avoiding the need for loops or for Matlab to check the
%structure of the matrix on its own.
%
%August 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<3||isempty(matType))
        matType=0;
    end
    
    switch(matType)
        case 0
            C=cholSemiDef(M,'lower');
        case 1
            C=M;
        otherwise
            error('Invalid matrix type specified');
    end

    opts.LT=true;
    y=linsolve(C,x,opts);
    %y=C\x;
    dist=sum(y.*y,1);
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
