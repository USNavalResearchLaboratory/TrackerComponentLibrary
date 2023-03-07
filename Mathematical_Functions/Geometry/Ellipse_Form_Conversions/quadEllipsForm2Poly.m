function coeffs=quadEllipsForm2Poly(A,x0,c)
%%QUADELLIPSFORM2POLY Given an ellipse (in 2D) or an ellipsoid (in >2D)
%          given by the equation such that a point x is in or on the
%          ellipsoid if (x-x0)'*A*(x-x0)<=c, convert it into a set of
%          coefficients for a quadratic multivariate polynomial such
%          that if a point is on or in the ellipsoid, then
%          f(x1,x2,x3,...)<=0, where f is the polynomial and the x's are
%          all of the variables. These coefficients are returned in an
%          ordering such that one can use polyValMultiDim to evaluate the
%          polynomial.
%
%INPUTS: A A numDimXnumDim real matrix that specifies the shape and
%          rotation of the ellipsoid. This doesn't actually have to be
%          symmetric or positive definite. numDim>=1.
%       x0 The numDimX1 real center of the ellipsoid.
%        c The real, scalar right-hand side of the ellipsoid equation,
%          which affects the scale of the ellipsoid.
%
%OUTPUTS: coeffs A 3X3X...X3 hypermatrix of the coefficients for the
%               multivariate quadratic polynomial. These are arranged such
%               that coeffs(a1,a2,a3...an) corresponds to the coefficient
%               of an x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1) term.
%               Thus, the number of indices coeffs takes is equal to the
%               dimensionality of x (not counting singleton dimensions at
%               the end of coeffs).
%
%The solution was obtained by just multiplying out (x-x0)'*A*(x-x0)-c and
%grouping terms. The inverse of this function is polyEllipsForm2Quad, which
%returns a symmetric A.
%
%EXAMPLE:
%In this example, we generate random values for A, x0, and c and change the
%ellipse format. For a random test point xTest, we then verify that
%(xTest-x0)'*A*(xTest-x0)-c produces the same value as using coeffs and
%xTest in polyValMultiDim. The relative error is shown and it is on the
%order of finite precision limitations.
% numDim=3;
% A=randn(numDim,numDim);
% x0=randn(numDim,1);
% c=randn(1);
% coeffs=quadEllipsForm2Poly(A,x0,c);
% xTest=randn(numDim,1);
% val1=(xTest-x0)'*A*(xTest-x0)-c;
% val2=polyValMultiDim(coeffs,xTest);
% RelDiff=(val1-val2)/val1
%
%January 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Make A symmetric. This doesn't change the value (x-x0)'*A*(x-x0). 
A=(A+A')/2;

numDim=size(A,1);

if(numDim==1)
    coeffs=zeros(3,1);
    coeffs(1)=x0^2*A-c;
    coeffs(2)=-2*x0*A;
    coeffs(3)=A;

elseif(numDim==2)
    coeffs=zeros(3,3);

    %The constant term.
    coeffs(1,1)=x0'*A*x0-c;

    %The linear terms.
    linTerms=-2*A*x0;%This works if A=A' --which we ensured above.
    coeffs(2,1)=linTerms(1);
    coeffs(1,2)=linTerms(2);

    %The squared terms.
    coeffs(3,1)=A(1,1);
    coeffs(1,3)=A(2,2);

    %The cross term.
    coeffs(2,2)=2*A(1,2);
elseif(numDim==3)
    coeffs=zeros(3,3,3);

    %The constant term.
    coeffs(1,1,1)=x0'*A*x0-c;

    %The linear terms.
    linTerms=-2*A*x0;
    coeffs(2,1,1)=linTerms(1);
    coeffs(1,2,1)=linTerms(2);
    coeffs(1,1,2)=linTerms(3);

    %The squared terms.
    coeffs(3,1,1)=A(1,1);
    coeffs(1,3,1)=A(2,2);
    coeffs(1,1,3)=A(3,3);

    %The cross terms.
    coeffs(2,2,1)=2*A(1,2);
    coeffs(2,1,2)=2*A(1,3);
    coeffs(1,2,2)=2*A(2,3);
else
    dims=3*ones(1,numDim);
    coeffs=zeros(dims);

    %The constant term.
    coeffs(1)=x0'*A*x0-c;

    idx=repmat({1},[1,numDim]);

    %The linear terms.
    linTerms=-2*A*x0;
    for curDim=1:numDim
        idx{curDim}=2;
        coeffs(idx{:})=linTerms(curDim);
        idx{curDim}=1;
    end

    %The squared terms.
    for curDim=1:numDim
        idx{curDim}=3;
        coeffs(idx{:})=A(curDim,curDim);
        idx{curDim}=1;
    end

    %The linear terms
    for curDim1=1:(numDim-1)
        idx{curDim1}=2;
        for curDim2=(curDim1+1):numDim
            idx{curDim2}=2;
            coeffs(idx{:})=2*A(curDim1,curDim2);
            idx{curDim2}=1;
        end
        idx{curDim1}=1;
    end
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
