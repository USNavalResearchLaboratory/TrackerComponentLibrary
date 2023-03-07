function [A,x0,c]=polyEllipsForm2Quad(coeffs)
%%POLYELLIPSFORM2QUAD Given a set of coefficients defining an ellipsoid in
%       multivariate polynomial form such that a point is in or on the
%       ellipsoid if f(x1,x2,x3,...)<=0, obtain A, b, and c such that a
%       point x is inside of the ellipsoid if (x-x0)'*A*(x-x0)<=c.
%
%INPUTS: coeffs A 3X3X...X3 hypermatrix of the coefficients for the
%               multivariate quadratic polynomial. These are arranged such
%               that coeffs(a1,a2,a3...an) corresponds to the coefficient
%               of an x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1) term.
%               Thus, the number of indices coeffs takes is equal to the
%               dimensionality of x (not counting singleton dimensions at
%               the end of coeffs). This is the format used in
%               polyValMultiDim.
%
%OUTPUTS: A A numDimXnumDim symmetric matrix.
%        x0 A numDimX1 vector.
%         c A scalar constant.
%
%The solution was obtained by just multiplying out (x-x0)'*A*(x-x0)-c and
%seeing how terms relate. The inverse of this function is
%quadEllipsForm2Poly.
%
%EXAMPLE 1:
%This converts a set of A (symmetric), x0, and c into coefficients and then
%uses this function to retrieve the original values. The relavtive error
%between the original and retrieved values is on the order of finite
%precision limits.
% numDim=2;
% A=randn(numDim,numDim);
% A=(A+A')/2;%Make it is symmetric.
% x0=randn(numDim,1);
% c=randn(1);
% coeffs=quadEllipsForm2Poly(A,x0,c);
% [A1,x01,c1]=polyEllipsForm2Quad(coeffs);
% v1=[A(:);x0(:);c(:)];
% v2=[A1(:);x01(:);c1(:)];
% RelErr=max(abs((v2-v1)./v1))
%
%EXAMPLE 2:
%This is essentially the same as the first example, except, we scale the
%coefficents  by an arbitrary amount. That ends up producing A different A
%and c value. However, if we transforms the equations so that c=1 (side A
%by c and c by c), then we see that again the relative error between the
%solutions is within finite precision limits.
% numDim=2;
% A=randn(numDim,numDim);
% A=(A+A')/2;%Make it is symmetric.
% x0=randn(numDim,1);
% c=randn(1);
% coeffs=quadEllipsForm2Poly(A,x0,c);
% coeffs=coeffs*10;%Scale the coefficients by an arbitrary amount.
% [A1,x01,c1]=polyEllipsForm2Quad(coeffs);
% v1=[A(:)./c(:);x0(:)];
% v2=[A1(:)./c1(:);x01(:)];
% RelErr=max(abs((v2-v1)./v1))
%
%January 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

sizeVals=size(coeffs);

numDim=length(sizeVals);

if(numDim==2)
    if(sizeVals(2)==1)
        %Account for a trailing singleton dimension.
        numDim=1;
    end
end

switch(numDim)
    case 1
        c3=coeffs(3);
        c2=coeffs(2);
        c1=coeffs(1);

        A=c3;
        x0=-c2/(2*c3);
        c=A*x0^2-c1;
    case 2
        A=zeros(2,2);
        A(1,1)=coeffs(3,1);
        A(2,2)=coeffs(1,3);
        A(1,2)=coeffs(2,2)/2;
        A(2,1)=A(1,2);

        x0=(-2*A)\[coeffs(2,1);coeffs(1,2)];
        c=x0'*A*x0-coeffs(1,1);
    case 3
        A=zeros(3,3);
        A(1,1)=coeffs(3,1,1);
        A(2,2)=coeffs(1,3,1);
        A(3,3)=coeffs(1,1,3);
        A(1,2)=coeffs(2,2,1)/2;
        A(2,1)=A(1,2);
        A(1,3)=coeffs(2,1,2)/2;
        A(3,1)=A(1,3);
        A(2,3)=coeffs(1,2,2)/2;
        A(3,2)=A(2,3);

        x0=(-2*A)\[coeffs(2,1,1);coeffs(1,2,1);coeffs(1,1,2)];
        c=x0'*A*x0-coeffs(1,1,1);
    otherwise
        A=zeros(numDim,numDim);
        idx=repmat({1},[1,numDim]);

        %The squared terms.
        for curDim=1:numDim
            idx{curDim}=3;
            A(curDim,curDim)=coeffs(idx{:});
            idx{curDim}=1;
        end

        %For the cross terms
        for curDim1=1:(numDim-1)
            idx{curDim1}=2;
            for curDim2=(curDim1+1):numDim
                idx{curDim2}=2;
                A(curDim1,curDim2)=coeffs(idx{:})/2;
                A(curDim2,curDim1)=A(curDim1,curDim2);
                idx{curDim2}=1;
            end
            idx{curDim1}=1;
        end
        
        %Collect terms to solve for x0.
        b=zeros(numDim,1);
        for curDim=1:numDim
            idx{curDim}=2;
            b(curDim)=coeffs(idx{:});
            idx{curDim}=1;
        end
        x0=(-2*A)\b;
        c=x0'*A*x0-coeffs(1);
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
