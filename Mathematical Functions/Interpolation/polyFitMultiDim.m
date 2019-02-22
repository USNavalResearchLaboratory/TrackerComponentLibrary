function [coeffs,termMat]=polyFitMultiDim(varVals,funVals,degree)
%%POLYFITMULTIDIM Obtain a least-squares multivariate polynomial fit of a
%                 specified order to a given set of data.
%
%INPUTS: varVals This is a numDimXnumPoints set oif dtapoints where the
%                function was evaluated.
%        funVals This is a numPointsX1 or 1XnumPoints vector of the values
%                of the functions at the specified data points.
%        degree The degree of the desired interpolating polynomial. The
%               degree must be less than the total number of monomial
%               terms, which is
%               (1+degree)*binomial(degree+numDim,numDim-1)/numDim.
%
%OUTPUTS: coeffs A hypermatrix of the coefficients for the multivariate
%               polynomial that can be used in the polyValMultiDim
%               function. These are arranged such that
%               coeffs(a1,a2,a3...an) corresponds to the coefficient of an
%               x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1) term.  Thus, the
%               number of indices coeffs takes is equal to the
%               dimensionality of x (not counting singleton dimensions at
%               the end of coeffs). Note that this ordering is the reverse
%               that used in the 1D polyval function that is built into
%               Matlab. The number of elements for each index in coeffs is
%               the maximum order of that dimension +1.
%       termMat An (n+1)XnumTerms matrix such that
%               termMat(:,i)=[c,a1,a2,...,an] is a monomial term where c is
%               the value of of the monomial coefficient and the monomial
%               is x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1).
%
%Multivariate polynomials consist of a sum of monomial terms of the form
%c*x1^(a1)*x2^(a2)*x3^(a3)...xn^(an)
%Here, the variables x are known as are the values z of the polynomial
%evaluated at those points. Thus, the unknowns are the c terms, the
%coefficients of the monomials. This function evaluates the monomial terms
%and build a matrix with all of the monomial term values for all of the
%equations. The solution to the coefficient values is thus the solution to
%a linear system of equations. The 
%
%EXAMPLE:
%Here, we have a two-dimensional function evaluated at a number of random
%points and we want to fit a polynomial to it. In this problem, we do not
%add noise to the function.
% points=linspace(-1,1,10);
% [x,y]=meshgrid(points,points);
% z=2*x.^3+2*x.*y-3*y.^3+4*x.^2.*y-5*y.^2.*x.^3;
% varVals=[x(:)';y(:)'];
% funVals=z(:);
% degree=5;
% [~,termMat]=polyFitMultiDim(varVals,funVals,degree);
% %One will get coefficients that are almost exact. Given that we know that
% %the coefficients are integers, we can eliminate all small coefficients by
% %rounding.
% termMat=round(termMat);
% %Get rid of terms that are numerically zero.
% sel=termMat(1,:)~=0;
% termMat=termMat(:,sel);
% terms2String(termMat,{'x','y'})
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(varVals,1);
numEq=size(varVals,2);

%The total number of monomial terms for the multivariate polynomial. This
%is the sum of binomial(degCur+numDim-1,numDim-1) from degCur=0 to degree,
%which is the number of compositions of the degree degCur into n parts.
totalNumMonomials=(1+degree)*binomial(degree+numDim,numDim-1)/numDim;

if(totalNumMonomials>numEq)
   error('The total number of monomial coefficients to find is greater than the number of data points provided.') 
end

termMat=zeros(numDim+1,totalNumMonomials);

V=zeros(numEq,totalNumMonomials);

%x^0,y^0, etc. value.
V(:,1)=1;
curTerm=2;
for curDegree=1:degree
    numMonomials=binomial(curDegree+numDim-1,numDim-1);
    
    for j=0:(numMonomials-1)
        curMonomial=unrankTComposition(j,numDim,curDegree+numDim,true)-1;
        termMat(2:end,curTerm)=curMonomial;
        for curEq=1:numEq
            V(curEq,curTerm)=prod(varVals(:,curEq).^curMonomial);
        end
        curTerm=curTerm+1;
    end
end

termMat(1,:)=V\funVals(:);
%Get rid of terms that are numerically zero.
sel=termMat(1,:)~=0;
termMat=termMat(:,sel);
coeffs=terms2MultiDimPolyMat(termMat);

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
