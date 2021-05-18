function [a,c]=HermiteInterpPoly(x,y)
%%HERMITEINTERPPOL Given a scalar function y(x) evaluated at certain
%                  values of x (control points) along with an arbitrary
%                  number of derivatives of y(x) at those points, this
%                  function returns the coefficients of the Hermite
%                  interpolating polynomial. If a large number of control
%                  points is used, then numerical precision problems can
%                  lead to bad results. In general, exceeding about 15
%                  control points can lead to a loss of precision.
%
%INPUTS: x An NpX1 or 1XNp vector of values at which the function y and its
%          derivative are given. The values are assumed to be given in
%          increasing order.
%        y A numMatchXNp matrix of values of the function y(x) and its
%          derivatives evaluated at values of the points in x. y(row,col)
%          corresponds to the row-1 derivative of y evaluated at x(col).
%
%OUTPUTS: a Coefficients of the interpolating polynomial of y in Newton
%           form as a ROW vector (1XnumCoeff). The interpolating polynomial
%           can be evaluated at a desired point using the function
%           polyValNewton with a and the c point vector. By stacking
%           multiple row vectors, one can use the function polyValNewton to
%           evaluate multiple polynomials at the same time (i.e for vector
%           valued functions).
%         c Control points for the interpolating polynomial of y in Newton
%           form as a 1X(numCoeff-1) row vector.
%
%Note that the function call
%aPower=convertPolynomialForm('Newton','PowerSeries',a,c);
%can convert a and c into power series (for use with the polyval function).
%However, a loss of precision will occur and some interpolation problems
%that would work very well in Newton's form might produce very bad results.
%
%The interpolating polynomial matches all values of y and its derivatives
%at the control points in x. Newton's form of a polynomial evaluated at the
%point z is
%y(z)=a(1)+sum_{k=1}^{n-1}a(k+1)(z-c(1))*(z-c(2))*...*(z-c(k))
%
%The algorithm is an implementation of the general divided difference in
%[1]. In Section 2.1.3 of [1], it is mentioned that the error in the Newton
%polynomial solution can be reduced by ordering the coefficients in a
%certain manner. As long as the number of points to match is kept small,
%that is not necessary.
%
%If that is difficult to understand, one might want to start with Chapter
%3.4 of [2], which describes the basic algorithm only for matching the
%position and the first derivative.
%
%Often it is best to perform interpolation just between a few points at a
%time to avoid  wild fluctuations of the fitted polynomial that can result
%when many points are used and finite precision errors will increase.
%
%REFERENCES:
%[1] Hermite interpolation algorithm described in Chapter 2.1.5 of J. Stoer
%    and R. Bulirsch, Introduction to Numerical Analysis, 3rd ed. New York:
%    Springer, 2002.
%[2] R. L. Burden and J. D. Faires, Numerical Analysis, 9th ed. Boston, MA:
%    Brooks/Cole, 2011.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    y=y';

    Np=length(x);
    %The number of moments to match from position to velocity,
    %acceleration, etc.
    numMatch=size(y,2);
    numCoeff=numMatch*Np;
    
    %To make indexing simple, z is built as described in the book.
    z=dupEls(x,numMatch);
    
    %This will hold the coefficients for the polynomial in Newton's form.
    %We will later change this to the expanded form that it returned.
    fCoeffs=zeros(numCoeff,1);
        
    %Fill in the entries for the 0th divided difference level.
    DDList=dupEls(y(:,1),numMatch);
    fCoeffs(1)=DDList(1);%The first coefficient is just the first point.
    
    %Calculate the entries for the each subsequent divided difference 
    %level. This overwrites the old values. When z(curEl+level)==z(curEl),
    %the appropriate given derivative must be inserted. curDoty and
    %shouldIncDoty are used to keep track of when the variable should be
    %incremented.
    for level=1:(numCoeff-1)
        curDoty=1;
        shouldIncDoty=false;
        for curEl=1:(numCoeff-level)
            if(z(curEl+level)==z(curEl))
                if(shouldIncDoty==true)
                    curDoty=curDoty+1;
                    shouldIncDoty=false;
                end
                %If they are equal, then the appropriate derivative should be
                %put in.
                DDList(curEl)=y(curDoty,level+1)/factorial(level);
            else%Otherwise, perform differencing.
                DDList(curEl)=(DDList(curEl+1)-DDList(curEl))/(z(curEl+level)-z(curEl));
                shouldIncDoty=true;
            end
        end
        fCoeffs(level+1)=DDList(1);%Extract the coefficient.
    end
    
    a=fCoeffs';
    %The last z value is not used in Newton's form of the interpolating
    %polynomial.
    c=z(1:(end-1));
    c=c(:).';
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
