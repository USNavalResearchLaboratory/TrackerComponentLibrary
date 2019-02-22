function [out1,out2,out3]=LagrangeInterpPoly(x,y)
%%LAGRANGEINTERPPOLY  Obtain either the Lagrange interpolating polynomial
%                     coefficients in x that match the given (x,y) pairs
%                     along with its first and second derivatives, if two
%                     inputs are provided, or obtain matrices Li, dLi, and
%                     d2Li for a set of points in x such that Li*yArb, for
%                     some arbitrary yArb provides the coefficients of an
%                     interpolating polynomial matching the points (x,yArb)
%                     and dLi*yArb and d2Li*yArd are the coefficients of
%                     interpolating polynomials for the first and second
%                     derivatives.
%
%INPUTS: x  An NX1 set of N-1D points at which the interpolating
%           polynomial should match given points in the second variable.
%        y  An optional NX1 parameter. If included, out1 is the
%           coefficients of a Lagrange interpolating polynomial matching
%           (x,y) and out2 is the coefficients for the derivatives of the
%           polynomial. Otherwise, if matrices are desired as outputs, this
%           parameter should be omitted.
%
%OUTPUTS: out1 If y is provided, this is the NX1 set of coefficients for
%              the Lagrange interpolating polynomial matching the values
%              y at the points x. If y is omitted, this is the NXN matrix
%              Li such that Li*y for an arbitrary y will provide the
%              Lagrange interpolating polynomial.
%         out2 If y is provided, this it the set of coefficients for the
%              derivative of the Langrange interpolating polynomial
%              matching values of y at the points x. If N-1=0, then this is
%              just zero. If y is omitted, this is the (N-1)XN matrix dLi
%              such that dLi*y for an arbitrary y will provide the
%              coefficients of the derivative of the Lagrange interpolating
%              polynomial. If N-1=0, then this is a 1XN vector of zeros.
%         out3 If y is provided, this is the set of coefficients 
%              for the second derivative of the Langrange interpolating
%              polynomial matching values of y at the points x. If N-2=0,
%              then this is just zero. If y is omitted, this is the (N-2)XN
%              matrix d2Li such that d2Li*y for an arbitrary y will provide
%              the coefficients of the derivative of the Lagrange
%              interpolating polynomial. If N-2=0, then this is a 1XN
%              vector of zeros.
%
%Lagrange Interpolating polynomials are described in Chapter 3.1 of [1].
%The interpolating polynomial at a point xp matching the values in y at the
%points x is given by
%P(x)=sum_{i=1}^N L_i(x) *y(i)
%where
%L_i(x)=prod_{k=1,k\neq i}^N(x-xp(k))/(xp(i)-xp(k))
%The derivative of the interpolating polynomial is the derivative of the
%above expression. That is,
%dP(x)=sum_{i=1}^N dL_i(x) *y(i)
%where
%dL_i(x)=\sum_{k=1,k\neq i}^N(1/(xp(i)-xp(k)))*\prod_{j=1,j\neq i,j\neq k}^N (x-xp(j))/(xp(i)-xp(j))
%Similarly, the second derivative is
%d2P(x)=sum_{i=1}^N d2L_i(x) *y(i)
%where
%d2L_i(x)=\sum_{k=1,k\neq i}^N(1/(xp(i)-xp(k)))*\sum_{n=1,n\neq i,n\neq k}^N(1/(xp(i)-xp(n)))*\prod_{j=1,j\neq i,j\neq k,j\neq n}^N (x-xp(j))/(xp(i)-xp(j))
%
%If one wishes to evaluate the interpolating polynomial at another point z,
%it can be done as
% interpValZ=polyval(LagrangeInterpPoly(x,y),z);
%On the other hand, if matrices are returned, one could use
% Li=LagrangeInterpPoly(x);
% interpPoly=Li*y;
% interpValZ=polyval(interpPoly,z);
%
%REFERENCES:
%[1] R. L. Burden and J. D. Faires, Numerical Analysis, 9th ed. Boston,
%    MA: Brooks/Cole, 2011.
%
%April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=length(x);

%Li(:,i) is the set of coefficients (in x) for the ith Lagrange polynomial.
%The Lagrange polynomials in x must be multiplied by the points y to get
%the interpolating polynomial that matches the points y for a given x.
Li=zeros(N,N);
for curOrder=1:N
    sel=((1:N)~=curOrder);
    %The poly function computes the coefficients for the polynomial with the
    %specified roots
    numPoly=poly(x(sel))';
    %polyval evaluates the polynomial at the specified point.
    Li(:,curOrder)=numPoly./polyval(numPoly,x(curOrder));
end

if(nargin>1)
    out1=Li*y;
else
    out1=Li;
end

if(nargout>1&&nargin<2)
    %dLi(:,i) is the set of coefficients (in x) that corresponds to the
    %derivative of the ith Lagrange polynomial. After multiplying by y, one
    %obtains the derivative of the interpolating polynomial evaluated at y.
    
    if(N-1<=0)
        out2=zeros(1,N);
        out3=zeros(1,N);
        return;
    end
    
    dLi=zeros(N-1,N);
    for curOrder=1:N
        for k=1:N
            if(k~=curOrder)
                sel=((1:N)~=curOrder)&((1:N)~=k);
                numPoly=poly(x(sel))';
                dLi(:,curOrder)=dLi(:,curOrder)+1/(x(curOrder)-x(k))*numPoly./polyval(numPoly,x(curOrder));
            end
        end
    end
    out2=dLi;
    
    if(nargout>2)
        if(N-2<=0)
            out2=zeros(1,N);
            out3=zeros(1,N);
            return;
        end
        
        %d2Li(:,i) is the set of coefficients (in x) that corresponds to
        %the second derivative of the ith Lagrange polynomial. After 
        %multiplying by y, one obtains the derivative of the interpolating
        %polynomial evaluated at y.
        d2Li=zeros(N-2,N);
        for curOrder=1:N
            for k=1:N
                if(k~=curOrder)
                    for n=1:N
                        if(n~=k&&n~=curOrder)
                            sel=((1:N)~=curOrder)&((1:N)~=k)&((1:N)~=n);
                            numPoly=poly(x(sel))';

                            d2Li(:,curOrder)=d2Li(:,curOrder)+(1/(x(curOrder)-x(k)))*(1/(x(curOrder)-x(n)))*numPoly/polyval(numPoly,x(curOrder));
                        end
                    end
                end
            end
        end
        out3=d2Li;
    end
elseif(nargout>1&&nargin>1)
    %There is no need to compute the dLi matrix if it does not have to be
    %returned. Just directly differentiate the coefficients.
    out2=polyder(out1);
    
    if(nargout>2)
        %If the second derivative equation is desired, then just
        %differentiate the coefficients again.
        out3=polyder(out2);
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
