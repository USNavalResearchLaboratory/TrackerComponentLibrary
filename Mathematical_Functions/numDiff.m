function J=numDiff(x,f,fDim,N,epsilon)
%%NUMDIFF Numerical differentiation using central difference formulae based
%         on Lagrange interpolating polynomials. The algorithm is made to
%         alert (and not crash) in the event that f fails. Specifically, if
%         f returns an empty matrix or it any of the components of f are
%         NaNs, then the numDiff function will terminate early, returning
%         an empty matrix to indicate failure.          
%
%INPUTS: x The xDimX1 vector or scalar point at which the derivative of the
%          (possibly vector) function is desired.
%        f The scalar or vector function that is to be differentiated. The
%          function f must take x as its parameter and its output should be
%          a scalar or a column vector.
%     fDim The dimensionality of the output of f.
%        N An integer >=1 specifying the order of the derivative
%          approximation. N+1 is the order of the error terms. Thus, N=1
%          means that the error terms scale as O(N^2). If N is omitted
%          or an empty matrix is passed, then N=2 is assumed. Values for
%          N=1 through 8 are explicitly coded in. For values 9 and above,
%          the coefficients of the derivative of the Lagrange interpolating
%          polynomial are solved.
%  epsilon A scalar or xDimX1 vector quantity specifying the finite step
%          size used for numerical differentiation. If a scalar value is
%          given, that value is used for differentiating with respect to
%          elements of x. If an xDimX1 value is given, then the
%          corresponding element of epsilon is used to differentiate each
%          element of x. If epsilon is omitted or an empty matrix is
%          passed, then epsilon=max(1e-5*x,1e-7); is used.
%
%OUTPUTS: J An fDimXxDim Jacobian matrix. Each column is the derivative
%           vector of f with respect to the corresponding element of x. If
%           at any point the function f returned a NaN or an empty matrix,
%           J will be an empty matrix.
%
%Central-difference numerical differentiation is discussed in terms of
%Lagrange interpolating polynomials in Chapter 4.1 of [1]. The Lagrange
%interpolating polynomials themselves are discussed in Chapter 3 of the
%book. In general, the coefficients of the interpolating polynomial for
%central difference numerical differentiation come from Equation 4.2,
%which expresses them in terms of the first derivative of a Lagrange
%interpolating polynomial.
%
%REFERENCES:
%[1] R. L. Burden and J. D. Faires, Numerical Analysis, 9th ed. Boston, MA:
%    Brooks/ Cole, 2011.
%
%January 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(N))
    N=2;
end

xDim=size(x,1);
if(nargin<5||isempty(epsilon))
    %If epsilon is not specified, then use some ad-hoc default value.
    epsilon=max(1e-5*abs(x),1e-7);
end
if(isscalar(epsilon))
   epsilon=repmat(epsilon,[xDim,1]); 
end

switch(N)
    case 1
        a=1;
        d=2;
    case 2
        a=[8;-1];
        d=12;
    case 3
        a=[45;-9;1];
        d=60;
    case 4
        a=[672;-168;32;-3];
        d=840;
    case 5
        a=[2100;-600;150;-25;2];
        d=2520;
    case 6
        a=[23760;-7425;2200;-495;72;-5];
        d=27720;
    case 7
        a=[315315;-105105;35035;-9555;1911;-245;15];
        d=360360;
    case 8
        a=[640640;-224224;81536;-25480;6272;-1120;128;-7];
        d=720720;
    otherwise
        [~,dLi]=LagrangeInterpPoly(-N:1:N);
        a=-dLi(end,N:-1:1)';
        d=1;
end

numP=length(a);

J=zeros(fDim,xDim);
for curEl=1:xDim
    epsCur=epsilon(curEl);
    
    for curP=1:numP
        xP=x;
        xP(curEl)=xP(curEl)+curP*epsCur;
        fxP=f(xP);
        if(isempty(fxP)||any(isnan(fxP)))
            J=[];
            return;
        end
        J(:,curEl)=J(:,curEl)+a(curP)*fxP;
        xP=x;
        xP(curEl)=xP(curEl)-curP*epsCur;
        fxP=f(xP);
        if(isempty(fxP)||any(isnan(fxP)))
            J=[];
            return;
        end
        J(:,curEl)=J(:,curEl)-a(curP)*fxP;
    end
    J(:,curEl)=J(:,curEl)/(d*epsCur);
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
